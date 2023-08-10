pub mod ncbi_compress {
    use std::ops::AddAssign;
    use std::path::Path;
    use std::{borrow::BorrowMut, fs};

    use bio::io::fasta;
    use rayon::prelude::*;
    use sourmash::encodings::HashFunctions;
    use sourmash::errors::SourmashError;
    use sourmash::signature::SigsTrait;
    use sourmash::sketch::minhash::KmerMinHash;
    use tempdir::TempDir;
    use trie_rs::TrieBuilder;

    use crate::logging::logging;
    use crate::minhashtree::minhashtree::{MinHashTree, MinHashTreeFunctionality};
    use crate::minhashtree_w_logging::minhashtree_w_logging::{MinHashTreeWithLogging, MinHashTreeWithLoggingFunctionality};
    use crate::trie_store::trie_store::TrieStoreBuilder;
 
    pub fn containment(needle: &KmerMinHash, haystack: &KmerMinHash) -> Result<f64, SourmashError> {
        let (intersect_size, _) = needle.intersection_size(haystack)?;
        Ok(intersect_size as f64 / needle.mins().len() as f64)
    }

    fn remove_accession_version(accession: &str) -> &str {
        accession.splitn(2, |c| c == '.').next().unwrap()
    }

    pub fn split_accessions_by_taxid<P: AsRef<Path> + std::fmt::Debug, Q: AsRef<Path> + std::fmt::Debug>(
        input_fasta_path: P,
        mapping_file_path: Vec<Q>,
        taxids_to_drop: &Vec<u64>,
        temp_file_output_dir: &str

    ) -> TempDir {
        // create a temp dir containing one file per taxid that input fasta accessions are sorted into
        // based on taxid in the mapping files (input fasta does not have taxid in header)
        log::info!("Creating accession to taxid mapping");
        let taxid_dir = TempDir::new_in(temp_file_output_dir, "accessions_by_taxid").unwrap();
        let reader = fasta::Reader::from_file(&input_fasta_path).unwrap();
        // Build a trie of the accessions in the input fasta
        let mut builder = TrieBuilder::new(); // is this supposed to be TrieStoreBuilder?
        reader.records().enumerate().for_each(|(i, result)| {
            let record = result.unwrap();
            let accession_id = record.id().split_whitespace().next().unwrap();
            let accession_no_version = remove_accession_version(accession_id);
            builder.push(accession_no_version);
            if i % 10_000 == 0 {
                log::info!("  Processed {} accessions", i);
            }
        });
        log::info!(" Started building accession trie");
        let accessions_trie = builder.build();
        log::info!(" Finished building accession trie");

        // Build a trie of the accessions in the mapping files
        let mut builder = TrieStoreBuilder::new();
        mapping_file_path.iter().for_each(|mapping_file_path| {
            log::info!(" Processing mapping file {:?}", mapping_file_path);
            let reader = csv::ReaderBuilder::new()
                .delimiter(b'\t')
                .from_path(mapping_file_path)
                .unwrap();
            let mut added = 0;
            reader.into_records().enumerate().for_each(|(i, result)| {
                if i % 10_000 == 0 {
                    log::info!("  Processed {} mappings, added {}", i, added);
                }

                let record = result.unwrap();
                let accession = &record[0];
                let accession_no_version = remove_accession_version(accession);

                // Only output mappings if the accession is in the source fasta file
                if !accessions_trie.exact_match(accession_no_version) {
                    return;
                }

                // If using the prot.accession2taxid.FULL file
                let taxid = if record.len() < 3 {
                    // The taxid will be at index 1
                    record[1].parse::<u64>().unwrap()
                } else {
                    // Otherwise there is a versionless accession ID at index 0 and the taxid is at index 2
                    record[2].parse::<u64>().unwrap()
                };

                if !taxids_to_drop.contains(&taxid) {
                    added += 1;
                    builder.push(accession_no_version, taxid);
                }
            });
            log::info!(
                " Finished Processing mapping file {:?}, added {} mappings",
                mapping_file_path,
                added
            );
        });
        log::info!(" Started building accession to taxid trie");
        let accession_to_taxid = builder.build();
        log::info!("Finished building accession to taxid trie");

        log::info!("Splitting accessions by taxid");
        let reader = fasta::Reader::from_file(&input_fasta_path).unwrap();
        // Split the input fasta accessions into one file per taxid
        for (i, record) in reader.records().enumerate() {
            if i % 10_000 == 0 {
                log::info!("  Split {} accessions", i);
            }
            let record = record.unwrap();
            let accession_id = record.id().split_whitespace().next().unwrap();
            let accession_no_version = remove_accession_version(accession_id);
            let taxid = if let Some(taxid) = accession_to_taxid.get(accession_no_version) {
                taxid
            } else {
                continue;
            };

            let file_path = taxid_dir.path().join(format!("{}.fasta", taxid));
            let file = fs::OpenOptions::new()
                .create(true)
                .append(true)
                .open(file_path)
                .unwrap();
            let mut writer = fasta::Writer::new(file);
            writer.write_record(&record).unwrap();
        }
        log::info!("Finished splitting accessions by taxid");
        taxid_dir
    }

    pub fn fasta_compress_taxid_w_logging<P: AsRef<Path> + std::fmt::Debug>(
        input_fasta_path: P,
        writer: &mut fasta::Writer<std::fs::File>,
        scaled: u64,
        k: u32,
        seed: u64,
        similarity_threshold: f64,
        chunk_size: usize,
        branch_factor: usize,
        accession_count: &mut u64,
        unique_accession_count: &mut u64,
        logging_contained_in_tree_fn: &str,
        logging_contained_in_chunk_fn: &str,
    ) {
        // take in a fasta file and output a fasta file with only unique accessions (based on similarity threshold)
        let reader = fasta::Reader::from_file(&input_fasta_path).unwrap();
        let mut tree = MinHashTreeWithLogging::new(branch_factor);

        let mut records_iter = reader.records();
        let mut chunk = records_iter
            .borrow_mut()
            .take(chunk_size)
            .collect::<Vec<_>>();
        while chunk.len() > 0 {
            let chunk_signatures = chunk
                .par_iter()
                .filter_map(|r| {
                    let record = r.as_ref().unwrap();
                    let mut hash =
                        KmerMinHash::new(scaled, k, HashFunctions::murmur64_DNA, seed, false, 0);
                    hash.add_sequence(record.seq(), true).unwrap();
                    // Run an initial similarity check here against the full tree, this is slow so we can parallelize it
                    let contained_in_tree = tree.contains(&hash, similarity_threshold);
                    match contained_in_tree {
                        Ok(None) => {
                            // If the tree doesn't contain the hash, we need to insert it
                            Some((hash, record))
                        }
                        Ok(Some(found_accessions)) => {
                            // log when tree already contains hash
                            let accession_id = record.id().split_whitespace().next().unwrap();
                            let found_accession_ids: Vec<String> = found_accessions.iter().map(|(accession_id, _)| accession_id.to_string()).collect::<Vec<_>>();
                            let found_accession_containments: Vec<String> = found_accessions.iter().map(|(_, containment)| containment.to_string()).collect::<Vec<_>>();
                            logging::write_to_file(
                                // log discarded, retained, containment
                                vec![accession_id, &found_accession_ids.join(", "), &found_accession_containments.join(", ")],
                                logging_contained_in_tree_fn,
                            );
                            None
                        }
                        Err(e) => {
                            // If there was an error, we need to log it and continue
                            log::error!("Error checking similarity: {}", e);
                            None
                        }
                    }
                })
                .collect::<Vec<_>>();

            let mut tmp: Vec<(KmerMinHash, &str)> = Vec::with_capacity(chunk_signatures.len() / 2); // initialize with a guess of what you think the size will be
            for (hash, record) in chunk_signatures {
                let accession_id = record.id().split_whitespace().next().unwrap();
                accession_count.add_assign(1); // += 1, used for logging
                // Perform a faster similarity check over just this chunk because we may have similarities within a chunk
                let similar_seqs = tmp
                    .par_iter()
                    .filter_map(|(other, accession_id)| {
                        let containment_value = containment(&hash, &other).unwrap();
                        if containment_value >= similarity_threshold {
                            Some((accession_id.to_string(), containment_value.to_string()))
                        } else {
                            None
                        }
                    }).collect::<Vec<_>>();

                let similar = !similar_seqs.is_empty();
                if !similar {
                    unique_accession_count.add_assign(1);
                    tmp.push((hash, accession_id));
                    writer.write_record(record).unwrap();

                    if *unique_accession_count % 10_000 == 0 {
                        log::info!(
                            "Processed {} accessions, {} unique",
                            accession_count,
                            unique_accession_count
                        );
                    }
                } else {
                    let similar_accession_ids: Vec<String> = similar_seqs.iter().map(|(accession_id, _)| accession_id.to_string()).collect::<Vec<_>>();
                    let similar_accession_containments: Vec<String> = similar_seqs.iter().map(|(_, containment)| containment.to_string()).collect::<Vec<_>>();
                    logging::write_to_file(
                        vec![accession_id, &similar_accession_ids.join(","), &similar_accession_containments.join(",")],
                        logging_contained_in_chunk_fn,
                    );
                }
            }
            for (hash, accession_id) in tmp {
                // only insert accession_id into tree if logging is enabled

                tree.insert(hash, accession_id).unwrap();
            }
            chunk = records_iter
                .borrow_mut()
                .take(chunk_size)
                .collect::<Vec<_>>();
        }
    }

    fn fasta_compress_taxid<P: AsRef<Path> + std::fmt::Debug>(
        input_fasta_path: P,
        writer: &mut fasta::Writer<std::fs::File>,
        scaled: u64,
        k: u32,
        seed: u64,
        similarity_threshold: f64,
        chunk_size: usize,
        branch_factor: usize,
        accession_count: &mut u64,
        unique_accession_count: &mut u64,
    ) {
        let reader = fasta::Reader::from_file(&input_fasta_path).unwrap();
        let mut tree = MinHashTree::new(branch_factor);
    
        let mut records_iter = reader.records();
        let mut chunk = records_iter
            .borrow_mut()
            .take(chunk_size)
            .collect::<Vec<_>>();
        while chunk.len() > 0 {
            let chunk_signatures = chunk
                .par_iter()
                .filter_map(|r| {
                    let record = r.as_ref().unwrap();
                    let mut hash =
                        KmerMinHash::new(scaled, k, HashFunctions::murmur64_DNA, seed, false, 0);
                    hash.add_sequence(record.seq(), true).unwrap();
                    // Run an initial similarity check here against the full tree, this is slow so we can parallelize it
                    if tree.contains(&hash, similarity_threshold).unwrap() {
                        None
                    } else {
                        Some((record, hash))
                    }
                })
                .collect::<Vec<_>>();
    
            let mut tmp = Vec::with_capacity(chunk_signatures.len() / 2);
            for (record, hash) in chunk_signatures {
                accession_count.add_assign(1);
                // Perform a faster similarity check over just this chunk because we may have similarities within a chunk
                let similar = tmp
                    .par_iter()
                    .any(|other| containment(&hash, other).unwrap() >= similarity_threshold);
    
                if !similar {
                    unique_accession_count.add_assign(1);
                    tmp.push(hash);
                    writer.write_record(record).unwrap();
    
                    if *unique_accession_count % 10_000 == 0 {
                        log::info!(
                            "Processed {} accessions, {} unique",
                            accession_count,
                            unique_accession_count
                        );
                    }
                }
            }
            for hash in tmp {
                tree.insert(hash).unwrap();
            }
            chunk = records_iter
                .borrow_mut()
                .take(chunk_size)
                .collect::<Vec<_>>();
        }
    }

    pub fn fasta_compress<P: AsRef<Path> + std::fmt::Debug>(
        input_fasta_path: P,
        accession_mapping_files: Vec<P>,
        output_fasta_path: P,
        taxids_to_drop: Vec<u64>,
        scaled: u64,
        k: u32,
        seed: u64,
        similarity_threshold: f64,
        chunk_size: usize,
        branch_factor: usize,
        temp_file_output_dir: &str,
        skip_split_by_taxid: bool,
        enable_sequence_retention_logging: bool,
        logging_contained_in_tree_fn: &str,
        logging_contained_in_chunk_fn: &str

    ) {
        let mut accession_count = 0;
        let mut unique_accession_count = 0;
        let mut writer = fasta::Writer::to_file(output_fasta_path).unwrap();

        if skip_split_by_taxid {
            log::info!("Skipping splitting accessions by taxid");
            if enable_sequence_retention_logging {
                fasta_compress_taxid_w_logging(
                    input_fasta_path,
                    &mut writer,
                    scaled,
                    k,
                    seed,
                    similarity_threshold,
                    chunk_size,
                    branch_factor,
                    &mut accession_count,
                    &mut unique_accession_count,
                    logging_contained_in_tree_fn,
                    logging_contained_in_chunk_fn,
                );
            } else {
                fasta_compress_taxid(
                    input_fasta_path,
                    &mut writer,
                    scaled,
                    k,
                    seed,
                    similarity_threshold,
                    chunk_size,
                    branch_factor,
                    &mut accession_count,
                    &mut unique_accession_count,
                );
            }
            if accession_count % 10_000 == 0 {
                log::info!(
                    "  Compressed {} accessions, {} uniqe accessions",
                    accession_count,
                    unique_accession_count
                );
            }
        } else {
            log::info!("Splitting accessions by taxid");
            let taxid_dir =
                split_accessions_by_taxid(
                    &input_fasta_path,
                    accession_mapping_files,
                    &taxids_to_drop,
                    temp_file_output_dir
                );
            log::info!("Finished splitting accessions by taxid");
            log::info!("Starting compression by taxid");
            for (i, entry) in fs::read_dir(taxid_dir.path()).unwrap().enumerate() {
                let entry = entry.unwrap();
                let path = entry.path();
                let input_fasta_path = path.to_str().unwrap();
                if enable_sequence_retention_logging {
                    fasta_compress_taxid_w_logging(
                        input_fasta_path,
                        &mut writer,
                        scaled,
                        k,
                        seed,
                        similarity_threshold,
                        chunk_size,
                        branch_factor,
                        &mut accession_count,
                        &mut unique_accession_count,
                        logging_contained_in_tree_fn,
                        logging_contained_in_chunk_fn,
                    );
                } else {
                    fasta_compress_taxid(
                        input_fasta_path,
                        &mut writer,
                        scaled,
                        k,
                        seed,
                        similarity_threshold,
                        chunk_size,
                        branch_factor,
                        &mut accession_count,
                        &mut unique_accession_count,
                    );
                }
                if i % 10_000 == 0 {
                    log::info!(
                        "  Compressed {} taxids, {} accessions, {} uniqe accessions",
                        i,
                        accession_count,
                        unique_accession_count
                    );
                }
            }
            taxid_dir.close().unwrap();
            log::info!(
                "Finished compression by taxid, {} accessions, {} uniqe accessions",
                accession_count,
                unique_accession_count
            );
        }
    }
}
