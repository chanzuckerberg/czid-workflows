pub mod ncbi_compress {
    use std::collections::HashMap;
    use std::ops::AddAssign;
    use std::path::Path;
    use std::{borrow::BorrowMut, fs};
    use std::sync::Mutex;

    use bio::io::fasta;
    use lazy_static::lazy_static;
    use rand::distributions::Alphanumeric;
    use rand::Rng;
    use rayon::prelude::*;
    use sourmash::encodings::HashFunctions;
    use sourmash::errors::SourmashError;
    use sourmash::signature::SigsTrait;
    use sourmash::sketch::minhash::KmerMinHash;

    use crate::logging::logging;
    use crate::minhashtree_w_logging::minhashtree_w_logging::{
        MinHashTreeWithLogging, MinHashTreeWithLoggingFunctionality,
    };

    pub fn containment(needle: &KmerMinHash, haystack: &KmerMinHash) -> Result<f64, SourmashError> {
        let (intersect_size, _) = needle.intersection_size(haystack)?;
        Ok(intersect_size as f64 / needle.mins().len() as f64)
    }

    fn remove_accession_version(accession: &str) -> &str {
        accession.splitn(2, |c| c == '.').next().unwrap()
    }

    // Mutex-protected HashMap for file handles
    // lazy_static! {
    //     static ref FILE_HANDLES: Mutex<HashMap<String, std::fs::File>> = Mutex::new(HashMap::new());
    // }

    // pub fn split_accessions_by_taxid(
    //     input_fasta_path: &str,
    //     mapping_file_path: Vec<String>,
    //     output_dir: &str,
    // ) -> std::path::PathBuf {
    //     // create a temp dir containing one file per taxid that input fasta accessions are sorted into
    //     // based on taxid in the mapping files (input fasta does not have taxid in header)

    //     log::info!("Splitting accessions by taxid");
    //     // Use a random string for the rocksdb directoru name
    //     // A tempdir can be removed by the operating system before we are done with it which
    //     // break the rocksdb. We need a random name so we can call this function multiple
    //     // times from different threads during testing.
    //     let rng = rand::thread_rng();
    //     let dir: String = rng
    //         .sample_iter(&Alphanumeric)
    //         .take(20)
    //         .map(char::from)
    //         .collect();
    //     let accession_to_taxid = rocksdb::DB::open_default(&dir).unwrap();
    //     fs::create_dir_all(&output_dir).expect("Error creating output directory");

    //     log::info!("Creating accession to taxid db");

    //     let taxid_path = Path::new(output_dir);
    //     log::info!("Creating taxid dir {:?}", taxid_path);
    //     let reader = fasta::Reader::from_file(&input_fasta_path).unwrap();
    //     let records: Vec<_> = reader
    //         .records()
    //         .enumerate()
    //         .collect();
    //     // Build a trie of the accessions in the input fasta
    //     records.par_iter().for_each(|(i, result)| {
    //         let record = result.as_ref().unwrap();
    //         let accession_id = record.id().split_whitespace().next().unwrap();
    //         let accession_no_version = remove_accession_version(accession_id);
    //         accession_to_taxid.put(accession_no_version, b"").unwrap();
    //         if i % 1_000_000 == 0 {
    //             log::info!("  Processed {} accessions", i);
    //         }
    //     });
    //     log::info!(" Finished loading accessions");

    //     mapping_file_path.par_iter().for_each(|mapping_file_path| {
    //         log::info!(" Processing mapping file {:?}", mapping_file_path);
    //         let reader = csv::ReaderBuilder::new()
    //             .delimiter(b'\t')
    //             .from_path(mapping_file_path)
    //             .unwrap();
    //         let mut added = 0;
    //         reader.into_records().enumerate().for_each(|(_i, result)| {
    //             // if i % 10_000 == 0 {
    //             //     log::info!("  Processed {} mappings, added {}", i, added);
    //             // }

    //             let record = result.unwrap();
    //             let accession = &record[0];
    //             let accession_no_version = remove_accession_version(accession);

    //             // Only output mappings if the accession is in the source fasta file
    //             if accession_to_taxid
    //                 .get(accession_no_version)
    //                 .unwrap()
    //                 .is_none()
    //             {
    //                 return;
    //             }

    //             // If using the prot.accession2taxid.FULL file
    //             let taxid = if record.len() < 3 {
    //                 // The taxid will be at index 1
    //                 record[1].parse::<u64>().unwrap()
    //             } else {
    //                 // Otherwise there is a versionless accession ID at index 0 and the taxid is at index 2
    //                 record[2].parse::<u64>().unwrap()
    //             };
    //             added += 1;
    //             // convert taxid to u64_vec
    //             accession_to_taxid
    //                 .put(accession_no_version, taxid.to_be_bytes())
    //                 .unwrap();
    //         });
    //         log::info!(
    //             " Finished Processing mapping file {:?}, added {} mappings",
    //             mapping_file_path,
    //             added
    //         );
    //     });
    //     log::info!("Finished building accession to taxid db");

    //     log::info!("Splitting accessions by taxid");
    //     let reader = fasta::Reader::from_file(&input_fasta_path).unwrap();
    //     // Split the input fasta accessions into one file per taxid
    //     let records: Vec<_> = reader
    //         .records()
    //         .enumerate()
    //         .collect();
    //     records.par_iter().for_each(|(i, record)| {
    //         if i % 10_000 == 0 {
    //             log::info!("  Split {} accessions", i);
    //         }
    //         let record = record.as_ref().unwrap();
    //         let accession_id = record.id().split_whitespace().next().unwrap();
    //         let accession_no_version = remove_accession_version(accession_id);
    //         let taxid = if let Some(taxid) = accession_to_taxid.get(accession_no_version).unwrap() {
    //             if taxid.len() == 0 {
    //                 0 // no taxid found
    //             } else {
    //                 u64::from_be_bytes(taxid.as_slice().try_into().unwrap())
    //             }
    //         } else {
    //             0 // no taxid found
    //         };
    //         let mut handles = FILE_HANDLES.lock().unwrap(); // Lock the mutex here
    //         // Check if the HashMap size exceeds the threshold
    //         if handles.len() >= 2000 { // 65k is the max number of open files
    //             if let Some(key_to_remove) = handles.keys().next().cloned() {
    //                 if let Some(file) = handles.remove(&key_to_remove) {
    //                     drop(file); // Close the file
    //                 }
    //             }
    //         }

    //         let file_path = format!("{}/{}.fasta", output_dir, taxid);
    //         let file = handles.entry(file_path.clone()).or_insert_with(|| {
    //             fs::OpenOptions::new()
    //                 .write(true)
    //                 .append(true)
    //                 .create(true)
    //                 .open(file_path)
    //                 .expect("Error opening fasta file")
    //         });

    //         let mut writer = fasta::Writer::new(file);
    //         writer.write_record(&record).unwrap();
    //     });
    //     fs::remove_dir_all(dir).expect("Error removing db");
    //     log::info!("Finished splitting accessions by taxid");
    //     taxid_path.to_path_buf()
    // }

    pub fn split_accessions_by_taxid(
        input_fasta_path: &str,
        mapping_file_path: Vec<String>,
        output_dir: &str,
    ) -> std::path::PathBuf {
        // create a temp dir containing one file per taxid that input fasta accessions are sorted into
        // based on taxid in the mapping files (input fasta does not have taxid in header)

        log::info!("Splitting accessions by taxid");
        // Use a random string for the rocksdb directoru name
        // A tempdir can be removed by the operating system before we are done with it which
        // break the rocksdb. We need a random name so we can call this function multiple
        // times from different threads during testing.
        let rng = rand::thread_rng();
        let dir: String = rng
            .sample_iter(&Alphanumeric)
            .take(20)
            .map(char::from)
            .collect();
        let accession_to_taxid = rocksdb::DB::open_default(&dir).unwrap();
        fs::create_dir_all(&output_dir).expect("Error creating output directory");

        log::info!("Creating accession to taxid db");

        //let taxid_dir = TempDir::new_in(temp_file_output_dir, "accessions_by_taxid").unwrap();
        // let taxid_path_str = format!("{}/accessions_by_taxid", temp_file_output_dir);
        let taxid_path = Path::new(output_dir);
        log::info!("Creating taxid dir {:?}", taxid_path);
        let reader = fasta::Reader::from_file(&input_fasta_path).unwrap();
        // Build a trie of the accessions in the input fasta
        reader.records().enumerate().for_each(|(i, result)| {
            let record = result.unwrap();
            let accession_id = record.id().split_whitespace().next().unwrap();
            let accession_no_version = remove_accession_version(accession_id);
            accession_to_taxid.put(accession_no_version, b"").unwrap();
            if i % 1_000_000 == 0 {
                log::info!("  Processed {} accessions", i);
            }
        });
        log::info!(" Finished loading accessions");

        mapping_file_path.par_iter().for_each(|mapping_file_path| {
            log::info!(" Processing mapping file {:?}", mapping_file_path);
            let reader = csv::ReaderBuilder::new()
                .delimiter(b'\t')
                .from_path(mapping_file_path)
                .unwrap();
            let mut added = 0;
            reader.into_records().enumerate().for_each(|(_i, result)| {
                // if i % 10_000 == 0 {
                //     log::info!("  Processed {} mappings, added {}", i, added);
                // }

                let record = result.unwrap();
                let accession = &record[0];
                let accession_no_version = remove_accession_version(accession);

                // Only output mappings if the accession is in the source fasta file
                if accession_to_taxid
                    .get(accession_no_version)
                    .unwrap()
                    .is_none()
                {
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
                added += 1;
                // convert taxid to u64_vec
                accession_to_taxid
                    .put(accession_no_version, taxid.to_be_bytes())
                    .unwrap();
            });
            log::info!(
                " Finished Processing mapping file {:?}, added {} mappings",
                mapping_file_path,
                added
            );
        });
        log::info!("Finished building accession to taxid db");

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
            let taxid = if let Some(taxid) = accession_to_taxid.get(accession_no_version).unwrap() {
                if taxid.len() == 0 {
                    0 // no taxid found
                } else {
                    u64::from_be_bytes(taxid.as_slice().try_into().unwrap())
                }
            } else {
                0 // no taxid found
            };
            let file_path = taxid_path.join(format!("{}.fasta", taxid));
            let file = fs::OpenOptions::new()
                .create(true)
                .append(true)
                .open(file_path)
                .unwrap();
            let mut writer = fasta::Writer::new(file);
            writer.write_record(&record).unwrap();
        }
        fs::remove_dir_all(dir).expect("Error removing db");
        log::info!("Finished splitting accessions by taxid");
        taxid_path.to_path_buf()
    }

    pub fn fasta_compress_taxid_w_logging(
        input_fasta_path: &str,
        writer: &mut fasta::Writer<std::fs::File>,
        scaled: u64,
        k: u32,
        seed: u64,
        similarity_threshold: f64,
        chunk_size: usize,
        branch_factor: usize,
        is_protein_fasta: bool,
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
                    let mut hash;
                    if is_protein_fasta {
                        hash = KmerMinHash::new(
                            scaled,
                            k,
                            HashFunctions::murmur64_protein,
                            seed,
                            false,
                            0,
                        );
                        hash.add_protein(record.seq()).unwrap();
                    } else {
                        hash = KmerMinHash::new(
                            scaled,
                            k,
                            HashFunctions::murmur64_DNA,
                            seed,
                            false,
                            0,
                        );
                        hash.add_sequence(record.seq(), true).unwrap();
                    }
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
                            let found_accession_ids: Vec<String> = found_accessions
                                .iter()
                                .map(|(accession_id, _)| accession_id.to_string())
                                .collect::<Vec<_>>();
                            let found_accession_containments: Vec<String> = found_accessions
                                .iter()
                                .map(|(_, containment)| containment.to_string())
                                .collect::<Vec<_>>();
                            logging::write_to_file(
                                // log discarded, retained, containment
                                vec![
                                    accession_id,
                                    &found_accession_ids.join(", "),
                                    &found_accession_containments.join(", "),
                                ],
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
                    })
                    .collect::<Vec<_>>();

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
                    let similar_accession_ids: Vec<String> = similar_seqs
                        .iter()
                        .map(|(accession_id, _)| accession_id.to_string())
                        .collect::<Vec<_>>();
                    let similar_accession_containments: Vec<String> = similar_seqs
                        .iter()
                        .map(|(_, containment)| containment.to_string())
                        .collect::<Vec<_>>();
                    logging::write_to_file(
                        vec![
                            accession_id,
                            &similar_accession_ids.join(","),
                            &similar_accession_containments.join(","),
                        ],
                        logging_contained_in_chunk_fn,
                    );
                }
            }
            for (hash, accession_id) in tmp {
                tree.insert(hash, accession_id).unwrap();
            }
            chunk = records_iter
                .borrow_mut()
                .take(chunk_size)
                .collect::<Vec<_>>();
        }
    }

    // struct KMerMinHashWrapper {
    //     sketch: KmerMinHash,
    // }

    // impl BloomSetTreeable for KMerMinHashWrapper {
    //     fn hashes(&self) -> Vec<u64> {
    //         self.sketch.mins()
    //     }

    //     fn containment(&self, haystack: &Self) -> f64 {
    //         let (common, _) = self.sketch.intersection_size(&haystack.sketch).unwrap();
    //         common as f64 / usize::max(1, self.sketch.size()) as f64
    //     }
    // }
    //

    pub fn fasta_compress_taxid(
        input_fasta_path: &str,
        writer: &mut fasta::Writer<std::fs::File>,
        scaled: u64,
        k: u32,
        seed: u64,
        similarity_threshold: f64,
        chunk_size: usize,
        _branch_factor: usize,
        is_protein_fasta: bool,
        accession_count: &mut u64,
        unique_accession_count: &mut u64,
    ) {
        let reader = fasta::Reader::from_file(&input_fasta_path).unwrap();
        // let mut tree: BloomSetTree<_, 4096> = BloomSetTree::new(branch_factor);
        let mut sketches: Vec<KmerMinHash> = Vec::new();
        let mut records_iter = reader.records();

        // Initialize a temporary vector to store the unique items from each chunk
        let mut unique_in_chunk: Vec<(KmerMinHash, fasta::Record)> = Vec::with_capacity(chunk_size);
        let mut unique_in_tree_and_chunk: Vec<(KmerMinHash, &fasta::Record)> = Vec::with_capacity(chunk_size);

        loop {
            let chunk = records_iter
                .borrow_mut()
                .take(chunk_size)
                .collect::<Vec<_>>();

            accession_count.add_assign(chunk.len() as u64);

            if chunk.len() == 0 {
                break;
            }

            // create signatures for each record in the chunk
            let chunk_signatures = chunk
            .par_iter()
            .map(|r| {
                let record = r.as_ref().unwrap();
                let mut hash;
                if is_protein_fasta {
                    hash = KmerMinHash::new(
                        scaled,
                        k,
                        HashFunctions::murmur64_protein,
                        seed,
                        false,
                        0,
                    );
                    hash.add_protein(record.seq()).unwrap();
                    (hash, record.clone())
                } else {
                    hash = KmerMinHash::new(
                        scaled,
                        k,
                        HashFunctions::murmur64_DNA,
                        seed,
                        false,
                        0,
                    );
                    hash.add_sequence(record.seq(), true).unwrap();
                    (hash, record.clone())
                }
            }).collect::<Vec<_>>();

            // we need to make sure records within the chunk arn't similar to each other before
            // we check them against the larger tree
            for (hash, record) in chunk_signatures {
                let similar = unique_in_chunk
                    .par_iter()
                    .any(|(other, _record)| containment(&hash, &other).unwrap() >= similarity_threshold);

                if !similar {
                    unique_in_chunk.push((hash, record));
                }
            }

            // Check if any of unique hashes in the chunk have similarity to the ones already in the tree
            let mut unique_in_tree_and_chunk = unique_in_chunk
                .par_iter()
                .filter_map(|(hash, record)| {
                    // Take a chunk of accessions and return the ones that are not already in sketches

                    // Do a parallel search over the sketches to see if any are similar don't
                    // return the item if you find a similar item. We are implementing this
                    // manually instead of using the linear index built into sourmash because
                    // we want to give up after finding the first hit while the linear index
                    // will return them all.
                    if sketches
                        .par_iter()
                        .find_any(|other| {
                            containment(&hash, other).unwrap() >= similarity_threshold
                        })
                        .is_some()
                    {
                        None
                    } else {
                        Some((hash.clone(), record))
                    }
                })
                .collect::<Vec<_>>();

            for (hash, record) in unique_in_tree_and_chunk.iter() {
                unique_accession_count.add_assign(1);
                if *unique_accession_count % 1_000_000 == 0 {
                    log::info!(
                        "Processed {} accessions, {} unique",
                        accession_count,
                        unique_accession_count
                    );
                }

                writer.write_record(record).unwrap();
                sketches.push(hash.clone());
            }

            // need to clear the vectors in this order because of the borrow checker
            unique_in_tree_and_chunk.clear();
            unique_in_chunk.clear();

            // for hash in tmp {
            //     tree.insert(hash);
            // }

        }
    }
}

#[cfg(test)]
mod tests {
    use tempfile::tempdir;

    use crate::ncbi_compress::ncbi_compress;
    use crate::util::util;
    use std::fs;

    #[test]
    fn test_split_accessions_by_taxid() {
        let input_fasta_path = "test_data/fasta_tools/inputs/nt";
        let mapping_files_directory =
            "test_data/ncbi_compress/split_accessions_by_taxid/inputs/accession2taxid";
        let input_mapping_file_paths = vec![
            format!("{}/{}", mapping_files_directory, "nucl_gb.accession2taxid"),
            format!("{}/{}", mapping_files_directory, "nucl_wgs.accession2taxid"),
            format!("{}/{}", mapping_files_directory, "pdb.accession2taxid"),
            format!(
                "{}/{}",
                mapping_files_directory, "prot.accession2taxid.FULL"
            ),
        ];

        let truth_output_dir = "test_data/ncbi_compress/split_accessions_by_taxid/truth_outputs";
        let test_output_dir = tempdir().unwrap();
        let test_dir_path_str = test_output_dir.path().to_str().unwrap();

        ncbi_compress::split_accessions_by_taxid(
            input_fasta_path,
            input_mapping_file_paths,
            test_dir_path_str,
        );

        let entries = fs::read_dir(test_dir_path_str)
        .expect("Failed to read directory");

        let entries: Vec<_> = entries
        .filter_map(Result::ok)
        .collect();

        // Assert that the directory is not empty
        assert!(!entries.is_empty(), "Directory is empty");

        // for entry in fs::read_dir(test_dir_path_str).unwrap() {
        for entry in fs::read_dir(truth_output_dir).unwrap() {
            // get test file path
            let entry = entry.unwrap();
            let truth_file_path = entry.path();
            let truth_file_name = truth_file_path.file_name().unwrap().to_str().unwrap();

            // get the test file path
            let test_file_path = format!("{}/{}", test_dir_path_str, truth_file_name);
            util::compare_fasta_records_from_files(&test_file_path, &truth_file_path.to_str().unwrap());
        }
    }

    #[test]
    fn test_conversion() {
        let taxid: u64 = 123;

        let dir = tempdir().unwrap();
        let db = rocksdb::DB::open_default(&dir).unwrap();
        db.put(b"test", taxid.to_be_bytes()).unwrap();

        let result = db.get(b"test").unwrap().unwrap();
        let second_taxid = u64::from_be_bytes(result.as_slice().try_into().unwrap());

        assert_eq!(taxid, second_taxid);
    }
}
