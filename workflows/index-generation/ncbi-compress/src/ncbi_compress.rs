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
    use trie_rs::{Trie, TrieBuilder};

    use crate::logging::logging;

        /// A trie that stores u64 values
    struct TrieStore {
        trie: Trie<u8>,
    }

    impl TrieStore {
        pub fn get(&self, key: &str) -> Option<u64> {
            self.trie
                .predictive_search(key.as_bytes())
                .first()
                .map(|bytes| {
                    let mut bytes = bytes.to_vec();
                    let value_bytes = bytes.split_off(bytes.len() - 8);
                    u64::from_be_bytes(value_bytes.try_into().unwrap())
                })
        }
    }

    struct TrieStoreBuilder {
        builder: TrieBuilder<u8>,
    }

    impl TrieStoreBuilder {
        pub fn new() -> Self {
            TrieStoreBuilder {
                builder: TrieBuilder::new(),
            }
        }

        pub fn push(&mut self, key: &str, value: u64) {
            let mut key = key.as_bytes().to_vec();
            key.extend_from_slice(&value.to_be_bytes());
            self.builder.push(key);
        }

        pub fn build(self) -> TrieStore {
            let trie = self.builder.build();
            TrieStore { trie }
        }
    }

    struct MinHashTreeNode {
        own: KmerMinHash,
        children_aggregate: KmerMinHash,
        accession: String,
    }

    struct MinHashTree {
        branching_factor: usize,
        nodes: Vec<MinHashTreeNode>,
    }

    impl MinHashTree {
        fn new(branching_factor: usize) -> Self {
            MinHashTree {
                branching_factor,
                nodes: Vec::new(),
            }
        }

        fn parent_idx(&self, node: usize) -> Option<usize> {
            if node == 0 {
                None
            } else {
                Some((node - 1) / self.branching_factor)
            }
        }

        fn child_idxes(&self, node: usize) -> Vec<usize> {
            let first = self.branching_factor * node + 0 + 1;
            let last = self.branching_factor * node + self.branching_factor + 1;
            (first..last.min(self.nodes.len() - 1)).collect()
        }

        fn merge_to_parent(
            &mut self,
            parent_idx: usize,
            child_idx: usize,
        ) -> Result<(), SourmashError> {
            let (left, right) = self.nodes.split_at_mut(child_idx);
            left[parent_idx]
                .children_aggregate
                .merge(&right[0].children_aggregate)
        }

        pub fn insert(&mut self, hash: KmerMinHash, accession: &str) -> Result<(), SourmashError> {
            let node = MinHashTreeNode {
                own: hash.clone(),
                children_aggregate: hash.clone(),
                accession: accession.to_string(),
            };
            let mut current_idx = self.nodes.len();
            self.nodes.push(node);

            while let Some(parent_idx) = self.parent_idx(current_idx) {
                // no need to aggregate to the root, it would just contain everything thus providing no information
                if parent_idx == 0 {
                    break;
                }

                self.merge_to_parent(parent_idx, current_idx)?;
                current_idx = parent_idx;
            }
            Ok(())
        }

        pub fn contains(
            &self,
            hash: &KmerMinHash,
            similarity_threshold: f64,
        ) -> Result<Option<Vec<String>>, SourmashError> {
            if self.nodes.is_empty() {
                return Ok(None);
            }

            let mut to_visit = vec![0];
            while !to_visit.is_empty() {
                // Check if any of the nodes in the to_visit list are similar enough
                let found_accession = to_visit.par_iter().filter_map(|node_idx| {
                    let node = self.nodes.get(*node_idx).unwrap();
                    if containment(hash, &node.own).unwrap() >= similarity_threshold {
                        Some(node.accession.clone())
                    } else {
                        None
                    }
                }).collect::<Vec<_>>();
                // If we found a similar node, we can stop searching
                if !found_accession.is_empty() {
                    return Ok(Some(found_accession));
                }
                // Otherwise, we need to search the children of the nodes in the to_visit list
                to_visit = to_visit
                    .par_iter()
                    .flat_map(|node_idx| {
                        let node = self.nodes.get(*node_idx).unwrap();
                        // If the children are similar enough, we need to search them
                        if containment(hash, &node.children_aggregate).unwrap() >= similarity_threshold
                        {
                            self.child_idxes(*node_idx)
                        } else {
                            vec![]
                        }
                    })
                    .collect();
            }
            // we didn't find any similar nodes in the tree
            Ok(None)
        }
    }

    fn containment(needle: &KmerMinHash, haystack: &KmerMinHash) -> Result<f64, SourmashError> {
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
    ) -> TempDir {
        // create a temp dir containing one file per taxid that input fasta accessions are sorted into
        // based on taxid in the mapping files (input fasta does not have taxid in header)
        log::info!("Creating accession to taxid mapping");
        let taxid_dir = TempDir::new("accessions_by_taxid").unwrap();
        let reader = fasta::Reader::from_file(&input_fasta_path).unwrap();
        // Build a trie of the accessions in the input fasta
        let mut builder = TrieBuilder::new();
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

    pub fn fasta_compress_taxid<P: AsRef<Path> + std::fmt::Debug>(
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
        // take in a fasta file and output a fasta file with only unique accessions (based on similarity threshold)
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
                    let contained_in_tree = tree.contains(&hash, similarity_threshold);
                    match contained_in_tree {
                        Ok(None) => {
                            // If the tree doesn't contain the hash, we need to insert it
                            Some((hash, record))
                        }
                        Ok(Some(found_accessions)) => {
                            // log when tree already contains hash
                            let accession_id = record.id().split_whitespace().next().unwrap();
                            logging::write_to_file(format!("accession_id {} found similar sequence in the tree and is being represented by : {}", accession_id, found_accessions.join(",")));
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

            let mut tmp: Vec<(KmerMinHash, &str)> = Vec::with_capacity(chunk_signatures.len() / 2); // why is this /2?
            for (hash, record) in chunk_signatures {
                let accession_id = record.id().split_whitespace().next().unwrap();
                accession_count.add_assign(1); // += 1, used for logging
                // Perform a faster similarity check over just this chunk because we may have similarities within a chunk
                let similar_seqs = tmp
                    .par_iter()
                    .filter_map(|(other, accession_id)| {
                        if containment(&hash, &other).unwrap() >= similarity_threshold {
                            Some(accession_id.to_string())
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
                    logging::write_to_file(
                        format!("accession_id {} is being representing by {} \n",accession_id, similar_seqs.join(", "))
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
    ) {
        log::info!("Splitting accessions by taxid");
        let taxid_dir =
            split_accessions_by_taxid(&input_fasta_path, accession_mapping_files, &taxids_to_drop);
        log::info!("Finished splitting accessions by taxid");
        let mut writer = fasta::Writer::to_file(output_fasta_path).unwrap();

        log::info!("Starting compression by taxid");
        let mut accession_count = 0;
        let mut unique_accession_count = 0;
        for (i, entry) in fs::read_dir(taxid_dir.path()).unwrap().enumerate() {
            let entry = entry.unwrap();
            let path = entry.path();
            let input_fasta_path = path.to_str().unwrap();
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
