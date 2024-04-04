pub mod ncbi_compress {
    use std::ops::AddAssign;
    use std::path::Path;
    use std::{borrow::BorrowMut, fs};

    use bio::io::fasta;
    use rand::distributions::Alphanumeric;
    use rand::Rng;
    use rayon::prelude::*;
    use sourmash::encodings::HashFunctions;
    use sourmash::errors::SourmashError;
    use sourmash::signature::SigsTrait;
    use sourmash::sketch::minhash::KmerMinHash;

    pub fn containment(needle: &KmerMinHash, haystack: &KmerMinHash) -> Result<f64, SourmashError> {
        let (intersect_size, _) = needle.intersection_size(haystack)?;
        Ok(intersect_size as f64 / needle.mins().len() as f64)
    }

    pub fn remove_accession_version(accession: &str) -> &str {
        accession.splitn(2, |c| c == '.').next().unwrap()
    }

    pub fn create_accession_to_taxid_db(
        dir: &str,
        input_fasta_path: &str,
        mapping_file_path: Vec<String>,
    ) -> rocksdb::DB {
        log::info!("Creating accession to taxid db");

        let accession_to_taxid = rocksdb::DB::open_default(&dir).unwrap();

        let reader = fasta::Reader::from_file(&input_fasta_path).unwrap();
        // Build a trie of the accessions in the input fasta
        reader
            .records()
            .enumerate()
            .par_bridge()
            .for_each(|(i, result)| {
                let record = result.as_ref().unwrap();
                let accession_id = record.id().split_whitespace().next().unwrap();
                let accession_no_version = remove_accession_version(accession_id);
                // RocksDB supports concurrent reads and writes so this is safe
                accession_to_taxid.put(accession_no_version, b"").unwrap();
                if i % 1_000_000 == 0 {
                    log::info!("  Processed {} accessions", i);
                }
            });
        log::info!(" Finished loading accessions");
        mapping_file_path.par_iter().for_each(|mapping_file_path| {
            let reader = csv::ReaderBuilder::new()
                .delimiter(b'\t')
                .from_path(mapping_file_path)
                .unwrap();
            let mut added = 0;
            reader.into_records().enumerate().for_each(|(_i, result)| {
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
        accession_to_taxid
    }

    pub fn split_accessions_by_taxid(
        input_fasta_path: &str,
        mapping_file_path: Vec<String>,
        output_dir: &str,
    ) -> std::path::PathBuf {
        log::info!("Splitting accessions by taxid");
        log::info!("Creating taxid dir {:?}", output_dir);
        fs::create_dir_all(&output_dir).expect("Error creating output directory");
        // Use a random string for the rocksdb directory name
        // A tempdir can be removed by the operating system before we are done with it which
        // break the rocksdb. We need a random name so we can call this function multiple
        // times from different threads during testing.
        let rng = rand::thread_rng();
        let dir: String = rng
            .sample_iter(&Alphanumeric)
            .take(20)
            .map(char::from)
            .collect();
        let accession_to_taxid =
            create_accession_to_taxid_db(&dir, input_fasta_path, mapping_file_path);
        let outpath = write_accessions_to_taxid(input_fasta_path, &accession_to_taxid, output_dir);
        fs::remove_dir_all(dir).expect("Error removing db");
        outpath
    }

    pub fn write_accessions_to_taxid(
        input_fasta_path: &str,
        accession_to_taxid: &rocksdb::DB,
        output_dir: &str,
    ) -> std::path::PathBuf {
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
            let file_path = format!("{}/{}.fasta", output_dir, taxid);
            let file = fs::OpenOptions::new()
                .create(true)
                .append(true)
                .open(file_path)
                .unwrap();
            let mut writer = fasta::Writer::new(file);
            writer.write_record(&record).unwrap();
        }
        log::info!("Finished splitting accessions by taxid");
        Path::new(output_dir).to_path_buf()
    }

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
                })
                .collect::<Vec<_>>();

            // we need to make sure records within the chunk arn't similar to each other before
            // we check them against the larger tree
            for (hash, record) in chunk_signatures {
                let similar = unique_in_chunk.par_iter().any(|(other, _record)| {
                    containment(&hash, &other).unwrap() >= similarity_threshold
                });

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
        }
    }
}

#[cfg(test)]
mod tests {
    use bio::io::fasta;
    use sourmash::signature::SigsTrait;
    use tempfile::tempdir;

    use crate::fasta_tools::fasta_tools;
    use crate::ncbi_compress::ncbi_compress;
    use crate::util::util;
    use std::fs;

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

    #[test]
    fn test_remove_accession_version() {
        let accession = "NC_000913.3";
        let accession_no_version = ncbi_compress::remove_accession_version(accession);
        assert_eq!(accession_no_version, "NC_000913");
    }

    #[test]
    fn test_containment() {
        let ksize = 31;
        let mut haystack_hash = sourmash::sketch::minhash::KmerMinHash::new(
            1,
            ksize.try_into().unwrap(),
            sourmash::encodings::HashFunctions::murmur64_DNA,
            62,
            false,
            0,
        );
        let mut needle_hash = sourmash::sketch::minhash::KmerMinHash::new(
            1,
            ksize.try_into().unwrap(),
            sourmash::encodings::HashFunctions::murmur64_DNA,
            62,
            false,
            0,
        );

        let needle = util::create_random_sequence(ksize, 25);
        let desired_containment = 0.6;
        let haystack_w_containment =
            util::create_sequence_w_containment(&needle, ksize, 25, desired_containment);
        haystack_hash
            .add_sequence(haystack_w_containment.as_bytes(), true)
            .unwrap();
        needle_hash.add_sequence(needle.as_bytes(), true).unwrap();

        let result = ncbi_compress::containment(&needle_hash, &haystack_hash).unwrap();
        // Assert that value is greater than 0.58 and less than 0.62
        assert!(
            result > (desired_containment - 0.02) && result < desired_containment + 0.02,
            "{}",
            format!("Containment is not between 0.58 and 0.62: {}", result)
        );
    }

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

        let entries = fs::read_dir(test_dir_path_str).expect("Failed to read directory");

        let entries: Vec<_> = entries.filter_map(Result::ok).collect();

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
            util::compare_fasta_records_from_files(
                &test_file_path,
                &truth_file_path.to_str().unwrap(),
            );
        }
    }

    #[test]
    fn test_create_accession_to_taxid_db() {
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

        let dir = tempdir().unwrap();
        let db = ncbi_compress::create_accession_to_taxid_db(
            dir.path().to_str().unwrap(),
            input_fasta_path,
            input_mapping_file_paths,
        );

        let result = db.get(b"X52703").unwrap().unwrap();
        let taxid = u64::from_be_bytes(result.as_slice().try_into().unwrap());
        assert_eq!(taxid, 9771);

        let result = db.get(b"Z14040").unwrap().unwrap();
        let taxid = u64::from_be_bytes(result.as_slice().try_into().unwrap());
        assert_eq!(taxid, 9913);
    }

    #[test]
    fn test_write_accessions_to_taxid() {
        // a little redundant compared to the split_accessions_by_taxid test
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

        let dir = tempdir().unwrap();
        let db = ncbi_compress::create_accession_to_taxid_db(
            dir.path().to_str().unwrap(),
            input_fasta_path,
            input_mapping_file_paths,
        );

        let output_dir = tempdir().unwrap();
        let output_dir_path_str = output_dir.path().to_str().unwrap();
        let outpath =
            ncbi_compress::write_accessions_to_taxid(input_fasta_path, &db, output_dir_path_str);

        let entries = fs::read_dir(output_dir_path_str).expect("Failed to read directory");

        let entries: Vec<_> = entries.filter_map(Result::ok).collect();

        // Assert that the directory is not empty
        assert!(!entries.is_empty(), "Directory is empty");

        for entry in fs::read_dir(outpath).unwrap() {
            // get test file path
            let entry = entry.unwrap();
            let test_file_path = entry.path();
            let test_file_name = test_file_path.file_name().unwrap().to_str().unwrap();

            // get the truth file path
            let truth_file_path = format!("{}/{}", truth_output_dir, test_file_name);
            util::compare_fasta_records_from_files(
                &test_file_path.to_str().unwrap(),
                &truth_file_path,
            );
        }
    }

    #[test]
    fn test_fasta_compress_taxid() {
        let ksize = 31;
        // create some sequences with containment between them
        let seq_shortest = util::create_random_sequence(ksize, 50);
        let seq_short = util::create_sequence_w_containment(&seq_shortest, ksize, 100, 0.65);
        let seq_medium = util::create_sequence_w_containment(&seq_short, ksize, 250, 0.65);
        let seq_long = util::create_sequence_w_containment(&seq_medium, ksize, 500, 0.65);
        let seq_longest = util::create_sequence_w_containment(&seq_long, ksize, 1000, 0.65);

        // create some "unique sequences" from the above sequences
        let unique_1 = util::create_sequence_w_containment(&seq_shortest, ksize, 100, 0.2);
        let unique_2 = util::create_sequence_w_containment(&seq_long, ksize, 600, 0.3);

        // organize the sequences into what we expect to be kept and what we expect to be discarded after compression
        let records_that_should_be_discarded = vec![
            fasta::Record::with_attrs("seq_long", None, seq_long.as_bytes()),
            fasta::Record::with_attrs("seq_medium", None, seq_medium.as_bytes()),
            fasta::Record::with_attrs("seq_short", None, seq_short.as_bytes()),
            fasta::Record::with_attrs("seq_shortest", None, seq_shortest.as_bytes()),
        ];
        let records_that_should_be_kept = vec![
            fasta::Record::with_attrs("seq_longest", None, seq_longest.as_bytes()),
            fasta::Record::with_attrs("unique_1", None, unique_1.as_bytes()),
            fasta::Record::with_attrs("unique_2", None, unique_2.as_bytes()),
        ];

        // combine the records into a single vector and write to file
        let combined_records: Vec<fasta::Record> = records_that_should_be_discarded
            .iter() // Borrow each element of the first vector.
            .cloned() // Clone each element (needed because fasta::Record likely doesn't implement Copy).
            .chain(records_that_should_be_kept.iter().cloned()) // Do the same for the second vector.
            .collect(); // Collect into a new Vec<fasta::Record>.

        // write test records to file
        let test_file = tempfile::NamedTempFile::new().unwrap();
        let test_file_path = test_file.path().to_str().unwrap();
        let mut writer = fasta::Writer::to_file(test_file_path).unwrap();
        for record in combined_records {
            writer.write_record(&record).unwrap();
        }

        // sort the fasta
        let test_sorted_file = tempfile::NamedTempFile::new().unwrap();
        let test_sorted_path = test_sorted_file.path().to_str().unwrap();
        let _ = fasta_tools::sort_fasta_by_sequence_length(test_file_path, test_sorted_path);

        // compress the fasta
        let output_fasta_path = tempfile::NamedTempFile::new().unwrap();
        let output_fasta_path = output_fasta_path.path().to_str().unwrap();
        let mut writer = fasta::Writer::to_file(output_fasta_path).unwrap();
        let scaled = 1;
        let k = 31;
        let seed = 62;
        let similarity_threshold = 0.6;
        let chunk_size = 1000;
        let branch_factor = 100;
        let is_protein_fasta = false;
        let mut accession_count = 0;
        let mut unique_accession_count = 0;

        ncbi_compress::fasta_compress_taxid(
            test_sorted_path,
            &mut writer,
            scaled,
            k,
            seed,
            similarity_threshold,
            chunk_size,
            branch_factor,
            is_protein_fasta,
            &mut accession_count,
            &mut unique_accession_count,
        );
        drop(writer);

        // get the records after compression
        let output_fasta_records = util::create_fasta_records_from_file(output_fasta_path);

        // convert the records to a format that can be compared and logged easily
        let actual_output_fasta_records = output_fasta_records
            .into_iter()
            .map(|record| (record.id().to_owned(), record.seq().to_owned()))
            .collect::<Vec<_>>();

        let mut expected_records_tuples = vec![];
        for record in &records_that_should_be_kept {
            expected_records_tuples.push((record.id().to_owned(), record.seq().to_owned()));
        }

        let mut records_that_should_be_discarded_tuples = vec![];
        for record in &records_that_should_be_discarded {
            records_that_should_be_discarded_tuples
                .push((record.id().to_owned(), record.seq().to_owned()));
        }

        // verify that the records that should be kept are in the output fasta
        for record in expected_records_tuples {
            assert!(
                actual_output_fasta_records.contains(&record),
                "Record {:?} is missing from the output fasta",
                record.0
            );
        }

        // verify that the records that should be discarded are not in the output fasta
        for record in &records_that_should_be_discarded_tuples {
            assert!(
                !actual_output_fasta_records.contains(&record),
                "Record {:?} should not be in the output fasta",
                record.0
            );
        }
    }
}
