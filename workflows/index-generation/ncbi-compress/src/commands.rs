// Wrappers around common actions / groups of actions for the CLI

pub mod commands {
    use std::fs;

    use bio::io::fasta;

    use crate::fasta_tools::fasta_tools;
    use crate::ncbi_compress::ncbi_compress;

    pub fn count_fasta_reads(input_fasta_path: &str) -> usize {
        let reader = fasta::Reader::from_file(&input_fasta_path).unwrap();
        let records_iter = reader.records();
        return records_iter.count();
    }

    pub fn split_fasta_into_chunks(
        input_fasta_path: &str,
        output_dir: &str,
        total_sequence_count: &usize,
        maximum_records_per_file: &usize,
        taxid: &str,
    ) -> Result<(), String> {
        fs::create_dir_all(&output_dir).expect("Error creating output directory");
        let reader = fasta::Reader::from_file(&input_fasta_path).unwrap();
        let records_iter = reader.records();

        let total_splits = (*total_sequence_count as f32 / *maximum_records_per_file as f32).ceil();
        let sequence_count_per_split = (*total_sequence_count as f32 / total_splits as f32).ceil();
        println!("total_splits: {}", total_splits);
        println!("sequence_count_per_split: {}", sequence_count_per_split);

        let mut seq_count = 0;
        let mut split_count = 1;
        let mut writer =
            fasta::Writer::to_file(format!("{}/{}_{}.fa", output_dir, taxid, split_count)).unwrap();
        for (_i, record) in records_iter.enumerate() {
            let record = record.unwrap();
            writer.write_record(&record).unwrap();
            seq_count += 1;
            if seq_count as f32 == sequence_count_per_split {
                seq_count = 0;
                if split_count as f32 == total_splits {
                    break;
                }
                split_count += 1;
                writer =
                    fasta::Writer::to_file(format!("{}/{}_{}.fa", output_dir, taxid, split_count))
                        .unwrap();
            }
        }

        Ok(())
    }

    pub fn fasta_compress_from_taxid_dir(
        input_taxid_dir: &str,
        output_fasta_path: &str,
        scaled: u64,
        k: u32,
        seed: u64,
        similarity_threshold: f64,
        chunk_size: usize,
        branch_factor: usize,
        is_protein_fasta: bool,
        split_apart_taxid_dir: &str,
    ) {
        log::info!("Starting compression by taxid");
        let mut accession_count = 0;
        let mut unique_accession_count = 0;

        let mut writer = fasta::Writer::to_file(output_fasta_path).unwrap();
        for (_i, entry) in fs::read_dir(input_taxid_dir).unwrap().enumerate() {
            let entry = entry.unwrap();
            let path = entry.path();
            let input_fasta_path = path.to_str().unwrap();
            let taxid = path
                .file_name()
                .unwrap()
                .to_str()
                .unwrap()
                .split(".")
                .collect::<Vec<&str>>()[0];
            let reads_count = count_fasta_reads(input_fasta_path);
            if reads_count >= 9500000 {
                // back of the envelope calculation for how many 50,000 character reads we can store in 488GB of RAM
                log::info!("Breaking apart taxid {} into smaller chunks", taxid);
                let input_taxid_dir = format!("{}/{}_split", split_apart_taxid_dir, taxid);
                split_fasta_into_chunks(
                    &input_fasta_path,
                    &input_taxid_dir,
                    &reads_count,
                    &3,
                    taxid,
                )
                .expect("error splitting fasta into chunks");
                log::info!(
                    "Finished breaking apart taxid {} into smaller chunks",
                    taxid
                );
                // recursively call fasta_compress_from_taxid_dir on the new directory containing the smaller chunks
                for (_i, entry) in fs::read_dir(input_taxid_dir).unwrap().enumerate() {
                    let entry = entry.unwrap();
                    let path = entry.path();
                    let split_taxid_fasta_path = path.to_str().unwrap();
                    ncbi_compress::fasta_compress_taxid(
                        split_taxid_fasta_path,
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
                }
            } else {
                ncbi_compress::fasta_compress_taxid(
                    input_fasta_path,
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
            }
            log::info!("Finished compressing taxid {}", input_fasta_path);
        }
    }

    pub fn fasta_compress_from_fasta_skip_split_by_taxid(
        input_fasta_path: &str,
        output_fasta_path: &str,
        scaled: u64,
        k: u32,
        seed: u64,
        similarity_threshold: f64,
        chunk_size: usize,
        branch_factor: usize,
        is_protein_fasta: bool,
    ) {
        let mut accession_count = 0;
        let mut unique_accession_count = 0;
        let mut writer = fasta::Writer::to_file(output_fasta_path).unwrap();
        ncbi_compress::fasta_compress_taxid(
            input_fasta_path,
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
    }

    pub fn fasta_compress_end_to_end(
        input_fasta_path: &str,
        accession_mapping_files: Vec<String>,
        output_fasta_path: &str,
        temp_file_output_dir: &str,
        scaled: u64,
        k: u32,
        seed: u64,
        similarity_threshold: f64,
        chunk_size: usize,
        branch_factor: usize,
        is_protein_fasta: bool,
    ) {
        // compress from starting fasta file, splitting by taxid

        // create a temp dir containing one file per taxid that input fasta accessions are sorted into
        fs::create_dir_all(&temp_file_output_dir).expect("Error creating output directory");
        let temp_taxid_dir_str = format!("{}/accessions_by_taxid", temp_file_output_dir);
        ncbi_compress::split_accessions_by_taxid(
            input_fasta_path,
            accession_mapping_files,
            &temp_taxid_dir_str,
        );

        let sorted_taxid_dir = format!("{}/sorted_taxid_dir", temp_file_output_dir);
        fasta_tools::sort_taxid_dir_by_sequence_length(&temp_taxid_dir_str, &sorted_taxid_dir);

        let mut accession_count = 0;
        let mut unique_accession_count = 0;
        let mut writer = fasta::Writer::to_file(output_fasta_path).unwrap();
        for (_i, entry) in fs::read_dir(sorted_taxid_dir).unwrap().enumerate() {
            let entry = entry.unwrap();
            let path = entry.path();
            let input_fasta_path = path.to_str().unwrap();
            ncbi_compress::fasta_compress_taxid(
                input_fasta_path,
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
        }
        // shuffle output_fasta_path
        let _ = fasta_tools::shuffle_fasta_by_sequence_index(
            output_fasta_path,
            &format!("shuffled_{}", output_fasta_path),
        );
    }
}

#[cfg(test)]
mod tests {
    use tempfile::tempdir;

    use crate::commands::commands;
    use crate::util::util;

    #[test]
    fn test_fasta_compress_from_fasta_skip_split_by_taxid() {
        let input_fasta_path = "test_data/fasta_tools/inputs/nt";
        let truth_fasta_path = "test_data/commands/common_truth_output/nt_out.fa";
        let test_directory = tempdir().unwrap();
        let temp_dir_path_str = test_directory.path().to_str().unwrap();
        let test_fasta_path = format!("{}/nt.fa", temp_dir_path_str);

        let scaled = 100;
        let k = 31;
        let seed = 42;
        let similarity_threshold = 0.5;
        let chunk_size = 100;
        let branch_factor = 1000;
        let is_protein_fasta = false;

        commands::fasta_compress_from_fasta_skip_split_by_taxid(
            input_fasta_path,
            &test_fasta_path,
            scaled,
            k,
            seed,
            similarity_threshold,
            chunk_size,
            branch_factor,
            is_protein_fasta,
        );
        util::compare_fasta_records_from_files(&truth_fasta_path, &test_fasta_path);
    }

    #[test]
    fn test_fasta_compress_from_fasta_end_to_end() {
        let input_fasta_path = "test_data/fasta_tools/inputs/nt";
        let truth_fasta_path = "test_data/commands/common_truth_output/nt_out.fa";
        let test_directory = tempdir().unwrap();
        let temp_dir_path_str = test_directory.path().to_str().unwrap();
        let test_fasta_path = format!("{}/nt.fa", temp_dir_path_str);

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

        let scaled = 100;
        let k = 31;
        let seed = 42;
        let similarity_threshold = 0.5;
        let chunk_size = 100;
        let branch_factor = 1000;
        let is_protein_fasta = false;

        commands::fasta_compress_end_to_end(
            input_fasta_path,
            input_mapping_file_paths,
            &test_fasta_path,
            temp_dir_path_str,
            scaled,
            k,
            seed,
            similarity_threshold,
            chunk_size,
            branch_factor,
            is_protein_fasta,
        );
        util::compare_fasta_records_from_files(&truth_fasta_path, &test_fasta_path);
    }

    #[test]
    fn test_fasta_compress_from_taxid_dir() {
        let input_taxid_dir = "test_data/commands/fasta_compress_from_taxid_dir/inputs";
        let truth_fasta_path = "test_data/commands/common_truth_output/nt_out.fa";

        let test_directory = tempdir().unwrap();
        let temp_dir_path_str = test_directory.path().to_str().unwrap();
        let output_fasta_path = format!("{}/nt.fa", temp_dir_path_str);

        let scaled = 100;
        let k = 31;
        let seed = 42;
        let similarity_threshold = 0.5;
        let chunk_size = 100;
        let branch_factor = 1000;
        let is_protein_fasta = false;

        commands::fasta_compress_from_taxid_dir(
            input_taxid_dir,
            &output_fasta_path,
            scaled,
            k,
            seed,
            similarity_threshold,
            chunk_size,
            branch_factor,
            is_protein_fasta,
            temp_dir_path_str,
        );

        util::compare_fasta_records_from_files(truth_fasta_path, &output_fasta_path);
    }

    #[test]
    fn test_split_fasta_into_chunks() {
        use bio::io::fasta;

        use crate::commands::commands::count_fasta_reads;
        use crate::commands::commands::split_fasta_into_chunks;

        let taxid_path_name = tempfile::NamedTempFile::new().unwrap();
        let taxid_path_name_str = taxid_path_name.path().to_str().unwrap();
        let mut writer = fasta::Writer::to_file(&taxid_path_name_str).unwrap();

        let records = vec![
            fasta::Record::with_attrs("1", None, b"ACGT"),
            fasta::Record::with_attrs("2", None, b"ACGT"),
            fasta::Record::with_attrs("3", None, b"ACGT"),
            fasta::Record::with_attrs("4", None, b"ACGT"),
            fasta::Record::with_attrs("5", None, b"ACGT"),
            fasta::Record::with_attrs("6", None, b"ACGT"),
            fasta::Record::with_attrs("7", None, b"ACGT"),
        ];

        for rec in records {
            writer.write_record(&rec).expect("error writing to file");
        }
        writer.flush().expect("error flushing file");

        let reads_count = count_fasta_reads(&taxid_path_name_str);
        assert_eq!(reads_count, 7);

        let test_directory = tempdir().unwrap();
        let temp_dir_path_str = test_directory.path().to_str().unwrap();

        let maximum_records_per_file = 5;
        split_fasta_into_chunks(
            &taxid_path_name_str,
            &temp_dir_path_str,
            &reads_count,
            &maximum_records_per_file,
            "123",
        )
        .unwrap();

        let mut chunk_1_reader = util::create_fasta_records_from_file(
            format!("{}/123_1.fa", temp_dir_path_str).as_str(),
        );
        let mut chunk_2_reader = util::create_fasta_records_from_file(
            format!("{}/123_2.fa", temp_dir_path_str).as_str(),
        );

        let mut expected_chunk_1_records = vec![
            fasta::Record::with_attrs("1", None, b"ACGT"),
            fasta::Record::with_attrs("2", None, b"ACGT"),
            fasta::Record::with_attrs("3", None, b"ACGT"),
            fasta::Record::with_attrs("4", None, b"ACGT"),
        ];

        let mut expected_chunk_2_records = vec![
            fasta::Record::with_attrs("5", None, b"ACGT"),
            fasta::Record::with_attrs("6", None, b"ACGT"),
            fasta::Record::with_attrs("7", None, b"ACGT"),
        ];

        assert_eq!(expected_chunk_1_records.sort(), chunk_1_reader.sort());
        assert_eq!(expected_chunk_2_records.sort(), chunk_2_reader.sort());
    }
}
