pub mod commands {
    use std::fs;

    use tempdir::TempDir;
    use bio::io::fasta;

    use crate::ncbi_compress::ncbi_compress;


    fn fasta_compress_w_logging_option (
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
        enable_sequence_retention_logging: bool,
        logging_contained_in_tree_fn: &str,
        logging_contained_in_chunk_fn: &str,
    ) {
        if enable_sequence_retention_logging {
            ncbi_compress::fasta_compress_taxid_w_logging(
                input_fasta_path,
                writer,
                scaled,
                k,
                seed,
                similarity_threshold,
                chunk_size,
                branch_factor,
                is_protein_fasta,
                accession_count,
                unique_accession_count,
                logging_contained_in_tree_fn,
                logging_contained_in_chunk_fn,
            );
        } else {
            ncbi_compress::fasta_compress_taxid(
                input_fasta_path,
                writer,
                scaled,
                k,
                seed,
                similarity_threshold,
                chunk_size,
                branch_factor,
                is_protein_fasta,
                accession_count,
                unique_accession_count,
            );
        }
    }

    pub fn fasta_compress_from_taxid_dir (
        input_taxid_dir: &str,
        output_fasta_path: &str,
        scaled: u64,
        k: u32,
        seed: u64,
        similarity_threshold: f64,
        chunk_size: usize,
        branch_factor: usize,
        is_protein_fasta: bool,
        enable_sequence_retention_logging: bool,
        logging_contained_in_tree_fn: &str,
        logging_contained_in_chunk_fn: &str,
    ) {
        log::info!("Starting compression by taxid");
        let mut accession_count = 0;
        let mut unique_accession_count = 0;
        let mut writer = fasta::Writer::to_file(output_fasta_path).unwrap();
        for (i, entry) in fs::read_dir(input_taxid_dir).unwrap().enumerate() {
            let entry = entry.unwrap();
            let path = entry.path();
            let input_fasta_path = path.to_str().unwrap();
            fasta_compress_w_logging_option(
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
                enable_sequence_retention_logging,
                logging_contained_in_tree_fn,
                logging_contained_in_chunk_fn,
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
    }

    pub fn fasta_compress_from_fasta_skip_split_by_taxid (
        input_fasta_path: &str,
        output_fasta_path: &str,
        scaled: u64,
        k: u32,
        seed: u64,
        similarity_threshold: f64,
        chunk_size: usize,
        branch_factor: usize,
        is_protein_fasta: bool,
        enable_sequence_retention_logging: bool,
        logging_contained_in_tree_fn: &str,
        logging_contained_in_chunk_fn: &str,
    ) {
        let mut accession_count = 0;
        let mut unique_accession_count = 0;
        let mut writer = fasta::Writer::to_file(output_fasta_path).unwrap();
        fasta_compress_w_logging_option(
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
            enable_sequence_retention_logging,
            logging_contained_in_tree_fn,
            logging_contained_in_chunk_fn,
        );
    }

    pub fn fasta_compress_end_to_end (
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
        enable_sequence_retention_logging: bool,
        logging_contained_in_tree_fn: &str,
        logging_contained_in_chunk_fn: &str,
    ) {
        // compress from starting fasta file, splitting by taxid

        // create a temp dir containing one file per taxid that input fasta accessions are sorted into
        fs::create_dir_all(&temp_file_output_dir).expect("Error creating output directory");
        let temp_taxid_dir = TempDir::new_in(temp_file_output_dir, "accessions_by_taxid").expect("Error creating tempdir");
        let temp_taxid_dir_str = temp_taxid_dir.path().to_str().unwrap();

        // let taxid_dir_str = format!("{}/accessions_by_taxid", temp_file_output_dir);
        let taxid_dir = ncbi_compress::split_accessions_by_taxid(
            input_fasta_path,
            accession_mapping_files,
            &temp_taxid_dir_str
        );

        let mut accession_count = 0;
        let mut unique_accession_count = 0;
        let mut writer = fasta::Writer::to_file(output_fasta_path).unwrap();
        for (i, entry) in fs::read_dir(taxid_dir).unwrap().enumerate() {
            let entry = entry.unwrap();
            let path = entry.path();
            let input_fasta_path = path.to_str().unwrap();
            fasta_compress_w_logging_option(
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
                enable_sequence_retention_logging,
                logging_contained_in_tree_fn,
                logging_contained_in_chunk_fn,
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
    }

}


// testing starts here
use tempfile::tempdir;

use crate::util::util;

#[test]
fn test_fasta_compress_from_fasta_skip_split_by_taxid() {
    let input_fasta_path = "test_data/fasta_tools/inputs/nt";
    let truth_fasta_path = "test_data/commands/fasta_compress_from_fasta_skip_split_by_taxid/truth-ouputs/nt_out.fa";
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
    let enable_sequence_retention_logging = false;
    let logging_contained_in_tree_fn = "";
    let logging_contained_in_chunk_fn = "";

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
        enable_sequence_retention_logging,
        logging_contained_in_tree_fn,
        logging_contained_in_chunk_fn,
    );

    assert!(util::are_files_equal(truth_fasta_path, &test_fasta_path));
}

#[test]
fn test_fasta_compress_from_fasta_end_to_end() {
    let input_fasta_path = "test_data/fasta_tools/inputs/nt";
    let truth_fasta_path = "test_data/commands/common_truth_output/nt_out.fa";
    let test_directory = tempdir().unwrap();
    let temp_dir_path_str = test_directory.path().to_str().unwrap();
    let test_fasta_path = format!("{}/nt.fa", temp_dir_path_str);

    let mapping_files_directory = "test_data/ncbi_compress/split_accessions_by_taxid/inputs/accession2taxid";
    let input_mapping_file_paths = vec![
        format!("{}/{}", mapping_files_directory, "nucl_gb.accession2taxid"),
        format!("{}/{}", mapping_files_directory, "nucl_wgs.accession2taxid"),
        format!("{}/{}", mapping_files_directory, "pdb.accession2taxid"),
        format!("{}/{}", mapping_files_directory, "prot.accession2taxid.FULL"),
    ];

    let scaled = 100;
    let k = 31;
    let seed = 42;
    let similarity_threshold = 0.5;
    let chunk_size = 100;
    let branch_factor = 1000;
    let is_protein_fasta = false;
    let enable_sequence_retention_logging = false;
    let logging_contained_in_tree_fn = "";
    let logging_contained_in_chunk_fn = "";

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
        enable_sequence_retention_logging,
        logging_contained_in_tree_fn,
        logging_contained_in_chunk_fn,
    );

    assert!(util::are_files_equal(truth_fasta_path, &test_fasta_path));
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
    let enable_sequence_retention_logging = false;
    let logging_contained_in_tree_fn = "";
    let logging_contained_in_chunk_fn = "";

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
        enable_sequence_retention_logging,
        logging_contained_in_tree_fn,
        logging_contained_in_chunk_fn,
    );

    assert!(util::are_files_equal(truth_fasta_path, &output_fasta_path));
}