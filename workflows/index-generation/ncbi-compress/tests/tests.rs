use std::path::Path;
use std::fs;

// use tempfile::NamedTempFile;

use ::ncbi_compress::ncbi_compress::ncbi_compress;
use ::ncbi_compress::util::util;

#[test]
fn test_split_accessions_by_taxis() {
    let pathogens = ["chkv", "streptococcus", "rhinovirus"];

    for pathogen in pathogens .iter() {
        let mapping_files = vec![
            Path::new("tests/test_data/accession2taxid/nucl_gb.accession2taxid.subset"),
            Path::new("tests/test_data/accession2taxid/nucl_wgs.accession2taxid.subset"),
            Path::new("tests/test_data/accession2taxid/pdb.accession2taxid.subset"),
        ];
        let sorted_seq = format!("tests/test_data/simulated_seqs/all_simulated_seqs_{}_sorted_subsample.fasta", pathogen);
        let taxids_to_drop: Vec<u64> = vec![9606];
        let taxid_dir = ncbi_compress::split_accessions_by_taxid(
            sorted_seq,
            mapping_files,
            &taxids_to_drop,
        );
        for (i, entry) in fs::read_dir(taxid_dir.path()).unwrap().enumerate() {
            let entry = entry.unwrap();
            let path = entry.path();
            let input_fasta_path = path.to_str().unwrap();
            let expected = format!("tests/test_data/expected_split_accessions_by_taxid/test_file_{}_{}.txt", i, pathogen);
            assert!(util::are_files_equal(input_fasta_path, &expected))
        }
    }
}

// #[test]
// fn test_fasta_compress_taxid() {
//     let input_fasta_path = Path::new("tests/test_data/simulated_seqs/all_simulated_seqs_chkv_sorted.fasta");
//     let expected_fasta_path = "tests/test_data/expected_fasta_compress_taxid_chkv.fasta";
//     let mut temp_file = NamedTempFile::new().unwrap();
//     let temp_file_path = temp_file.path();
//     let mut writer = fasta::Writer::to_file(temp_file_path).unwrap();
//     let scaled = 1000;
//     let k = 31;
//     let seed = 42;
//     let similarity_threshold = 0.6;
//     let chunk_size = 1000;
//     let branch_factor = 1000;
//     let mut accession_count = 0;
//     let mut unique_accession_count = 0;
//     ncbi_compress::fasta_compress_taxid(
//         input_fasta_path,
//         &mut writer,
//         scaled,
//         k,
//         seed,
//         similarity_threshold,
//         chunk_size,
//         branch_factor,
//         &mut accession_count,
//         &mut unique_accession_count,
//     );

//     assert!(util::are_files_equal(temp_file_path.file_name().unwrap().to_str().unwrap(), &expected_fasta_path));
// }
