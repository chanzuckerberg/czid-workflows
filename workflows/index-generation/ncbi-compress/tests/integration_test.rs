use std::path::Path;
use std::fs;

use tempfile::NamedTempFile;

use ::ncbi_compress::ncbi_compress::ncbi_compress::fasta_compress;
use ::ncbi_compress::util::util;

#[test]
fn test_fasta_compress() {
    assert!(true);
    let pathogens = ["chkv", "streptococcus", "rhinovirus"];

    for pathogen in pathogens .iter() {

        let mapping_files = vec![
            Path::new("tests/test_data/accession2taxid/nucl_gb.accession2taxid.subset"),
            Path::new("tests/test_data/accession2taxid/nucl_wgs.accession2taxid.subset"),
            Path::new("tests/test_data/accession2taxid/pdb.accession2taxid.subset"),
        ];

        let sorted_seq = format!("tests/test_data/simulated_seqs/all_simulated_seqs_{}_sorted_subsample.fasta", pathogen);
        let expected_compressed = format!("tests/test_data/expected_compression_results/nt_compressed_0.6_{}.fa", pathogen);

        let mut temp_file = NamedTempFile::new().unwrap();
        let temp_file_path = temp_file.path();

        fasta_compress(
            Path::new(&sorted_seq),
            mapping_files,
            Path::new("tests/test_data/test_output.fasta"),
            vec![9606],
            1000,
            31,
            42,
            0.6,
            1000,
            1000
        );

        assert!(util::are_files_equal("tests/test_data/test_output.fasta", &expected_compressed));
    }

}
