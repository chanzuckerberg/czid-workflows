use std::path::Path;
use std::io::Read;

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
        let temp_file = NamedTempFile::new().unwrap();
        let temp_file_path = temp_file.path();

        fasta_compress(
            Path::new(&sorted_seq),
            mapping_files,
            temp_file_path,
            vec![9606],
            1000,
            31,
            42,
            0.6,
            1000,
            1000,
            true,
        );

        let mut file_content = String::new();
        let mut file = std::fs::File::open(&temp_file_path).expect("Failed to open file");
        file.read_to_string(&mut file_content)
            .expect("Failed to read from file");

        let expected_contents = util::read_contents(&expected_compressed);
        assert!(file_content == expected_contents);
    }

}
