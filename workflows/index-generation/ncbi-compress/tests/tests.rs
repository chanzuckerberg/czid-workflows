use std::path::Path;
use std::fs;

use sourmash::sketch::minhash::KmerMinHash;
use sourmash::encodings::HashFunctions;
use sourmash::signature::SigsTrait;

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

#[test]
fn test_containment() {
    // TODO: fix this so that it works (create_sequenxe_w_containment not working as expected)

    let branch_factor = 1000;
    let scaled = 1000;
    let k: u32 = 31;
    let seed = 42;
    let similarity_threshold = 0.6;

    let mut tree = ncbi_compress::MinHashTree::new(branch_factor);
    let seq = util::create_random_sequence(k as usize, 1000);
    let seq_60_containment = util::create_sequence_w_containment(k as usize, &seq, 0.6);
    let seq_50_containment = util::create_sequence_w_containment(k as usize, &seq, 0.5);
    let seq_40_containment = util::create_sequence_w_containment(k as usize, &seq, 0.4);

    let seq_chunks = util::split_string_into_chunks(&seq, 1000);
    let seq_60_chunks = util::split_string_into_chunks(&seq_60_containment, 1000);
    // check that seq_chunks and seq_60_chunks have 60% chunk overlap
    let mut overlap_count = 0;
    for chunk in seq_chunks.iter() {
        if seq_60_chunks.contains(&chunk) {
            overlap_count += 1;
        }
    }


    let mut original_hash = KmerMinHash::new(scaled, k, HashFunctions::murmur64_DNA, seed, false, 0);
    original_hash.add_sequence(seq.as_bytes(), true).unwrap();

    let mut hash_40 = KmerMinHash::new(scaled, k, HashFunctions::murmur64_DNA, seed, false, 0);
    hash_40.add_sequence(seq_40_containment.as_bytes(), true).unwrap();

    let cont = ncbi_compress::containment(&original_hash, &hash_40).unwrap();

    println!("cont: {}", cont);

    // // create tree of initial sequences
    // for seq in seqs.iter() {
    //     let mut hash = KmerMinHash::new(scaled, k, HashFunctions::murmur64_DNA, seed, false, 0);
    //     hash.add_sequence(seq.as_bytes(), true).unwrap();
    //     println!("hash mins {:?}", hash.mins());
    //     let _ = tree.insert(hash, seq);
    // }

    // for seq in seqs2.iter() {
    //     let mut hash = KmerMinHash::new(scaled, k, HashFunctions::murmur64_DNA, seed, false, 0);
    //     hash.add_sequence(seq.as_bytes(), true).unwrap();
    //     let contains = tree.contains(&hash, similarity_threshold);
    //     println!("hash mins {:?}", hash.mins());
    //     match contains {
    //         Ok(Some(found_accessions)) => {
    //             let found_accession_ids: Vec<String> = found_accessions.iter().map(|(accession_id, _)| accession_id.to_string()).collect::<Vec<_>>();
    //             let found_accession_containments: Vec<String> = found_accessions.iter().map(|(_, containment)| containment.to_string()).collect::<Vec<_>>();
    //             println!("{} Found: {:?} Containments: {:?}", seq, found_accession_ids, found_accession_containments);
    //         },
    //         Ok(None) => println!("{} Not found ", seq),
    //         Err(e) => println!("Error: {}", e),
    //     }

        // assert!(tree.contains(seq));
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
