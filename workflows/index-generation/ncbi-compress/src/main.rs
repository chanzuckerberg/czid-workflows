use clap::{Arg, Command};

use ncbi_compress::commands::commands::{
    fasta_compress_end_to_end, fasta_compress_from_fasta_skip_split_by_taxid,
    fasta_compress_from_taxid_dir,
};
use ncbi_compress::fasta_tools::fasta_tools::{
    sort_fasta_by_sequence_length,
    count_accessions_by_taxid,
    sort_taxid_dir_by_sequence_length,
    shuffle_fasta_by_sequence_index
};
use ncbi_compress::ncbi_compress::ncbi_compress::split_accessions_by_taxid;

use ncbi_compress::logging::logging;

pub fn main() {
    logging::init_stdout_logging();
    let matches = Command::new("ncbi-compress")
        .version("1.0")
        .author("Your Name")
        .about("Does awesome things with sequences")
        .subcommand(
            Command::new("sort-fasta-by-sequence-length")
                .about("sort_fasta_by_sequence_length")
                .arg(
                    Arg::new("input_fasta")
                        .help("Input file to be sorted")
                        .long("input-fasta")
                        .required(true),
                )
                .arg(
                    Arg::new("output_fasta")
                        .help("Output file for the sorted fasta")
                        .long("output")
                        .required(true),
                ),
        )
        .subcommand(
            Command::new("sort-taxid-dir-by-sequence-length")
                .about("sort_taxid_dir_by_sequence_length")
                .arg(
                    Arg::new("input_taxid_dir")
                        .help("Input directory of taxid fasta files")
                        .long("input-taxid-dir")
                        .required(true),
                )
                .arg(
                    Arg::new("output_taxid_dir")
                        .help("Output directory for the sorted taxid fasta files")
                        .long("output-taxid-dir")
                        .required(true),
                ),
        )
        .subcommand(
            Command::new("shuffle-fasta")
                .about("shuffle-fasta")
                .arg(
                    Arg::new("input_fasta")
                        .help("Input fasta file to be shuffled")
                        .long("input-fasta")
                        .required(true),
                )
                .arg(
                    Arg::new("output_fasta")
                        .help("Output fasta file for the shuffled fasta")
                        .long("output-fasta")
                        .required(true),
                ),
        )
        .subcommand(
            Command::new("break-into-individual-taxids-only")
                .about("break_into_individual_taxids_only")
                .arg(
                    Arg::new("input_fasta")
                        .help("Input file to be ordered")
                        .long("input-fasta")
                        .required(true),
                )
                .arg(
                    Arg::new("accession_mapping_files")
                        .help("Accession mapping files")
                        .long("accession-mapping-files")
                        .required(true)
                        .num_args(1..=4),
                ) // Allows the flag to appear multiple times
                .arg(
                    Arg::new("output_dir")
                        .help("Output directory for the split fasta files")
                        .long("output-dir")
                        .required(true),
                ),
        )
        .subcommand(
            Command::new("count-accessions-by-taxid")
                .about("count-accessions-by-taxid")
                .arg(
                    Arg::new("input-taxid-dir")
                        .help("input directory for fasta files broken into taxids")
                        .long("input-taxid-dir")
                        .required(true),
                )
                .arg(
                    Arg::new("output-summary-path")
                        .help("output-summary-path")
                        .long("output-summary-path")
                        .required(true)
                )
        )
        .subcommand(
            Command::new("fasta-compress-from-fasta-skip-split-by-taxid")
                .about("fasta_compress_from_fasta_skip_split_by_taxid")
                .arg(
                    Arg::new("input_fasta")
                        .help("Input file to be ordered")
                        .long("input-fasta")
                        .required(true),
                )
                .arg(
                    Arg::new("output_fasta")
                        .help("Output file for the ordered fasta")
                        .long("output-fasta")
                        .required(true),
                )
                .arg(
                    Arg::new("scaled")
                        .help("Scaled value for the minhash")
                        .long("scaled")
                        .value_parser(clap::value_parser!(u64))
                        .required(true),
                )
                .arg(
                    Arg::new("k")
                        .help("Kmer size for the minhash")
                        .long("k")
                        .value_parser(clap::value_parser!(u32))
                        .required(true),
                )
                .arg(
                    Arg::new("seed")
                        .help("Seed for the minhash")
                        .long("seed")
                        .value_parser(clap::value_parser!(u64))
                        .default_value("42"),
                )
                .arg(
                    Arg::new("similarity_threshold")
                        .help("Similarity threshold for the minhash")
                        .long("similarity-threshold")
                        .value_parser(clap::value_parser!(f64))
                        .required(true),
                )
                .arg(
                    Arg::new("chunk_size")
                        .help("Chunk size for the parallelization")
                        .long("chunk-size")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("10000"),
                )
                .arg(
                    Arg::new("branch_factor")
                        .help("Branching factor for the tree")
                        .long("branch-factor")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("64"),
                )
                .arg(
                    Arg::new("is_protein_fasta")
                        .help("Is protein fasta")
                        .long("is-protein-fasta")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("enable_sequence_retention_logging")
                        .help("Enable sequence retention logging")
                        .long("enable-sequence-retention-logging")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("logging_contained_in_tree_fn")
                        .help("Logging file for containment in the tree")
                        .long("logging-contained-in-tree-fn"),
                )
                .arg(
                    Arg::new("logging_contained_in_chunk_fn")
                        .help("Logging file for containment in the chunk")
                        .long("logging-contained-in-chunk-fn"),
                ),
        )
        .subcommand(
            Command::new("fasta-compress-from-taxid-dir")
                .about("fasta_compress_from_taxid_dir")
                .arg(
                    Arg::new("input_fasta_dir")
                        .help("Input directory of taxid fasta files")
                        .long("input-fasta-dir")
                        .required(true),
                )
                .arg(
                    Arg::new("output_fasta")
                        .help("Output file for the ordered fasta")
                        .long("output-fasta")
                        .required(true),
                )
                .arg(
                    Arg::new("scaled")
                        .help("Scaled value for the minhash")
                        .long("scaled")
                        .value_parser(clap::value_parser!(u64))
                        .required(true),
                )
                .arg(
                    Arg::new("k")
                        .help("Kmer size for the minhash")
                        .long("k")
                        .value_parser(clap::value_parser!(u32))
                        .required(true),
                )
                .arg(
                    Arg::new("seed")
                        .help("Seed for the minhash")
                        .long("seed")
                        .value_parser(clap::value_parser!(u64))
                        .default_value("42"),
                )
                .arg(
                    Arg::new("similarity_threshold")
                        .help("Similarity threshold for the minhash")
                        .long("similarity-threshold")
                        .value_parser(clap::value_parser!(f64))
                        .required(true),
                )
                .arg(
                    Arg::new("chunk_size")
                        .help("Chunk size for the parallelization")
                        .long("chunk-size")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("10000"),
                )
                .arg(
                    Arg::new("branch_factor")
                        .help("Branching factor for the tree")
                        .long("branch-factor")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("64"),
                )
                .arg(
                    Arg::new("split_apart_taxid_dir_name")
                        .help("split_apart_taxid_dir_name")
                        .long("split-apart-taxid-dir-name")
                        .default_value("split_taxid_dir"),
                )
                .arg(
                    Arg::new("is_protein_fasta")
                        .help("Is protein fasta")
                        .long("is-protein-fasta")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("enable_sequence_retention_logging")
                        .help("Enable sequence retention logging")
                        .long("enable-sequence-retention-logging")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("logging_contained_in_tree_fn")
                        .help("Logging file for containment in the tree")
                        .long("logging-contained-in-tree-fn"),
                )
                .arg(
                    Arg::new("logging_contained_in_chunk_fn")
                        .help("Logging file for containment in the chunk")
                        .long("logging-contained-in-chunk-fn"),
                ),
        )
        .subcommand(
            Command::new("fasta-compress-end-to-end")
                .about("fasta_compress_end_to_end")
                .arg(
                    Arg::new("input_fasta")
                        .help("Input file to be ordered")
                        .long("input-fasta")
                        .required(true),
                )
                .arg(
                    Arg::new("accession_mapping_files")
                        .help("Accession mapping files")
                        .long("accession-mapping-files")
                        .num_args(1..=4) // Allows the flag to appear multiple times
                        .required(true),
                )
                .arg(
                    Arg::new("temp_file_output_dir")
                        .help("temp_file_output_dir")
                        .long("temp-file-output-dir")
                        .required(true),
                )
                .arg(
                    Arg::new("output_fasta")
                        .help("Output file for the ordered fasta")
                        .long("output-fasta")
                        .required(true),
                )
                .arg(
                    Arg::new("scaled")
                        .help("Scaled value for the minhash")
                        .long("scaled")
                        .value_parser(clap::value_parser!(u64))
                        .required(true),
                )
                .arg(
                    Arg::new("k")
                        .help("Kmer size for the minhash")
                        .long("k")
                        .value_parser(clap::value_parser!(u32))
                        .required(true),
                )
                .arg(
                    Arg::new("seed")
                        .help("Seed for the minhash")
                        .long("seed")
                        .value_parser(clap::value_parser!(u64))
                        .default_value("42"),
                )
                .arg(
                    Arg::new("similarity_threshold")
                        .help("Similarity threshold for the minhash")
                        .long("similarity-threshold")
                        .value_parser(clap::value_parser!(f64))
                        .required(true),
                )
                .arg(
                    Arg::new("chunk_size")
                        .help("Chunk size for the parallelization")
                        .long("chunk-size")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("10000"),
                )
                .arg(
                    Arg::new("branch_factor")
                        .help("Branching factor for the tree")
                        .long("branch-factor")
                        .value_parser(clap::value_parser!(usize))
                        .default_value("64"),
                )
                .arg(
                    Arg::new("is_protein_fasta")
                        .help("Is protein fasta")
                        .long("is-protein-fasta")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("enable_sequence_retention_logging")
                        .help("Enable sequence retention logging")
                        .long("enable-sequence-retention-logging")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("logging_contained_in_tree_fn")
                        .help("Logging file for containment in the tree")
                        .long("logging-contained-in-tree-fn"),
                )
                .arg(
                    Arg::new("logging_contained_in_chunk_fn")
                        .help("Logging file for containment in the chunk")
                        .long("logging-contained-in-chunk-fn"),
                ),
        )
        .get_matches();

    match matches.subcommand() {
        Some(("sort-fasta-by-sequence-length", sub_m)) => {
            let input_fasta = sub_m.get_one::<String>("input_fasta").unwrap();
            let output_fasta = sub_m.get_one::<String>("output_fasta").unwrap();
            sort_fasta_by_sequence_length(input_fasta, output_fasta).unwrap();
        }
        Some(("sort-taxid-dir-by-sequence-length", sub_m)) => {
            let input_taxid_dir = sub_m.get_one::<String>("input_taxid_dir").unwrap();
            let output_taxid_dir = sub_m.get_one::<String>("output_taxid_dir").unwrap();
            sort_taxid_dir_by_sequence_length(&input_taxid_dir, &output_taxid_dir);
        }
        Some(("shuffle-fasta", sub_m)) => {
            let input_fasta = sub_m.get_one::<String>("input_fasta").unwrap();
            let output_fasta = sub_m.get_one::<String>("output_fasta").unwrap();
            shuffle_fasta_by_sequence_index(&input_fasta, &output_fasta);
        }
        Some(("break-into-individual-taxids-only", sub_m)) => {
            let input_fasta = sub_m.get_one::<String>("input_fasta").unwrap();
            let accession_mapping_files: Vec<String> = sub_m
                .get_many("accession_mapping_files")
                .expect("Error getting accession mapping files")
                .cloned()
                .collect();
            let output_dir = sub_m.get_one::<String>("output_dir").unwrap();
            split_accessions_by_taxid(input_fasta, accession_mapping_files, output_dir);
        }
        Some(("count-accessions-by-taxid", sub_m)) => {
            let input_taxid_dir = sub_m.get_one::<String>("input-taxid-dir").unwrap();
            let output_summary_path = sub_m.get_one::<String>("output-summary-path").unwrap();
            count_accessions_by_taxid(input_taxid_dir, output_summary_path).unwrap();
        }
        Some(("fasta-compress-from-fasta-skip-split-by-taxid", sub_m)) => {
            let input_fasta = sub_m.get_one::<String>("input_fasta").unwrap();
            let output_fasta = sub_m.get_one::<String>("output_fasta").unwrap();
            let scaled = sub_m.get_one("scaled").unwrap();
            let k = sub_m.get_one("k").unwrap();
            let seed = sub_m.get_one("seed").unwrap();
            let similarity_threshold = sub_m.get_one("similarity_threshold").unwrap();
            let chunk_size = sub_m.get_one("chunk_size").unwrap();
            let branch_factor = sub_m.get_one("branch_factor").unwrap();
            let is_protein_fasta = sub_m.get_flag("is_protein_fasta");
            let enable_sequence_retention_logging =
                sub_m.get_flag("enable_sequence_retention_logging");

            let mut logging_contained_in_tree_fn = "";
            let mut logging_contained_in_chunk_fn = "";

            if enable_sequence_retention_logging {
                // log discarded, retained, containment
                logging_contained_in_tree_fn = sub_m
                    .get_one::<String>("logging_contained_in_tree_fn")
                    .unwrap();
                logging_contained_in_chunk_fn = sub_m
                    .get_one::<String>("logging_contained_in_chunk_fn")
                    .unwrap();
                logging::initialize_tsv(
                    logging_contained_in_tree_fn,
                    vec!["discarded", "retained", "containment"],
                );
                logging::initialize_tsv(
                    logging_contained_in_chunk_fn,
                    vec!["discarded", "retained", "containment"],
                );
            }

            fasta_compress_from_fasta_skip_split_by_taxid(
                input_fasta,
                output_fasta,
                *scaled,
                *k,
                *seed,
                *similarity_threshold,
                *chunk_size,
                *branch_factor,
                is_protein_fasta,
                enable_sequence_retention_logging,
                logging_contained_in_tree_fn,
                logging_contained_in_chunk_fn,
            );
        }
        Some(("fasta-compress-from-taxid-dir", sub_m)) => {
            let input_fasta_dir = sub_m.get_one::<String>("input_fasta_dir").unwrap();
            let output_fasta = sub_m.get_one::<String>("output_fasta").unwrap();
            let scaled = sub_m.get_one("scaled").unwrap();
            let k = sub_m.get_one("k").unwrap();
            let seed = sub_m.get_one("seed").unwrap();
            let similarity_threshold = sub_m.get_one("similarity_threshold").unwrap();
            let split_apart_taxid_dir_name = sub_m
                .get_one::<String>("split_apart_taxid_dir_name")
                .unwrap();
            let enable_sequence_retention_logging =
                sub_m.get_flag("enable_sequence_retention_logging");
            let chunk_size = sub_m.get_one("chunk_size").unwrap();
            let branch_factor = sub_m.get_one("branch_factor").unwrap();
            let is_protein_fasta = sub_m.get_one("is_protein_fasta").unwrap();

            let mut logging_contained_in_tree_fn = "";
            let mut logging_contained_in_chunk_fn = "";

            if enable_sequence_retention_logging {
                // log discarded, retained, containment
                logging_contained_in_tree_fn = sub_m
                    .get_one::<String>("logging_contained_in_tree_fn")
                    .unwrap();
                logging_contained_in_chunk_fn = sub_m
                    .get_one::<String>("logging_contained_in_chunk_fn")
                    .unwrap();
                logging::initialize_tsv(
                    logging_contained_in_tree_fn,
                    vec!["discarded", "retained", "containment"],
                );
                logging::initialize_tsv(
                    logging_contained_in_chunk_fn,
                    vec!["discarded", "retained", "containment"],
                );
            }

            fasta_compress_from_taxid_dir(
                input_fasta_dir,
                output_fasta,
                *scaled,
                *k,
                *seed,
                *similarity_threshold,
                *chunk_size,
                *branch_factor,
                *is_protein_fasta,
                split_apart_taxid_dir_name,
                enable_sequence_retention_logging,
                logging_contained_in_tree_fn,
                logging_contained_in_chunk_fn,
            );
        }
        Some(("fasta-compress-end-to-end", sub_m)) => {
            println!("fasta_compress_end_to_end");
            let input_fasta = sub_m.get_one::<String>("input_fasta").unwrap();
            println!("input_fasta: {}", input_fasta);
            let accession_mapping_files: Vec<String> = sub_m
                .get_many("accession_mapping_files")
                .expect("Error getting accession mapping files")
                .cloned()
                .collect();
            let temp_file_output_dir = sub_m.get_one::<String>("temp_file_output_dir").unwrap();
            let output_fasta = sub_m.get_one::<String>("output_fasta").unwrap();
            let scaled = sub_m.get_one("scaled").unwrap();
            let k = sub_m.get_one("k").unwrap();
            let seed = sub_m.get_one("seed").unwrap();
            let similarity_threshold = sub_m.get_one("similarity_threshold").unwrap();
            let enable_sequence_retention_logging =
                sub_m.get_flag("enable_sequence_retention_logging");
            let chunk_size = sub_m.get_one("chunk_size").unwrap();
            let branch_factor = sub_m.get_one("branch_factor").unwrap();
            let is_protein_fasta = sub_m.get_one("is_protein_fasta").unwrap();

            let mut logging_contained_in_tree_fn = "";
            let mut logging_contained_in_chunk_fn = "";
            if enable_sequence_retention_logging {
                // log discarded, retained, containment
                logging_contained_in_tree_fn = sub_m
                    .get_one::<String>("logging_contained_in_tree_fn")
                    .unwrap();
                logging_contained_in_chunk_fn = sub_m
                    .get_one::<String>("logging_contained_in_chunk_fn")
                    .unwrap();
                logging::initialize_tsv(
                    logging_contained_in_tree_fn,
                    vec!["discarded", "retained", "containment"],
                );
                logging::initialize_tsv(
                    logging_contained_in_chunk_fn,
                    vec!["discarded", "retained", "containment"],
                );
            }

            fasta_compress_end_to_end(
                input_fasta,
                accession_mapping_files,
                output_fasta,
                temp_file_output_dir,
                *scaled,
                *k,
                *seed,
                *similarity_threshold,
                *chunk_size,
                *branch_factor,
                *is_protein_fasta,
                enable_sequence_retention_logging,
                logging_contained_in_tree_fn,
                logging_contained_in_chunk_fn,
            );
        }
        _ => unreachable!("Command not found"),
    }
}

// #[test]
