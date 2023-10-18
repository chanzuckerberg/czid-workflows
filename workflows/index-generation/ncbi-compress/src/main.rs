use clap::{Parser, ArgGroup};

use ncbi_compress::ncbi_compress::ncbi_compress::fasta_compress;
use ncbi_compress::logging::logging;




#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]

#[clap(group = ArgGroup::new("break-apart-by-taxid").required(true))]
#[clap(group = ArgGroup::new("compress").required(false))]
struct Args {
    /// Path to the input fasta file
    #[arg(short, long, required = true, group = "compress", group="break-apart-by-taxid")]
    input_fasta: String,

    // Path to directory of input fasta files broken up by taxids
    #[arg(short, long, required = true, group = "compress")]
    input_fasta_dir: String,

    /// Path to the output fasta file
    #[arg(short, long, required = true, group = "compress")]
    output_fasta: String,

    /// Path to the accession to taxid csv file
    /// At least one required if not skipping split by taxid
    #[arg(short, long, group = "compress", group = "break-apart-by-taxid")]
    accession_mapping_files: Vec<String>,

    /// Taxids to drop from the output
    #[arg(long, group = "compress", group = "break-apart-by-taxid")]
    taxids_to_drop: Vec<u64>,

    /// Scaled value for the minhash
    /// (default: 1000)
    #[arg(short, long, default_value = "1000", group="compress")]
    scaled: u64,

    /// Kmer size for the minhash
    /// (default: 31)
    #[arg(short, long, default_value = "31", group="compress")]
    k: u32,

    /// Seed for the minhash
    /// (default: 42)
    #[arg(long, default_value = "42", group="compress")]
    seed: u64,

    /// Similarity threshold for the minhash
    /// (default: 0.6)
    /// (must be between 0 and 1)
    #[arg(short = 't', long, default_value = "0.6", group="compress")]
    similarity_threshold: f64,

    /// Chunk size for the parallelization
    /// (default: 1000)
    #[arg(short, long, default_value = "10000", group`="compress")]
    chunk_size: usize,

    /// Branching factor for the tree
    /// (default: 1000)
    #[arg(short, long, default_value = "1000", group="compress")]
    branch_factor: usize,

    /// Skip splitting the input fasta by taxid
    /// (default: false)
    #[arg(long, default_value = "false", group="compress")]
    skip_split_by_taxid: bool,

    /// Logging file for containment in the tree
    /// (default: logging_contained_in_tree.tsv)
    #[arg(long, default_value = "logging_contained_in_tree.tsv", group="compress")]
    logging_contained_in_tree_fn: String,

    /// Logging file for containment in the chunk
    /// (default: logging_contained_in_chunk.tsv)
    #[arg(long, default_value = "logging_contained_in_chunk.tsv", group="compress")]
    logging_contained_in_chunk_fn: String,

    // Enable logging
    #[arg(long, default_value = "false", group="compress")]
    enable_sequence_retention_logging: bool,

    // temp file directory for the split fasta
    #[arg(long, default_value = "/mnt", group="compress", group="break-apart-by-taxid")]
    temp_file_output_dir: String,

    // is protein fasta
    #[arg(long, default_value = "false", group="compress")]
    is_protein_fasta: bool,

    // break into individual taxids
    #[arg(long, default_value = "false")]
    break_into_individual_taxids_only: bool,

}


pub fn main() {

    logging::init_stdout_logging();

    let args = Args::parse();

    if args.enable_sequence_retention_logging {
        // log discarded, retained, containment
        logging::initialize_tsv(&args.logging_contained_in_tree_fn, vec!["discarded", "retained", "containment"]);
        logging::initialize_tsv(&args.logging_contained_in_chunk_fn, vec!["discarded", "retained", "containment"]);
    }

    fasta_compress(
        args.input_fasta,
        args.accession_mapping_files,
        args.output_fasta,
        args.taxids_to_drop,
        args.scaled,
        args.k,
        args.seed,
        args.similarity_threshold,
        args.chunk_size,
        args.branch_factor,
        &args.temp_file_output_dir,
        args.skip_split_by_taxid,
        args.is_protein_fasta,
        args.enable_sequence_retention_logging,
        &args.logging_contained_in_tree_fn,
        &args.logging_contained_in_chunk_fn,
        args.break_into_individual_taxids_only,
    );

}
