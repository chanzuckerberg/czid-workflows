use std::borrow::BorrowMut;
use std::io::Write;
use std::path::Path;

use bio::io::fasta;
use chrono::Local;
use clap::Parser;
use env_logger::Builder;
use log::LevelFilter;
use rayon::prelude::*;
use sourmash::encodings::HashFunctions;
use sourmash::errors::SourmashError;
use sourmash::signature::SigsTrait;
use sourmash::sketch::minhash::KmerMinHash;

fn containment(needle: &KmerMinHash, haystack: &KmerMinHash) -> Result<f64, SourmashError> {
    let (intersect_size, _) = needle.intersection_size(haystack)?;
    Ok(intersect_size as f64 / needle.mins().len() as f64)
}

struct MinHashTreeNode {
    own: KmerMinHash,
    children_aggregate: KmerMinHash,
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

    fn insert(&mut self, hash: KmerMinHash) -> Result<(), SourmashError> {
        let node = MinHashTreeNode {
            own: hash.clone(),
            children_aggregate: hash.clone(),
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

    fn contains(
        &self,
        hash: &KmerMinHash,
        similarity_threshold: f64,
    ) -> Result<bool, SourmashError> {
        if self.nodes.is_empty() {
            return Ok(false);
        }

        let mut to_visit = vec![0];
        while !to_visit.is_empty() {
            let found = to_visit.par_iter().any(|node_idx| {
                let node = self.nodes.get(*node_idx).unwrap();
                containment(hash, &node.own).unwrap() >= similarity_threshold
            });

            if found {
                return Ok(true);
            }

            to_visit = to_visit
                .par_iter()
                .flat_map(|node_idx| {
                    let node = self.nodes.get(*node_idx).unwrap();
                    if containment(hash, &node.children_aggregate).unwrap() >= similarity_threshold
                    {
                        self.child_idxes(*node_idx)
                    } else {
                        vec![]
                    }
                })
                .collect();
        }
        Ok(false)
    }
}

fn fasta_compress<P: AsRef<Path> + std::fmt::Debug>(
    input_fasta_path: P,
    output_fasta_path: P,
    scaled: u64,
    k: u32,
    seed: u64,
    similarity_threshold: f64,
    chunk_size: usize,
    branch_factor: usize,
) {
    let reader = fasta::Reader::from_file(&input_fasta_path).unwrap();
    let mut writer = fasta::Writer::to_file(output_fasta_path).unwrap();

    let mut tree = MinHashTree::new(branch_factor);
    let mut unique_accessions: i64 = 0;

    let mut records_iter = reader.records().enumerate();
    let mut chunk = records_iter
        .borrow_mut()
        .take(chunk_size)
        .collect::<Vec<_>>();
    while chunk.len() > 0 {
        let chunk_signatures = chunk
            .par_iter()
            .filter_map(|(i, r)| {
                let record = r.as_ref().unwrap();
                let mut hash =
                    KmerMinHash::new(scaled, k, HashFunctions::murmur64_DNA, seed, false, 0);
                hash.add_sequence(record.seq(), true).unwrap();
                // Run an initial similarity check here against the full tree, this is slow so we can parallelize it
                if tree.contains(&hash, similarity_threshold).unwrap() {
                    None
                } else {
                    Some((*i, record, hash))
                }
            })
            .collect::<Vec<_>>();

        let mut tmp = Vec::with_capacity(chunk_signatures.len() / 2);
        for (i, record, hash) in chunk_signatures {
            // Perform a faster similarity check over just this chunk because we may have similarities within a chunk
            let similar = tmp
                .par_iter()
                .any(|other| containment(&hash, other).unwrap() >= similarity_threshold);

            if !similar {
                unique_accessions += 1;
                tmp.push(hash);
                writer.write_record(record).unwrap();

                if unique_accessions % 10_000 == 0 {
                    log::info!(
                        "Processed {} accessions, {} unique",
                        i + 1,
                        unique_accessions
                    );
                }
            }
        }
        for hash in tmp {
            tree.insert(hash).unwrap();
        }
        chunk = records_iter
            .borrow_mut()
            .take(chunk_size)
            .collect::<Vec<_>>();
    }
}

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Path to the input fasta file
    #[arg(short, long)]
    input_fasta: String,

    /// Path to the output fasta file
    #[arg(short, long)]
    output_fasta: String,

    /// Scaled value for the minhash
    /// (default: 1000)
    #[arg(short, long, default_value = "1000")]
    scaled: u64,

    /// Kmer size for the minhash
    /// (default: 31)
    #[arg(short, long, default_value = "31")]
    k: u32,

    /// Seed for the minhash
    /// (default: 42)
    #[arg(short, long, default_value = "42")]
    seed: u64,

    /// Similarity threshold for the minhash
    /// (default: 0.6)
    /// (must be between 0 and 1)
    #[arg(short, long, default_value = "0.6")]
    similarity_threshold: f64,

    /// Chunk size for the parallelization
    /// (default: 1000)
    #[arg(short, long, default_value = "10000")]
    chunk_size: usize,

    /// Branching factor for the tree
    /// (default: 1000)
    #[arg(short, long, default_value = "1000")]
    branch_factor: usize,
}

fn main() {
    Builder::new()
        .format(|buf, record| {
            writeln!(
                buf,
                "{} {}: {}",
                record.level(),
                Local::now().format("%Y-%m-%d %H:%M:%S%.3f"),
                record.args()
            )
        })
        .filter(None, LevelFilter::Info)
        .init();

    let args = Args::parse();

    fasta_compress(
        args.input_fasta,
        args.output_fasta,
        args.scaled,
        args.k,
        args.seed,
        args.similarity_threshold,
        args.chunk_size,
        args.branch_factor,
    );
}
