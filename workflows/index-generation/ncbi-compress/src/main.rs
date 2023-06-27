use std::io::Write;
use std::ops::AddAssign;
use std::path::Path;
use std::{borrow::BorrowMut, fs};

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
use tempdir::TempDir;
use trie_rs::{Trie, TrieBuilder};

/// A trie that stores u64 values
struct TrieStore {
    trie: Trie<u8>,
}

impl TrieStore {
    pub fn get(&self, key: &str) -> Option<u64> {
        self.trie
            .common_prefix_search(key.as_bytes())
            .first()
            .map(|bytes| {
                let mut bytes = bytes.to_vec();
                let value_bytes = bytes.split_off(bytes.len() - 8);
                u64::from_be_bytes(value_bytes.try_into().unwrap())
            })
    }
}

struct TrieStoreBuilder {
    builder: TrieBuilder<u8>,
}

impl TrieStoreBuilder {
    pub fn new() -> Self {
        TrieStoreBuilder {
            builder: TrieBuilder::new(),
        }
    }

    pub fn push(&mut self, key: &str, value: u64) {
        let mut key = key.as_bytes().to_vec();
        key.extend_from_slice(&value.to_be_bytes());
        self.builder.push(key);
    }

    pub fn build(self) -> TrieStore {
        let trie = self.builder.build();
        TrieStore { trie }
    }
}

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

    pub fn insert(&mut self, hash: KmerMinHash) -> Result<(), SourmashError> {
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

    pub fn contains(
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

fn split_accessions_by_taxid<P: AsRef<Path> + std::fmt::Debug, Q: AsRef<Path> + std::fmt::Debug>(
    input_fasta_path: P,
    mapping_file_path: Vec<Q>,
    taxids_to_drop: &Vec<u64>,
) -> TempDir {
    log::info!("Creating accession to taxid mapping");
    let taxid_dir = TempDir::new("accessions_by_taxid").unwrap();
    let reader = fasta::Reader::from_file(&input_fasta_path).unwrap();
    let mut builder = TrieBuilder::new();
    reader.records().enumerate().for_each(|(i, result)| {
        let record = result.unwrap();
        let accession_id = record.id().split_whitespace().next().unwrap();
        builder.push(accession_id);
        if i % 10_000 == 0 {
            log::info!("  Processed {} accessions", i);
        }
    });
    log::info!(" Started building accession trie");
    let accessions_trie = builder.build();
    log::info!(" Finished building accession trie");

    let mut builder = TrieStoreBuilder::new();
    mapping_file_path.iter().for_each(|mapping_file_path| {
        log::info!(" Processing mapping file {:?}", mapping_file_path);
        let reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_path(mapping_file_path)
            .unwrap();
        let mut added = 0;
        reader.into_records().enumerate().for_each(|(i, result)| {
            if i % 10_000 == 0 {
                log::info!("  Processed {} mappings, added {}", i, added);
            }

            let record = result.unwrap();
            let accession = record[0].as_bytes();
            let accession_no_version = accession.splitn(2, |b| *b == b'.').next().unwrap();

            // Only output mappings if the accession is in the source files

            // If using the prot.accession2taxid.FULL file
            let (accession, taxid) =
                if record.len() < 3 && accessions_trie.exact_match(accession_no_version) {
                    // Remove the version number and the taxid will be at index 1
                    (
                        std::str::from_utf8(accession_no_version).unwrap(),
                        record[1].parse::<u64>().unwrap(),
                    )
                } else if accessions_trie.exact_match(accession) {
                    // Otherwise there is a versionless accession ID at index 0 and the taxid is at index 2
                    (
                        std::str::from_utf8(accession).unwrap(),
                        record[2].parse::<u64>().unwrap(),
                    )
                } else {
                    return;
                };

            if !taxids_to_drop.contains(&taxid) {
                added += 1;
                builder.push(accession, taxid);
            }
        });
        log::info!(
            " Finished Processing mapping file {:?}, added {} mappings",
            mapping_file_path,
            added
        );
    });
    log::info!(" Started building accession to taxid trie");
    let accession_to_taxid = builder.build();
    log::info!("Finished building accession to taxid trie");

    log::info!("Splitting accessions by taxid");
    let reader = fasta::Reader::from_file(&input_fasta_path).unwrap();
    for (i, record) in reader.records().enumerate() {
        let record = record.unwrap();
        let accession_id = record.id().split_whitespace().next().unwrap();
        let taxid = if let Some(taxid) = accession_to_taxid.get(accession_id) {
            taxid
        } else {
            continue;
        };

        let file_path = taxid_dir.path().join(format!("{}.fasta", taxid));
        let file = fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(file_path)
            .unwrap();
        let mut writer = fasta::Writer::new(file);
        writer.write_record(&record).unwrap();

        if i % 10_000 == 0 {
            log::info!("  Split {} accessions", i);
        }
    }
    log::info!("Finished splitting accessions by taxid");

    taxid_dir
}

fn fasta_compress_taxid<P: AsRef<Path> + std::fmt::Debug>(
    input_fasta_path: P,
    writer: &mut fasta::Writer<std::fs::File>,
    scaled: u64,
    k: u32,
    seed: u64,
    similarity_threshold: f64,
    chunk_size: usize,
    branch_factor: usize,
    accession_count: &mut u64,
    unique_accession_count: &mut u64,
) {
    let reader = fasta::Reader::from_file(&input_fasta_path).unwrap();
    let mut tree = MinHashTree::new(branch_factor);

    let mut records_iter = reader.records();
    let mut chunk = records_iter
        .borrow_mut()
        .take(chunk_size)
        .collect::<Vec<_>>();
    while chunk.len() > 0 {
        let chunk_signatures = chunk
            .par_iter()
            .filter_map(|r| {
                let record = r.as_ref().unwrap();
                let mut hash =
                    KmerMinHash::new(scaled, k, HashFunctions::murmur64_DNA, seed, false, 0);
                hash.add_sequence(record.seq(), true).unwrap();
                // Run an initial similarity check here against the full tree, this is slow so we can parallelize it
                if tree.contains(&hash, similarity_threshold).unwrap() {
                    None
                } else {
                    Some((record, hash))
                }
            })
            .collect::<Vec<_>>();

        let mut tmp = Vec::with_capacity(chunk_signatures.len() / 2);
        for (record, hash) in chunk_signatures {
            accession_count.add_assign(1);
            // Perform a faster similarity check over just this chunk because we may have similarities within a chunk
            let similar = tmp
                .par_iter()
                .any(|other| containment(&hash, other).unwrap() >= similarity_threshold);

            if !similar {
                unique_accession_count.add_assign(1);
                tmp.push(hash);
                writer.write_record(record).unwrap();

                if *unique_accession_count % 10_000 == 0 {
                    log::info!(
                        "Processed {} accessions, {} unique",
                        accession_count,
                        unique_accession_count
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

fn fasta_compress<P: AsRef<Path> + std::fmt::Debug>(
    input_fasta_path: P,
    accession_mapping_files: Vec<P>,
    output_fasta_path: P,
    taxids_to_drop: Vec<u64>,
    scaled: u64,
    k: u32,
    seed: u64,
    similarity_threshold: f64,
    chunk_size: usize,
    branch_factor: usize,
) {
    log::info!("Splitting accessions by taxid");
    let taxid_dir =
        split_accessions_by_taxid(&input_fasta_path, accession_mapping_files, &taxids_to_drop);
    log::info!("Finished splitting accessions by taxid");
    let mut writer = fasta::Writer::to_file(output_fasta_path).unwrap();

    log::info!("Starting compression by taxid");
    let mut accession_count = 0;
    let mut unique_accession_count = 0;
    for (i, entry) in fs::read_dir(taxid_dir.path()).unwrap().enumerate() {
        let entry = entry.unwrap();
        let path = entry.path();
        let input_fasta_path = path.to_str().unwrap();
        fasta_compress_taxid(
            input_fasta_path,
            &mut writer,
            scaled,
            k,
            seed,
            similarity_threshold,
            chunk_size,
            branch_factor,
            &mut accession_count,
            &mut unique_accession_count,
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

    taxid_dir.close().unwrap();
    log::info!(
        "Finished compression by taxid, {} accessions, {} uniqe accessions",
        accession_count,
        unique_accession_count
    );
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

    /// Path to the accession to taxid csv file
    /// At least one required
    #[arg(short, long, required = true)]
    accession_mapping_files: Vec<String>,

    /// Taxids to drop from the output
    #[arg(long)]
    taxids_to_drop: Vec<u64>,

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
    #[arg(long, default_value = "42")]
    seed: u64,

    /// Similarity threshold for the minhash
    /// (default: 0.6)
    /// (must be between 0 and 1)
    #[arg(short = 't', long, default_value = "0.6")]
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
        args.accession_mapping_files,
        args.output_fasta,
        args.taxids_to_drop,
        args.scaled,
        args.k,
        args.seed,
        args.similarity_threshold,
        args.chunk_size,
        args.branch_factor,
    );
}
