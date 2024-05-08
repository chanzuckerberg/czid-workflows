# NCBI Compress

### install: 
1. [install rust and cargo](https://doc.rust-lang.org/cargo/getting-started/installation.html)
2. install ncbi-compress: `cargo build --release` (only necessary if you want a release vesion of the code that is more performant)

### test:
* to run all tests: `cargo test  -- --nocapture`
* to run only tests for one file: `cargo test util -- --nocapture`


### usage:

[Lucidchart](https://lucid.app/lucidchart/19dd8822-fec8-4c9b-b0cb-8199ea0c9ae4/edit?viewport_loc=-6806%2C-2082%2C13130%2C6583%2C0_0&invitationId=inv_275ac20b-9677-470f-8238-d10208629a4c) for algorithm overview with compute resources

#### sorting:

```
# individual fasta
target/release/ncbi-compress sort-fasta-by-sequence-length \
    --input-fasta nt \
    --output-fasta nt_sorted.fa

# directory of fastas that need to be sorted
target/release/ncbi-compress sort-taxid-dir-by-sequence-length \
    --input-taxid-dir  reads_by_taxid_nt \
    --output-taxid-dir sorted_reads_by_taxid_nt
```

#### shuffling
```
target/release/ncbi-compress shuffle-fasta \
    --input-fasta nt \
    --output-fasta nt_shuffled.fa
```

#### sort fasta into individual taxids
```
target/release/ncbi-compress break-into-individual-taxids-only \
    --input-fasta nt \
    --accession-mapping-files \
        ./mapping_files/nucl_gb.accession2taxid \
        ./mapping_files/nucl_wgs.accession2taxid \
    --output-dir reads_by_taxid_nt
```

#### count accessions by taxid
```
target/release/ncbi-compress count-accessions-by-taxid \
    --input-taxid-dir reads_by_taxid_nt \
    --output-summary-path output_counts.tsv
```

#### compress fasta or directory of fastas:
```
# compress from directory
target/release/ncbi-compress fasta-compress-from-taxid-dir \
    --input-fasta-dir sorted_reads_by_taxid_nt \
    --output-fasta nt_compressed_0.9.fa \
    --k 31 \
    --scaled 1000 \
    --similarity-threshold 0.9 \
    --split-apart-taxid-dir-name split_nt_2 # this is used for if we have any taxids with accessions over

# compress single fasta
target/release/ncbi-compress fasta-compress-from-fasta-skip-split-by-taxid \
    --input-fasta nt \
    --output-fasta nt_compressed_0.9.fa \
    --k 31 \
    --scaled 1000 \
    --similarity-threshold 0.9 \
    --split-apart-taxid-dir-name split_nt_2 # this is used for if we have any taxids with accessions over

# fasta compress from end to end (will sort into taxids, compress, and shuffle)
target/release/ncbi-compress fasta-compress-end-to-end \
    --input-fasta nt \
    --accession-mapping-files \
        ./mapping_files/nucl_gb.accession2taxid \
        ./mapping_files/nucl_wgs.accession2taxid \
    --output-fasta nt_compressed_0.9.fa \
    --k 31 \
    --scaled 1000 \
    --similarity-threshold 0.9 \
    --split-apart-taxid-dir-name split_nt_2 # this is used for if we have any taxids with accessions over
```
