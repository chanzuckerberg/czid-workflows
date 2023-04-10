#!/bin/bash

set -euxo pipefail

PATH_CONTIGS="contigs.fasta"
PATH_COMPREHENSIVE_REPORT="final_reports/comprehensive_AMR_metrics.tsv"
PATH_OUTPUT_SAM="test.sam"
PATH_OUTPUT_BAM="test.bam"

# Create index to enable querying fasta file
samtools faidx $PATH_CONTIGS

# Create SAM header with mock sequence sizes
echo -e "@HD\tVN:1.0\tSO:unsorted" > $PATH_OUTPUT_SAM
awk -F "\t" -v OFS="\t" '{
  contigName = $2
  geneId = $36

  # Ignore header, and lines with no info
  if(NR == 1 || contigName == "" || geneId == "")
    next;

  print "@SQ", "SN:"geneId, "LN:100"
}' $PATH_COMPREHENSIVE_REPORT | sort | uniq >> $PATH_OUTPUT_SAM

# Go through each line of the TSV and output a SAM record
awk -F "\t" -v OFS="\t" -v PATH_CONTIGS="$PATH_CONTIGS" '{
  contigName = $2
  geneId = $36

  # Ignore header, and lines with no info
  if(NR == 1 || contigName == "" || geneId == "")
    next;

  # Contig names here have an additional "_" followed by a number at the end,
  # but the contig name in contigs.fasta does not, so remove it.
  gsub(/_[^_]*$/, "", contigName);

  # Fetch contig sequence
  command = "samtools faidx -n 0 contigs.fasta \""contigName"\" | tail -n 1"
  command | getline contigSequence

  print contigName, 0, geneId, "1", 255, "*", "*", 0, 0, contigSequence, "*"
}' $PATH_COMPREHENSIVE_REPORT >> $PATH_OUTPUT_SAM

# Convert SAM to BAM and index (ignore warning about no CIGAR string; we don't need those)
samtools sort -o $PATH_OUTPUT_BAM $PATH_OUTPUT_SAM 2>/dev/null
samtools index $PATH_OUTPUT_BAM

# Test: samtools view -h test.bam "ARO:3000777|ID:153|Name:adeF|NCBI:CT025801.2" | samtools fasta
