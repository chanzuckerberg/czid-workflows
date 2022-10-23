version 1.1

# Build host genome indexes for host_filter.wdl (2022 version)
# - Bowtie2 (genome)
# - HISAT2 (genome + splice junctions)
# - kallisto (transcriptome)
# - minimap2 (used not in short-read-mngs host filtering, but rather the ONT equivalent)
# ERCC sequences are spiked-in to all three indexes. Lastly takes an array of other spike-ins for
# the Bowtie2 and HISAT2 indexes.
# Warning: HISAT2 requires huge RAM to build the spliced index (>200G for human).
#          But the index file size and aligner memory usage are relatively small.
workflow host_filter_indexing {
  input {
    String genome_name = "GRCh38_ERCC"

    # host genomic DNA
    File genome_fasta_gz = "ftp://ftp.ensembl.org/pub/release-94/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
    # host transcript models on the above genomic DNA (for HISAT2 spliced alignment)
    File? transcripts_gtf_gz = "ftp://ftp.ensembl.org/pub/release-94/gtf/homo_sapiens/Homo_sapiens.GRCh38.94.gtf.gz"
    # host transcript sequences (for kallisto)
    Array[File] transcripts_fasta_gz = ["ftp://ftp.ensembl.org/pub/release-94/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"]

    # ERCC sequences to spike in to the genome and transcript indexes
    File ERCC_fasta_gz

    # Additional FASTA file(s) to spike into the Bowtie2 & HISAT2 indexes (e.g. EBV, phiX)
    # Sequence names must be unique among all FASTAs!
    Array[File] other_fasta_gz = []

    String docker = "ghcr.io/chanzuckerberg/czid-workflows/czid-short-read-mngs-public:ff7af40"
  }

  call bowtie2_build {
    input:
    fasta_gz = flatten([[genome_fasta_gz, ERCC_fasta_gz], other_fasta_gz]),
    genome_name, docker
  }
  
  call hisat2_build {
    input:
    fasta_gz = flatten([[genome_fasta_gz, ERCC_fasta_gz], other_fasta_gz]),
    transcripts_gtf_gz, genome_name, docker
  }

  call kallisto_index {
    input:
    transcripts_fasta_gz = flatten([transcripts_fasta_gz, [ERCC_fasta_gz]]),
    genome_name, docker
  }

  call minimap2_index {
    input:
    fasta_gz = flatten([[genome_fasta_gz, ERCC_fasta_gz], other_fasta_gz]),
    genome_name, docker
  }
  
  output {
    File bowtie2_index_tar = bowtie2_build.index_tar
    File hisat2_index_tar = hisat2_build.index_tar
    File kallisto_idx = kallisto_index.idx
    File minimap2_mmi = minimap2_index.index_mmi

    # also output the input files, to facilitate archival/provenance
    File input_genome_fasta_gz = genome_fasta_gz
    File? input_transcripts_gtf_gz = transcripts_gtf_gz
    Array[File] input_transcripts_fasta_gz = transcripts_fasta_gz
    File input_ERCC_fasta_gz = ERCC_fasta_gz
    Array[File] input_other_fasta_gz = other_fasta_gz
  }
}

task bowtie2_build {
  input {
    Array[File] fasta_gz
    String genome_name
    Int seed = 42

    Int cpu = 16
    String docker
  }
  
  command <<<
    set -euxo pipefail
    TMPDIR=${TMPDIR:-/tmp}

    all_fasta="$TMPDIR/all.fasta"
    pigz -dc ~{sep(' ',fasta_gz)} > "$all_fasta"
    
    mkdir -p "$TMPDIR"'/bt2/~{genome_name}'
    >&2 bowtie2-build --seed ~{seed} --threads ~{cpu} "$all_fasta" "$TMPDIR"'/bt2/~{genome_name}/~{genome_name}'
    >&2 ls -lR "$TMPDIR/bt2"
    env -C "$TMPDIR/bt2" tar c . > '~{genome_name}.bowtie2.tar'
  >>>
  
  output {
      File index_tar = "~{genome_name}.bowtie2.tar"
  }
  
  runtime {
      docker: docker
      cpu: cpu
      memory: "~{cpu*2}GiB"
  }
}

task hisat2_build {
  input {
    Array[File] fasta_gz
    File? transcripts_gtf_gz
    String genome_name
    
    Int cpu = 32
    String docker
  }
  
  command <<<
    set -euxo pipefail
    TMPDIR=${TMPDIR:-/tmp}

    all_fasta="$TMPDIR/all.fasta"
    pigz -dc ~{sep(' ',fasta_gz)} > "$all_fasta"

    mkdir -p "$TMPDIR"'/hisat2/~{genome_name}'
    if [[ -n '~{transcripts_gtf_gz}' ]]; then
      # convert GTF per http://daehwankimlab.github.io/hisat2/howto/
      /hisat2/hisat2_extract_splice_sites.py <(pigz -dc '~{transcripts_gtf_gz}') > "$TMPDIR/genome.ss" & pid=$!
      /hisat2/hisat2_extract_exons.py <(pigz -dc '~{transcripts_gtf_gz}') > "$TMPDIR/genome.exon"
      wait $pid
      >&2 /hisat2/hisat2-build -p 16 \
        --exon "$TMPDIR/genome.exon" --ss "$TMPDIR/genome.ss" \
        "$all_fasta" "$TMPDIR"'/hisat2/~{genome_name}/~{genome_name}'
    else
      >&2 /hisat2/hisat2-build -p 16 "$all_fasta" "$TMPDIR"'/hisat2/~{genome_name}/~{genome_name}'
    fi
    >&2 ls -lR "$TMPDIR/hisat2"
    env -C "$TMPDIR/hisat2" tar c . > '~{genome_name}.hisat2.tar'
  >>>
  
  output {
    File index_tar = "~{genome_name}.hisat2.tar"
  }
  
  runtime {
    docker: docker
    cpu: cpu
    memory: "240G"
  }
}

task kallisto_index {
  input {
    Array[File] transcripts_fasta_gz
    String genome_name

    String docker
  }

  String idx_fn = "~{genome_name}.kallisto.idx"
  command <<<
    set -euxo pipefail
    /kallisto/kallisto index --index '~{idx_fn}' ~{sep(' ',transcripts_fasta_gz)}
    >&2 ls -l
  >>>
  
  output {
    File idx = idx_fn
  }
  
  runtime {
    docker: docker
    memory: "8GiB"
  }
}

task minimap2_index {
  input {
    Array[File] fasta_gz
    String genome_name
    String opts = ""

    String docker
  }

  command <<<
    set -euxo pipefail
    TMPDIR=${TMPDIR:-/tmp}

    all_fasta="$TMPDIR/all.fasta"
    pigz -dc ~{sep(' ',fasta_gz)} > "$all_fasta"

    >&2 minimap2 ~{opts} -d '~{genome_name}.mmi' "$all_fasta"
    >&2 ls -l
  >>>

  output {
      File index_mmi = "~{genome_name}.mmi"
  }

  runtime {
      docker: docker
      memory: "8GiB"
  }
}
