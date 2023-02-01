version 1.1

# Build host genome indexes for host_filter.wdl (2022 version)
# - Bowtie2 (genome)
# - HISAT2 (genome + splice junctions)
# - kallisto (transcriptome)
# - minimap2 (used not in short-read-mngs host filtering, but rather the ONT equivalent)
# - STAR (used in old version of short-read-mngs host filtering, kept temporarily so we can support both)
# ERCC sequences are spiked-in to all three indexes. Lastly takes an array of other spike-ins for
# the Bowtie2 and HISAT2 indexes.
# Warning: HISAT2 requires huge RAM to build the spliced index (>200G for human).
#          But the index file size and aligner memory usage are relatively small.
workflow host_filter_indexing {
  input {
    String genome_name

    # host genomic DNA
    File genome_fasta_gz
    # host transcript models on the above genomic DNA (for HISAT2 spliced alignment)
    File? transcripts_gtf_gz
    # host transcript sequences (for kallisto)
    Array[File] transcripts_fasta_gz = []

    # ERCC sequences to spike in to the genome and transcript indexes
    File ERCC_fasta_gz

    # Additional FASTA file(s) to spike into the Bowtie2 & HISAT2 indexes (e.g. EBV, phiX)
    # Sequence names must be unique among all FASTAs!
    Array[File] other_fasta_gz = []

    String docker
  }

  call ensure_gz as genome_fasta {
    # accommodate uncompressed genome_fasta_gz; this makes it more convenient to use some of our
    # existing host genome FASTAs which we archived without compression.
    input:
    maybe_gz = genome_fasta_gz,
    docker
  }

  call bowtie2_build {
    input:
    fasta_gz = flatten([[genome_fasta.gz, ERCC_fasta_gz], other_fasta_gz]),
    genome_name, docker
  }

  call hisat2_build {
    input:
    fasta_gz = flatten([[genome_fasta.gz, ERCC_fasta_gz], other_fasta_gz]),
    transcripts_gtf_gz, genome_name, docker
  }

  call kallisto_index {
    input:
    transcripts_fasta_gz = flatten([transcripts_fasta_gz, [ERCC_fasta_gz]]),
    genome_name, docker
  }

  call minimap2_index as minimap2_index_dna {
    input:
    fasta_gz = flatten([[genome_fasta.gz, ERCC_fasta_gz], other_fasta_gz]),
    nucleotide_type = "dna",
    genome_name, docker
  }

  call minimap2_index as minimap2_index_rna {
    input:
    fasta_gz = flatten([[genome_fasta.gz, ERCC_fasta_gz], other_fasta_gz]),
    nucleotide_type = "rna",
    genome_name, docker
  }

  call star_generate {
    input:
    fasta_gz = flatten([[genome_fasta.gz, ERCC_fasta_gz], other_fasta_gz]),
    transcripts_gtf_gz, genome_name, docker
  }

  output {
    File bowtie2_index_tar = bowtie2_build.index_tar
    File hisat2_index_tar = hisat2_build.index_tar
    File kallisto_idx = kallisto_index.idx
    File minimap2_dna_mmi = minimap2_index_dna.index_mmi
    File minimap2_rna_mmi = minimap2_index_rna.index_mmi
    File star_genome_tar = star_generate.star_genome_tar

    # also output the input files, to facilitate archival/provenance
    File original_genome_fasta_gz = genome_fasta.gz
    File? original_transcripts_gtf_gz = transcripts_gtf_gz
    Array[File] original_transcripts_fasta_gz = transcripts_fasta_gz
    File original_ERCC_fasta_gz = ERCC_fasta_gz
    Array[File] original_other_fasta_gz = other_fasta_gz
  }
}

task ensure_gz {
  input {
    File maybe_gz
    String docker
  }

  String name = basename(maybe_gz)

  command <<<
    set -euxo pipefail
    mkdir ans
    if gzip -t '~{maybe_gz}'; then
      cp '~{maybe_gz}' ans/
    else
      pigz -c -p 4 '~{maybe_gz}' > 'ans/~{name}.gz'
    fi
  >>>

  output {
      File gz = glob("ans/*")[0]
  }

  runtime {
      docker: docker
      cpu: 4
      memory: "4GiB"
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
    memory: "16GiB"
  }
}

task minimap2_index {
  input {
    Array[File] fasta_gz
    String genome_name
    String nucleotide_type

    String docker
  }

  command <<<
    set -euxo pipefail
    TMPDIR=${TMPDIR:-/tmp}

    all_fasta="$TMPDIR/all.fasta"
    pigz -dc ~{sep(' ',fasta_gz)} > "$all_fasta"

    if [ "~{nucleotide_type}" == "dna" ]; then
        >&2 minimap2 -x map-ont -d '~{genome_name}_{nucleotide_type}.mmi' "$all_fasta"
    else
        >&2 minimap2 -x splice -d '~{genome_name}_{nucleotide_type}.mmi' "$all_fasta"
    fi
    >&2 ls -l
  >>>

  output {
      File index_mmi = "~{genome_name}_{nucleotide_type}.mmi"
  }

  runtime {
      docker: docker
      memory: "32GiB"
  }
}

task star_generate {
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

    gtf_flag=""
    if [[ -n '~{transcripts_gtf_gz}' ]]; then
      transcripts_gtf="$TMPDIR/transcripts.gtf"
      pigz -dc '~{transcripts_gtf_gz}' > "$transcripts_gtf"
      gtf_flag = "--sjdbGTFfile \"$transcripts_gtf\""
    fi

    # Make directory for STAR genome
    STAR_GENOME="~{genome_name}_STAR_genome"
    # HACK: we used to support splitting star indexes into many parts, this made things slower
    # Here we generate the index as if it is in many parts, but there is only ever one part for
    # backwards compatibility
    mkdir -p "$STAR_GENOME/part-0"

    STAR \
      --sjdbGTFfile "~{transcripts_gtf_gz}" \
      --runThreadN ~{cpu} \
      --runMode genomeGenerate \
      --genomeFastaFiles "$all_fasta" \
      --limitGenomeGenerateRAM 64000000000 \
      --genomeDir "$STAR_GENOME/part-0" $gtf_flag

    # create a parts.txt file for backwards compatibility
    echo 1 > "$STAR_GENOME/parts.txt"

    # tar STAR genome
    tar cvf "$STAR_GENOME.tar" -C $(pwd) $STAR_GENOME
  >>>

  output {
      File star_genome_tar = "~{genome_name}_STAR_genome.tar"
  }

  runtime {
      docker: docker
      cpu: cpu
      memory: "64GiB"
  }
}
