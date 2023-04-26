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
    File? ERCC_fasta_gtf

    # Additional FASTA file(s) to spike into the Bowtie2 & HISAT2 indexes (e.g. EBV, phiX)
    # Sequence names must be unique among all FASTAs!
    Array[File] other_fasta_gz = []

    String docker_image_id
  }

  call ensure_gz as genome_fasta {
    # accommodate uncompressed genome_fasta_gz; this makes it more convenient to use some of our
    # existing host genome FASTAs which we archived without compression.
    input:
    maybe_gz = genome_fasta_gz,
    docker_image_id
  }

  call concatenate_and_unzip_fastas {
    input:
    fasta_gz = flatten([[genome_fasta.gz, ERCC_fasta_gz], other_fasta_gz]),
    docker_image_id,
  }

  call bowtie2_build {
    input:
    fasta = concatenate_and_unzip_fastas.fasta,
    genome_name,
    docker_image_id,
  }

  call hisat2_build {
    input:
    fasta = concatenate_and_unzip_fastas.fasta,
    transcripts_gtf_gz,
    genome_name,
    docker_image_id,
  }

  call kallisto_index {
    input:
    transcripts_fasta_gz = flatten([transcripts_fasta_gz, [ERCC_fasta_gz]]),
    genome_name,
    docker_image_id,
  }

  call minimap2_index as minimap2_index_dna {
    input:
    fasta = concatenate_and_unzip_fastas.fasta,
    nucleotide_type = "dna",
    genome_name,
    docker_image_id,
  }

  call minimap2_index as minimap2_index_rna {
    input:
    fasta = concatenate_and_unzip_fastas.fasta,
    nucleotide_type = "rna",
    genome_name,
    docker_image_id,
  }

  call star_generate {
    input:
    fasta = concatenate_and_unzip_fastas.fasta,
    ERCC_fasta_gtf, 
    transcripts_gtf_gz,
    genome_name,
    docker_image_id,
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
    String docker_image_id
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
      docker: docker_image_id
      cpu: 4
      memory: "4GiB"
  }
}

task concatenate_and_unzip_fastas {
  input {
    Array[File] fasta_gz
    String docker_image_id
  }

  command <<<
    pigz -dc ~{sep(' ',fasta_gz)} > "all.fasta"
  >>>

  output {
      File fasta = "all.fasta"
  }

  runtime {
      docker: docker_image_id
      cpu: 4
      memory: "4GiB"
  }
}

task bowtie2_build {
  input {
    File fasta
    String genome_name
    Int seed = 42

    Int cpu = 16
    String docker_image_id
  }

  command <<<
    set -euxo pipefail
    TMPDIR=${TMPDIR:-/tmp}

    mkdir -p "$TMPDIR"'/bt2/~{genome_name}'
    >&2 bowtie2-build --seed ~{seed} --threads ~{cpu} "~{fasta}" "$TMPDIR"'/bt2/~{genome_name}/~{genome_name}'
    >&2 ls -lR "$TMPDIR/bt2"
    ln -s "$TMPDIR"'/bt2/~{genome_name}' "$TMPDIR"'/bt2/~{genome_name}.bowtie2'
    env -C "$TMPDIR/bt2" tar c . > '~{genome_name}.bowtie2.tar'
  >>>

  output {
      File index_tar = "~{genome_name}.bowtie2.tar"
  }

  runtime {
      docker: docker_image_id
      cpu: cpu
      memory: "~{cpu*2}GiB"
  }
}

task hisat2_build {
  input {
    File fasta
    File? transcripts_gtf_gz
    String genome_name

    Int cpu = 32
    String docker_image_id
  }

  command <<<
    set -euxo pipefail
    TMPDIR=${TMPDIR:-/tmp}

    mkdir -p "$TMPDIR"'/hisat2/~{genome_name}'
    if [[ -n '~{transcripts_gtf_gz}' ]]; then
      # convert GTF per http://daehwankimlab.github.io/hisat2/howto/
      /hisat2/hisat2_extract_splice_sites.py <(pigz -dc '~{transcripts_gtf_gz}') > "$TMPDIR/genome.ss" & pid=$!
      /hisat2/hisat2_extract_exons.py <(pigz -dc '~{transcripts_gtf_gz}') > "$TMPDIR/genome.exon"
      wait $pid
      >&2 /hisat2/hisat2-build -p 16 \
        --exon "$TMPDIR/genome.exon" --ss "$TMPDIR/genome.ss" \
        "~{fasta}" "$TMPDIR"'/hisat2/~{genome_name}/~{genome_name}'
    else
      >&2 /hisat2/hisat2-build -p 16 "~{fasta}" "$TMPDIR"'/hisat2/~{genome_name}/~{genome_name}'
    fi
    >&2 ls -lR "$TMPDIR/hisat2"
    env -C "$TMPDIR/hisat2" tar c . > '~{genome_name}.hisat2.tar'
  >>>

  output {
    File index_tar = "~{genome_name}.hisat2.tar"
  }

  runtime {
    docker: docker_image_id
    cpu: cpu
    memory: "240G"
  }
}

task kallisto_index {
  input {
    Array[File] transcripts_fasta_gz
    String genome_name

    String docker_image_id
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
    docker: docker_image_id
    memory: "16GiB"
  }
}

task minimap2_index {
  input {
    File fasta
    String genome_name
    String nucleotide_type

    String docker_image_id
  }

  command <<<
    set -euxo pipefail
    TMPDIR=${TMPDIR:-/tmp}

    if [ "~{nucleotide_type}" == "dna" ]; then
        >&2 minimap2 -x map-ont -d '~{genome_name}_~{nucleotide_type}.mmi' "~{fasta}"
    else
        >&2 minimap2 -x splice -d '~{genome_name}_~{nucleotide_type}.mmi' "~{fasta}"
    fi
    >&2 ls -l
  >>>

  output {
      File index_mmi = "~{genome_name}_~{nucleotide_type}.mmi"
  }

  runtime {
      docker: docker_image_id
      memory: "32GiB"
  }
}

task star_generate {
  input {
    File fasta
    File? ERCC_fasta_gtf
    File? transcripts_gtf_gz
    String genome_name


    Int cpu = 32
    String docker_image_id
  }

  command <<<
    set -euxo pipefail
    TMPDIR=${TMPDIR:-/tmp}

    gtf_flag=""
    if [[ -n '~{transcripts_gtf_gz}' || -n '~{ERCC_fasta_gtf}' ]]; then
      transcripts_gtf="$TMPDIR/transcripts.gtf"
      gtf_flag="--sjdbGTFfile \"$transcripts_gtf\""
      if [[ -n '~{transcripts_gtf_gz}' ]]; then
        pigz -dc '~{transcripts_gtf_gz}' > "$transcripts_gtf"
      fi 
      if [[ -n '~{ERCC_fasta_gtf}' ]]; then
        cat '~{ERCC_fasta_gtf}' >> "$transcripts_gtf"
      fi 
    fi

    # Make directory for STAR genome
    STAR_GENOME="~{genome_name}_STAR_genome"
    # HACK: we used to support splitting star indexes into many parts, this made things slower
    # Here we generate the index as if it is in many parts, but there is only ever one part for
    # backwards compatibility
    mkdir -p "$STAR_GENOME/part-0"

    STAR \
      --runThreadN ~{cpu} \
      --runMode genomeGenerate \
      --genomeFastaFiles "~{fasta}" \
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
      docker: docker_image_id
      cpu: cpu
      memory: "64GiB"
  }
}
