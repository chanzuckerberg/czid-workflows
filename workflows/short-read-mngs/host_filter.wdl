version 1.0

# CZ ID short-read-mngs pipeline stage 1 (2022 version):
# - input validation & QC
# - host & human filtering
# - deduplication
# - subsampling
workflow czid_host_filter {
  input {
    File fastqs_0
    File? fastqs_1
    String nucleotide_type = "DNA"

    File adapter_fasta

    String host_genome
    File bowtie2_index_tar
    File hisat2_index_tar
    File kallisto_idx
    File? gtf_gz  # Ensembl GTF for host species

    File human_bowtie2_index_tar
    File human_hisat2_index_tar

    Int max_input_fragments
    Int max_subsample_fragments

    Int cpu = 16
    String docker_image_id

    # legacy idseq-dag inputs:
    String file_ext = "fastq"
    String s3_wd_uri
  }

  # Validate input reads (and truncate if very large)
  call RunValidateInput {
    input:
    reads1_fastq = fastqs_0,
    reads2_fastq = fastqs_1,
    file_ext = file_ext,
    max_input_fragments = max_input_fragments,
    docker_image_id = docker_image_id,
    s3_wd_uri = s3_wd_uri
  }

  # Adapter trimming and QC filtering
  call fastp_qc {
    input:
    valid_input1_fastq = RunValidateInput.valid_input1_fastq,
    valid_input2_fastq = RunValidateInput.valid_input2_fastq,
    adapter_fasta = adapter_fasta,
    docker_image_id = docker_image_id,
    cpu = cpu
  }

  # Quantify host transcripts and ERCC
  # NOTE: we run kallisto even if nucleotide_type == "DNA" in order to get ERCC read counts.
  # The transcript & gene abundances are ~meaningless in that case, of course. This isn't a big
  # wasted cost because kallisto is so fast.
  call kallisto {
    input:
    fastp1_fastq = fastp_qc.fastp1_fastq,
    fastp2_fastq = fastp_qc.fastp2_fastq,
    kallisto_idx = kallisto_idx,
    gtf_gz = gtf_gz,
    docker_image_id = docker_image_id,
    cpu = cpu
  }

  # Filter out host reads.
  # Two stages: bowtie2 --very-sensitive-local, followed by splice-aware HISAT2.
  call bowtie2_filter {
    input:
    fastp1_fastq = fastp_qc.fastp1_fastq,
    fastp2_fastq = fastp_qc.fastp2_fastq,
    index_tar = bowtie2_index_tar,
    docker_image_id = docker_image_id,
    cpu = cpu
  }

  call hisat2_filter {
    input:
    bowtie2_host_filtered1_fastq = bowtie2_filter.bowtie2_host_filtered1_fastq,
    bowtie2_host_filtered2_fastq = bowtie2_filter.bowtie2_host_filtered2_fastq,
    index_tar = hisat2_index_tar,
    docker_image_id = docker_image_id,
    cpu = cpu
  }

  # If paired-end, collect insert size metrics from unfiltered, host-aligned bowtie2 BAM.
  if (defined(fastqs_1)) {
    call collect_insert_size_metrics {
      input:
      bam = bowtie2_filter.bam,
      docker_image_id = docker_image_id
    }
  }

  # Additionally filter out human reads, if the host is non-human.
  if (host_genome != "human") {
    call bowtie2_human_filter {
      input:
      hisat2_host_filtered1_fastq = hisat2_filter.hisat2_host_filtered1_fastq,
      hisat2_host_filtered2_fastq = hisat2_filter.hisat2_host_filtered2_fastq,
      index_tar = human_bowtie2_index_tar,
      docker_image_id = docker_image_id,
      cpu = cpu
    }
    call hisat2_human_filter {
      input:
      bowtie2_human_filtered1_fastq = bowtie2_human_filter.bowtie2_human_filtered1_fastq,
      bowtie2_human_filtered2_fastq = bowtie2_human_filter.bowtie2_human_filtered2_fastq,
      index_tar = human_hisat2_index_tar,
      docker_image_id = docker_image_id,
      cpu = cpu
    }
  }

  # Collect effectively filtered reads from the previous conditional
  File hisat2_filtered1_fastq = select_first([hisat2_human_filter.hisat2_human_filtered1_fastq, hisat2_filter.hisat2_host_filtered1_fastq])
  File? hisat2_filtered2_fastq = if defined(hisat2_human_filter.hisat2_human_filtered2_fastq) then hisat2_human_filter.hisat2_human_filtered2_fastq
                                                                          else hisat2_filter.hisat2_host_filtered2_fastq

  # Deduplicate filtered reads using custom czid-dedup tool.
  # It retains one exemplar [pair] from each duplicate cluster, and produces mapping from exemplar
  # read name to cluster size.
  call RunCZIDDedup {
    input:
    hisat2_filtered1_fastq = hisat2_filtered1_fastq,
    hisat2_filtered2_fastq = hisat2_filtered2_fastq,
    docker_image_id = docker_image_id,
    s3_wd_uri = s3_wd_uri
  }

  # Subsample remaining reads.
  call RunSubsample {
    input:
    dedup1_fastq = RunCZIDDedup.dedup1_fastq,
    dedup2_fastq = RunCZIDDedup.dedup2_fastq,
    duplicate_cluster_sizes_tsv = RunCZIDDedup.duplicate_cluster_sizes_tsv,
    max_subsample_fragments = max_subsample_fragments,
    docker_image_id = docker_image_id,
    s3_wd_uri = s3_wd_uri
  }

  output {
    File input_read_count = RunValidateInput.reads_in_count
    File validate_input_out_valid_input1_fastq = RunValidateInput.valid_input1_fastq
    File? validate_input_out_valid_input2_fastq = RunValidateInput.valid_input2_fastq
    File validate_input_out_validate_input_summary_json = RunValidateInput.validate_input_summary_json
    File validate_input_out_count = RunValidateInput.reads_out_count

    File fastp_out_fastp1_fastq = fastp_qc.fastp1_fastq
    File? fastp_out_fastp2_fastq = fastp_qc.fastp2_fastq
    File fastp_out_count = fastp_qc.reads_out_count
    File fastp_html = fastp_qc.fastp_html
    File fastp_json = fastp_qc.fastp_json

    File kallisto_transcript_abundance_tsv = kallisto.transcript_abundance_tsv
    File kallisto_ERCC_counts_tsv = kallisto.ERCC_counts_tsv
    File? kallisto_gene_abundance_tsv = kallisto.gene_abundance_tsv

    File bowtie2_host_filtered1_fastq = bowtie2_filter.bowtie2_host_filtered1_fastq
    File? bowtie2_host_filtered2_fastq = bowtie2_filter.bowtie2_host_filtered2_fastq
    File bowtie2_host_filtered_out_count = bowtie2_filter.reads_out_count
    File bowtie2_host_filtered_bam = bowtie2_filter.bam
    File hisat2_host_filtered1_fastq = hisat2_filter.hisat2_host_filtered1_fastq
    File? hisat2_host_filtered2_fastq = hisat2_filter.hisat2_host_filtered2_fastq
    File hisat2_host_filtered_out_count = hisat2_filter.reads_out_count

    File? insert_size_metrics = collect_insert_size_metrics.insert_size_metrics
    File? insert_size_histogram = collect_insert_size_metrics.insert_size_histogram

    File? bowtie2_human_filtered1_fastq = bowtie2_human_filter.bowtie2_human_filtered1_fastq
    File? bowtie2_human_filtered2_fastq = bowtie2_human_filter.bowtie2_human_filtered2_fastq
    File? bowtie2_human_filtered_out_count = bowtie2_human_filter.reads_out_count
    File? hisat2_human_filtered1_fastq = hisat2_human_filter.hisat2_human_filtered1_fastq
    File? hisat2_human_filtered2_fastq = hisat2_human_filter.hisat2_human_filtered2_fastq
    File? hisat2_human_filtered_out_count = hisat2_human_filter.reads_out_count

    File czid_dedup_out_dedup1_fastq = RunCZIDDedup.dedup1_fastq
    File? czid_dedup_out_dedup2_fastq = RunCZIDDedup.dedup2_fastq
    File czid_dedup_out_duplicate_clusters_csv = RunCZIDDedup.duplicate_clusters_csv
    File czid_dedup_out_duplicate_cluster_sizes_tsv = RunCZIDDedup.duplicate_cluster_sizes_tsv
    File czid_dedup_out_count = RunCZIDDedup.reads_out_count

    File subsampled_out_subsampled_1_fa = RunSubsample.subsampled_1_fa
    File? subsampled_out_subsampled_2_fa = RunSubsample.subsampled_2_fa
    File? subsampled_out_subsampled_merged_fa = RunSubsample.subsampled_merged_fa
    File subsampled_out_count = RunSubsample.reads_out_count
  }
}

task RunValidateInput {
  input {
    File reads1_fastq
    File? reads2_fastq
    String file_ext

    Int max_input_fragments

    String docker_image_id
    String s3_wd_uri
  }
  Boolean paired = defined(reads2_fastq)
  command<<<
    set -euxo pipefail
    idseq-dag-run-step --workflow-name host_filter \
      --step-module idseq_dag.steps.run_validate_input \
      --step-class PipelineStepRunValidateInput \
      --step-name validate_input_out \
      --input-files '[["~{sep='","' select_all([reads1_fastq, reads2_fastq])}"]]' \
      --output-files '["validate_input_summary.json", ~{if paired then '"valid_input1.fastq", "valid_input2.fastq"' else '"valid_input1.fastq"'}]' \
      --output-dir-s3 '~{s3_wd_uri}' \
      --additional-files '{}' \
      --additional-attributes '{"truncate_fragments_to": ~{max_input_fragments}, "file_ext": "~{file_ext}"}'
  >>>
  output {
    String step_description_md = read_string("validate_input_out.description.md")
    File validate_input_summary_json = "validate_input_summary.json"
    File valid_input1_fastq = "valid_input1.fastq"
    File? valid_input2_fastq = "valid_input2.fastq"
    File reads_out_count = "validate_input_out.count"
    File reads_in_count = "fastqs.count"
  }
  runtime {
    docker: docker_image_id
    cpu: 4
    memory: "8G"
  }
}

task fastp_qc {
    # fastp all-in-one for
    # - adapter trimming
    # - quality filtering
    # - complexity filtering
  input {
    File valid_input1_fastq
    File? valid_input2_fastq
    File adapter_fasta

    # These default QC thresholds are loosely based on the pre-2022 pipeline using PriceSeq & LZW
    String fastp_options = "--dont_eval_duplication --length_required 35" +
                           " --qualified_quality_phred 17 --unqualified_percent_limit 15 --n_base_limit 15" +
                           " --sdust_complexity_filter --complexity_threshold 60"

    String docker_image_id
    Int cpu = 16
  }
  Boolean paired = defined(valid_input2_fastq)
  String fastp_invocation = "fastp"
        + " -i ${valid_input1_fastq} ${'-I ' + valid_input2_fastq}"
        + " -o fastp1.fastq ${if (paired) then '-O fastp2.fastq' else ''}"
        + " -w ${cpu} ${fastp_options}"
        + " --adapter_fasta ${adapter_fasta} ${if (paired) then '--detect_adapter_for_pe' else ''}"

  command<<<
    set -euxo pipefail
    ~{fastp_invocation}
    count="$(jq .read1_after_filtering.total_reads fastp.json)"
    if [ '~{paired}' == 'true' ]; then
        count=$((2 * count))
    fi
    jq --null-input --arg count "$count" '{"fastp_out":$count}' > fastp_out.count
    # TODO: extract insert size metrics from JSON, also render histogram?

    python3 - << 'EOF'
    import textwrap
    with open("fastp.description.md", "w") as outfile:
      print(textwrap.dedent("""
      **fastp read trimming & filtering**

      Processes the reads using [fastp](https://github.com/OpenGene/fastp):

      1. Trim adapters
      2. Quality score filter
      3. Non-called base (N) filter
      4. Length filter
      5. Complexity filter ([custom feature](https://github.com/mlin/fastp/tree/mlin/sdust)
         using the [SDUST algorithm](https://pubmed.ncbi.nlm.nih.gov/16796549/))

      fastp is run on the FASTQ file(s) from input validation:
      ```
      ~{fastp_invocation}
      ```

      fastp documentation can be found [here](https://github.com/OpenGene/fastp)
      """).strip(), file=outfile)
    EOF
  >>>
  output {
    String step_description_md = read_string("fastp.description.md")
    File fastp1_fastq = "fastp1.fastq"
    File? fastp2_fastq = "fastp2.fastq"
    File fastp_html = "fastp.html"
    File fastp_json = "fastp.json"
    File reads_out_count = "fastp_out.count"
  }
  runtime {
    docker: docker_image_id
    cpu: cpu
    memory: "~{cpu}G"
  }
}

task kallisto {
  input {
    File fastp1_fastq
    File? fastp2_fastq
    File kallisto_idx
    File? gtf_gz
    String kallisto_options = ""

    String docker_image_id
    Int cpu = 16
  }
  Boolean paired = defined(fastp2_fastq)
  # TODO: input fragment length parameters for non-paired-end (l = average, s = std dev)
  String kallisto_invocation = "/kallisto/kallisto quant"
      + " -i '${kallisto_idx}' -o $(pwd) --plaintext ${if (paired) then '' else '--single -l 200 -s 20'} ${kallisto_options} -t ${cpu}"
      + " '~{fastp1_fastq}'" + if (defined(fastp2_fastq)) then " '~{fastp2_fastq}'" else ""

  command <<<
    set -euxo pipefail

    # NOTE: kallisto exit code will be 1 if no reads pseudoalign, which we don't necessarily
    #       consider an error. Therefore decide success based on existence of run_info.json and
    #       abundance.tsv
    ~{kallisto_invocation} || true
    >&2 jq . run_info.json

    mv abundance.tsv reads_per_transcript.kallisto.tsv

    # extract ERCC counts
    echo -e "target_id\test_counts" > ERCC_counts.tsv
    grep ERCC- reads_per_transcript.kallisto.tsv | cut -f1,4 >> ERCC_counts.tsv

    # If we've been provided the GTF, then roll up the transcript abundance estimates by gene.
    if [[ -n '~{gtf_gz}' ]]; then
      python3 - reads_per_transcript.kallisto.tsv '~{gtf_gz}' << 'EOF'
    # Given kallisto output tsv based on index of Ensembl transcripts FASTA, and matching
    # Ensembl GTF, report the total est_counts and tpm for each gene (sum over all transcripts
    # of each gene).
    import sys
    import pandas as pd
    import gtfparse

    kallisto_df = pd.read_csv(sys.argv[1], sep="\t")

    gtf_df = gtfparse.read_gtf(sys.argv[2])
    tx_df = gtf_df[gtf_df["feature"] == "transcript"][
        ["transcript_id", "transcript_version", "gene_id"]
    ]
    # kallisto target_id is a versioned transcript ID e.g. "ENST00000390446.3", while the GTF
    # breaks out: transcript_id "ENST00000390446"; transcript_version "3";
    # synthesize a column with the versioned transcript ID for merging.
    tx_df = tx_df.assign(
        transcript_id_version=tx_df["transcript_id"] + "." + tx_df["transcript_version"]
    )

    merged_df = pd.merge(
        kallisto_df[["target_id", "est_counts", "tpm"]],
        tx_df[["transcript_id_version", "gene_id"]],
        left_on="target_id",
        right_on="transcript_id_version",
    )

    gene_abundance = merged_df.groupby("gene_id").sum(numeric_only=True)
    gene_abundance.to_csv("reads_per_gene.kallisto.tsv", sep="\t")
    EOF
    fi

    python3 - << 'EOF'
    import textwrap
    with open("kallisto.description.md", "w") as outfile:
      print(textwrap.dedent("""
      **kallisto RNA quantification**

      Quantifies host transcripts using [kallisto](https://pachterlab.github.io/kallisto/about).
      The host transcript sequences are sourced from Ensembl, along with
      [ERCC control sequences](https://www.nist.gov/programs-projects/external-rna-controls-consortium).
      Not all CZ ID host species have transcripts indexed; for those without, kallisto is run using ERCC
      sequences only.

      kallisto is run on the fastp-filtered FASTQ(s):

      ```
      ~{kallisto_invocation}
      ```

      kallisto documentation can be found [here](https://pachterlab.github.io/kallisto/manual), including
      details of the `transcript_abundance.tsv` output format.
      """).strip(), file=outfile)
    EOF
  >>>

  output {
    String step_description_md = read_string("kallisto.description.md")
    File transcript_abundance_tsv = "reads_per_transcript.kallisto.tsv"
    File ERCC_counts_tsv = "ERCC_counts.tsv"
    File? gene_abundance_tsv = "reads_per_gene.kallisto.tsv"
  }

  runtime {
    docker: docker_image_id
    cpu: cpu
    memory: "~{cpu}G"
  }
}

task bowtie2_filter {
  # Remove reads [pairs] with bowtie2 hits to the given index
  input {
    File fastp1_fastq
    File? fastp2_fastq

    # GENOME_NAME.bowtie2.tar should contain GENOME_NAME/GENOME_NAME.*.bt*
    File index_tar
    String bowtie2_options = "--very-sensitive-local"

    String docker_image_id
    Int cpu = 16
  }

  Boolean paired = defined(fastp2_fastq)
  String genome_name = basename(index_tar, ".bowtie2.tar")
  String bowtie2_invocation =
      "bowtie2 -x '/tmp/${genome_name}/${genome_name}' ${bowtie2_options} -p ${cpu}"
        + (if (paired) then " -1 '${fastp1_fastq}' -2 '${fastp2_fastq}'" else " -U '${fastp1_fastq}'")
        + " -q -S '/tmp/bowtie2.sam'"

  command <<<
    set -euxo pipefail

    tar xf '~{index_tar}' -C /tmp

    ~{bowtie2_invocation}

    # generate sort & compressed BAM file for archival
    samtools sort -n -o "bowtie2_host.bam" -@ 4 -T /tmp "/tmp/bowtie2.sam" & samtools_pid=$!

    # Extract reads [pairs] that did NOT map to the index
    if [[ '~{paired}' == 'true' ]]; then
        #    1 (read paired)
        #    4 (read unmapped)
        # +  8 (mate unmapped)
        # ----
        #   13
        samtools fastq -f 13 -1 'bowtie2_host_filtered1.fastq' -2 'bowtie2_host_filtered2.fastq' -0 /dev/null -s /dev/null /tmp/bowtie2.sam
    else
        samtools fastq -f 4 /tmp/bowtie2.sam > 'bowtie2_host_filtered1.fastq'
    fi

    count="$(cat bowtie2_host_filtered{1,2}.fastq | wc -l)"
    count=$((count / 4))
    jq --null-input --arg count "$count" '{"bowtie2_host_filtered_out":$count}' > 'bowtie2_host_filtered_out.count'

    python3 - << 'EOF'
    import textwrap
    with open("bowtie2.description.md", "w") as outfile:
      print(textwrap.dedent("""
      **bowtie2 host filtering**

      Filters out reads matching the host genome using
      [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml). Runs
      `bowtie2 ~{bowtie2_options}` using a precomputed index, then uses
      [samtools](http://www.htslib.org/) to keep reads *not* mapping to the host genome.

      Bowtie2 is run on the fastp-filtered FASTQ(s):

      ```
      ~{bowtie2_invocation}
      ```

      Then, non-mapping reads are selected using `samtools fastq -f ~{if (paired) then 13 else 4}`.

      Bowtie2 documentation can be found [here](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
      """).strip(), file=outfile)
    EOF

    wait $samtools_pid
  >>>

  output {
    String step_description_md = read_string("bowtie2.description.md")
    File bowtie2_host_filtered1_fastq = "bowtie2_host_filtered1.fastq"
    File? bowtie2_host_filtered2_fastq = "bowtie2_host_filtered2.fastq"
    File reads_out_count = "bowtie2_host_filtered_out.count"
    File bam = "bowtie2_host.bam"
  }
  runtime {
    docker: docker_image_id
    cpu: cpu
    memory: "~{cpu*2}G"
  }
}

task hisat2_filter {
  # Remove reads [pairs] with HISAT2 hits to the given index
  input {
    File bowtie2_host_filtered1_fastq
    File? bowtie2_host_filtered2_fastq

    # GENOME_NAME.hisat2.tar should contain GENOME_NAME/GENOME_NAME.*.ht2
    File index_tar
    String hisat2_options = ""

    String docker_image_id
    Int cpu = 16
  }

  Boolean paired = defined(bowtie2_host_filtered2_fastq)
  String genome_name = basename(index_tar, ".hisat2.tar")
  String hisat2_invocation =
      "/hisat2/hisat2 -x '/tmp/${genome_name}/${genome_name}' ${hisat2_options} -p ${cpu}"
        + (if (paired) then " -1 '${bowtie2_host_filtered1_fastq}' -2 '${bowtie2_host_filtered2_fastq}'" else " -U '${bowtie2_host_filtered2_fastq}'")
        + " -q -S /tmp/hisat2.sam"

  command <<<
    set -euxo pipefail

    tar xf '~{index_tar}' -C /tmp

    ~{hisat2_invocation}

    # Extract reads [pairs] that did NOT map to the index
    if [[ '~{paired}' == 'true' ]]; then
        #    1 (read paired)
        #    4 (read unmapped)
        # +  8 (mate unmapped)
        # ----
        #   13
        samtools fastq -f 13 -1 'hisat2_host_filtered1.fastq' -2 'hisat2_host_filtered2.fastq' -0 /dev/null -s /dev/null /tmp/hisat2.sam
    else
        samtools fastq -f 4 /tmp/hisat2.sam > 'hisat2_host_filtered1.fastq'
    fi

    count="$(cat hisat2_host_filtered{1,2}.fastq | wc -l)"
    count=$((count / 4))
    jq --null-input --arg count "$count" '{"hisat2_host_filtered_out":$count}' > 'hisat2_host_filtered_out.count'

    python3 - << 'EOF'
    import textwrap
    with open("hisat2.description.md", "w") as outfile:
      print(textwrap.dedent("""
      **HISAT2 host filtering**

      Filters out reads matching the host genome using
      [HISAT2](http://daehwankimlab.github.io/hisat2/). Runs `hisat2` using a precomputed index,
      then uses [samtools](http://www.htslib.org/) to keep reads *not* mapping to the
      host genome.

      HISAT2 complements Bowtie2 with a different algorithm that also models potential RNA splice
      junctions (if CZ ID indexes transcript models for the host).

      HISAT2 is run on the bowtie2-filtered FASTQ(s):

      ```
      ~{hisat2_invocation}
      ```

      Then, non-mapping reads are selected using `samtools fastq -f ~{if (paired) then 13 else 4}`.

      HISAT2 documentation can be found [here](http://daehwankimlab.github.io/hisat2/)
      """).strip(), file=outfile)
    EOF
  >>>

  output {
    String step_description_md = read_string("hisat2.description.md")
    File hisat2_host_filtered1_fastq = "hisat2_host_filtered1.fastq"
    File? hisat2_host_filtered2_fastq = "hisat2_host_filtered2.fastq"
    File reads_out_count = "hisat2_host_filtered_out.count"
  }
  runtime {
    docker: docker_image_id
    cpu: cpu
    memory: "~{cpu*2}G"
  }
}

###################################################################################################
### NOTE: bowtie2_human_filter and hisat2_human_filter are roughly copy/paste of the _host_filter
###       tasks above. We'd much prefer to consolidate them, but the webapp pipeline visualization
###       isn't yet able to handle WDL tasks used multiple times with dynamic output filenames.
###################################################################################################

task bowtie2_human_filter {
  # Remove reads [pairs] with bowtie2 hits to the given index
  input {
    File hisat2_host_filtered1_fastq
    File? hisat2_host_filtered2_fastq

    # GENOME_NAME.bowtie2.tar should contain GENOME_NAME/GENOME_NAME.*.bt*
    File index_tar
    String bowtie2_options = "--very-sensitive-local"

    String docker_image_id
    Int cpu = 16
  }

  Boolean paired = defined(hisat2_host_filtered2_fastq)
  String genome_name = basename(index_tar, ".bowtie2.tar")
  String bowtie2_invocation =
      "bowtie2 -x '/tmp/${genome_name}/${genome_name}' ${bowtie2_options} -p ${cpu}"
        + (if (paired) then " -1 '${hisat2_host_filtered1_fastq}' -2 '${hisat2_host_filtered2_fastq}'" else " -U '${hisat2_host_filtered1_fastq}'")
        + " -q -S '/tmp/bowtie2.sam'"

  command <<<
    set -euxo pipefail

    tar xf '~{index_tar}' -C /tmp

    ~{bowtie2_invocation}

    # generate sort & compressed BAM file for archival
    samtools sort -n -o "bowtie2_human.bam" -@ 4 -T /tmp "/tmp/bowtie2.sam" & samtools_pid=$!

    # Extract reads [pairs] that did NOT map to the index
    if [[ '~{paired}' == 'true' ]]; then
        #    1 (read paired)
        #    4 (read unmapped)
        # +  8 (mate unmapped)
        # ----
        #   13
        samtools fastq -f 13 -1 'bowtie2_human_filtered1.fastq' -2 'bowtie2_human_filtered2.fastq' -0 /dev/null -s /dev/null /tmp/bowtie2.sam
    else
        samtools fastq -f 4 /tmp/bowtie2.sam > 'bowtie2_human_filtered1.fastq'
    fi

    count="$(cat bowtie2_human_filtered{1,2}.fastq | wc -l)"
    count=$((count / 4))
    jq --null-input --arg count "$count" '{"bowtie2_human_filtered_out":$count}' > 'bowtie2_human_filtered_out.count'

    python3 - << 'EOF'
    import textwrap
    with open("bowtie2.description.md", "w") as outfile:
      print(textwrap.dedent("""
      **bowtie2 human filtering**

      Filters out reads matching the human genome using
      [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml). This is similar to the
      host filtering task, but CZ ID also filters non-human samples against human genome indexes
      to alleviate any potential data privacy concerns.
      """).strip(), file=outfile)
    EOF

    wait $samtools_pid
  >>>

  output {
    String step_description_md = read_string("bowtie2.description.md")
    File bowtie2_human_filtered1_fastq = "bowtie2_human_filtered1.fastq"
    File? bowtie2_human_filtered2_fastq = "bowtie2_human_filtered2.fastq"
    File reads_out_count = "bowtie2_human_filtered_out.count"
    File bam = "bowtie2_human.bam"
  }
  runtime {
    docker: docker_image_id
    cpu: cpu
    memory: "~{cpu*2}G"
  }
}

task hisat2_human_filter {
  # Remove reads [pairs] with HISAT2 hits to the given index
  input {
    File bowtie2_human_filtered1_fastq
    File? bowtie2_human_filtered2_fastq

    # GENOME_NAME.hisat2.tar should contain GENOME_NAME/GENOME_NAME.*.ht2
    File index_tar
    String hisat2_options = ""

    String docker_image_id
    Int cpu = 16
  }

  Boolean paired = defined(bowtie2_human_filtered2_fastq)
  String genome_name = basename(index_tar, ".hisat2.tar")
  String hisat2_invocation =
      "/hisat2/hisat2 -x '/tmp/${genome_name}/${genome_name}' ${hisat2_options} -p ${cpu}"
        + (if (paired) then " -1 '${bowtie2_human_filtered1_fastq}' -2 '${bowtie2_human_filtered2_fastq}'" else " -U '${bowtie2_human_filtered2_fastq}'")
        + " -q -S /tmp/hisat2.sam"

  command <<<
    set -euxo pipefail

    tar xf '~{index_tar}' -C /tmp

    ~{hisat2_invocation}

    # Extract reads [pairs] that did NOT map to the index
    if [[ '~{paired}' == 'true' ]]; then
        #    1 (read paired)
        #    4 (read unmapped)
        # +  8 (mate unmapped)
        # ----
        #   13
        samtools fastq -f 13 -1 'hisat2_human_filtered1.fastq' -2 'hisat2_human_filtered2.fastq' -0 /dev/null -s /dev/null /tmp/hisat2.sam
    else
        samtools fastq -f 4 /tmp/hisat2.sam > 'hisat2_human_filtered1.fastq'
    fi

    count="$(cat hisat2_human_filtered{1,2}.fastq | wc -l)"
    count=$((count / 4))
    jq --null-input --arg count "$count" '{"hisat2_human_filtered_out":$count}' > 'hisat2_human_filtered_out.count'

    python3 - << 'EOF'
    import textwrap
    with open("hisat2.description.md", "w") as outfile:
      print(textwrap.dedent("""
      **HISAT2 human filtering**

      Filters out reads matching the human genome using
      [HISAT2](http://daehwankimlab.github.io/hisat2/). This is similar to the host filtering
      task, but CZ ID also filters non-human samples against human genome indexes to alleviate any
      potential data privacy concerns.
      """).strip(), file=outfile)
    EOF
  >>>

  output {
    String step_description_md = read_string("hisat2.description.md")
    File hisat2_human_filtered1_fastq = "hisat2_human_filtered1.fastq"
    File? hisat2_human_filtered2_fastq = "hisat2_human_filtered2.fastq"
    File reads_out_count = "hisat2_human_filtered_out.count"
  }
  runtime {
    docker: docker_image_id
    cpu: cpu
    memory: "~{cpu*2}G"
  }
}

task collect_insert_size_metrics {
  input {
    File bam
    String docker_image_id
  }

  command <<<
    picard CollectInsertSizeMetrics 'I=~{bam}' O=picard_insert_metrics.txt H=insert_size_histogram.pdf
    python3 - << 'EOF'
    import textwrap
    with open("collect_insert_size_metrics.description.md", "w") as outfile:
      print(textwrap.dedent("""
      **Picard CollectInsertSizeMetrics**

      This step computes insert size metrics for Paired End samples. These metrics are computed by
      the Broad Institute's Picard toolkit.

      Picard is run on the output BAM file obtained from running Bowtie2 on the host genome:

      ```
      picard CollectInsertSizeMetrics 'I=~{bam}' O=picard_insert_metrics.txt H=insert_size_histogram.pdf
      ```

      Picard documentation can be found [here](https://gatk.broadinstitute.org/hc/en-us/articles/360037055772-CollectInsertSizeMetrics-Picard-)
      """).strip(), file=outfile)
    EOF
  >>>

  output {
    String step_description_md = read_string("collect_insert_size_metrics.description.md")
    # If no reads mapped to the host, then picard exits "successfully" without creating these files.
    File? insert_size_metrics = "picard_insert_metrics.txt"
    File? insert_size_histogram = "insert_size_histogram.pdf"
  }

  runtime {
    docker: docker_image_id
    cpu: 1
    memory: "8G"
  }
}

task RunCZIDDedup {
  input {
    File hisat2_filtered1_fastq
    File? hisat2_filtered2_fastq
    String docker_image_id
    String s3_wd_uri
  }
  Boolean paired = defined(hisat2_filtered2_fastq)
  command<<<
    set -euxo pipefail

    >&2 idseq-dag-run-step --workflow-name host_filter \
      --step-module idseq_dag.steps.run_czid_dedup \
      --step-class PipelineStepRunCZIDDedup \
      --step-name czid_dedup_out \
      --input-files '[["~{sep='","' select_all([hisat2_filtered1_fastq, hisat2_filtered2_fastq])}"]]' \
      --output-files '[~{if paired then '"dedup1.fastq","dedup2.fastq"' else '"dedup1.fastq"'}, "clusters.csv", "duplicate_cluster_sizes.tsv"]' \
      --output-dir-s3 '~{s3_wd_uri}' \
      --additional-files '{}' \
      --additional-attributes '{}'
  >>>
  output {
    String step_description_md = read_string("czid_dedup_out.description.md")
    File dedup1_fastq = "dedup1.fastq"
    File? dedup2_fastq = "dedup2.fastq"
    File duplicate_clusters_csv = "clusters.csv"
    File duplicate_cluster_sizes_tsv = "duplicate_cluster_sizes.tsv"
    File reads_out_count = "czid_dedup_out.count"
  }
  runtime {
    docker: docker_image_id
    cpu: 4
    memory: "16G"
  }
}

task RunSubsample {
  input {
    File dedup1_fastq
    File? dedup2_fastq
    File duplicate_cluster_sizes_tsv
    Int max_subsample_fragments

    String docker_image_id
    String s3_wd_uri
  }
  Boolean paired = defined(dedup2_fastq)
  command<<<
  set -euxo pipefail
  TMPDIR="${TMPDIR:-/tmp}"

  # Convert FASTQs to FASTAs: the idseq-dag subsampling tool inputs and outputs FASTAs, and
  # downstream pipeline stages consume the FASTAs.
  seqtk seq -a '~{dedup1_fastq}' > "$TMPDIR/reads1.fasta" & pid=$!
  fastas="\"$TMPDIR/reads1.fasta\""
  if [[ '~{paired}' == 'true' ]]; then
    seqtk seq -a '~{dedup2_fastq}' > "$TMPDIR/reads2.fasta"
    wait $pid
    # also generate merged FASTA. `seqtk mergepe` interleaves the reads but doesn't append /1 /2 to
    # the names, so we add an awk kludge to do that.
    seqtk mergepe "$TMPDIR/reads1.fasta" "$TMPDIR/reads2.fasta" | awk '
        BEGIN {
          name = "";
        }
        /^>.*/ {
          if ($0 != name) {
            name = $0;
            printf("%s/1\n", $0);
          } else {
            printf("%s/2\n", $0);
          }
        }
        ! /^>.*/ { print; }
      ' > "$TMPDIR/reads_merged.fasta"
    fastas="\"$TMPDIR/reads1.fasta\",\"$TMPDIR/reads2.fasta\",\"$TMPDIR/reads_merged.fasta\""
  else
    wait $pid
  fi

  # subsample FASTAs
  idseq-dag-run-step --workflow-name host_filter \
    --step-module idseq_dag.steps.run_subsample \
    --step-class PipelineStepRunSubsample \
    --step-name subsampled_out \
    --input-files '[['"$fastas"'], ["~{duplicate_cluster_sizes_tsv}"]]' \
    --output-files '[~{if paired then '"subsampled_1.fa", "subsampled_2.fa", "subsampled_merged.fa"' else '"subsampled_1.fa"'}]' \
    --output-dir-s3 '~{s3_wd_uri}' \
    --additional-files '{}' \
    --additional-attributes '{"max_fragments": ~{max_subsample_fragments}}'
  >>>
  output {
    String step_description_md = read_string("subsampled_out.description.md")
    File subsampled_1_fa = "subsampled_1.fa"
    File? subsampled_2_fa = "subsampled_2.fa"
    File? subsampled_merged_fa = "subsampled_merged.fa"
    File reads_out_count = "subsampled_out.count"
  }
  runtime {
    docker: docker_image_id
    cpu: 4
    memory: "8G"
  }
}
