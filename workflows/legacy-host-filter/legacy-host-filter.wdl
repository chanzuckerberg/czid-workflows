version 1.0

task RunValidateInput {
  input {
    String docker_image_id
    String s3_wd_uri
    Array[File] fastqs
    Int max_input_fragments
    String file_ext
  }
  command<<<
  set -euxo pipefail
  idseq-dag-run-step --workflow-name host_filter \
    --step-module idseq_dag.steps.run_validate_input \
    --step-class PipelineStepRunValidateInput \
    --step-name validate_input_out \
    --input-files '[["~{sep='","' fastqs}"]]' \
    --output-files '["validate_input_summary.json", ~{if length(fastqs) == 2 then '"valid_input1.fastq", "valid_input2.fastq"' else '"valid_input1.fastq"'}]' \
    --output-dir-s3 '~{s3_wd_uri}' \
    --additional-files '{}' \
    --additional-attributes '{"truncate_fragments_to": ~{max_input_fragments}, "file_ext": "~{file_ext}"}'
  >>>
  output {
    String step_description_md = read_string("validate_input_out.description.md")
    File validate_input_summary_json = "validate_input_summary.json"
    File valid_input1_fastq = "valid_input1.fastq"
    File? valid_input2_fastq = "valid_input2.fastq"
    File? output_read_count = "validate_input_out.count"
    File? input_read_count = "fastqs.count"
  }
  runtime {
    docker: docker_image_id
  }
}

task RunStar {
  input {
    String docker_image_id
    String s3_wd_uri
    File validate_input_summary_json
    Array[File] valid_input_fastq
    File star_genome
    String nucleotide_type
    String host_genome
    String genome_dir = "STAR_genome/part-0/"
  }
  command<<<
  # TODO(Ryan): remove when status upload is not dependent on idseq-dag see: https://app.shortcut.com/idseq/story/163323
  # this comment is for the miniwdl plugin uploader to parse:
  # --step-name star_out 
  set -euxo pipefail

  python3 <<CODE
  """ save description to file """
  from idseq_utils.save_descriptions import star_description
  from idseq_utils.exceptions import print_exceptions

  def main():
    description = star_description("~{nucleotide_type}")
    with open("star_out.description.md", "w+") as f:
      f.write(description)

  print_exceptions(main)
  CODE

  mkdir "STAR_genome"
  tar xf "~{star_genome}" -C "STAR_genome" --strip-components 1
  # Set Parameters
  SAMMODE="None"
  SAMTYPE="None"

  # Currently we always use 'GeneCounts', 
  QUANTMODE="~{if nucleotide_type == 'RNA' && host_genome == 'human' then 'TranscriptomeSAM GeneCounts' else 'GeneCounts'}"
  if [[ "~{length(valid_input_fastq)}" -eq "2" ]] && [[ "~{host_genome}" == "human" ]]; then
    SAMMODE="NoQS"
    SAMTYPE="BAM Unsorted"
  fi

  if [[ $(jq '."500-10000"' "~{validate_input_summary_json}") -gt "1" ]] || [[ $(jq '."10000+"' "~{validate_input_summary_json}") -gt "1" ]]; then 
    STARlong \
    --outFilterMultimapNmax 99999 \
    --outFilterScoreMinOverLread 0.5 \
    --outFilterMatchNminOverLread 0.5 \
    --outReadsUnmapped Fastx \
    --outFilterMismatchNmax 999 \
    --clip3pNbases 0 \
    --runThreadN "$(nproc --all)" \
    --genomeDir "~{genome_dir}" \
    --readFilesIn "~{sep='" "' valid_input_fastq}" \
    --seedSearchStartLmax 20 \
    --seedPerReadNmax 100000 \
    --seedPerWindowNmax 1000 \
    --alignTranscriptsPerReadNmax 100000 \
    --outSAMmode $SAMMODE \
    --outSAMtype $SAMTYPE \
    --quantMode $QUANTMODE
  else
    STAR --outFilterMultimapNmax 99999 \
    --outFilterScoreMinOverLread 0.5 \
    --outFilterMatchNminOverLread 0.5 \
    --outReadsUnmapped Fastx \
    --outFilterMismatchNmax 999 \
    --outSAMmode $SAMMODE \
    --outSAMtype $SAMTYPE \
    --clip3pNbases 0 \
    --runThreadN "$(nproc --all)" \
    --limitOutSJcollapsed 2000000 \
    --runRNGseed 777 \
    --genomeDir "~{genome_dir}" \
    --quantMode $QUANTMODE \
    --readFilesIn "~{sep='" "' valid_input_fastq}" 
  fi

  if [ -f "Aligned.toTranscriptome.out.bam" ]; then 
    mv "Aligned.toTranscriptome.out.bam" "Aligned.out.bam"
  fi

  python3 <<CODE
  """ sync pairs of files, sort by entry id, count reads """
  import idseq_utils.sync_pairs as sp
  import shutil
  import glob
  from idseq_utils.exceptions import print_exceptions, BrokenReadPairError
  
  def main():
    unmapped = sorted(glob.glob("Unmapped.out.mate*"))
    output_files, too_discrepant = sp.sync_pairs(unmapped)
    if too_discrepant:
        raise BrokenReadPairError("Paired input files were detected to be out of order. Both pairs must have the same read order")
    for unmapped_file in output_files:
        sp.sort_fastx_by_entry_id(unmapped_file)

    for ind, unmapped_file in enumerate(output_files):
        shutil.move(unmapped_file, f"unmapped{ind+1}.fastq")

  print_exceptions(main)
  CODE

  if [ -f "Aligned.out.bam" ]; then 
    picard CollectInsertSizeMetrics I=Aligned.out.bam O=picard_insert_metrics.txt H=insert_size_histogram.pdf || echo WARNING: picard failed
  fi 

  python3 <<CODE
  """ count reads """
  import idseq_utils.count_reads as cr
  import glob
  from idseq_utils.exceptions import print_exceptions

  def main():
    input_files = sorted(glob.glob("unmapped*.fastq"))
    cr.main("star_out", input_files)

  print_exceptions(main)
  CODE

  if [ -f "ReadsPerGene.out.tab" ]; then 
    mv ReadsPerGene.out.tab reads_per_gene.star.tab
  fi

  STAR --version > star_out_version.txt
  rm "~{genome_dir}"/SAindex # the star genome is pretty big (1.5G)
  rm "~{genome_dir}"/Genome 
  >>>
  output {
    String step_description_md = read_string("star_out.description.md")
    File unmapped1_fastq = "unmapped1.fastq"
    File output_log_file = "Log.final.out"
    File? unmapped2_fastq = "unmapped2.fastq"
    File? aligned_file = "Aligned.out.bam"
    File? output_read_count = "star_out.count"
    File? output_gene_file = "reads_per_gene.star.tab"
    File? output_metrics_file = "picard_insert_metrics.txt"
    File? output_histogram_file = "insert_size_histogram.pdf"
    File? version = "star_out_version.txt"
  }
  runtime {
    docker: docker_image_id
  }
}

task RunTrimmomatic {
  input {
    String docker_image_id
    String s3_wd_uri
    Array[File] unmapped_fastq
    File adapter_fasta
  }
  command<<<
  set -euxo pipefail
  idseq-dag-run-step --workflow-name host_filter \
    --step-module idseq_dag.steps.run_trimmomatic \
    --step-class PipelineStepRunTrimmomatic \
    --step-name trimmomatic_out \
    --input-files '[["~{sep='","' unmapped_fastq}"]]' \
    --output-files '[~{if length(unmapped_fastq) == 2 then '"trimmomatic1.fastq", "trimmomatic2.fastq"' else '"trimmomatic1.fastq"'}]' \
    --output-dir-s3 '~{s3_wd_uri}' \
    --additional-files '{"adapter_fasta": "~{adapter_fasta}"}' \
    --additional-attributes '{}'
  java -jar /usr/local/bin/trimmomatic-0.38.jar -version > trimmomatic_version.txt
  
  >>>
  output {
    String step_description_md = read_string("trimmomatic_out.description.md")
    File trimmomatic1_fastq = "trimmomatic1.fastq"
    File? trimmomatic2_fastq = "trimmomatic2.fastq"
    File? output_read_count = "trimmomatic_out.count"
    File? version = "trimmomatic_version.txt"
  }
  runtime {
    docker: docker_image_id
  }
}

task RunPriceSeq {
  input {
    String docker_image_id
    String s3_wd_uri
    Array[File] trimmomatic_fastq
  }
  command<<<
  set -euxo pipefail
  idseq-dag-run-step --workflow-name host_filter \
    --step-module idseq_dag.steps.run_priceseq \
    --step-class PipelineStepRunPriceSeq \
    --step-name priceseq_out \
    --input-files '[["~{sep='","' trimmomatic_fastq}"]]' \
    --output-files '[~{if length(trimmomatic_fastq) == 2 then '"priceseq1.fa", "priceseq2.fa"' else '"priceseq1.fa"'}]' \
    --output-dir-s3 '~{s3_wd_uri}' \
    --additional-files '{}' \
    --additional-attributes '{}'
  PriceSeqFilter 2> /dev/null | head -n1 > priceseq_version.txt
  >>>
  output {
    String step_description_md = read_string("priceseq_out.description.md")
    File priceseq1_fa = "priceseq1.fa"
    File? priceseq2_fa = "priceseq2.fa"
    File? output_read_count = "priceseq_out.count"
    File? version = "priceseq_version.txt"
  }
  runtime {
    docker: docker_image_id
  }
}

task RunCZIDDedup {
  input {
    String docker_image_id
    String s3_wd_uri
    Array[File] priceseq_fa
  }
  command<<<
  set -euxo pipefail
  idseq-dag-run-step --workflow-name host_filter \
    --step-module idseq_dag.steps.run_czid_dedup \
    --step-class PipelineStepRunCZIDDedup \
    --step-name czid_dedup_out \
    --input-files '[["~{sep='","' priceseq_fa}"]]' \
    --output-files '[~{if length(priceseq_fa) == 2 then '"dedup1.fa", "dedup2.fa"' else '"dedup1.fa"'}, "clusters.csv", "duplicate_cluster_sizes.tsv"]' \
    --output-dir-s3 '~{s3_wd_uri}' \
    --additional-files '{}' \
    --additional-attributes '{}'
  czid-dedup --version > czid_dedup_version.txt
  >>>
  output {
    String step_description_md = read_string("czid_dedup_out.description.md")
    File dedup1_fa = "dedup1.fa"
    File? dedup2_fa = "dedup2.fa"
    File duplicate_clusters_csv = "clusters.csv"
    File duplicate_cluster_sizes_tsv = "duplicate_cluster_sizes.tsv"
    File? output_read_count = "czid_dedup_out.count"
    File? version = "czid_dedup_version.txt"
  }
  runtime {
    docker: docker_image_id
  }
}

task RunLZW {
  input {
    String docker_image_id
    String s3_wd_uri
    Array[File] dedup_fa
    File duplicate_clusters_csv
    File duplicate_cluster_sizes_tsv
  }
  command<<<
  set -euxo pipefail
  idseq-dag-run-step --workflow-name host_filter \
    --step-module idseq_dag.steps.run_lzw \
    --step-class PipelineStepRunLZW \
    --step-name lzw_out \
    --input-files '[["~{sep='","' dedup_fa}", "~{duplicate_clusters_csv}", "~{duplicate_cluster_sizes_tsv}"]]' \
    --output-files '[~{if length(dedup_fa) == 2 then '"lzw1.fa", "lzw2.fa"' else '"lzw1.fa"'}]' \
    --output-dir-s3 '~{s3_wd_uri}' \
    --additional-files '{}' \
    --additional-attributes '{"thresholds": [0.45, 0.42], "threshold_readlength": 150}'
  >>>
  output {
    String step_description_md = read_string("lzw_out.description.md")
    File lzw1_fa = "lzw1.fa"
    File? lzw2_fa = "lzw2.fa"
    File? output_read_count = "lzw_out.count"
  }
  runtime {
    docker: docker_image_id
  }
}

task RunBowtie2_bowtie2_out {
  input {
    String docker_image_id
    String s3_wd_uri
    Array[File] lzw_fa
    Array[File] dedup_fa
    File duplicate_clusters_csv
    File duplicate_cluster_sizes_tsv
    File bowtie2_genome
  }
  command<<<
  set -euxo pipefail
  idseq-dag-run-step --workflow-name host_filter \
    --step-module idseq_dag.steps.run_bowtie2 \
    --step-class PipelineStepRunBowtie2 \
    --step-name bowtie2_out \
    --input-files '[["~{sep='","' lzw_fa}"], ["~{sep='","' dedup_fa}", "~{duplicate_clusters_csv}", "~{duplicate_cluster_sizes_tsv}"]]' \
    --output-files '[~{if length(lzw_fa) == 2 then '"bowtie2_1.fa", "bowtie2_2.fa", "bowtie2_merged.fa"' else '"bowtie2_1.fa"'}]' \
    --output-dir-s3 '~{s3_wd_uri}' \
    --additional-files '{"bowtie2_genome": "~{bowtie2_genome}"}' \
    --additional-attributes '{"output_sam_file": "bowtie2.sam"}'
  bowtie2 --version > bowtie2_version.txt
  >>>
  output {
    String step_description_md = read_string("bowtie2_out.description.md")
    File bowtie2_1_fa = "bowtie2_1.fa"
    File? bowtie2_2_fa = "bowtie2_2.fa"
    File? bowtie2_merged_fa = "bowtie2_merged.fa"
    File? output_read_count = "bowtie2_out.count"
    File? version = "bowtie2_version.txt"
  }
  runtime {
    docker: docker_image_id
  }
}

task RunSubsample {
  input {
    String docker_image_id
    String s3_wd_uri
    Array[File] bowtie2_fa
    Array[File] dedup_fa
    File duplicate_clusters_csv
    File duplicate_cluster_sizes_tsv
    Int max_subsample_fragments
  }
  command<<<
  set -euxo pipefail
  idseq-dag-run-step --workflow-name host_filter \
    --step-module idseq_dag.steps.run_subsample \
    --step-class PipelineStepRunSubsample \
    --step-name subsampled_out \
    --input-files '[["~{sep='","' bowtie2_fa}"], ["~{sep='","' dedup_fa}", "~{duplicate_clusters_csv}", "~{duplicate_cluster_sizes_tsv}"]]' \
    --output-files '[~{if length(dedup_fa) == 2 then '"subsampled_1.fa", "subsampled_2.fa", "subsampled_merged.fa"' else '"subsampled_1.fa"'}]' \
    --output-dir-s3 '~{s3_wd_uri}' \
    --additional-files '{}' \
    --additional-attributes '{"max_fragments": ~{max_subsample_fragments}}'
  >>>
  output {
    String step_description_md = read_string("subsampled_out.description.md")
    File subsampled_1_fa = "subsampled_1.fa"
    File? subsampled_2_fa = "subsampled_2.fa"
    File? subsampled_merged_fa = "subsampled_merged.fa"
    File? output_read_count = "subsampled_out.count"
  }
  runtime {
    docker: docker_image_id
  }
}

task RunStarDownstream {
  input {
    String docker_image_id
    String s3_wd_uri
    Array[File] subsampled_fa
    File validate_input_summary_json
    Array[File] valid_input_fastq
    Array[File] dedup_fa
    File duplicate_clusters_csv
    File duplicate_cluster_sizes_tsv
    File human_star_genome
  }
  command<<<
  set -euxo pipefail
  idseq-dag-run-step --workflow-name host_filter \
    --step-module idseq_dag.steps.run_star_downstream \
    --step-class PipelineStepRunStarDownstream \
    --step-name star_human_out \
    --input-files '[["~{sep='","' subsampled_fa}"], ["~{validate_input_summary_json}", "~{sep='","' valid_input_fastq}"], ["~{sep='","' dedup_fa}", "~{duplicate_clusters_csv}", "~{duplicate_cluster_sizes_tsv}"]]' \
    --output-files '[~{if length(dedup_fa) == 2 then '"unmapped_human_1.fa", "unmapped_human_2.fa"' else '"unmapped_human_1.fa"'}]' \
    --output-dir-s3 '~{s3_wd_uri}' \
    --additional-files '{"star_genome": "~{human_star_genome}"}' \
    --additional-attributes '{}'
  STAR --version > star_human_version.txt
  >>>
  output {
    String step_description_md = read_string("star_human_out.description.md")
    File unmapped_human_1_fa = "unmapped_human_1.fa"
    File? unmapped_human_2_fa = "unmapped_human_2.fa"
    File? output_read_count = "star_human_out.count"
    File? version = "star_human_version.txt"
  }
  runtime {
    docker: docker_image_id
  }
}

task RunBowtie2_bowtie2_human_out {
  input {
    String docker_image_id
    String s3_wd_uri
    Array[File] unmapped_human_fa
    Array[File] dedup_fa
    File duplicate_clusters_csv
    File duplicate_cluster_sizes_tsv
    File human_bowtie2_genome
  }
  command<<<
  set -euxo pipefail
  idseq-dag-run-step --workflow-name host_filter \
    --step-module idseq_dag.steps.run_bowtie2 \
    --step-class PipelineStepRunBowtie2 \
    --step-name bowtie2_human_out \
    --input-files '[["~{sep='","' unmapped_human_fa}"], ["~{sep='","' dedup_fa}", "~{duplicate_clusters_csv}", "~{duplicate_cluster_sizes_tsv}"]]' \
    --output-files '[~{if length(dedup_fa) == 2 then '"bowtie2_human_1.fa", "bowtie2_human_2.fa", "bowtie2_human_merged.fa"' else '"bowtie2_human_1.fa"'}]' \
    --output-dir-s3 '~{s3_wd_uri}' \
    --additional-files '{"bowtie2_genome": "~{human_bowtie2_genome}"}' \
    --additional-attributes '{"output_sam_file": "bowtie2_human.sam"}'
  bowtie2 --version > bowtie2_human_version.txt
  >>>
  output {
    String step_description_md = read_string("bowtie2_human_out.description.md")
    File bowtie2_human_1_fa = "bowtie2_human_1.fa"
    File? bowtie2_human_2_fa = "bowtie2_human_2.fa"
    File? bowtie2_human_merged_fa = "bowtie2_human_merged.fa"
    File? output_read_count = "bowtie2_human_out.count"
    File? version = "bowtie2_human_version.txt"
  }
  runtime {
    docker: docker_image_id
  }
}

task RunGsnapFilter {
  input {
    String docker_image_id
    String s3_wd_uri
    Array[File] subsampled_fa
    Array[File] dedup_fa
    File duplicate_clusters_csv
    File duplicate_cluster_sizes_tsv
    File gsnap_genome
  }
  command<<<
  set -euxo pipefail
  idseq-dag-run-step --workflow-name host_filter \
    --step-module idseq_dag.steps.run_gsnap_filter \
    --step-class PipelineStepRunGsnapFilter \
    --step-name gsnap_filter_out \
    --input-files '[["~{sep='","' subsampled_fa}"], ["~{sep='","' dedup_fa}", "~{duplicate_clusters_csv}", "~{duplicate_cluster_sizes_tsv}"]]' \
    --output-files '[~{if length(dedup_fa) == 2 then '"gsnap_filter_1.fa", "gsnap_filter_2.fa", "gsnap_filter_merged.fa"' else '"gsnap_filter_1.fa"'}]' \
    --output-dir-s3 '~{s3_wd_uri}' \
    --additional-files '{"gsnap_genome": "~{gsnap_genome}"}' \
    --additional-attributes '{"output_sam_file": "gsnap_filter.sam"}'
  gsnap --version > gsnap_filter_version.txt
  >>>
  output {
    String step_description_md = read_string("gsnap_filter_out.description.md")
    File gsnap_filter_1_fa = "gsnap_filter_1.fa"
    File? gsnap_filter_2_fa = "gsnap_filter_2.fa"
    File? gsnap_filter_merged_fa = "gsnap_filter_merged.fa"
    File? output_read_count = "gsnap_filter_out.count"
    File? version = "gsnap_filter_version.txt"
  }
  runtime {
    docker: docker_image_id
  }
}


workflow czid_host_filter {
  input {
    String docker_image_id
    String s3_wd_uri
    File fastqs_0
    File? fastqs_1
    String file_ext
    String nucleotide_type
    String host_genome
    File adapter_fasta
    File star_genome
    File bowtie2_genome
    File gsnap_genome = "s3://czid-public-references/host_filter/human/2018-02-15-utc-1518652800-unixtime__2018-02-15-utc-1518652800-unixtime/hg38_pantro5_k16.tar"
    String human_star_genome
    String human_bowtie2_genome
    Int max_input_fragments
    Int max_subsample_fragments
  }

  call RunValidateInput {
    input:
      docker_image_id = docker_image_id,
      s3_wd_uri = s3_wd_uri,
      fastqs = select_all([fastqs_0, fastqs_1]),
      file_ext = file_ext,
      max_input_fragments = max_input_fragments
  }

  call RunStar {
    input:
      docker_image_id = docker_image_id,
      s3_wd_uri = s3_wd_uri,
      validate_input_summary_json = RunValidateInput.validate_input_summary_json,
      valid_input_fastq = select_all([RunValidateInput.valid_input1_fastq, RunValidateInput.valid_input2_fastq]),
      star_genome = star_genome,
      nucleotide_type = nucleotide_type,
      host_genome = host_genome
  }

  call RunTrimmomatic {
    input:
      docker_image_id = docker_image_id,
      s3_wd_uri = s3_wd_uri,
      unmapped_fastq = select_all([RunStar.unmapped1_fastq, RunStar.unmapped2_fastq]),
      adapter_fasta = adapter_fasta
  }

  call RunPriceSeq {
    input:
      docker_image_id = docker_image_id,
      s3_wd_uri = s3_wd_uri,
      trimmomatic_fastq = select_all([RunTrimmomatic.trimmomatic1_fastq, RunTrimmomatic.trimmomatic2_fastq])
  }

  call RunCZIDDedup {
    input:
      docker_image_id = docker_image_id,
      s3_wd_uri = s3_wd_uri,
      priceseq_fa = select_all([RunPriceSeq.priceseq1_fa, RunPriceSeq.priceseq2_fa])
  }

  call RunLZW {
    input:
      docker_image_id = docker_image_id,
      s3_wd_uri = s3_wd_uri,
      dedup_fa = select_all([RunCZIDDedup.dedup1_fa, RunCZIDDedup.dedup2_fa]),
      duplicate_clusters_csv = RunCZIDDedup.duplicate_clusters_csv,
      duplicate_cluster_sizes_tsv = RunCZIDDedup.duplicate_cluster_sizes_tsv
  }

  call RunBowtie2_bowtie2_out {
    input:
      docker_image_id = docker_image_id,
      s3_wd_uri = s3_wd_uri,
      lzw_fa = select_all([RunLZW.lzw1_fa, RunLZW.lzw2_fa]),
      dedup_fa = select_all([RunCZIDDedup.dedup1_fa, RunCZIDDedup.dedup2_fa]),
      duplicate_clusters_csv = RunCZIDDedup.duplicate_clusters_csv,
      duplicate_cluster_sizes_tsv = RunCZIDDedup.duplicate_cluster_sizes_tsv,
      bowtie2_genome = bowtie2_genome
  }

  call RunSubsample {
    input:
      docker_image_id = docker_image_id,
      s3_wd_uri = s3_wd_uri,
      bowtie2_fa = select_all([RunBowtie2_bowtie2_out.bowtie2_1_fa, RunBowtie2_bowtie2_out.bowtie2_2_fa, RunBowtie2_bowtie2_out.bowtie2_merged_fa]),
      dedup_fa = select_all([RunCZIDDedup.dedup1_fa, RunCZIDDedup.dedup2_fa]),
      duplicate_clusters_csv = RunCZIDDedup.duplicate_clusters_csv,
      duplicate_cluster_sizes_tsv = RunCZIDDedup.duplicate_cluster_sizes_tsv,
      max_subsample_fragments = max_subsample_fragments
  }

  if (host_genome != "human") {
    call RunStarDownstream {
      input:
        docker_image_id = docker_image_id,
        s3_wd_uri = s3_wd_uri,
        subsampled_fa = select_all([RunSubsample.subsampled_1_fa, RunSubsample.subsampled_2_fa, RunSubsample.subsampled_merged_fa]),
        validate_input_summary_json = RunValidateInput.validate_input_summary_json,
        valid_input_fastq = select_all([RunValidateInput.valid_input1_fastq, RunValidateInput.valid_input2_fastq]),
        dedup_fa = select_all([RunCZIDDedup.dedup1_fa, RunCZIDDedup.dedup2_fa]),
        duplicate_clusters_csv = RunCZIDDedup.duplicate_clusters_csv,
        duplicate_cluster_sizes_tsv = RunCZIDDedup.duplicate_cluster_sizes_tsv,
        human_star_genome = human_star_genome
    }

    call RunBowtie2_bowtie2_human_out {
      input:
        docker_image_id = docker_image_id,
        s3_wd_uri = s3_wd_uri,
        unmapped_human_fa = select_all([RunStarDownstream.unmapped_human_1_fa, RunStarDownstream.unmapped_human_2_fa]),
        dedup_fa = select_all([RunCZIDDedup.dedup1_fa, RunCZIDDedup.dedup2_fa]),
        duplicate_clusters_csv = RunCZIDDedup.duplicate_clusters_csv,
        duplicate_cluster_sizes_tsv = RunCZIDDedup.duplicate_cluster_sizes_tsv,
        human_bowtie2_genome = human_bowtie2_genome
    }
  }

  Array[File] gsnap_filter_input = if (host_genome == "human")
    then select_all([RunSubsample.subsampled_1_fa, RunSubsample.subsampled_2_fa, RunSubsample.subsampled_merged_fa])
    else select_all([RunBowtie2_bowtie2_human_out.bowtie2_human_1_fa, RunBowtie2_bowtie2_human_out.bowtie2_human_2_fa, RunBowtie2_bowtie2_human_out.bowtie2_human_merged_fa])

  call RunGsnapFilter {
    input:
      docker_image_id = docker_image_id,
      s3_wd_uri = s3_wd_uri,
      subsampled_fa = gsnap_filter_input,
      dedup_fa = select_all([RunCZIDDedup.dedup1_fa, RunCZIDDedup.dedup2_fa]),
      duplicate_clusters_csv = RunCZIDDedup.duplicate_clusters_csv,
      duplicate_cluster_sizes_tsv = RunCZIDDedup.duplicate_cluster_sizes_tsv,
      gsnap_genome = gsnap_genome
  }

  output {
    File validate_input_out_validate_input_summary_json = RunValidateInput.validate_input_summary_json
    File? validate_input_out_count = RunValidateInput.output_read_count
    File star_out_unmapped1_fastq = RunStar.unmapped1_fastq
    File? star_out_unmapped2_fastq = RunStar.unmapped2_fastq
    File? star_out_log_file = RunStar.output_log_file
    File? star_out_count = RunStar.output_read_count
    File? star_version = RunStar.version
    File trimmomatic_out_trimmomatic1_fastq = RunTrimmomatic.trimmomatic1_fastq
    File? trimmomatic_out_trimmomatic2_fastq = RunTrimmomatic.trimmomatic2_fastq
    File? trimmomatic_out_count = RunTrimmomatic.output_read_count
    File? trimmomatic_version = RunTrimmomatic.version
    File priceseq_out_priceseq1_fa = RunPriceSeq.priceseq1_fa
    File? priceseq_out_priceseq2_fa = RunPriceSeq.priceseq2_fa
    File? priceseq_out_count = RunPriceSeq.output_read_count
    File? priceseq_version = RunPriceSeq.version
    File czid_dedup_out_dedup1_fa = RunCZIDDedup.dedup1_fa
    File? czid_dedup_out_dedup2_fa = RunCZIDDedup.dedup2_fa
    File czid_dedup_out_duplicate_clusters_csv = RunCZIDDedup.duplicate_clusters_csv
    File czid_dedup_out_duplicate_cluster_sizes_tsv = RunCZIDDedup.duplicate_cluster_sizes_tsv
    File? czid_dedup_out_count = RunCZIDDedup.output_read_count
    File? czid_dedup_version = RunCZIDDedup.version
    File lzw_out_lzw1_fa = RunLZW.lzw1_fa
    File? lzw_out_lzw2_fa = RunLZW.lzw2_fa
    File? lzw_out_count = RunLZW.output_read_count
    File bowtie2_out_bowtie2_1_fa = RunBowtie2_bowtie2_out.bowtie2_1_fa
    File? bowtie2_out_bowtie2_2_fa = RunBowtie2_bowtie2_out.bowtie2_2_fa
    File? bowtie2_out_bowtie2_merged_fa = RunBowtie2_bowtie2_out.bowtie2_merged_fa
    File? bowtie2_out_count = RunBowtie2_bowtie2_out.output_read_count
    File? bowtie2_version = RunBowtie2_bowtie2_out.version
    File subsampled_out_subsampled_1_fa = RunSubsample.subsampled_1_fa
    File? subsampled_out_subsampled_2_fa = RunSubsample.subsampled_2_fa
    File? subsampled_out_subsampled_merged_fa = RunSubsample.subsampled_merged_fa
    File? subsampled_out_count = RunSubsample.output_read_count
    File? star_human_out_unmapped_human_1_fa = RunStarDownstream.unmapped_human_1_fa
    File? star_human_out_unmapped_human_2_fa = RunStarDownstream.unmapped_human_2_fa
    File? star_human_out_count = RunStarDownstream.output_read_count
    File? star_human_version = RunStarDownstream.version
    File? bowtie2_human_out_bowtie2_human_1_fa = RunBowtie2_bowtie2_human_out.bowtie2_human_1_fa
    File? bowtie2_human_out_bowtie2_human_2_fa = RunBowtie2_bowtie2_human_out.bowtie2_human_2_fa
    File? bowtie2_human_out_bowtie2_human_merged_fa = RunBowtie2_bowtie2_human_out.bowtie2_human_merged_fa
    File? bowtie2_human_out_count = RunBowtie2_bowtie2_human_out.output_read_count
    File? bowtie2_human_version = RunBowtie2_bowtie2_human_out.version
    File gsnap_filter_out_gsnap_filter_1_fa = RunGsnapFilter.gsnap_filter_1_fa
    File? gsnap_filter_out_gsnap_filter_2_fa = RunGsnapFilter.gsnap_filter_2_fa
    File? gsnap_filter_out_gsnap_filter_merged_fa = RunGsnapFilter.gsnap_filter_merged_fa
    File? gsnap_filter_out_count = RunGsnapFilter.output_read_count
    File? gsnap_filter_version = RunGsnapFilter.version
    File? input_read_count = RunValidateInput.input_read_count
    File? output_gene_file = RunStar.output_gene_file
    File? output_metrics_file = RunStar.output_metrics_file
    File? output_histogram_file = RunStar.output_histogram_file
  }
}
