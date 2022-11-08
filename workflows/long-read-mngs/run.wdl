version 1.0

task RunValidateInput {
    input {
        File input_fastq
        String docker_image_id
    }

    command <<<
        set -euxo pipefail
        filter_count "~{input_fastq}" original "No reads provided"
        fastp -i "~{input_fastq}" -o sample_validated.fastq
        filter_count sample_validated.fastq validated "No reads remaining after input validation"
    >>>

    output {
        File validated_output = "sample_validated.fastq"
        File raw_reads = "original_reads.count"
        File raw_bases = "original_bases.count"
        File validated_reads = "validated_reads.count"
        File validated_bases = "validated_bases.count"
    }

    runtime {
        docker: docker_image_id
    }
}

task RunQualityFilter {
    input {
        File input_fastq
        String docker_image_id
    }

    command <<<
        set -euxo pipefail
        fastp -i "~{input_fastq}" --qualified_quality_phred 9 --length_required 100 --low_complexity_filter --complexity_threshold 30 --dont_eval_duplication -o sample_quality_filtered.fastq
        filter_count sample_quality_filtered.fastq quality_filtered "No reads remaining after quality filtering"
    >>>

    output {
        File fastp_output = "sample_quality_filtered.fastq"
        File quality_filtered_reads = "quality_filtered_reads.count"
        File quality_filtered_bases = "quality_filtered_bases.count"
    }

    runtime {
        docker: docker_image_id
    }
}

task RunHostFilter {
    input {
        File input_fastq
        String library_type
        File minimap_host_db
        String docker_image_id
    }

    command <<<
        set -euxo pipefail
        # run minimap2 against host genome
        if [[ "~{library_type}" == "RNA" ]]
        then
            minimap2 -t $(nproc) -ax splice "~{minimap_host_db}" "~{input_fastq}" -o sample.hostfiltered.sam --split-prefix temp_name
        else # assuming DNA
            minimap2 -t $(nproc) -ax map-ont "~{minimap_host_db}" "~{input_fastq}" -o sample.hostfiltered.sam --split-prefix temp_name
        fi
        samtools fastq -n -f 4 sample.hostfiltered.sam > sample.hostfiltered.fastq
        samtools view -S -b sample.hostfiltered.sam > sample.hostfiltered.bam

        filter_count sample.hostfiltered.fastq host_filtered "No reads remaining after host filtering"
    >>>

    output {
        File host_filter_sam = "sample.hostfiltered.bam"
        File host_filter_fastq = "sample.hostfiltered.fastq"
        File host_filtered_reads = "host_filtered_reads.count"
        File host_filtered_bases = "host_filtered_bases.count"
    }

    runtime {
        docker: docker_image_id
    }
}

task RunHumanFilter {
    input {
        File input_fastq
        String library_type
        File minimap_human_db
        String docker_image_id
    }

    command <<<
        set -euxo pipefail
        # run minimap2 against human genome
        if [[ "~{library_type}" == "RNA" ]]
        then
            echo "DEBUG: inside library_type == RNA, running minimap2 -ax splice" >> output.txt
            minimap2 -t $(nproc) -ax splice "~{minimap_human_db}" "~{input_fastq}" -o sample.humanfiltered.sam -t 63 --split-prefix temp_name
        else # assuming DNA
            echo "DEBUG: inside library_type == DNA, running minimap2 -ax map-ont" >> output.txt
            minimap2 -t $(nproc) -ax map-ont "~{minimap_human_db}" "~{input_fastq}" -o sample.humanfiltered.sam -t 63 --split-prefix temp_name
        fi
        # extract the unmapped reads for downstream processing
        echo "about to run samtools fastq" >> output.txt
        samtools fastq -n -f 4 sample.humanfiltered.sam > sample.humanfiltered.fastq
        samtools view -S -b sample.humanfiltered.sam > sample.humanfiltered.bam

        filter_count sample.humanfiltered.fastq human_filtered "No reads remaining after human read filtering"
    >>>

    output {
        File human_filter_sam = "sample.humanfiltered.bam"
        File human_filter_fastq = "sample.humanfiltered.fastq"
        File human_filtered_reads = "human_filtered_reads.count"
        File human_filtered_bases = "human_filtered_bases.count"
    }

    runtime {
        docker: docker_image_id
    }
}


task RunSubsampling {
    input {
        File input_fastq
        Int subsample_depth
        String docker_image_id
    }

    command <<<
        set -euxo pipefail
        head -"~{subsample_depth}" "~{input_fastq}" > sample.subsampled.fastq

        # We should always have reads after subsampling, but adding for consistency with other steps
        filter_count sample.subsampled.fastq subsampled "No reads remaining after subsampling"
    >>>

    output {
        File subsampled_fastq = "sample.subsampled.fastq"
        File subsampled_reads = "subsampled_reads.count"
        File subsampled_bases = "subsampled_bases.count"
    }

    runtime {
        docker: docker_image_id
    }
}

task RunAssembly {
    input {
        File input_fastq
        String guppy_basecaller_setting
        Int polishing_iterations
        String docker_image_id
    }

    command <<<
        set -euxo pipefail
        flye_setting="--nano-raw"
        if [[ "~{guppy_basecaller_setting}" == "super" ]]
        then
            flye_setting="--nano-hq"
        fi

        # run flye to assembly contigs
        flye --threads $(nproc) --meta $flye_setting "~{input_fastq}" --out-dir temp_flye_out --iterations "~{polishing_iterations}"

        # ERROR HANDLING - assembly somethings fails (due to low coverage) and is then missing...
        #                  ... the temp_flye_out/assembly.fasta file
        if [ -f temp_flye_out/assembly.fasta ]
        then
            cat temp_flye_out/assembly.fasta > sample.assembled_reads.fasta
        else
            #just copy original .fastq to .fasta
            seqtk seq -a "~{input_fastq}" > sample.assembled_reads.fasta 
        fi

        zip -r temp_flye_out.zip temp_flye_out
    >>>

    output {
        File assembled_fasta = "sample.assembled_reads.fasta"
        File temp_assembly_dir = "temp_flye_out.zip"
    }

    runtime {
        docker: docker_image_id
    }
}

task RunReadsToContigs {
    input {
        File input_fastq
        File assembled_reads
        String docker_image_id
    }

    command <<<
        set -euxo pipefail
        # use minimap2 to align reads back to contigs
        minimap2 -ax map-ont "~{assembled_reads}" "~{input_fastq}" -o sample.reads_to_contigs.sam -t 15
        samtools view -b sample.reads_to_contigs.sam | samtools sort > sample.reads_to_contigs.bam
        samtools index sample.reads_to_contigs.bam sample.reads_to_contigs.bam.bai

        # extract reads that didn't map to contigs as non-contig reads 
        samtools fastq -n -f 4 sample.reads_to_contigs.sam > sample.non_contigs.fastq

        # convert non-contigs.fastq file to .fasta file
        seqtk seq -a sample.non_contigs.fastq > sample.non_contigs.fasta
        cat sample.reads_to_contigs.sam | grep -v "^@" | cut -f1,3,10 > reads_to_contigs.tsv
    >>>

    output {
        File reads_to_contigs_sam = "sample.reads_to_contigs.sam"
        File reads_to_contigs_bam = "sample.reads_to_contigs.bam"
        File reads_to_contigs_bai = "sample.reads_to_contigs.bam.bai"
        File reads_to_contigs_tsv = "reads_to_contigs.tsv"
        File non_contigs_fastq = "sample.non_contigs.fastq"
        File non_contigs_fasta = "sample.non_contigs.fasta"
    }

    runtime {
        docker: docker_image_id
    }
}

task PrepareNTAlignmentInputs {
    input {
        File non_contig_reads_fa
        File assembled_reads_fa
        String docker_image_id
    }

    command <<<
        set -euxo pipefail
        # combine contigs and non-contig reads into a single file
        cat "~{assembled_reads_fa}" "~{non_contig_reads_fa}" > sample.all_sequences_to_align_full.fasta

        # remove duplicates from full .fasta (important when assembly failed)
        seqkit rmdup -s < sample.all_sequences_to_align_full.fasta > sample.all_sequences_to_align.fasta
    >>>

    output {
        File all_sequences_to_align = "sample.all_sequences_to_align.fasta"
    }

    runtime {
        docker: docker_image_id
    }
}

task RunNTAlignment {
    input {
        File all_sequences_to_align
        String? db_path
        String minimap2_args 
        Boolean run_locally = false
        File? local_minimap2_index 
        String prefix
        # only required for remote alignment
        String? s3_wd_uri
        String docker_image_id
    }

    command <<<
        set -euxo pipefail
        if [[ "~{run_locally}" == true ]]; then
          minimap2-scatter ~{minimap2_args} "~{local_minimap2_index}" "~{all_sequences_to_align}" > "~{prefix}.paf"
        else
          python3 <<CODE
        import os
        from idseq_utils.batch_run_helpers import run_alignment

        run_alignment(
            input_dir="~{s3_wd_uri}",
            db_path="~{db_path}",
            result_path="gsnap.paf",
            aligner="minimap2",
            aligner_args="~{minimap2_args}",
            queries=["~{all_sequences_to_align}"],
        )
        CODE
        fi
        python3 /usr/local/lib/python3.10/dist-packages/idseq_utils/paf2blast6.py gsnap.paf
        mv *frompaf.m8 "gsnap.m8" # TODO: rewrite paf2blast6.py to output in this format
    >>>

    output {
        File out_paf = "gsnap.paf"
        File out_m8 = "gsnap.m8"
    }

    runtime {
        docker: docker_image_id
    }
}

# TODO: (tmorse) fuse me
task RunCallHitsNT {
    input {
        File m8_file
        File lineage_db
        File taxon_blacklist
        File deuterostome_db
        File accession2taxid
        Int min_read_length = 36
        String docker_image_id
        String count_type = "NT"
    }

    command <<<
        set -euxo pipefail
        python3 <<CODE
        from idseq_dag.util.m8 import call_hits_m8, generate_taxon_count_json_from_m8
        call_hits_m8(
            input_m8="~{m8_file}",
            lineage_map_path="~{lineage_db}",
            accession2taxid_dict_path="~{accession2taxid}",
            output_m8="gsnap.deduped.m8",
            output_summary="gsnap.hitsummary.tab",
            min_alignment_length=~{min_read_length},
            deuterostome_path="~{deuterostome_db}",
            taxon_whitelist_path=None,
            taxon_blacklist_path="~{taxon_blacklist}",
        )
        generate_taxon_count_json_from_m8(
            blastn_6_path="gsnap.deduped.m8",
            hit_level_path="gsnap.hitsummary.tab",
            count_type="~{count_type}",
            lineage_map_path="~{lineage_db}",
            deuterostome_path="~{deuterostome_db}",
            taxon_whitelist_path=None,
            taxon_blacklist_path="~{taxon_blacklist}",
            duplicate_cluster_sizes_path=None,
            output_json_file="gsnap_counts_with_dcr.json",
        )
        CODE
    >>>

    output {
        File deduped_out_m8 = "gsnap.deduped.m8"
        File hitsummary = "gsnap.hitsummary.tab"
        File counts_json = "gsnap_counts_with_dcr.json"
    }

    runtime {
        docker: docker_image_id
    }
}

task RunNRAlignment {
    input {
        File assembled_reads_fa
        String? db_path
        String diamond_args 
        Boolean run_locally = false
        File? local_diamond_index 
        # only required for remote alignment
        String? s3_wd_uri
        String docker_image_id
    }

    command <<<
        set -euxo pipefail  
        if [[ "~{run_locally}" == true ]]; then 
          diamond makedb --in "~{local_diamond_index}" -d reference
          diamond blastx -d reference -q "~{assembled_reads_fa}" -o "diamond.m8" "--~{diamond_args}"
        else
          python3 <<CODE
        import os 
        from idseq_utils.batch_run_helpers import run_alignment

        run_alignment(
            input_dir="~{s3_wd_uri}",
            db_path="~{db_path}",
            result_path="diamond.m8",
            aligner="diamond",
            aligner_args="~{diamond_args}",
            queries=["~{assembled_reads_fa}"],
        )
        CODE
        diamond --version > diamond_version.txt
        fi
    >>>

    output {
        File out_m8 = "diamond.m8"
    }

    runtime {
        docker: docker_image_id
    }
}

# TODO: (tmorse) fuse me
task RunCallHitsNR {
    input {
        File m8_file
        File lineage_db
        File taxon_blacklist
        File deuterostome_db
        File accession2taxid
        Int min_read_length = 36
        String docker_image_id
        String count_type = "NR"
    }

    command <<<
        set -euxo pipefail
        python3 <<CODE
        from idseq_dag.util.m8 import call_hits_m8, generate_taxon_count_json_from_m8
        call_hits_m8(
            input_m8="~{m8_file}",
            lineage_map_path="~{lineage_db}",
            accession2taxid_dict_path="~{accession2taxid}",
            output_m8="rapsearch2.deduped.m8",
            output_summary="rapsearch2.hitsummary.tab",
            min_alignment_length=~{min_read_length},
            deuterostome_path="~{deuterostome_db}",
            taxon_whitelist_path=None,
            taxon_blacklist_path="~{taxon_blacklist}",
        )
        generate_taxon_count_json_from_m8(
            blastn_6_path="rapsearch2.deduped.m8",
            hit_level_path="rapsearch2.hitsummary.tab",
            count_type="~{count_type}",
            lineage_map_path="~{lineage_db}",
            deuterostome_path="~{deuterostome_db}",
            taxon_whitelist_path=None,
            taxon_blacklist_path="~{taxon_blacklist}",
            duplicate_cluster_sizes_path=None,
            output_json_file="rapsearch2_counts_with_dcr.json",
        )
        CODE
    >>>

    output {
        File deduped_out_m8 = "rapsearch2.deduped.m8"
        File hitsummary = "rapsearch2.hitsummary.tab"
        File counts_json = "rapsearch2_counts_with_dcr.json"
    }

    runtime {
        docker: docker_image_id
    }
}

# TODO: (tmorse) fuse me
task ReassignM8NT {
    input {
        File m8
        File reads_to_contigs_tsv
        String docker_image_id
    }

    command <<<
        set -euxo pipefail  
        python3 /usr/local/bin/reassign_m8.py \
            --m8-filepath "~{m8}" \
            --reads-to-contigs-filepath "~{reads_to_contigs_tsv}" \
            --output-filepath "m8_reassigned_nt.tab" \
    >>>

    output {
        File m8_reassigned = "m8_reassigned_nt.tab"
    }

    runtime {
        docker: docker_image_id
    }
}

# TODO: (tmorse) fuse me
task ReassignM8NR {
    input {
        File m8
        File reads_to_contigs_tsv
        String docker_image_id
    }

    command <<<
        set -euxo pipefail  
        python3 /usr/local/bin/reassign_m8.py \
            --m8-filepath "~{m8}" \
            --reads-to-contigs-filepath "~{reads_to_contigs_tsv}" \
            --output-filepath "m8_reassigned_nr.tab" \
    >>>

    output {
        File m8_reassigned = "m8_reassigned_nr.tab"
    }

    runtime {
        docker: docker_image_id
    }
}

# TODO: (tmorse) fuse me
task TallyHitsNT {
    input {
        File reads_fastq
        File m8
        File hitsummary
        File reads_to_contigs_tsv
        String docker_image_id
    }

    command <<<
        set -euxo pipefail  
        python3 /usr/local/bin/tally_counts.py \
            --reads-fastq-filepath "~{reads_fastq}" \
            --m8-filepath "~{m8}" \
            --hitsummary-filepath "~{hitsummary}" \
            --reads-to-contigs-filepath "~{reads_to_contigs_tsv}" \
            --output-filepath "tallied_hits_nt.csv" \
    >>>

    output {
        File tallied_hits = "tallied_hits_nt.csv"
    }

    runtime {
        docker: docker_image_id
    }
}

# TODO: (tmorse) fuse me
task TallyHitsNR {
    input {
        File reads_fastq
        File m8
        File hitsummary
        File reads_to_contigs_tsv
        String docker_image_id
    }

    command <<<
        set -euxo pipefail  
        python3 /usr/local/bin/tally_counts.py \
            --reads-fastq-filepath "~{reads_fastq}" \
            --m8-filepath "~{m8}" \
            --hitsummary-filepath "~{hitsummary}" \
            --reads-to-contigs-filepath "~{reads_to_contigs_tsv}" \
            --output-filepath "tallied_hits_nr.csv" \
    >>>

    output {
        File tallied_hits = "tallied_hits_nr.csv"
    }

    runtime {
        docker: docker_image_id
    }
}

task UnmappedReads {
    input {
        File input_file
        File hitsummary_nt
        File hitsummary_nr
        File reads_to_contigs_tsv
        String docker_image_id
    }
    command <<<
        set -euxo pipefail
        cut -f1 "~{hitsummary_nt}" "~{hitsummary_nr}" | sort | uniq > mapped.txt
        cut -f1,2 "~{reads_to_contigs_tsv}" > all_reads.txt
        python3 << CODE
        with open("mapped.txt", "r") as m:
            mapped_hits = set(m.read().splitlines())

        unmapped_reads = set()
        with open("all_reads.txt", "r") as ar:
            for line in ar:
                read, contig = line.strip().split("\t")
                if read not in mapped_hits and contig not in mapped_hits:
                    unmapped_reads.add(read)

        with open("unmapped_reads.txt", "w") as ur:
            for read in unmapped_reads:
                ur.write(read)
                ur.write("\n")
        CODE
        seqkit grep -f "unmapped_reads.txt" "~{input_file}" > unmapped_reads.fastq
    >>>

    output {
        File unmapped_reads = "unmapped_reads.fastq"
    }

    runtime {
        docker: docker_image_id
    }
}

task GenerateAnnotatedFasta {
    input {
        File pre_alignment_fasta
        File nt_m8
        File nr_m8
        String docker_image_id
    }

    command <<<
        set -euxo pipefail
        python3 <<CODE
        from idseq_dag.steps.generate_annotated_fasta import generate_annotated_fasta

        generate_annotated_fasta(
            pre_alignment_fa_path = "~{pre_alignment_fasta}",
            nt_m8_path = "~{nt_m8}",
            nr_m8_path = "~{nr_m8}",
            annotated_fasta_path = "refined_annotated_merged.fa",
            unidentified_fasta_path = "refined_unidentified.fa",
        )
        CODE
    >>>

    output {
        # String step_description_md = read_string("refined_annotated_out.description.md")
        File assembly_refined_annotated_merged_fa = "refined_annotated_merged.fa"
        File assembly_refined_unidentified_fa = "refined_unidentified.fa"
    }

  runtime {
    docker: docker_image_id
  }
}

task GenerateTaxidFasta {
    input {
        File annotated_merged_fa
        File nt_hitsummary_tab
        File nr_hitsummary_tab
        File lineage_db
        String docker_image_id
    }

    command <<<
        set -euxo pipefail
        python3 <<CODE
        from idseq_dag.steps.generate_taxid_fasta import generate_taxid_fasta

        generate_taxid_fasta(
            "~{annotated_merged_fa}",
            "~{nt_hitsummary_tab}",
            "~{nr_hitsummary_tab}",
            "~{lineage_db}",
            "refined_taxid_annot.fasta",
        )
        CODE
    >>>

    output {
        # String step_description_md = read_string("refined_taxid_fasta_out.description.md")
        File assembly_refined_taxid_annot_fasta = "refined_taxid_annot.fasta"
        File? output_read_count = "refined_taxid_fasta_out.count"
    }

    runtime {
        docker: docker_image_id
    }
}

task GenerateTaxidLocator {
    input {
        File assembly_refined_taxid_annot_fasta
        String docker_image_id
    }

    command <<<
        set -euxo pipefail
        idseq-dag-run-step --workflow-name postprocess \
            --step-module idseq_dag.steps.generate_taxid_locator \
            --step-class PipelineStepGenerateTaxidLocator \
            --step-name refined_taxid_locator_out \
            --input-files '[["~{assembly_refined_taxid_annot_fasta}"]]' \
            --output-files '["assembly/refined_taxid_annot_sorted_nt.fasta", "assembly/refined_taxid_locations_nt.json", "assembly/refined_taxid_annot_sorted_nr.fasta", "assembly/refined_taxid_locations_nr.json", "assembly/refined_taxid_annot_sorted_genus_nt.fasta", "assembly/refined_taxid_locations_genus_nt.json", "assembly/refined_taxid_annot_sorted_genus_nr.fasta", "assembly/refined_taxid_locations_genus_nr.json", "assembly/refined_taxid_annot_sorted_family_nt.fasta", "assembly/refined_taxid_locations_family_nt.json", "assembly/refined_taxid_annot_sorted_family_nr.fasta", "assembly/refined_taxid_locations_family_nr.json", "assembly/refined_taxid_locations_combined.json"]' \
            --output-dir-s3 '' \
            --additional-files '{}' \
            --additional-attributes '{}'
    >>>

    output {
        String step_description_md = read_string("refined_taxid_locator_out.description.md")
        File assembly_refined_taxid_annot_sorted_nt_fasta = "assembly/refined_taxid_annot_sorted_nt.fasta"
        File assembly_refined_taxid_locations_nt_json = "assembly/refined_taxid_locations_nt.json"
        File assembly_refined_taxid_annot_sorted_nr_fasta = "assembly/refined_taxid_annot_sorted_nr.fasta"
        File assembly_refined_taxid_locations_nr_json = "assembly/refined_taxid_locations_nr.json"
        File assembly_refined_taxid_annot_sorted_genus_nt_fasta = "assembly/refined_taxid_annot_sorted_genus_nt.fasta"
        File assembly_refined_taxid_locations_genus_nt_json = "assembly/refined_taxid_locations_genus_nt.json"
        File assembly_refined_taxid_annot_sorted_genus_nr_fasta = "assembly/refined_taxid_annot_sorted_genus_nr.fasta"
        File assembly_refined_taxid_locations_genus_nr_json = "assembly/refined_taxid_locations_genus_nr.json"
        File assembly_refined_taxid_annot_sorted_family_nt_fasta = "assembly/refined_taxid_annot_sorted_family_nt.fasta"
        File assembly_refined_taxid_locations_family_nt_json = "assembly/refined_taxid_locations_family_nt.json"
        File assembly_refined_taxid_annot_sorted_family_nr_fasta = "assembly/refined_taxid_annot_sorted_family_nr.fasta"
        File assembly_refined_taxid_locations_family_nr_json = "assembly/refined_taxid_locations_family_nr.json"
        File assembly_refined_taxid_locations_combined_json = "assembly/refined_taxid_locations_combined.json"
        File? output_read_count = "refined_taxid_locator_out.count"
    }

    runtime {
        docker: docker_image_id
    }
}

workflow czid_long_read_mngs {
    input {
        String docker_image_id
        # this is required for remote alignment
        String? s3_wd_uri

        File input_fastq

        String library_type = "RNA"
        String guppy_basecaller_setting = "hac" # fast, hac, super

        Int subsample_depth = 400000 # should be 4x the number of reads desired

        File minimap_host_db
        File minimap_human_db

        Int polishing_iterations = 1

        File? minimap2_local_db_path
        File? diamond_local_db_path

        File lineage_db
        File accession2taxid_db
        File taxon_blacklist
        Int min_read_length = 36
        File deuterostome_db

        String? minimap2_db
        String minimap2_args = "-cx asm20 --secondary=yes"
        String minimap2_prefix = "gsnap"

        String? diamond_db
        String diamond_args = "long-reads"
    }

    call RunValidateInput {
        input:
            input_fastq = input_fastq,
            docker_image_id = docker_image_id
    }

    call RunQualityFilter {
        input:
            input_fastq = RunValidateInput.validated_output,
            docker_image_id = docker_image_id
    }

    call RunHostFilter {
        input:
            input_fastq = RunQualityFilter.fastp_output,
            library_type = library_type,
            minimap_host_db = minimap_host_db,
            docker_image_id = docker_image_id
    }

    call RunHumanFilter {
        input:
            input_fastq = RunHostFilter.host_filter_fastq,
            library_type = library_type,
            minimap_human_db = minimap_human_db,
            docker_image_id = docker_image_id
    }

    call RunSubsampling {
        input:
            input_fastq = RunHumanFilter.human_filter_fastq,
            subsample_depth = subsample_depth,
            docker_image_id = docker_image_id
    }

    call RunAssembly {
        input:
            input_fastq = RunSubsampling.subsampled_fastq,
            guppy_basecaller_setting = guppy_basecaller_setting,
            polishing_iterations = polishing_iterations,
            docker_image_id = docker_image_id
    }

    call RunReadsToContigs {
        input:
            input_fastq = RunSubsampling.subsampled_fastq,
            assembled_reads = RunAssembly.assembled_fasta,
            docker_image_id = docker_image_id
    }


    call PrepareNTAlignmentInputs {
        input:
            non_contig_reads_fa = RunReadsToContigs.non_contigs_fasta,
            assembled_reads_fa = RunAssembly.assembled_fasta,
            docker_image_id = docker_image_id
    }

    call RunNTAlignment {
        input:
            all_sequences_to_align = PrepareNTAlignmentInputs.all_sequences_to_align,
            s3_wd_uri = s3_wd_uri,

            db_path = minimap2_db,
            minimap2_args = minimap2_args,
            run_locally = defined(minimap2_local_db_path),
            local_minimap2_index = minimap2_local_db_path,
            prefix= minimap2_prefix,
            docker_image_id = docker_image_id
    }

    call RunCallHitsNT { 
        input:
            m8_file = RunNTAlignment.out_m8,
            lineage_db = lineage_db,
            taxon_blacklist = taxon_blacklist,
            deuterostome_db = deuterostome_db,
            accession2taxid = accession2taxid_db,
            min_read_length = min_read_length,
            count_type = library_type,
            docker_image_id = docker_image_id,
    }

    call RunNRAlignment {
        input:
            assembled_reads_fa=RunAssembly.assembled_fasta,
            db_path=diamond_db,
            diamond_args=diamond_args,
            run_locally=defined(diamond_local_db_path),
            local_diamond_index=diamond_local_db_path,
            s3_wd_uri=s3_wd_uri,
            docker_image_id=docker_image_id
    }

    call RunCallHitsNR { 
      input:
        m8_file = RunNRAlignment.out_m8,
        lineage_db = lineage_db,
        taxon_blacklist = taxon_blacklist,
        deuterostome_db = deuterostome_db,
        accession2taxid = accession2taxid_db,
        min_read_length = min_read_length,
        count_type = library_type,
        docker_image_id = docker_image_id,
    }

    call TallyHitsNT {
      input:
        reads_fastq = RunSubsampling.subsampled_fastq,
        m8 = RunCallHitsNT.deduped_out_m8,
        hitsummary = RunCallHitsNT.hitsummary,
        reads_to_contigs_tsv = RunReadsToContigs.reads_to_contigs_tsv,
        docker_image_id = docker_image_id,
    }

    call UnmappedReads {
      input:
        input_file = RunSubsampling.subsampled_fastq,
        hitsummary_nt = RunCallHitsNT.hitsummary,
        hitsummary_nr = RunCallHitsNR.hitsummary,
        reads_to_contigs_tsv = RunReadsToContigs.reads_to_contigs_tsv,
        docker_image_id = docker_image_id
    }

    call TallyHitsNR {
      input:
        reads_fastq = RunSubsampling.subsampled_fastq,
        m8 = RunCallHitsNR.deduped_out_m8,
        hitsummary = RunCallHitsNR.hitsummary,
        reads_to_contigs_tsv = RunReadsToContigs.reads_to_contigs_tsv,
        docker_image_id = docker_image_id,
    }

    call ReassignM8NT {
      input:
        m8 = RunCallHitsNT.deduped_out_m8,
        reads_to_contigs_tsv = RunReadsToContigs.reads_to_contigs_tsv,
        docker_image_id = docker_image_id,
    }

    call ReassignM8NR {
      input:
        m8 = RunCallHitsNR.deduped_out_m8,
        reads_to_contigs_tsv = RunReadsToContigs.reads_to_contigs_tsv,
        docker_image_id = docker_image_id,
    }

    call GenerateAnnotatedFasta {
      input:
        pre_alignment_fasta = RunSubsampling.subsampled_fastq,
        nt_m8 = ReassignM8NT.m8_reassigned,
        nr_m8 = ReassignM8NR.m8_reassigned,
        docker_image_id = docker_image_id,
    }

    call GenerateTaxidFasta {
      input:
        annotated_merged_fa = GenerateAnnotatedFasta.assembly_refined_annotated_merged_fa,
        nt_hitsummary_tab = RunCallHitsNT.hitsummary,
        nr_hitsummary_tab = RunCallHitsNR.hitsummary,
        lineage_db = lineage_db,
        docker_image_id = docker_image_id,
    }

    call GenerateTaxidLocator {
        input:
            assembly_refined_taxid_annot_fasta = GenerateTaxidFasta.assembly_refined_taxid_annot_fasta,
            docker_image_id = docker_image_id
    }

    output {
        File nt_deduped_out_m8 = RunCallHitsNT.deduped_out_m8
        File nt_hitsummary = RunCallHitsNT.hitsummary
        File nt_counts_json = RunCallHitsNT.counts_json
        File nr_deduped_out_m8 = RunCallHitsNR.deduped_out_m8
        File nr_hitsummary = RunCallHitsNR.hitsummary
        File nr_counts_json = RunCallHitsNR.counts_json
        File nt_tallied_hits = TallyHitsNT.tallied_hits
        File nr_tallied_hits = TallyHitsNR.tallied_hits
        File unmapped_reads = UnmappedReads.unmapped_reads
        File refined_taxid_fasta_out_assembly_refined_taxid_annot_fasta = GenerateTaxidFasta.assembly_refined_taxid_annot_fasta
        File? refined_taxid_fasta_out_count = GenerateTaxidFasta.output_read_count
        File refined_taxid_locator_out_assembly_refined_taxid_annot_sorted_nt_fasta = GenerateTaxidLocator.assembly_refined_taxid_annot_sorted_nt_fasta
        File refined_taxid_locator_out_assembly_refined_taxid_locations_nt_json = GenerateTaxidLocator.assembly_refined_taxid_locations_nt_json
        File refined_taxid_locator_out_assembly_refined_taxid_annot_sorted_nr_fasta = GenerateTaxidLocator.assembly_refined_taxid_annot_sorted_nr_fasta
        File refined_taxid_locator_out_assembly_refined_taxid_locations_nr_json = GenerateTaxidLocator.assembly_refined_taxid_locations_nr_json
        File refined_taxid_locator_out_assembly_refined_taxid_annot_sorted_genus_nt_fasta = GenerateTaxidLocator.assembly_refined_taxid_annot_sorted_genus_nt_fasta
        File refined_taxid_locator_out_assembly_refined_taxid_locations_genus_nt_json = GenerateTaxidLocator.assembly_refined_taxid_locations_genus_nt_json
        File refined_taxid_locator_out_assembly_refined_taxid_annot_sorted_genus_nr_fasta = GenerateTaxidLocator.assembly_refined_taxid_annot_sorted_genus_nr_fasta
        File refined_taxid_locator_out_assembly_refined_taxid_locations_genus_nr_json = GenerateTaxidLocator.assembly_refined_taxid_locations_genus_nr_json
        File refined_taxid_locator_out_assembly_refined_taxid_annot_sorted_family_nt_fasta = GenerateTaxidLocator.assembly_refined_taxid_annot_sorted_family_nt_fasta
        File refined_taxid_locator_out_assembly_refined_taxid_locations_family_nt_json = GenerateTaxidLocator.assembly_refined_taxid_locations_family_nt_json
        File refined_taxid_locator_out_assembly_refined_taxid_annot_sorted_family_nr_fasta = GenerateTaxidLocator.assembly_refined_taxid_annot_sorted_family_nr_fasta
        File refined_taxid_locator_out_assembly_refined_taxid_locations_family_nr_json = GenerateTaxidLocator.assembly_refined_taxid_locations_family_nr_json
        File refined_taxid_locator_out_assembly_refined_taxid_locations_combined_json = GenerateTaxidLocator.assembly_refined_taxid_locations_combined_json
        File? refined_taxid_locator_out_count = GenerateTaxidLocator.output_read_count
    }
}
