# IDseq Consensus Genome workflow
# Based on original work at:
# - CZ Biohub SARS-CoV-2 pipeline, https://github.com/czbiohub/sc2-illumina-pipeline
# - ARTIC Oxford Nanopore MinION SARS-CoV-2 SOP, https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html
# With enhancements and additional modules by the CZI Infectious Disease team

version 1.1

workflow consensus_genome {
    input {
        # Required parameters
        File fastqs_0
        File? fastqs_1

        Int max_reads = 75000000

        String docker_image_id
        File ercc_fasta = "s3://idseq-public-references/consensus-genome/ercc_sequences.fasta"
        File kraken2_db_tar_gz  # TODO: make this optional; only required if filter_reads == true, even for Illumina
        File primer_bed = "s3://idseq-public-references/consensus-genome/artic_v3_primers.bed" # Only required for Illumina

        File? ref_fasta # Only required for Illumina (ONT SC2 reference is built into ARTIC); takes precedence over ref_accession_id
        String? ref_accession_id # Only required for Illumina; has no effect if ref_fasta is set

        File ref_host
        String technology # Input sequencing technology ("Illumina" or "ONT"); ONT only works with SC2 samples (SC2 reference is built into ARTIC)

        # Sample name: include in tags and files
        String sample

        # Optional prefix to add to the output filenames
        String prefix = ""

        # ONT-specific inputs
        File primer_schemes = "s3://idseq-public-references/consensus-genome/artic-primer-schemes_v2.tar.gz"
        String primer_set = "nCoV-2019/V3"
        # filters in accordance with recommended parameters in ARTIC SARS-CoV-2 bioinformatics protocol are...
        # ...intended to remove obviously chimeric reads.
        Boolean apply_length_filter = true # Set to False for Clear Labs samples

        # set default min_length to 350 unless midnight primers are used
        Int min_length = if primer_set == "nCoV-2019/V1200" then 250 else 350
        # set default max_length to 1500 unless midnight primers are used
        Int max_length = if primer_set == "nCoV-2019/V1200" then 1500 else 700
        # normalise: default is set to 1000 to avoid spurious indels observed in validation
        Int normalise  = 1000
        # medaka_model: default is selected to support current Clear Labs workflow
        String medaka_model = "r941_min_high_g360"
        String vadr_options = "-s -r --nomisc --mkey sarscov2 --lowsim5term 2 --lowsim3term 2 --fstlowthr 0.0 --alt_fail lowscore,fsthicnf,fstlocnf --noseqnamemax"
        File vadr_model = "s3://idseq-public-references/consensus-genome/vadr-models-sarscov2-1.2-2.tar.gz"

        # Illumina-specific parameters
        # Step parameters
        Boolean filter_reads = true
        Boolean trim_adapters = true

        Float ivarFreqThreshold = 0.75
        Int   ivarQualThreshold  = 20
        Int   minDepth          = if "~{primer_bed}" == "s3://idseq-public-references/consensus-genome/na_primers.bed" then 5 else 10

        # If no_reads_quast is true, quast runs without considering the raw reads (only considering the reference genome and the consensus.fa).
        # This reduces the number of informative metrics that quast provides, but speeds up the step since quast is faster when it doesn't consider raw reads.
        # (Not expected to be used in idseq production)
        String no_reads_quast = false

        # assumes about 20 mutations between 2 random samples
        # (this is an overestimate to increase sensitivity)
        Float bcftoolsCallTheta = 0.0006

        # Dummy values - required by SFN interface
        String s3_wd_uri = ""
    }

    call ValidateInput {
        input:
            prefix = prefix,
            fastqs = select_all([fastqs_0, fastqs_1]),
            technology = technology,
            max_reads = max_reads,
            docker_image_id = docker_image_id
    }

    if (ref_accession_id != None) {
        call FetchSequenceByAccessionId {
            input: accession_id = select_first([ref_accession_id]),
            docker_image_id = docker_image_id
        }
    }

    if (technology == "ONT" && apply_length_filter) {
        call ApplyLengthFilter {
            input:
                prefix = prefix,
                fastqs = ValidateInput.validated_fastqs,
                min_length = min_length,
                max_length = max_length,
                docker_image_id = docker_image_id
        }    
    }

    call RemoveHost {
        input:
            prefix = prefix,
            fastqs = select_first([ApplyLengthFilter.filtered_fastqs, ValidateInput.validated_fastqs]),
            ref_host = ref_host,
            technology = technology,
            docker_image_id = docker_image_id
    }

    if (technology == "Illumina") {
        call QuantifyERCCs {
            input:
                prefix = prefix,
                fastqs = RemoveHost.host_removed_fastqs,
                ercc_fasta = ercc_fasta,
                docker_image_id = docker_image_id
        }    
        if (filter_reads) {
            call FilterReads {
                input:
                    prefix = prefix,
                    fastqs = RemoveHost.host_removed_fastqs,
                    ref_fasta = select_first([ref_fasta, FetchSequenceByAccessionId.sequence_fa]),
                    kraken2_db_tar_gz = kraken2_db_tar_gz,
                    docker_image_id = docker_image_id
            }
        }
        if (trim_adapters) {
            call TrimReads {
                input:
                    fastqs = select_first([FilterReads.filtered_fastqs, RemoveHost.host_removed_fastqs]),
                    docker_image_id = docker_image_id
            }
        }
        call AlignReads {
            input:
                prefix = prefix,
                sample = sample,
                # use trimReads output if we ran it; otherwise fall back to FilterReads output if we ran it; 
                # otherwise fall back to RemoveHost output
                fastqs = select_first([TrimReads.trimmed_fastqs, FilterReads.filtered_fastqs, RemoveHost.host_removed_fastqs]),
                ref_fasta = select_first([ref_fasta, FetchSequenceByAccessionId.sequence_fa]),
                docker_image_id = docker_image_id
        }
        call TrimPrimers {
            input:
                prefix = prefix,
                alignments = AlignReads.alignments,
                primer_bed = primer_bed,
                docker_image_id = docker_image_id
        }
        call MakeConsensus {
            input:
                prefix = prefix,
                sample = sample,
                bam = TrimPrimers.trimmed_bam_ch,
                ivarFreqThreshold = ivarFreqThreshold,
                minDepth = minDepth,
                ivarQualThreshold = ivarQualThreshold,
                docker_image_id = docker_image_id
        }
        # this step does not rely on outputs of QUAST, so we can move it here to avoid complex logic
        call CallVariants {
            input:
                prefix = prefix,
                call_variants_bam = TrimPrimers.trimmed_bam_ch,
                ref_fasta = select_first([ref_fasta, FetchSequenceByAccessionId.sequence_fa]),
                bcftoolsCallTheta = bcftoolsCallTheta,
                ivarQualThreshold = ivarQualThreshold,
                minDepth = minDepth,
                docker_image_id = docker_image_id
        }
        call RealignConsensus {
            input:
                prefix = prefix,
                sample = sample,
                ref_fasta = select_first([ref_fasta, FetchSequenceByAccessionId.sequence_fa]),
                consensus = MakeConsensus.consensus_fa,
                docker_image_id = docker_image_id
        }
    }

    if (technology == "ONT"){
        call RunMinion {
            input:
                prefix = prefix,
                sample = sample,
                fastqs = RemoveHost.host_removed_fastqs,
                primer_schemes = primer_schemes,
                normalise = normalise,
                medaka_model = medaka_model,
                primer_set = primer_set,
                docker_image_id = docker_image_id
        }
    }

    call Quast {
        input:
            prefix = prefix,
            assembly = select_first([MakeConsensus.consensus_fa, RunMinion.consensus_fa]),
            bam = select_first([TrimPrimers.trimmed_bam_ch, RunMinion.primertrimmedbam]),
            # use trimReads output if we ran it; otherwise fall back to FilterReads output, or RemoveHost for ONT
            fastqs = select_first([TrimReads.trimmed_fastqs, FilterReads.filtered_fastqs, RemoveHost.host_removed_fastqs]),
            ref_fasta = select_first([ref_fasta, FetchSequenceByAccessionId.sequence_fa]), # FIXME: (AK) primer_schemes/nCoV-2019.reference.fasta
            no_reads_quast = no_reads_quast,
            technology = technology,
            primer_schemes = primer_schemes, # Only required for ONT; contains reference genome for ARTIC
            primer_set = primer_set,
            docker_image_id = docker_image_id
    }

    call ComputeStats {
        input:
            prefix = prefix,
            sample = sample,
            cleaned_bam = select_first([TrimPrimers.trimmed_bam_ch, RunMinion.primertrimmedbam]),
            assembly = select_first([MakeConsensus.consensus_fa, RunMinion.consensus_fa]),
            ercc_stats = QuantifyERCCs.ercc_out,       # does not exist - NO ERCC results for ONT, this argument must be optional
            vcf = select_first([CallVariants.variants_ch, RunMinion.vcf_pass]),
            fastqs = RemoveHost.host_removed_fastqs, # select_all([fastqs_0, fastqs_1]), # FIXME: (AK) verify the correct value for this
            ref_host = ref_host,                       # FIXME: (AK) primer_schemes/nCoV-2019.reference.fasta?
            technology = technology,
            docker_image_id = docker_image_id
    }

    # TODO: generalize VADR to run on any coronavirus reference or any viral reference with a VADR model available
    if (ref_accession_id == None) {
        call Vadr {
            input:
                prefix = prefix,
                assembly = select_first([MakeConsensus.consensus_fa, RunMinion.consensus_fa]),
                vadr_options = vadr_options,
                vadr_model = vadr_model,
                docker_image_id = docker_image_id
        }
    }

    call ZipOutputs {
        input:
            prefix = prefix,
            outputFiles = select_all(flatten([
                RemoveHost.host_removed_fastqs,
                select_all([
                    MakeConsensus.consensus_fa,
                    RunMinion.consensus_fa,
                    ComputeStats.depths_fig,
                    TrimPrimers.trimmed_bam_ch,
                    RunMinion.primertrimmedbam,
                    TrimPrimers.trimmed_bam_bai,
                    RunMinion.primertrimmedbai,
                    Quast.quast_txt,
                    Quast.quast_tsv,
                    AlignReads.alignments,
                    RunMinion.alignedbam,
                    QuantifyERCCs.ercc_out,            # No ERCC results for ONT
                    QuantifyERCCs.ercc_out,            # No ERCC results for ONT
                    ComputeStats.output_stats,
                    ComputeStats.sam_depths,
                    CallVariants.variants_ch,
                    RunMinion.vcf,
                    RealignConsensus.muscle_output,
                    RunMinion.muscle_output,
                    Vadr.vadr_quality,                 # Optional (VADR only runs on default (coronavirus) reference)
                    Vadr.vadr_alerts,                  # Optional (VADR only runs on default (coronavirus) reference)
                    Vadr.vadr_errors                   # Optional (only present if VADR ran and exited with an error)
                ])
            ])),
            docker_image_id = docker_image_id
    }

    output {
        Array[File] remove_host_out_host_removed_fastqs = RemoveHost.host_removed_fastqs 
        File? quantify_erccs_out_ercc_out = QuantifyERCCs.ercc_out # does not exist for ONT
        Array[File]+? filter_reads_out_filtered_fastqs = FilterReads.filtered_fastqs # does not exist for ONT
        Array[File]+? trim_reads_out_trimmed_fastqs = TrimReads.trimmed_fastqs # does not exist for ONT
        File? align_reads_out_alignments = select_first([AlignReads.alignments, RunMinion.alignedbam])
        File? trim_primers_out_trimmed_bam_ch = select_first([TrimPrimers.trimmed_bam_ch, RunMinion.primertrimmedbam])
        File? trim_primers_out_trimmed_bam_bai = select_first([TrimPrimers.trimmed_bam_bai, RunMinion.primertrimmedbai])
        File? make_consensus_out_consensus_fa = select_first([MakeConsensus.consensus_fa, RunMinion.consensus_fa])
        File? realign_consensus_fa = select_first([RealignConsensus.muscle_output, RunMinion.muscle_output])
        File? quast_out_quast_txt = Quast.quast_txt
        File? quast_out_quast_tsv = Quast.quast_tsv
        File? call_variants_out_variants_ch = select_first([CallVariants.variants_ch, RunMinion.vcf_pass])
        File? compute_stats_out_depths_fig = ComputeStats.depths_fig
        File? compute_stats_out_output_stats = ComputeStats.output_stats
        File? compute_stats_out_sam_depths = ComputeStats.sam_depths
        File? vadr_quality_out = Vadr.vadr_quality  # Optional (VADR only runs on default (coronavirus) reference)
        File? vadr_alerts_out = Vadr.vadr_alerts    # Optional (VADR only runs on default (coronavirus) reference)
        File? vadr_errors = Vadr.vadr_errors        # Optional (only present if VADR ran and exited with an error)
        File? minion_log = RunMinion.log
        File zip_outputs_out_output_zip = ZipOutputs.output_zip
    }
}

task ValidateInput{
    input {
        String prefix
        Array[File]+ fastqs
        String technology
        Int max_reads

        String docker_image_id
    }

    command <<<
        set -uxo pipefail 
        function raise_error {
            set +x
            export error=$1 cause=$2
            jq -nc ".wdl_error_message=true | .error=env.error | .cause=env.cause" > /dev/stderr
            exit 1 
        }
        if [[ "~{technology}" == "ONT" ]] && [[ "~{length(fastqs)}" -gt 1 ]]; then
            # ONT pipeline should only have one input
            raise_error InvalidInputFileError "Multiple fastqs provided for ONT"
        fi 

        counter=1
        for fastq in ~{sep=' ' fastqs}; do 
            # limit max # of reads to max_reads
            seqkit head -n "~{max_reads}" $fastq -o "~{prefix}validated_$counter.fastq.gz" 2> read_error.txt
            if [[ -s read_error.txt ]]; then 
                # Checks if seqkit can parse input files
                raise_error InvalidFileFormatError "Error parsing one of the input files: ""$(cat read_error.txt)"
            fi 
           ((counter++))
        done
        set -e
        seqkit stats "~{prefix}"validated*fastq.gz -T > input_stats.tsv
        if grep  -q "FASTA" <<< $(cut -f 2 input_stats.tsv ); then 
            # Input files cannot be in FASTA format
            raise_error InvalidInputFileError "One or more of the input files is in FASTA format"
        fi 
        if [[ "~{technology}" == "Illumina" ]]; then 
            # check if any of the input files has max length > 500bp
            MAXLEN=$(cut -f 8 input_stats.tsv | tail -n "~{length(fastqs)}" | sort -n | tail -n 1)
            if [[ $MAXLEN > 500 ]]; then 
                raise_error InvalidInputFileError "Read longer than 500bp for Illumina"
            fi 
        fi
    >>>

    output {
        Array[File]+ validated_fastqs = glob("~{prefix}validated*.fastq.gz")
        File? input_stats = "input_stats.tsv"
    }

    runtime {
        docker: docker_image_id
    }
}

task FetchSequenceByAccessionId {
    input {
        String accession_id
        String docker_image_id
    }

    command <<<
        function incrementAccession { 
            ( [[ $1 =~ ([A-Z0-9_]*)\.([0-9]+) ]] && echo "${BASH_REMATCH[1]}.$(( ${BASH_REMATCH[2]} + 1 ))");
        }
        # Try fetching accession id. If not found, try incrementing the version. 
        ({ taxoniq get-from-s3 --accession-id "~{accession_id}"; } || \
        { taxoniq get-from-s3 --accession-id $(incrementAccession "~{accession_id}"); } \
        || exit 4; ) > sequence.fa;
        if [[ $? == 4 ]]; then
            export error=AccessionIdNotFound cause="Accession ID ~{accession_id} not found in the index"
            jq -nc ".wdl_error_message=true | .error=env.error | .cause=env.cause" > /dev/stderr
            exit 4
        fi
        exit $?
    >>>

    output {
        File sequence_fa = "sequence.fa"
    }

    runtime {
        docker: docker_image_id
    }
}

task ApplyLengthFilter {
    input {
        String prefix
        Array[File]+ fastqs
        Int min_length
        Int max_length
        String docker_image_id
    }

    command <<<
        artic guppyplex --min-length ~{min_length} --max-length ~{max_length} --directory $(dirname "~{fastqs[0]}") --output filtered.fastq
    >>>

    output {
        Array[File]+ filtered_fastqs = ["filtered.fastq"]
    }

    runtime {
        docker: docker_image_id
    }
}

task RemoveHost {
    input {
        String prefix
        Array[File]+ fastqs  
        File ref_host
        String technology
        String docker_image_id
    }

    command <<<
        set -euxo pipefail

        export CORES=`nproc --all`
        if [[ "~{length(fastqs)}" == 1 ]]; then
            if [[ "~{technology}" == "Illumina" ]]; then
                minimap2 -t $CORES -ax sr ~{ref_host} ~{fastqs[0]} | \
                samtools view --no-PG -@ $CORES -b -f 4 | \
                samtools fastq -@ $CORES -0 "~{prefix}no_host_1.fq.gz" -n -c 6 -
            else # if technology == ONT
                minimap2 -t $CORES -ax map-ont ~{ref_host} ~{fastqs[0]} | \
                samtools view --no-PG -@ $CORES -b -f 4 | \
                samtools fastq -@ $CORES -0 "~{prefix}no_host_1.fq.gz" -n -c 6 -
            fi
        else
            minimap2 -t $CORES -ax sr ~{ref_host} ~{sep=' ' fastqs} | \
            samtools view --no-PG -@ $CORES -b -f 4 | \
            samtools fastq -@ $CORES -1 "~{prefix}no_host_1.fq.gz" -2 "~{prefix}no_host_2.fq.gz" -0 /dev/null -s /dev/null -n -c 6 -
        fi

        if [[ -z $(gzip -cd "~{prefix}no_host_1.fq.gz" | head -c1) ]]; then
            set +x
            export error=InsufficientReadsError cause="No reads after RemoveHost"
            jq -nc ".wdl_error_message=true | .error=env.error | .cause=env.cause" > /dev/stderr
            exit 1
        fi
    >>>

    output {
        Array[File] host_removed_fastqs = glob("~{prefix}no_host_*.fq.gz")
    }

    runtime {
        docker: docker_image_id
    }
}

task QuantifyERCCs {
    input {
        String prefix
        Array[File]+ fastqs
        File ercc_fasta

        String docker_image_id
    }

    command <<<
        set -euxo pipefail

        minimap2 -ax sr ~{ercc_fasta} ~{sep=' ' fastqs} | samtools view --no-PG -bo ercc_mapped.bam
        samtools stats ercc_mapped.bam > "~{prefix}ercc_stats.txt"
    >>>

    output {
        File ercc_out = "~{prefix}ercc_stats.txt"
    }

    runtime {
        docker: docker_image_id
    }
}

task FilterReads {
    input {
        String prefix
        # SARS-CoV-2 default
        String taxid = "2697049"
        Array[File]+ fastqs

        File ref_fasta
        File kraken2_db_tar_gz

        String docker_image_id
    }

    command <<<
        _no_reads_error() {
            set +x
            export error=InsufficientReadsError cause="No reads after FilterReads"
            jq -nc ".wdl_error_message=true | .error=env.error | .cause=env.cause" > /dev/stderr
            exit 1
        }

        set -euxo pipefail

        export TMPDIR=${TMPDIR:-/tmp}
        export CORES=`nproc --all`

        minimap2 -ax sr -t $CORES "~{ref_fasta}" ~{sep=' ' fastqs} \
            | samtools sort -@ $CORES -n -O bam -o "${TMPDIR}/mapped.bam"

        if [[ "~{length(fastqs)}" == 1 ]]; then
            samtools fastq -@ $CORES -G 12 -0 "${TMPDIR}/paired1.fq.gz" -n -c 6 "${TMPDIR}/mapped.bam"
        else
            samtools fastq -@ $CORES -G 12 -1 "${TMPDIR}/paired1.fq.gz" -2 "${TMPDIR}/paired2.fq.gz" \
                -0 /dev/null -s /dev/null -n -c 6 "${TMPDIR}/mapped.bam"
        fi

        paired1size=$(stat --printf="%s" "${TMPDIR}/paired1.fq.gz")
        if (( paired1size > 28 )); then
            if [[ "~{length(fastqs)}" == 1 ]]; then
                KRAKEN_ARGS="${TMPDIR}/paired1.fq.gz"
                KRAKEN_OUTPUT_ARG="${TMPDIR}/~{prefix}classified_1.fq"
            else
                KRAKEN_ARGS="--paired ${TMPDIR}/paired1.fq.gz ${TMPDIR}/paired2.fq.gz"
                KRAKEN_OUTPUT_ARG="${TMPDIR}/~{prefix}classified#.fq"
            fi

            mkdir "${TMPDIR}/kraken_db"
            tar -xv --use-compress-program=pigz -C "${TMPDIR}/kraken_db" -f "~{kraken2_db_tar_gz}"
            kraken2 --db ${TMPDIR}/kraken_db/* \
                --threads $CORES \
                --report ${TMPDIR}/~{prefix}kraken2_report.txt \
                --classified-out $KRAKEN_OUTPUT_ARG \
                --output - \
                --memory-mapping --gzip-compressed \
                $KRAKEN_ARGS

            grep --no-group-separator -A3 "kraken:taxid|~{taxid}" \
                "${TMPDIR}/~{prefix}classified_1.fq" \
                > "${TMPDIR}/~{prefix}filtered_1.fq" || _no_reads_error
            [[ "~{length(fastqs)}" == 1 ]] || grep --no-group-separator -A3 "kraken:taxid|~{taxid}" \
                "${TMPDIR}/~{prefix}classified_2.fq" \
                > "${TMPDIR}/~{prefix}filtered_2.fq" || _no_reads_error
            bgzip -@ $CORES -c "${TMPDIR}/~{prefix}filtered_1.fq" > "~{prefix}filtered_1.fq.gz"
            [[ "~{length(fastqs)}" == 1 ]] || bgzip -@ $CORES -c "${TMPDIR}/~{prefix}filtered_2.fq" > "~{prefix}filtered_2.fq.gz"
        else
            mv "${TMPDIR}/paired1.fq.gz" "~{prefix}filtered_1.fq.gz"
            [[ "~{length(fastqs)}" == 1 ]] || mv "${TMPDIR}/paired2.fq.gz" "~{prefix}filtered_2.fq.gz"
        fi

        if [[ -z $(gzip -cd "~{prefix}filtered_1.fq.gz" | head -c1) ]]; then
            _no_reads_error
        fi
    >>>

    output {
        Array[File]+ filtered_fastqs = glob("~{prefix}filtered_*.fq.gz")
    }

    runtime {
        docker: docker_image_id
    }
}

task TrimReads {
    input {
        Array[File]+ fastqs

        String docker_image_id
    }

    command <<<
        set -euxo pipefail

        BASENAME="trim_reads"
        if [[ "~{length(fastqs)}" == 1 ]]; then
            trim_galore --gzip --fastqc "~{fastqs[0]}" --basename $BASENAME
            mv "${BASENAME}_trimmed.fq.gz" "${BASENAME}_val_1.fq.gz"
        else
            trim_galore --gzip --fastqc --paired --basename $BASENAME ~{sep=' ' fastqs}
        fi

        if [[ -z $(gzip -cd "${BASENAME}_val_1.fq.gz" | head -c1) ]]; then
            set +x
            export error=InsufficientReadsError cause="No reads after TrimReads"
            jq -nc ".wdl_error_message=true | .error=env.error | .cause=env.cause" > /dev/stderr
            exit 1
        fi
    >>>

    output {
        Array[File]+ trimmed_fastqs = glob("*_val_[12].fq.gz")
    }

    runtime {
        docker: docker_image_id
    }
}

task AlignReads {
    # TODO: process errors: no reads left (unlikely)
    input {
        String prefix
        String sample
        Array[File]+ fastqs
        File ref_fasta

        String docker_image_id
    }

    command <<<
        set -euxo pipefail

        export CORES=`nproc --all`

        # Sample id included in the bam files
        minimap2 -ax sr -t $CORES -R '@RG\tID:~{sample}\tSM:~{sample}' "~{ref_fasta}" ~{sep=' ' fastqs} \
            | samtools sort -@ $CORES -O bam -o "~{prefix}aligned_reads.bam"
    >>>

    output {
        File alignments = "~{prefix}aligned_reads.bam"
    }

    runtime {
        docker: docker_image_id
    }
}

task TrimPrimers {
    input {
        String prefix
        File alignments
        File primer_bed

        String docker_image_id

        Int samQualThreshold = 20
    }

    command <<<
        set -euxo pipefail
        samtools view -F4 -q "~{samQualThreshold}" -o ivar.bam "~{alignments}"
        samtools index ivar.bam
        # The SNAP protocol may result in primer position offsets due to polymerases adding additional bases 
        # (https://github.com/andersen-lab/ivar/pull/88)... if the primer bed file given is for the SNAP protocol, 
        # then add offset allowance of 5 bp in accordance with ivar developer recommendation
        if [[ "$(basename '~{primer_bed}')" == "snap_primers.bed" ]]; then
            primerOffset=5
        elif [[ "$(basename '~{primer_bed}')" == "artic_v3_short_275_primers.bed" ]]; then
            primerOffset=2
        else
            primerOffset=0
        fi
        ivar trim -x $primerOffset -e -i ivar.bam -b "~{primer_bed}" -p ivar.out
        samtools sort -O bam -o "~{prefix}primertrimmed.bam" ivar.out.bam
        samtools index "~{prefix}primertrimmed.bam"
    >>>

    output {
        File trimmed_bam_ch = "~{prefix}primertrimmed.bam"
        File trimmed_bam_bai = "~{prefix}primertrimmed.bam.bai"
    }

    runtime {
        docker: docker_image_id
    }

}

task MakeConsensus {
    input {
        String prefix
        String sample
        File bam

        Float ivarFreqThreshold
        Int minDepth
        Int ivarQualThreshold

        String docker_image_id
    }

    command <<<
        set -euxo pipefail

        samtools index "~{bam}"
        samtools mpileup -A -d 0 -Q0 "~{bam}" | ivar consensus -q "~{ivarQualThreshold}" -t "~{ivarFreqThreshold}" -m "~{minDepth}" -n N -p "~{prefix}primertrimmed.consensus"
        echo ">""~{sample}" > "~{prefix}consensus.fa"
        seqtk seq -l 50 "~{prefix}primertrimmed.consensus.fa" | tail -n +2 >> "~{prefix}consensus.fa"

        # One-line file means just the fasta header with no reads
        if [[ $(wc -l "~{prefix}consensus.fa" | cut -d' ' -f1) == 1 ]]; then
            set +x
            export error=InsufficientReadsError cause="No reads after MakeConsensus"
            jq -nc ".wdl_error_message=true | .error=env.error | .cause=env.cause" > /dev/stderr
            exit 1
        fi
    >>>

    output {
        File consensus_fa = "~{prefix}consensus.fa"
    }

    runtime {
        docker: docker_image_id
    }
}

task CallVariants {
    input {
        String prefix
        File call_variants_bam  # same as primertrimmed_bam produced by trimPrimers
        File ref_fasta
        Int ivarQualThreshold
        Float bcftoolsCallTheta
        Int minDepth

        String docker_image_id
    }

    command <<<
        set -euxo pipefail

        # NOTE: we use samtools instead of bcftools mpileup because bcftools 1.9 ignores -d0
        samtools mpileup -aa -u -Q "~{ivarQualThreshold}" -d 100000000 -L 100000000 -t AD -f "~{ref_fasta}" "~{call_variants_bam}" | bcftools call --ploidy 1 -m -P "~{bcftoolsCallTheta}" -v - | bcftools view -i 'DP>=~{minDepth}' > "~{prefix}variants.vcf"
        bgzip "~{prefix}variants.vcf"
        tabix "~{prefix}variants.vcf.gz"
        bcftools stats "~{prefix}variants.vcf.gz" > "~{prefix}bcftools_stats.txt"
    >>>

    output {
        File variants_ch = "~{prefix}variants.vcf.gz"
        File bcftools_stats_ch = "~{prefix}bcftools_stats.txt"
    }

    runtime {
        docker: docker_image_id
    }
}

task RealignConsensus {
    input {
        String prefix
        String sample
        File ref_fasta
        File consensus
        String docker_image_id
    }

    command <<<
        # MUSCLE accepts a fasta file containing all the sequences that are to be aligned (in this case we want
        # to align the reference and the consensus genome) and outputs a multiple sequence alignment file. 
        # Some documentation here: https://www.drive5.com/muscle/manual/basic_alignment.html
        cat "~{consensus}" "~{ref_fasta}" > "~{sample}.muscle.in.fasta" 
        muscle -in "~{sample}.muscle.in.fasta" -out "~{sample}.muscle.out.fasta"

    >>>

    output{
        File muscle_output = "~{sample}.muscle.out.fasta"
    }

    runtime {
        docker: docker_image_id
    }
}

task RunMinion {
    input {
        String prefix
        String sample
        Array[File]+ fastqs
        File primer_schemes
        String primer_set
        Int normalise
        String medaka_model
        String docker_image_id
    }

    command <<<
        set -euxo pipefail
        export CORES=`nproc --all`

        tar -xzf "~{primer_schemes}"

        # TODO: upgrade to artic 1.3.0 when released (https://github.com/artic-network/fieldbioinformatics/pull/70)
        artic minion --medaka --no-longshot --normalise "~{normalise}" --threads 4 --scheme-directory primer_schemes --read-file ~{sep=' ' fastqs} --medaka-model "~{medaka_model}" "~{primer_set}" "~{sample}"
        # the .bam file doesn't seem to be sorted when it comes out, so explicitly sorting it here because a
        # ...sorted .bam is necessary for ComputeStats step downstream
        samtools sort "~{sample}.primertrimmed.rg.sorted.bam" > "~{sample}.primertrimmed.rg.resorted.bam"
        mv "~{sample}.primertrimmed.rg.resorted.bam" "~{sample}.primertrimmed.rg.sorted.bam"
        samtools index "~{sample}.primertrimmed.rg.sorted.bam" # to create "~{sample}.primertrimmed.rg.sorted.bai"
        gunzip "~{sample}.pass.vcf.gz"
    >>>

    output {
        File primertrimmedbam = "~{sample}.primertrimmed.rg.sorted.bam"
        File primertrimmedbai = "~{sample}.primertrimmed.rg.sorted.bam.bai"
        File alignedbam = "~{sample}.sorted.bam"
        File vcf_pass = "~{sample}.pass.vcf"
        File vcf = "~{sample}.merged.vcf"
        File consensus_fa = "~{sample}.consensus.fasta"
        File muscle_output = "~{sample}.muscle.out.fasta"
        File log = "~{sample}.minion.log.txt"
    }

    runtime {
        docker: docker_image_id
    }
}


task Quast {
    input {
        String prefix
        File assembly   # same as consensus_fa
        File bam
        Array[File]+ fastqs
        File? ref_fasta
        File? primer_schemes
        String? primer_set
        String no_reads_quast
        String technology
        Int threads = 4

        String docker_image_id
    }

    command <<<
        set -euxo pipefail

        export CORES=`nproc --all`
        a=`cat "~{assembly}" | wc -l`
        cp "~{bam}" .
        export BAM=$(basename "~{bam}")
        cp "~{assembly}" .
        export ASSEMBLY=$(basename "~{assembly}")

        if [[ $a -ne 0 ]]; then
            if [[ ~{no_reads_quast} = true ]]; then
                quast.py --min-contig 0 -o quast -r "~{ref_fasta}" -t $CORES --ref-bam "~{bam}" "~{assembly}"
            else
                if [[ "~{technology}" == "Illumina" ]]; then
                    if [[ "~{length(fastqs)}" == 1 ]]; then
                        quast.py --min-contig 0 -o quast -r "~{ref_fasta}" -t $CORES --ref-bam "$BAM" "$ASSEMBLY" --single "~{fastqs[0]}"
                    else
                        quast.py --min-contig 0 -o quast -r "~{ref_fasta}" -t $CORES --ref-bam "$BAM" "$ASSEMBLY" -1 ~{sep=' -2 ' fastqs}
                    fi
                else # technology == "ONT"
                    # Currently, the ONT branch only supports the ARTIC SARS-CoV-2 SOP, which bundles its own reference genome.
                    # The ref_fasta parameter is ignored and the bundled genome reference from ARTIC primer_schemes is used instead.
                    tar -xzf ~{primer_schemes}
                    quast.py --min-contig 0 -o quast -r "primer_schemes/~{primer_set}/nCoV-2019.reference.fasta" -t $CORES --ref-bam "$BAM" "$ASSEMBLY" --nanopore "~{fastqs[0]}"
                fi
            fi
        else
            mkdir quast
            echo "quast folder is empty" > "quast/report.txt"
        fi
    >>>

    output {
        Array[File] quast_dir = glob("quast/*")
        File quast_txt = "quast/report.txt"
        File? quast_tsv = "quast/report.tsv"
    }

    runtime {
        docker: docker_image_id
    }
}

task ComputeStats {
    input {
        String prefix
        String sample
        File cleaned_bam
        File assembly
        File? ercc_stats  # optional, will only exist for Illumina runs
        File vcf

        Array[File]+ fastqs
        File ref_host

        String technology

        String docker_image_id
    }

    command <<<
        set -euxo pipefail

        samtools index "~{cleaned_bam}"
        samtools stats "~{cleaned_bam}" > "~{prefix}samtools_stats.txt"
        samtools depth -aa -d 0 "~{cleaned_bam}" | awk '{print $3}' > "~{prefix}samtools_depth.txt"

        python3 <<CODE
        import argparse
        import collections
        import gzip
        import json
        import re
        import pysam
        from Bio import SeqIO
        import numpy as np
        from matplotlib import pyplot as plt
        import seaborn as sns

        stats = {"sample_name": "~{sample}"}

        depths = open("~{prefix}samtools_depth.txt").read().splitlines()
        if depths:
            depths = np.array([int(d) for d in depths])
        else:
            raise Exception("Insufficient coverage to proceed with CG analysis")

        stats["depth_avg"] = depths.mean()
        stats["depth_q.25"] = np.quantile(depths, .25)
        stats["depth_q.5"] = np.quantile(depths, .5)
        stats["depth_q.75"] = np.quantile(depths, .75)
        stats["depth_frac_above_10x"] = (depths >= 10).mean()
        stats["depth_frac_above_25x"] = (depths >= 30).mean()
        stats["depth_frac_above_50x"] = (depths >= 30).mean()
        stats["depth_frac_above_100x"] = (depths >= 100).mean()

        ax = sns.lineplot(np.arange(1, len(depths)+1), depths)
        ax.set_title("~{sample}")
        ax.set(xlabel="position", ylabel="depth")
        plt.yscale("symlog")
        plt.savefig("~{prefix}depths.png")

        seq, = SeqIO.parse("~{assembly}", "fasta")
        stats["allele_counts"] = dict(collections.Counter(str(seq.seq)))

        fq_list=list(filter(None, ["~{sep='\",\"' fastqs}"]))
        try:
            fq_lines = 0
            for fq_file in fq_list:
                with gzip.open(fq_file, 'r') as f:
                    for line in f: fq_lines += 1
        except OSError:
            fq_lines = 0
            for fq_file in fq_list:
                with open(fq_file, 'r') as f:
                    for line in f: fq_lines += 1

        stats["total_reads"] = int(int(fq_lines) / 4)

        with open("~{prefix}samtools_stats.txt") as f:
            sam_stats_re = re.compile(r"SN\s+([^\s].*):\s+(\d+)")
            for line in f:
                matched = sam_stats_re.match(line)
                if matched:
                    if matched.group(1) == "reads mapped":
                        stats["mapped_reads"] = int(matched.group(2))
                    elif matched.group(1) == "reads mapped and paired":
                        stats["mapped_paired"] = int(matched.group(2))
                    elif matched.group(1) == "inward oriented pairs":
                        stats["paired_inward"] = int(matched.group(2)) * 2
                    elif matched.group(1) == "outward oriented pairs":
                        stats["paired_outward"] = int(matched.group(2)) * 2
                    elif matched.group(1) == "pairs with other orientation":
                        stats["paired_other_orientation"] = int(matched.group(2)) * 2

        if "~{technology}" == "Illumina":
            with open("~{ercc_stats}") as f:
                ercc_stats_re = re.compile(r"SN\s+([^\s].*):\s+(\d+)")
                for line in f:
                    matched = ercc_stats_re.match(line)
                    if matched:
                        if matched.group(1) == "reads mapped":
                            stats["ercc_mapped_reads"] = int(matched.group(2))
                        elif matched.group(1) == "reads mapped and paired":
                            stats["ercc_mapped_paired"] = int(matched.group(2))

        def countVCF(vcf_file, snpcol, mnpcol, indelcol, statsdict):
            vcf = pysam.VariantFile(vcf_file)
            statsdict[snpcol] = 0
            statsdict[mnpcol] = 0
            statsdict[indelcol] = 0
            for rec in vcf.fetch():
                allele_lens = set([len(a) for a in [rec.ref] + list(rec.alts)])
                if len(allele_lens) > 1:
                    statsdict[indelcol] += 1
                else:
                    l, = allele_lens
                    if l == 1:
                        statsdict[snpcol] += 1
                    else:
                        statsdict[mnpcol] += 1
            return statsdict

        stats = {**stats, **countVCF("~{vcf}", 'ref_snps', 'ref_mnps', 'ref_indels', stats)}

        allele_counts = stats["allele_counts"]
        stats["n_actg"] = sum(v for k, v in allele_counts.items() if k in "ACTGU")
        if "N" in allele_counts.keys():
            stats["n_missing"] = allele_counts["N"]
        else:
            stats["n_missing"] = 0
        if "-" in allele_counts.keys():
            stats["n_gap"] = allele_counts["-"]
        else:
            stats["n_gap"] = 0
        stats["n_ambiguous"] = sum(v for k, v in allele_counts.items() if k not in "ACTGUN-")

        with open("~{prefix}stats.json", "w") as f:
            json.dump(stats, f, indent=2)

        CODE
    >>>

    output {
        File sam_depths = "~{prefix}samtools_depth.txt"
        File sam_stats = "~{prefix}samtools_stats.txt"
        File depths_fig = "~{prefix}depths.png"
        File output_stats = "~{prefix}stats.json"
    }

    runtime {
        docker: docker_image_id
    }
}

# NOTE: if we add this step, we need to make it conditional on whether or not the pipeline is running for SARS-CoV-2.
# Expanding to other viruses will require downloading the full set of VADR models
task Vadr {
    # Based on original work at https://github.com/AndrewLangvt/genomic_analyses/blob/v0.4.5/tasks/task_ncbi.wdl
    # Requires coronavirus VADR models from https://ftp.ncbi.nlm.nih.gov/pub/nawrocki/vadr-models/coronaviridae/CURRENT/
    input {
        String prefix
        File assembly
        String vadr_options
        File vadr_model
        String docker_image_id
    }

    command <<<
        set -e
        source /etc/profile
        mkdir -p /usr/local/share/vadr/models
        tar xzvf "~{vadr_model}" -C /usr/local/share/vadr/models --strip-components 1
        # find available RAM
        RAM_MB=$(free -m | head -2 | tail -1 | awk '{print $2}')
        
        {
            # run VADR
            v-annotate.pl ~{vadr_options} --mxsize $RAM_MB "~{assembly}" "vadr-output"        
        } || {
            # in validation, some samples fail with errors: 
            # ... ERROR in cmalign_run(), cmalign failed in a bad way...
            # ... ERROR, at least one sequence name exceeds the maximum GenBank allowed length of 50...
            # we want to capture VADR errors in outputs but these should not cause the workflow to fail entirely
            grep "ERROR" vadr-output/vadr-output.vadr.log > vadr_error.txt
        }
    >>>

    output {
        File? vadr_errors = "vadr_error.txt"
        File? vadr_quality = "vadr-output/vadr-output.vadr.sqc"
        File? vadr_alerts = "vadr-output/vadr-output.vadr.alt.list"
    }

    runtime {
        docker: docker_image_id
    }
}

task ZipOutputs {
    input {
        String prefix
        Array[File] outputFiles
        String docker_image_id
    }

    command <<<
        set -euxo pipefail

        export TMPDIR=${TMPDIR:-/tmp}

        mkdir ${TMPDIR}/outputs
        cp ~{sep=' ' outputFiles} ${TMPDIR}/outputs/
        zip -r -j ~{prefix}outputs.zip ${TMPDIR}/outputs/
    >>>

    output {
        File output_zip = "~{prefix}outputs.zip"
    }

    runtime {
        docker: docker_image_id
    }
}
