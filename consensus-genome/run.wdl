# The following pipeline was initially based on previous work at: https://github.com/czbiohub/sc2-illumina-pipeline
# workflow version: consensus-genomes-1.5.0
version 1.0

workflow consensus_genome {
    input {
        # Required parameters
        File fastqs_0
        File? fastqs_1

        String docker_image_id
        File ercc_fasta
        File kraken2_db_tar_gz
        File primer_bed
        File ref_fasta
        File ref_host

        # Sample name : include in tags and files
        String sample

        # Optional prefix to add to the output filenames
        String prefix = ""

        # Step parameters

        # all of the following seem to be parameters in params.* of the nextflow pipeline;
        # could they be extracted to a separate params file? put in .json input to wdl?
        Boolean filter_reads = true
        Boolean trim_adapters = true

        Float ivarFreqThreshold = 0.9
        Int   ivarQualTreshold  = 20
        Int   minDepth          = 10

        String no_reads_quast = false

        # assumes about 20 mutations between 2 random samples
        # (this is an overestimate to increase sensitivity)
        Float bcftoolsCallTheta = 0.0006

        # Dummy values - required by SFN interface
        String s3_wd_uri = ""
    }

    call RemoveHost {
        input:
            prefix = prefix,
            fastqs = select_all([fastqs_0, fastqs_1]),
            ref_host = ref_host,
            docker_image_id = docker_image_id
    }

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
                ref_fasta = ref_fasta,
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
            ref_fasta = ref_fasta,
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
            ivarQualThreshold = ivarQualTreshold,
            docker_image_id = docker_image_id
    }

    call Quast {
        input:
            prefix = prefix,
            assembly = MakeConsensus.consensus_fa,
            bam = TrimPrimers.trimmed_bam_ch,
            # use trimReads output if we ran it; otherwise fall back to FilterReads output
            fastqs = select_first([TrimReads.trimmed_fastqs, FilterReads.filtered_fastqs]),
            ref_fasta = ref_fasta,
            no_reads_quast = no_reads_quast,
            docker_image_id = docker_image_id
    }

    call CallVariants {
        input:
            prefix = prefix,
            call_variants_bam = TrimPrimers.trimmed_bam_ch,
            ref_fasta = ref_fasta,
            bcftoolsCallTheta = bcftoolsCallTheta,
            minDepth = minDepth,
            docker_image_id = docker_image_id
    }

    call ComputeStats {
        input:
            prefix = prefix,
            sample = sample,
            cleaned_bam = TrimPrimers.trimmed_bam_ch,
            assembly = MakeConsensus.consensus_fa,
            ercc_stats = QuantifyERCCs.ercc_out,
            vcf = CallVariants.variants_ch,
            fastqs = select_all([fastqs_0, fastqs_1]),
            ref_host = ref_host,
            docker_image_id = docker_image_id
    }

    call ZipOutputs {
        input:
            prefix = prefix,
            outputFiles = select_all(flatten([
                RemoveHost.host_removed_fastqs,
                select_all([
                    MakeConsensus.consensus_fa,
                    ComputeStats.depths_fig,
                    TrimPrimers.trimmed_bam_ch,
                    TrimPrimers.trimmed_bam_bai,
                    Quast.quast_txt,
                    Quast.quast_tsv,
                    AlignReads.alignments,
                    QuantifyERCCs.ercc_out,
                    QuantifyERCCs.ercc_out,
                    ComputeStats.output_stats,
                    ComputeStats.sam_depths,
                    CallVariants.variants_ch,
                ])
            ])),
            docker_image_id = docker_image_id
    }

    output {
        Array[File] remove_host_out_host_removed_fastqs = RemoveHost.host_removed_fastqs
        File quantify_erccs_out_ercc_out = QuantifyERCCs.ercc_out
        Array[File]+? filter_reads_out_filtered_fastqs = FilterReads.filtered_fastqs
        Array[File]+? trim_reads_out_trimmed_fastqs = TrimReads.trimmed_fastqs
        File? align_reads_out_alignments = AlignReads.alignments
        File? trim_primers_out_trimmed_bam_ch = TrimPrimers.trimmed_bam_ch
        File? trim_primers_out_trimmed_bam_bai = TrimPrimers.trimmed_bam_bai
        File? make_consensus_out_consensus_fa = MakeConsensus.consensus_fa
        File? quast_out_quast_txt = Quast.quast_txt
        File? quast_out_quast_tsv = Quast.quast_tsv
        File? call_variants_out_variants_ch = CallVariants.variants_ch
        File? compute_stats_out_depths_fig = ComputeStats.depths_fig
        File? compute_stats_out_output_stats = ComputeStats.output_stats
        File? compute_stats_out_sam_depths = ComputeStats.sam_depths
        File zip_outputs_out_output_zip = ZipOutputs.output_zip
    }
}

# TODO: task to validate input

task RemoveHost {
    input {
        String prefix
        Array[File]+ fastqs
        File ref_host

        String docker_image_id
    }

    command <<<
        set -euxo pipefail

        export CORES=`nproc --all`
        if [[ "~{length(fastqs)}" == 1 ]]; then
            minimap2 -t $CORES -ax sr ~{ref_host} ~{fastqs[0]} | \
            samtools view -@ $CORES -b -f 4 | \
            samtools fastq -@ $CORES -0 "~{prefix}no_host_1.fq.gz" -n -c 6 -
        else
            minimap2 -t $CORES -ax sr ~{ref_host} ~{sep=' ' fastqs} | \
            samtools view -@ $CORES -b -f 4 | \
            samtools fastq -@ $CORES -1 "~{prefix}no_host_1.fq.gz" -2 "~{prefix}no_host_2.fq.gz" -0 /dev/null -s /dev/null -n -c 6 -
        fi

        if [ -z $(gzip -cd "~{prefix}no_host_1.fq.gz" | head -c1) ] && ([-z "~{prefix}no_host_2.fq.gz"] || [ -z $(gzip -cd "~{prefix}no_host_2.fq.gz" | head -c1) ]); then
            set +x
            >&2 echo "{\"wdl_error_message\": true, \"error\": \"InsufficientReadsError\", \"cause\": \"No reads after RemoveHost\"}"
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

        minimap2 -ax sr ~{ercc_fasta} ~{sep=' ' fastqs} | samtools view -bo ercc_mapped.bam
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
            >&2 echo "{\"wdl_error_message\": true, \"error\": \"InsufficientReadsError\", \"cause\": \"No reads after FilterReads\"}"
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
            tar -xzv -C "${TMPDIR}/kraken_db" -f "~{kraken2_db_tar_gz}"
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

        echo `gzip -cd "~{prefix}filtered_1.fq.gz" | head -c1`
        if [ -z $(gzip -cd "~{prefix}filtered_1.fq.gz" | head -c1) ] && ([[ "~{length(fastqs)}" == 1 ]] || [ -z $(gzip -cd "~{prefix}filtered_2.fq.gz" | head -c1) ]); then
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

        if [ -z $(gzip -cd "${BASENAME}_val_1.fq.gz" | head -c1) ] && ([[ "~{length(fastqs)}" == 1 ]] || [ -z $(gzip -cd "${BASENAME}_val_2.fq.gz" | head -c1) ]); then
            set +x
            >&2 echo "{\"wdl_error_message\": true, \"error\": \"InsufficientReadsError\", \"cause\": \"No reads after TrimReads\"}"
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
        ivar trim -e -i ivar.bam -b "~{primer_bed}" -p ivar.out
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
            >&2 echo "{\"wdl_error_message\": true, \"error\": \"InsufficientReadsError\", \"cause\": \"No reads after MakeConsensus\"}"
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

task Quast {
    input {
        String prefix
        File assembly   # same as consensus_fa
        File bam
        Array[File]+ fastqs
        File ref_fasta

        String no_reads_quast

        Int threads = 4

        String docker_image_id
    }

    command <<<
        set -euxo pipefail

        export CORES=`nproc --all`
        a=`cat "~{assembly}" | wc -l`

        if [ $a -ne 0 ]; then
            if [ ~{no_reads_quast} = true ]; then
                quast --min-contig 0 -o quast -r "~{ref_fasta}" -t $CORES --ref-bam "~{bam}" "~{assembly}"
            else
                if [[ "~{length(fastqs)}" == 1 ]]; then
                    quast --min-contig 0 -o quast -r "~{ref_fasta}" -t $CORES --ref-bam "~{bam}" "~{assembly}" --single "~{fastqs[0]}"
                else
                    quast --min-contig 0 -o quast -r "~{ref_fasta}" -t $CORES --ref-bam "~{bam}" "~{assembly}" -1 ~{sep=' -2 ' fastqs}
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

task CallVariants {
    input {
        String prefix
        File call_variants_bam  # same as primertrimmed_bam produced by trimPrimers
        File ref_fasta
        Float bcftoolsCallTheta
        Int minDepth

        String docker_image_id
    }

    command <<<
        set -euxo pipefail

        # NOTE: we use samtools instead of bcftools mpileup because bcftools 1.9 ignores -d0
        samtools mpileup -u -d 0 -t AD -f "~{ref_fasta}" "~{call_variants_bam}" | bcftools call --ploidy 1 -m -P "~{bcftoolsCallTheta}" -v - | bcftools view -i 'DP>=~{minDepth}' > "~{prefix}variants.vcf"
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

task ComputeStats {
    input {
        String prefix
        String sample
        File cleaned_bam
        File assembly
        File ercc_stats
        File vcf

        Array[File]+ fastqs
        File ref_host

        String docker_image_id
    }

    command <<<
        set -euxo pipefail

        samtools index "~{cleaned_bam}"
        samtools stats "~{cleaned_bam}" > "~{prefix}samtools_stats.txt"
        samtools depth -aa -d 0 "~{cleaned_bam}" | awk '{print $3}' > "~{prefix}samtools_depth.txt"

        python <<CODE

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
            depths = np.array([0]*pysam.AlignmentFile("~{cleaned_bam}", "rb").lengths([0]))

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
