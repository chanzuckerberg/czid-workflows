def star_description(collect_insert_size_metrics_for):
    outSAMmode = "--outSAMmode None"
    collect_insert_size_metrics_for = collect_insert_size_metrics_for.lower()
    if collect_insert_size_metrics_for == "dna":
        outSAMmode = """--outSAMtype BAM Unsorted
            --outSAMmode NoQS"""

    if collect_insert_size_metrics_for == "rna":
        outSAMmode = """--outSAMtype BAM Unsorted
            --outSAMmode NoQS
            --quantMode TranscriptomeSAM GeneCounts"""

    description = ""

    if collect_insert_size_metrics_for:
        description += """
            **Host Subtraction**
        """

    description += f"""
        Implements the step for Host Subtraction.

        The STAR aligner is used for rapid first-pass host filtration.
        Unmapped reads are passed to the subsequent step. The current implementation of STAR,
        will fail to remove host sequences that map to multiple regions, thus these are filtered
        out by a subsequent host filtration step using Bowtie2.

        Different parameters are required for alignment of short vs long reads using STAR.
        Therefore, based on the initial input validation, the appropriate parameters are selected.

        If short reads:
        ```
        STAR
        --outFilterMultimapNmax 99999
        --outFilterScoreMinOverLread 0.5
        --outFilterMatchNminOverLread 0.5
        --outReadsUnmapped Fastx
        --outFilterMismatchNmax 999
        {outSAMmode}
        --clip3pNbases 0
        --runThreadN {{cpus}}
        --genomeDir {{genome_dir}}"
        --readFilesIn {{input files}}
        ```

        If long reads (specifically if there are more than 1 reads with length greater than
        READ_LEN_CUTOFF_HIGH, as determined during input validation step):
        ```
        STARlong
        --outFilterMultimapNmax 99999
        --outFilterScoreMinOverLread 0.5
        --outFilterMatchNminOverLread 0.5
        --outReadsUnmapped Fastx
        --outFilterMismatchNmax 999
        {outSAMmode}
        --clip3pNbases 0
        --runThreadN {{cpus}}
        --genomeDir {{genome_dir}}
        --readFilesIn {{input files}}
        --seedSearchStartLmax 20
        --seedPerReadNmax 100000
        --seedPerWindowNmax 1000
        --alignTranscriptsPerReadNmax 100000
        ```

        STAR documentation can be found [here](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)
    """

    if collect_insert_size_metrics_for:
        description += """
            **Insert Metrics**

            This step also computes insert size metrics for Paired End samples from human hosts.
            These metrics are computed by the Broad Institute's Picard toolkit.

            Picard is run on the output BAM file obtained from running STAR:

            ```
            java -jar picard.jar CollectInsertSizeMetrics
                I={output bam file}
                O=picard_insert_metrics.txt
                H=insert_size_histogram.pdf
            ```

            Picard documentation can be found [here](https://broadinstitute.github.io/picard/)
        """

    return description


def main():
    import sys

    step_name = sys.argv[1]
    if step_name == "star_out":
        description = star_description(sys.argv[2])

    with open(f"{step_name}.description.md", "w+") as f:
        f.write(description)


if __name__ == "__main__":
    main()
