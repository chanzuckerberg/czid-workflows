def harvest_sample(sample, outputs_json, taxadb):
    ans = {}

    # collect read counts at various pipeline steps
    ans["paired"] = (
        outputs_json["czid_short_read_mngs.host_filter.fastp_out_fastp2_fastq"]
        is not None
    )
    ans["input_reads"] = read_output_jsonfile(outputs_json, "host_filter.input_read_count")[
        "fastqs"
    ]
    for step in [
        "validate_input",
        "fastp",
        "bowtie2_host_filtered",
        "hisat2_host_filtered",
        "czid_dedup",
        "subsampled",
    ]:
        ans[step + "_reads"] = read_output_jsonfile(
            outputs_json, "host_filter." + step + "_out_count"
        )[step + "_out"]

    # TODO: count reads mapped to ERCC controls (target name ERCC-*)
    # TODO: count reads with only poor E-value hits to either NR or NT (quasi-unmapped)
    ans["initial_unmapped_reads"] = read_output_jsonfile(
        outputs_json, "non_host_alignment.annotated_out_count"
    )["unidentified_fasta"]
    ans["refined_unmapped_reads"] = read_output_jsonfile(
        outputs_json, "postprocess.refined_annotated_out_count"
    )["unidentified_fasta"]

    contig_summary = read_output_jsonfile(
        outputs_json, "postprocess.contig_summary_out_assembly_combined_contig_summary_json"
    )
    contig_summary = {(elt["count_type"], elt["taxid"]): elt for elt in contig_summary}
    contigs_mapped = set()
    for elt in contig_summary.values():
        for key in elt["contig_counts"]:
            if key != "*":
                contigs_mapped.add(key)

    contig_lengths = load_contig_lengths(outputs_json)
    ans.update(contigs_stats(contig_lengths))
    for key, value in contigs_stats(contig_lengths, contigs_mapped).items():
        ans["mapped_" + key] = value

    # collect NR/NT taxon counts
    ans = {"counts": ans, "taxa": {}}
    ans["taxa"]["NT"] = harvest_sample_taxon_counts(
        sample, outputs_json, contig_summary, contig_lengths, "NT", taxadb
    )
    ans["taxa"]["NR"] = harvest_sample_taxon_counts(
        sample, outputs_json, contig_summary, contig_lengths, "NR", taxadb
    )

    # compute figures of merit vs. truth dataset
    ans["truth"] = read_truth_file(BENCHMARKS["samples"][sample]["truth"])

    total_reads = ans["counts"]["subsampled_reads"] * (2 if ans["counts"]["paired"] else 1)

    nt_aupr = truth_aupr(ans["taxa"]["NT"], ans["truth"], total_reads)
    ans["NT_aupr"] = nt_aupr["aupr"]
    ans["NT_precision"] = list(nt_aupr["precision"])
    ans["NT_recall"] = list(nt_aupr["recall"])
    ans["NT_l2_norm"] = truth_l2_norm(ans["taxa"]["NT"], ans["truth"], total_reads)

    nr_aupr = truth_aupr(ans["taxa"]["NR"], ans["truth"], total_reads)
    ans["NR_aupr"] = nr_aupr["aupr"]
    ans["NR_precision"] = list(nr_aupr["precision"])
    ans["NR_recall"] = list(nr_aupr["recall"])
    ans["NR_l2_norm"] = truth_l2_norm(ans["taxa"]["NR"], ans["truth"], total_reads)

    return ans