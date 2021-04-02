#!/usr/bin/env python3
"""
Harvest results of benchmarks.yml sample runs (by either run_local.py or run_dev.py) into JSON
input for short-read-mngs-benchmarks.ipynb

Uses taxadb.sqlite to resolve species names. That can be generated in a minute or two like so:
    taxadb download -o taxadb --type taxa
    taxadb create -i taxadb --dbname taxadb.sqlite
"""
import sys
import os
import argparse
import json
import itertools
from pathlib import Path
from contextlib import ExitStack
from urllib.parse import urlparse

import boto3
from botocore import UNSIGNED as BOTOCORE_UNSIGNED
from botocore.config import Config as BotocoreConfig
from taxadb.taxid import TaxID
import numpy as np

from _util import load_benchmarks_yml, adjusted_aupr

BENCHMARKS = load_benchmarks_yml()


def main():
    parser = argparse.ArgumentParser(
        sys.argv[0], description="harvest benchmark runs into importable JSON"
    )

    parser.add_argument(
        "outputs",
        metavar="sample=path/to/rundir",
        type=str,
        nargs="+",
        help="sample identifier & path to completed run directory",
    )
    parser.add_argument(
        "--taxadb", metavar="FILENAME", type=str, help="taxadb SQLite file, if available"
    )

    args = parser.parse_args(sys.argv[1:])
    harvest(**vars(args))


def harvest(outputs, taxadb):
    queue = []

    # process command line args
    i = 0
    while i < len(outputs):
        if outputs[i].endswith("="):
            assert i + 1 < len(outputs), f"missing path after {outputs[i]}"
            sample = outputs[i][:-1]
            rundir = outputs[i + 1]
            i += 2
        else:
            parts = outputs[i].split("=")
            assert len(parts) == 2, f"invalid SAMPLE_ID=RUNDIR pair: {outputs[i]}"
            sample = parts[0]
            rundir = parts[1]
            i += 1
        assert sample in BENCHMARKS["samples"], f"unknown sample {sample}"
        if not rundir.startswith("s3:"):
            assert (
                Path(rundir) / "outputs.json"
            ).is_file(), f"couldn't find outputs.json in {rundir}"
        queue.append((sample, rundir))

    if taxadb:
        taxadb = TaxID(dbtype="sqlite", dbname=taxadb)

    # harvest each supplied sample
    rslt = {}
    for sample, rundir in queue:
        assert sample not in rslt, f"repeated sample {sample}"
        outputs_json = read_outputs_json(rundir)
        rslt[sample] = harvest_sample(sample, outputs_json, taxadb)
        rslt[sample]["outputs"] = outputs_json

    print(json.dumps(rslt, indent=2))


def harvest_sample(sample, outputs_json, taxadb):
    ans = {}

    # collect read counts at various pipeline steps
    ans["paired"] = (
        outputs_json["idseq_short_read_mngs.host_filter.validate_input_out_valid_input2_fastq"]
        is not None
    )
    ans["input_reads"] = read_output_jsonfile(outputs_json, "host_filter.input_read_count")[
        "fastqs"
    ]
    for step in [
        "validate_input",
        "star",
        "trimmomatic",
        "priceseq",
        "idseq_dedup",
        "lzw",
        "bowtie2",
        "subsampled",
        "gsnap_filter",
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


def harvest_sample_taxon_counts(
    sample, outputs_json, contig_summary, contig_lengths, dbtype, taxadb
):
    assert dbtype in ("NR", "NT")

    # read in the taxon counts & contig summary JSON files
    taxon_counts = read_output_jsonfile(
        outputs_json,
        "postprocess.refined_taxon_count_out_assembly_refined_taxon_counts_with_dcr_json",
    )["pipeline_output"]["taxon_counts_attributes"]

    # for each species in the taxon counts (excluding genus/family for now)
    ans = {}
    for rslt in taxon_counts:
        if rslt["count_type"] == dbtype and rslt["tax_level"] == 1 and int(rslt["tax_id"]) >= 0:
            # lookup corresponding contig summary
            contig_summary_key = (dbtype, rslt["tax_id"])
            if contig_summary_key in contig_summary:
                rslt_contigs = contig_summary[contig_summary_key]["contig_counts"]
            else:
                rslt_contigs = {}
            assert rslt["tax_id"] not in ans

            # combine info
            info = {
                "tax_level": rslt["tax_level"],
                "reads": rslt["nonunique_count"],
                "reads_dedup": rslt["unique_count"],
                "avg_aln_len": rslt["alignment_length"],
            }
            info.update(contigs_stats(contig_lengths, (k for k in rslt_contigs if k != "*")))
            info["contigs_reads"] = sum(rslt_contigs[k] for k in rslt_contigs if k != "*")
            if taxadb:
                info["tax_name"] = taxadb.sci_name(int(rslt["tax_id"]))
            ans[rslt["tax_id"]] = info

    # sort by abundance
    return {
        k: v
        for k, v in sorted(
            ans.items(), key=lambda kv: (kv[1]["reads_dedup"], kv[1]["avg_aln_len"]), reverse=True
        )
    }


def read_outputs_json(rundir):
    if not rundir.startswith("s3:"):
        # local miniwdl run directory
        return read_json(str(Path(rundir) / "outputs.json"))

    # read SFN-WDL output JSONs from S3, reshape them to resemble the local outputs.json
    items = itertools.chain(
        *(
            read_json(os.path.join(rundir, stage + "_output.json")).items()
            for stage in (
                "host_filter",
                "non_host_alignment",
                "postprocess",
                "experimental",
            )
        )
    )
    ans = {}
    for key, value in items:
        ans["idseq_short_read_mngs." + key[6:]] = value
    return ans


def read_json(path):  # from either local or s3 path
    if not path.startswith("s3:"):
        with open(path) as infile:
            return json.load(infile)
    return json.load(s3object(path).get()["Body"])


def s3object(s3uri):
    parts = urlparse(s3uri)
    return boto3.resource("s3").Object(parts.netloc, parts.path.lstrip("/"))


def read_output_jsonfile(outputs_json, key):
    """
    From the WDL outputs dict, read the contents of an output JSON file
    """
    return read_json(outputs_json["idseq_short_read_mngs." + key])


def load_contig_lengths(outputs_json):
    """
    Generate dict contig id -> length
    """
    fasta = outputs_json["idseq_short_read_mngs.postprocess.assembly_out_assembly_contigs_fasta"]
    lengths = {}
    with ExitStack() as cleanup:
        if fasta.startswith("s3:"):
            lines = (line.decode("utf-8") for line in s3object(fasta).get()["Body"].iter_lines())
        else:
            lines = cleanup.enter_context(open(fasta))
        cur = None
        for line in lines:
            if line.startswith(">"):
                if cur is not None:
                    lengths[cur[0]] = cur[1]
                cur = (line[1:].strip(), 0)
            else:
                cur = (cur[0], cur[1] + len(line.strip()))
        if cur is not None:
            lengths[cur[0]] = cur[1]
    return lengths


def contigs_stats(contig_lengths, ids=None):
    lengths = []
    for id in ids if ids else contig_lengths.keys():
        lengths.append(contig_lengths[id])
    lengths.sort(reverse=True)
    total_nt = sum(lengths)
    cumlen = 0
    N50 = 0
    for i, length_i in enumerate(lengths):
        if cumlen + length_i >= total_nt / 2:
            N50 = length_i
            break
        cumlen += length_i
    return {"contigs": len(lengths), "contigs_nt": total_nt, "contigs_N50": N50}


def read_truth_file(path):
    taxon_abundances = {}
    s3url = urlparse(path)
    s3 = boto3.resource("s3", config=BotocoreConfig(signature_version=BOTOCORE_UNSIGNED))
    for line in (
        s3.Object(s3url.netloc, s3url.path.lstrip("/")).get()["Body"].read().decode().splitlines()
    ):
        taxid, _, abundance, rank, species_name = line.split("\t")[:5]
        # All current benchmark datasets have truth data on species level only,
        # so we use this as a simplifying assumption downstream
        assert rank == "species"
        taxon_abundances[taxid] = float(abundance)
    return taxon_abundances


"""
From Ye et al. 2019:
Precision is the proportion of true positive species identified in the sample divided by the number of total species
identified.
Recall is the proportion of true positive species divided by the number of distinct species actually in the sample.
A potential drawback of AUPR is that it is biased toward low-precision, high-recall classifiers.
Classifiers that do not recall all of the ground-truth taxa are penalized with zero AUPR from the highest achieved
recall to 100% recall. For classifiers that do reach 100% recall, additional false positive taxon calls do not
further penalize the AUPR score.

In addition to considering the number of correctly identified species, it is also important to evaluate how
accurately the abundance of each species or genera in the resulting classification reflects the abundance of each
species in the original biological sample ("ground truth"). This is especially critical for applications such as
microbiome sequencing studies, where changes in population composition can confer phenotypic effects. Abundance can
be considered either as the relative abundance of reads from each taxa ("raw") or by inferring abundance of the
number of individuals from each taxa by correcting read counts for genome size ("corrected"). Some programs
incorporate a correction for genome length into abundance estimates; this calculation can also be manually performed
by reweighting the read counts after classification. Here we use raw abundance profiles unless correction is
performed automatically by the software, as in the case of PathSeq and Bracken.

To evaluate the accuracy of abundance profiles, we can calculate the pairwise distances between ground-truth
abundances and normalized abundance counts for each identified taxon at a given taxonomic level (e.g., species or
genus). For this, we calculate the L2 distance for a given dataset’s classified output as the straight-line distance
between the observed and true abundance vectors. We can also use this measure to compare abundance profiles between
classifiers by instead computing L2 distances between classified abundances for pairs of classifiers. Abundance
profile distance is more sensitive to accurate quantification of the highly abundant taxa present in the sample
(Aitchison, 1982; Quinn et al., 2018). High numbers of very-low-abundance false positives will not dramatically
affect the measure because they comprise only a small portion of the total abundance. For this reason, using such a
measure alongside AUPR, which is highly sensitive to classifiers’ performance in correctly identifying low-abundance
taxa, allows comprehensive evaluation of classifier performance.

The L2 distance should be considered as a representation of the abundance profiles. Because metagenomic abundance
profiles are proportional data and not absolute data, it is important to remember that many common distance metrics
(including L2 distance) are not true mathematical metrics in proportional space. Generally, in proportional data
analysis, a common method is to normalize proportions by using the centered log-ratio transform to calculate
distances. However, the output of these metagenomic classifiers includes many low-abundance false positives, leading
to sparse zero counts for many taxa across the different reports. The log-transform of these zero counts is
undefined unless arbitrary pseudocounts are added to each taxa, which can negatively bias accurate classifiers
because false positive taxa will have added counts. Another commonly used metric to compare abundance profiles is
the UniFrac distance, which considers both the abundance proportion of component taxa as well as the evolutionary
distance for incorrectly called taxa. However, using this metric is complicated by the difficulty in assessing
evolutionary distance between microbial species’ whole genomes.
"""


def truth_aupr(classified_taxa, truth_taxa, total_reads):
    """
    Compute AUPR (and other metrics) of taxon classifications vs. truth data
    """
    missed_taxa = [tax_id for tax_id in truth_taxa if tax_id not in classified_taxa]
    correctness_labels = [1 if tax_id in truth_taxa else 0 for tax_id in classified_taxa]
    correctness_labels += [1] * len(missed_taxa)
    # Using raw abundances as proxies for confidence score per Ye2009 methodology
    # https://www.cell.com/cell/fulltext/S0092-8674(19)30775-5#fig2
    confidence_scores = [i["reads_dedup"] / total_reads for i in classified_taxa.values()]
    confidence_scores += [0] * len(missed_taxa)
    return adjusted_aupr(correctness_labels, confidence_scores, force_monotonic=False)


def truth_l2_norm(classified_taxa, truth_taxa, total_reads):
    """
    Measure accuracy of taxa relative abundance (L2 norm of vector difference from truth vector)
    """
    truth_sum = sum(truth_taxa.values())
    relative_abundances_diff = [
        truth_taxa[taxon] / truth_sum
        - classified_taxa.get(taxon, {}).get("reads_dedup", 1e-100) / total_reads
        for taxon in truth_taxa
    ]
    return np.linalg.norm(relative_abundances_diff, ord=2)


if __name__ == "__main__":
    main()
