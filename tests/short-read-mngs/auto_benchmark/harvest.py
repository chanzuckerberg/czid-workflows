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
import boto3
from pathlib import Path
from contextlib import ExitStack
from urllib.parse import urlparse
from _util import load_benchmarks_yml
from taxadb.taxid import TaxID

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
        if rundir.startswith("s3:"):
            assert os.environ.get(
                "AWS_PROFILE", False
            ), f"set environment AWS_PROFILE to read from {rundir}"
        else:
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


if __name__ == "__main__":
    main()
