import argparse
import boto3
import json
from urllib.parse import urlparse
import subprocess

s3 = boto3.resource("s3")


PAIRED = "czid_short_read_mngs.host_filter.validate_input_out_valid_input2_fastq"
TAXON_COUNTS = "czid_short_read_mngs.postprocess.refined_taxon_count_out_assembly_refined_taxon_counts_with_dcr_json"
CONTIGS_FASTA = "czid_short_read_mngs.postprocess.assembly_out_assembly_contigs_fasta"
CONTIGS_SUMMARY = "czid_short_read_mngs.postprocess.contig_summary_out_assembly_combined_contig_summary_json"


def s3object(s3uri):
    parts = urlparse(s3uri)
    return s3.Object(parts.netloc, parts.path.lstrip("/"))


def format_output_files(output_json_path, run):
    if output_json_path.startswith("s3://"):
        # when running short-read-mngs in the cloud, we have files for each stage
        # TODO: genericize this
        stages = ["host_filter", "non_host_alignment", "postprocess", "experimental"]
        output_json = {}
        for stage in stages:
            output_json.update(
                format_output_files_s3(
                    f"{output_json_path.rstrip('/')}/{stage}_output.json", run
                )
            )

    else:
        output_json = format_output_files_local(output_json_path, run)

    return {
        f"taxon_counts_{run}": output_json[TAXON_COUNTS],
        f"contig_summary_{run}": output_json[CONTIGS_SUMMARY],
        f"contig_fasta_{run}": output_json[CONTIGS_FASTA],
        f"step_counts_{run}": [
            v for k, v in output_json.items() if k.endswith("_count") and v
        ],
    }


def format_output_files_s3(output_json_path, run):
    output_file = json.load(s3object(output_json_path).get()["Body"])
    return {
        "czid_short_read_mngs." + k.lstrip("czid_"): v for k, v in output_file.items()
    }  # reformat keys based on what's expected from local runs


def format_output_files_local(output_json_path, run):
    with open(output_json_path) as f:
        output_json = json.load(f)
    return output_json


def run_miniwdl(input_json):
    benchmark_wdl_path = "workflows/benchmark/run.wdl"
    input_json["workflow_type"] = "short-read-mngs"
    input_json["docker_image_id"] = "czid-benchmark"
    cmd = [
        ".venv/bin/miniwdl",
        "run",
        benchmark_wdl_path,
        "--input",
        json.dumps(input_json),
        "--no-color",
        "--verbose",
    ]
    subprocess.run(cmd, check=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("run1")
    parser.add_argument("-r", "--ref", required=False)
    parser.add_argument("-t", "--truth", required=False)

    args = parser.parse_args()

    input_json = format_output_files(args.run1, run="run_1")
    if args.ref:
        ref_json = format_output_files(args.ref, run="run_2")
        input_json.update(ref_json)

    if args.truth:
        input_json["ground_truth"] = args.truth

    run_miniwdl(input_json)


if __name__ == "__main__":
    main()
