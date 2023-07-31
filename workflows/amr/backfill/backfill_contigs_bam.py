from collections import namedtuple
import csv
import json
import re
import sys
from typing import Tuple
from urllib.parse import urlparse


import boto3
import click

from create_indexed_contigs_bam import cli as index_contigs


# execution_array format:
# [
#     {
#         "final_summary": "s3://bucket/path/tsv",
#         "contigs_fasta": "s3://bucket/path/contigs.fasta",
#         "sorted_bam": "s3://bucket/path/sorted.bam",
#         "bai_index": "s3://bucket/path/sorted.bam.bai",
#         "id": 42
#     },
#     ...
# ]

FINAL_SUMMARY = "final_summary"
CONTIGS = "contigs_fasta"
BAM_FILE = "sorted_bam"
BAI_INDEX = "bai_index"
WORKFLOW_ID = "id"

ParsedS3Url = namedtuple("ParsedS3Url", ["bucket", "key"])

SUMMARY_FILE_LOC = "/tmp/final_summary.tsv"
CONTIGS_FILE_LOC = "/tmp/contigs.fasta"

OUTPUT_BAM = "./contig_amr_report.sorted.bam"
OUTPUT_BAI = "./contig_amr_report.sorted.bam.bai"


@click.command()
@click.argument("execution_array", type=click.File("r"), required=True)
@click.pass_context
def cli(ctx, execution_array: click.File):
    workflow_list = json.load(execution_array)
    log = []

    s3_client = boto3.client("s3")

    for entry in workflow_list:
        success = False
        workflow_id = entry[WORKFLOW_ID]
        try:
            s3_path_summary = parse_s3_url(entry[FINAL_SUMMARY])
            s3_path_contigs = parse_s3_url(entry[CONTIGS])

            # we have to write the files to disk for pysam compatibility
            s3_client.download_file(s3_path_summary.bucket, s3_path_summary.key, SUMMARY_FILE_LOC)
            s3_client.download_file(s3_path_contigs.bucket, s3_path_contigs.key, CONTIGS_FILE_LOC)

            # run the actual script for indexing contigs
            ctx.invoke(index_contigs, final_summary=SUMMARY_FILE_LOC, contigs=CONTIGS_FILE_LOC)

            # upload newly created BAM and BAI files
            s3_path_bam_file = parse_s3_url(entry[BAM_FILE])
            s3_path_bai_index = parse_s3_url(entry[BAI_INDEX])

            s3_client.upload_file(OUTPUT_BAM, s3_path_bam_file.bucket, s3_path_bam_file.key)
            s3_client.upload_file(OUTPUT_BAI, s3_path_bai_index.bucket, s3_path_bai_index.key)

            log.append({ "workflow_id": workflow_id, "success": True })
        except Exception as e:
            log.append({ "workflow_id": workflow_id, "success": False, "error": str(e)})
            continue
    with open("log.tsv", "w") as logfile:
        fieldnames = ["workflow_id", "success", "error"]
        writer = csv.DictWriter(logfile, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for item in log:
            writer.writerow(item)
    return


def parse_s3_url(s3_url: str) -> ParsedS3Url:
    parse = urlparse(s3_url, allow_fragments=False)
    return ParsedS3Url(bucket=parse.netloc, key=parse.path.lstrip("/"))


if __name__ == '__main__':
    cli()
