import csv
import re
import sys
from typing import Tuple
from urllib.parse import urlparse

import boto3
import click

OUTPUT_BAM = "contig_amr_report.sorted.bam"
OUTPUT_BAI = "contig_amr_report.sorted.bam.bai"

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
        except Exception as e:
            log.append("workflow_id": workflow_id, "success": False)
            continue
        finally:
            log.append("workflow_id": workflow_id, "success": True)
    
    return


def parse_s3_url(s3_url: str) -> Dict[str, str]:
    parse = urlparse(s3_url, allow_fragments=False)
    return {
        "bucket": parse.netloc,
        "key": parse.path.lstrip("/")
    }
