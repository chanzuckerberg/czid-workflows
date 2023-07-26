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
#         "amr.RunResultsPerSample.final_summary": "s3://bucket/path/tsv",
#         "amr.ZipOutputs.contigs": "s3://bucket/path/contigs.fasta",
#         "amr.tsvToSam.output_sorted": "s3://bucket/path/sorted.bam",
#         "amr.tsvToSam.output_sorted_bai": "s3://bucket/path/sorted.bam.bai",
#         "workflow_id": 42
#     },
#     {
#         "...": "..."
#     }
# ]

WORKFLOW_KEYS = {
    "final_summary": "amr.RunResultsPerSample.final_summary",
    "contigs_fasta": "amr.ZipOutputs.contigs",
    "sorted_bam": "amr.tsvToSam.output_sorted",
    "bai_index": "amr.tsvToSam.output_sorted_bai",
    "id": "workflow_id",
}

@click.command()
@click.argument("execution_array", type=click.File("r"), required=True)
@click.pass_context
def cli(ctx, execution_array: click.File):
    workflow_list = json.load(execution_array)
    log = []

    s3_client = boto3.client("s3")

    for workflow in workflow_list:
        success = False
        workflow_id = workflow["id"]
        try:
            old_bam_bais = []
            for key in ["sorted_bam", "bai_index"]:
                if key in workflow:
                    parsed_url = parse_s3_url(workflow[key])
                    old_bam_bais.append(parsed_url)
                else:
                    raise Exception
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
