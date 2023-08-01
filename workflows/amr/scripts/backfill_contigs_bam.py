from collections import namedtuple
import csv
import json
from urllib.parse import urlparse


import boto3
import click

from create_indexed_contigs_bam import cli as index_contigs


# backfill_data format:
# [
#     {
#         "final_summary": "s3://bucket/path/tsv",
#         "contigs": "s3://bucket/path/contigs.fasta",
#         "sorted_bam": "s3://bucket/path/sorted.bam",
#         "bai_index": "s3://bucket/path/sorted.bam.bai",
#         "id": 42
#     },
#     ...
# ]

FINAL_SUMMARY = "final_summary"
CONTIGS = "contigs"
BAM_FILE = "sorted_bam"
BAI_INDEX = "bai_index"
WORKFLOW_ID = "id"

ParsedS3Url = namedtuple("ParsedS3Url", ["Bucket", "Key"])

SUMMARY_FILE_LOC = "/tmp/final_summary.tsv"
CONTIGS_FILE_LOC = "/tmp/contigs.fasta"

OUTPUT_BAM = "./contig_amr_report.sorted.bam"
OUTPUT_BAI = "./contig_amr_report.sorted.bam.bai"

OUTPUT_LOG = "backfill_log.tsv"


@click.command()
@click.argument("backfill_data", type=click.File("r"), required=True)
@click.pass_context
def cli(ctx, backfill_data: click.File):
    workflow_list = json.load(backfill_data)
    log = []

    s3_client = boto3.client("s3")

    for entry in workflow_list:
        workflow_id = entry[WORKFLOW_ID]
        try:
            s3_path_summary = parse_s3_url(entry[FINAL_SUMMARY])
            s3_path_contigs = parse_s3_url(entry[CONTIGS])

            # we have to write the files to disk for pysam compatibility
            s3_client.download_file(s3_path_summary.Bucket, s3_path_summary.Key, SUMMARY_FILE_LOC)
            s3_client.download_file(s3_path_contigs.Bucket, s3_path_contigs.Key, CONTIGS_FILE_LOC)

            # run the actual script for indexing contigs
            ctx.invoke(index_contigs, contigs_file=CONTIGS_FILE_LOC, final_summary=SUMMARY_FILE_LOC)

            # make backup copies of old files
            # upload newly created BAM and BAI files
            s3_path_bam_file = parse_s3_url(entry[BAM_FILE])
            s3_path_bai_index = parse_s3_url(entry[BAI_INDEX])

            s3_client.copy(s3_path_bam_file._asdict(), s3_path_bam_file.Bucket, f"{s3_path_bam_file.Key}.pre-1-2-15")
            s3_client.copy(s3_path_bai_index._asdict(), s3_path_bai_index.Bucket, f"{s3_path_bai_index.Key}.pre-1-2-15")

            s3_client.upload_file(OUTPUT_BAM, s3_path_bam_file.Bucket, s3_path_bam_file.Key)
            s3_client.upload_file(OUTPUT_BAI, s3_path_bai_index.Bucket, s3_path_bai_index.Key)

            log.append({"workflow_id": workflow_id, "success": True})
        except Exception as e:
            log.append({"workflow_id": workflow_id, "success": False, "error": f"{type(e).__name__}: {str(e)}"})
            continue
    with open(OUTPUT_LOG, "w") as logfile:
        fieldnames = ["workflow_id", "success", "error"]
        writer = csv.DictWriter(logfile, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for item in log:
            writer.writerow(item)
    return


def parse_s3_url(s3_url: str) -> ParsedS3Url:
    parse = urlparse(s3_url, allow_fragments=False)
    return ParsedS3Url(Bucket=parse.netloc, Key=parse.path.lstrip("/"))


if __name__ == '__main__':
    cli()
