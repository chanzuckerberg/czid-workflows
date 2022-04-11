#!/usr/bin/env python3

import boto3
import botocore
import datetime
import os
import subprocess
import logging
import argparse

from collections import namedtuple
import concurrent.futures

S3_BASE_PREFIX = "ncbi-sources"

BUCKET = os.environ["BUCKET"]
REMOTE_SERVER = os.environ["NCBI_SERVER"]
FILES_TO_DOWNLOAD = [
    "blast/db/FASTA/nt.gz",
    "blast/db/FASTA/nt.gz.md5",
    "blast/db/FASTA/nr.gz",
    "blast/db/FASTA/nr.gz.md5",
    "pub/taxonomy/taxdump.tar.gz",
    "pub/taxonomy/taxdump.tar.gz.md5",
    "pub/taxonomy/accession2taxid/dead_nucl.accession2taxid.gz",
    "pub/taxonomy/accession2taxid/dead_nucl.accession2taxid.gz.md5",
    "pub/taxonomy/accession2taxid/dead_prot.accession2taxid.gz",
    "pub/taxonomy/accession2taxid/dead_prot.accession2taxid.gz.md5",
    "pub/taxonomy/accession2taxid/dead_wgs.accession2taxid.gz",
    "pub/taxonomy/accession2taxid/dead_wgs.accession2taxid.gz.md5",
    "pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz",
    "pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz.md5",
    "pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz",
    "pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz.md5",
    "pub/taxonomy/accession2taxid/pdb.accession2taxid.gz",
    "pub/taxonomy/accession2taxid/pdb.accession2taxid.gz.md5",
    "pub/taxonomy/accession2taxid/prot.accession2taxid.gz",
    "pub/taxonomy/accession2taxid/prot.accession2taxid.gz.md5",
    "pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz",
    "pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz.md5"
]

UNZIP_LIST = ["nt.gz", "nr.gz"]

TEMP_FOLDER = "ncbi_sources"

MAX_ATTEMPTS = 3
SEPARATE_FOLDER_PATTERNS = ["accession2taxid"]

# File namedtuple structure, for ease of use
# url = remote path the file lives on
# relative_path = the relative path to the file from the temp folder
File = namedtuple("File", ["url", "relative_path"])

s3 = boto3.resource("s3", endpoint_url=os.environ.get('ENDPOINT_URL'))

# Logger output format config. Example output:
# 2020/03/30 16:04:47 PDT INFO: Done file doesn't exist. Should run.
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(levelname)s: %(message)s',
    datefmt='%Y/%m/%d %H:%M:%S %Z'
    )
logger = logging.getLogger(__name__)


def main(date_tag_override=None):
    # Get current date and use the datestring as the name for our s3 subfolder
    date_tag = date_tag_override if date_tag_override else datetime.datetime.today().strftime("%Y-%m-%d")
    dated_subfolder = f"{S3_BASE_PREFIX}/{date_tag}"

    logger.info(msg=f"Using S3 bucket '{BUCKET}'")
    logger.info(msg=f"Using date tag '{date_tag}'")

    if done_file_exists(BUCKET, dated_subfolder):
        return

    files_list = consolidate_file_information(FILES_TO_DOWNLOAD, REMOTE_SERVER, SEPARATE_FOLDER_PATTERNS)

    multidownload(files_list, TEMP_FOLDER, MAX_ATTEMPTS)

    extracted_list = extract_files(UNZIP_LIST, TEMP_FOLDER)
    unzipped_files = consolidate_file_information(extracted_list, ".", SEPARATE_FOLDER_PATTERNS)
    files_list = list(filter(lambda f: f.relative_path not in UNZIP_LIST, files_list))
    files_list.extend(unzipped_files)

    dated_s3_path = f"{BUCKET}/{dated_subfolder}"
    multi_s3parcp_upload(files_list, TEMP_FOLDER, dated_s3_path)

    write_done_file(BUCKET, dated_subfolder)
    write_index_name(BUCKET, date_tag)

# --- main functions ---


def done_file_exists(BUCKET, dated_subfolder):
    try:
        # Don't run if the done file is there already
        s3.Object(BUCKET, f"{dated_subfolder}/done").load()
        logger.info(msg="Done file exists already. Skipping this run.")
        return True
    except botocore.exceptions.ClientError:
        logger.info(msg="Done file doesn't exist. Should run.")
        return False


def consolidate_file_information(file_list, remote_server, separate_folder_patterns):
    file_tuples = []
    for path in file_list:
        url = f"{remote_server}/{path}"
        filename = os.path.basename(path)
        relative_path = filename
        for pattern in separate_folder_patterns:
            if pattern in path:
                relative_path = f"{pattern}/{relative_path}"
                break
        relative_path
        file_tuples.append(File(url, relative_path))
    return file_tuples


def multidownload(file_list, dest_folder, max_attempts):
    needed_files = file_list
    for attempt in range(max_attempts):
        with concurrent.futures.ThreadPoolExecutor() as executor:
            future_to_url = {executor.submit(_download_file, url, dest_folder, relative_path): url
                             for (url, relative_path) in needed_files}
            for future in concurrent.futures.as_completed(future_to_url):
                url = future_to_url[future]
                try:
                    future.result()
                except Exception as exc:
                    filename = os.path.basename(url)
                    exc_type = type(exc).__name__
                    exc_msg = exc.args
                    logger.warning(
                        msg=f"Error downloading file {filename}.\nURL: {url}\n{exc_type}: {exc_msg}"
                        )
                else:
                    needed_files = list(filter(lambda file: file.url != url, needed_files))
                    logger.info(msg=f"{url} successfully downloaded.")

    if len(needed_files) > 0:
        errored_files = list(map(lambda f: f.relative_path, needed_files))
        format_error_list = "\n* ".join(errored_files)
        logger.error(msg="Maximum download attempts reached.")
        raise RuntimeError(f"Failed to download files:\n* {format_error_list}\nAborting run.")


def extract_files(unzip_list, src_folder):
    extracted_list = []
    for file_name in unzip_list:
        logger.info(msg=f"Extracting {file_name}...")
        unzip_file_name = os.path.splitext(file_name)[0]
        source = f"{src_folder}/{file_name}"
        destination = f"{src_folder}/{unzip_file_name}"
        # unzip with pigz; -f to force overwrite,
        # -d for decompression. output file is the same
        # name minus the extension, input file is deleted.
        cmd_args = ["pigz", "-f", "-d", source]
        logger.info(msg=f"Command: {' '.join(cmd_args)}")
        subprocess.check_call(cmd_args)
        extracted_list.append(destination)
    return extracted_list


def multi_s3parcp_upload(file_list, src_folder, s3_base_path):
    for file in file_list:
        logger.info(msg=f"Uploading {file.relative_path} to s3...")
        # compute proper file source and destination paths
        src_path = f"{src_folder}/{file.relative_path}"
        s3_destination = f"{s3_base_path}/{file.relative_path}"
        # upload to s3 using s3parcp with checksum
        cmd_args = ["s3parcp", "--checksum", src_path, f"s3://{s3_destination}"]
        logger.info(msg=f"Command: {' '.join(cmd_args)}")
        subprocess.check_call(cmd_args)


def write_done_file(BUCKET, dated_subfolder):
    logger.info(msg="Uploading done file...")
    timestamp = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M:%S %Z")
    s3.Object(BUCKET, f"{dated_subfolder}/done").put(Body=timestamp)


def write_index_name(BUCKET, date_tag):
    logger.info(msg="Uploading index name...")
    s3.Object(BUCKET, f"{S3_BASE_PREFIX}/index_name").put(Body=date_tag)
    logger.info(msg="Index sources done!")

# --- sub functions ---


def _download_file(url, dest_folder, relative_path):
    logger.info(msg=f"Downloading file at {url}...")
    subfolder = os.path.dirname(relative_path)
    output_folder = f"{dest_folder}/{subfolder}" if subfolder else dest_folder
    # -P for directory prefix, -cnv to continue a file if already
    # present and to turn verbose mode off
    cmd_args = ["wget", "-P", output_folder, "-cnv", url]
    logger.info(msg=f"Command: {' '.join(cmd_args)}")
    subprocess.check_call(cmd_args)


if __name__ == '__main__':
    # parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--date", help="manually assign a date tag")
    args = parser.parse_args()
    date = args.date
    if not date:
        date = os.environ.get('INDEX_NAME')
    main(date)
