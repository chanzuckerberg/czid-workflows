from typing import Iterable

from .go import Slice_uint64
from .s3quilt import DownloadChunks

import boto3


client = boto3.client('s3')


def download_chunks(bucket: str, key: str, output_file_path: str, starts: Iterable[int], lengths: Iterable[int]):
    DownloadChunks(client.meta.region_name, bucket, key, output_file_path, Slice_uint64(starts), Slice_uint64(lengths))
