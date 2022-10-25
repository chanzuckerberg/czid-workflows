import os
import sys
from typing import Iterable

sys.path.append(os.path.dirname(os.path.realpath(__file__)))

from .go import Slice_uint64
from .s3quilt import DownloadChunks


def download_chunks(bucket: str, key: str, output_file_path: str, starts: Iterable[int], lengths: Iterable[int]):
    DownloadChunks(bucket, key, output_file_path, Slice_uint64(starts), Slice_uint64(lengths))
