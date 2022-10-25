from typing import Iterable
from .go import Slice_uint64
from . import s3quilt


def download_chunks(bucket: str, key: str, output_file_path: str, starts: Iterable[int], lengths: Iterable[int]):
    s3quilt.DownloadChunks(bucket, key, output_file_path, Slice_uint64(starts), Slice_uint64(lengths))
