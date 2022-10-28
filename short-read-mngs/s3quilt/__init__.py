from typing import Iterable

from .go import Slice_uint64
from .s3quilt import DownloadChunks, DownloadChunksToFile

"""
A note on concurrency: The process is mostly IO bound and go uses green threads so you want high concurrency,
higher than the number of cores.
"""

def download_chunks(
    bucket: str,
    key: str,
    starts: Iterable[int],
    lengths: Iterable[int],
    concurrency: int=100,
):
    return DownloadChunks(bucket, key, Slice_uint64(starts), Slice_uint64(lengths), concurrency)


def download_chunks_to_file(
    bucket: str,
    key: str,
    filepath: str,
    starts: Iterable[int],
    lengths: Iterable[int],
    concurrency: int=100,
):
    return DownloadChunksToFile(
        bucket,
        key,
        filepath,
        Slice_uint64(starts),
        Slice_uint64(lengths),
        concurrency,
    )
