from typing import Iterable

from .go import Slice_uint64
from .s3quilt import DownloadChunks, DownloadChunksToFile


def download_chunks(bucket: str, key: str, starts: Iterable[int], lengths: Iterable[int]):
    return DownloadChunks(bucket, key, Slice_uint64(starts), Slice_uint64(lengths))


def download_chunks_to_file(bucket: str, key: str, filepath: str, starts: Iterable[int], lengths: Iterable[int]):
    return DownloadChunksToFile(bucket, key, filepath, Slice_uint64(starts), Slice_uint64(lengths))
