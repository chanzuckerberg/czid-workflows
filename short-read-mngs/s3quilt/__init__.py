from typing import Iterable

from .go import Slice_uint64
from .s3quilt import DownloadChunks


def download_chunks(bucket: str, key: str, starts: Iterable[int], lengths: Iterable[int]):
    return DownloadChunks(bucket, key, Slice_uint64(starts), Slice_uint64(lengths))
