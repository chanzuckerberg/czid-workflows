import shelve
from typing import Union
import marisa_trie


class WrappedTrie:
    def __init__(self, trie: marisa_trie.RecordTrie, stringify: bool):
        self.trie = trie
        self.stringify = stringify

    def _extract(self, raw):
        if self.stringify:
            head = [elem.decode() if isinstance(elem, bytes) else str(elem) for elem in raw[0]]
        else:
            head = list(raw[0])
        if len(head) == 1:
            return head[0]
        return head

    def __getitem__(self, key: str):
        return self._extract(self.trie[str(key)])

    def __contains__(self, key):
        return key in self.trie

    def __enter__(self):
        return self

    def __exit__(self, type, value, tb):
        pass

    def get(self, name, default=None):
        r = self.trie.get(str(name))
        if r:
            return self._extract(r)
        return default


def open_file_db_by_extension(db_path: str, fmt: Union[str, None] = None, stringify=True):
    if db_path.endswith(".db"):
        return shelve.open(db_path.replace('.db', ''), 'r')
    assert fmt, "fmt is required for loading marisa trie key value stores"
    return WrappedTrie(marisa_trie.RecordTrie(fmt).mmap(db_path), stringify)
