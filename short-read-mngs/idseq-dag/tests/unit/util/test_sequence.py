import unittest
from unittest.mock import patch, call

# Class under test
from idseq_dag.util.sequence import chunks

class TestSequence(unittest.TestCase):
    '''Tests for `util/sequence.py`'''

    def test_chunks_list(self):
        '''Chunking a list produces a sequence of arrays'''
        l = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        chunk_size = 3

        r = list(chunks(l, chunk_size))

        self.assertEqual(
            r,
            [
                [1, 2, 3],
                [4, 5, 6],
                [7, 8, 9],
                [10]
            ]
        )

    def test_chunks_str(self):
        '''Chunking a string produces a sequence of strings'''
        l = "not really a big string"
        chunk_size = 5

        r = list(chunks(l, chunk_size))

        self.assertEqual(
            r,
            [
                'not r',
                'eally',
                ' a bi',
                'g str',
                'ing'
            ]
        )

    def test_chunks_empty_list(self):
        '''Chunking an empty list produces an empty array'''
        l = []
        chunk_size = 4

        r = list(chunks(l, chunk_size))

        self.assertEqual(
            r,
            []
        )

    def test_chunks_small_list(self):
        '''Chunking a list smaller than chunk size produces a sequence containing a single array element'''
        l = [1, 2, 3]
        chunk_size = 4

        r = list(chunks(l, chunk_size))

        self.assertEqual(
            r,
            [
                [1, 2, 3]
            ]
        )

