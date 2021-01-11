import unittest
from tempfile import TemporaryFile

from idseq_dag.util.parsing import BlastnOutput6Reader, BlastnOutput6Writer
from idseq_dag.util.parsing import BlastnOutput6NTReader, BlastnOutput6NTWriter
from idseq_dag.util.parsing import BlastnOutput6NTRerankedReader, BlastnOutput6NTRerankedWriter
from idseq_dag.util.parsing import HitSummaryReader, HitSummaryWriter
from idseq_dag.util.parsing import HitSummaryMergedReader, HitSummaryMergedWriter

from tests.unit.unittest_helpers import relative_file_path

# These tests are fresh as of 2021-01-11, when we delete a lot of tests to make them work again KEEP THESE (please)


class TestBlastnOutput6Reader(unittest.TestCase):
    def test_read(self):
        blastn_input_6 = [
            "# a comment",  # a comment
            "1	MK468611.1	100.0	126	0	0	1	126	8433	8308	1.1e-74	290.5",
            "2	MK468611.1	90.0	126	0	0	1	126	8433	8308		290.5", # missing evalue
        ]
        rows = list(BlastnOutput6Reader(blastn_input_6))
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0]["pident"], 100.0)
        self.assertEqual(rows[1]["evalue"], "")

    def test_read_error_empty_line(self):
        blastn_input_6 = [""]
        with self.assertRaises(Exception):
            next(BlastnOutput6Reader(blastn_input_6))

    def test_read_error_too_many_columns(self):
        blastn_input_6 = [
            "# a comment",  # a comment
            "1	MK468611.1	100.0	126	0	0	1	126	8433	8308	1.1e-74	290.5	1",
        ]
        with self.assertRaises(Exception):
            next(BlastnOutput6Reader(blastn_input_6))
    
    def test_read_error_wrong_data_type(self):
        blastn_input_6 = [
            "# a comment",  # a comment
            "1	MK468611.1	not_num	126	0	0	1	126	8433	8308	1.1e-74	290.5",
        ]
        with self.assertRaises(Exception):
            next(BlastnOutput6Reader(blastn_input_6))

    def test_filtration(self):
        blastn_input_6 = [
            "# a comment",  # a comment
            "1	MK468611.1	100.0	126	0	0	1	126	8433	8308	1.1e-74	290.5",
            "2	MK468611.1	135.0	126	0	0	1	126	8433	8308	1.1e-74	290.5",  # pident too high
            "3 	MK468611.1	-0.25	126	0	0	1	126	8433	8308	1.1e-74	290.5",  # pident too low
            "4	MK468611.1	-0.25	126	0	0	1	126	8433	8308	NaN	290.5",  # NaN error
            "5	MK468611.1	-0.25	126	0	0	1	126	8433	8308	1	290.5",  # error too high
        ]

        rows = list(BlastnOutput6Reader(blastn_input_6, filter_invalid=True))
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["qseqid"], "1")


class TestBlastnOutput6Writer(unittest.TestCase):
    def test_write(self):
        with TemporaryFile("w") as f:
            writer = BlastnOutput6Writer(f)
            writer.writerow({
                "qseqid": "1",
                "sseqid": "2",
                "pident": 1.0,
                "length": 5,
                "mismatch": 10,
                "gapopen": 15,
                "qstart": 20,
                "qend": 25,
                "sstart": 30,
                "send": 35,
                "evalue": 0.1,
                "bitscore": 1.0,
            })
            writer.writerow({
                "qseqid": "2",
            })

    def test_write_error_bad_key(self):
        with TemporaryFile("w") as f:
            writer = BlastnOutput6Writer(f)
            with self.assertRaises(Exception):
                writer.writerow({"bad key": 1})



class TestBlastnOutput6NTReader(unittest.TestCase):
    def test_read(self):
        blastn_input_6 = [
            "# a comment",  # a comment
            "1	MK468611.1	100.0	126	0	0	1	126	8433	8308	1.1e-74	290.5	11	22",
            "2	MK468611.1	90.0	126	0	0	1	126	8433	8308		290.5	11	22", # missing evalue
        ]
        rows = list(BlastnOutput6NTReader(blastn_input_6))
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0]["qlen"], 11)
        self.assertEqual(rows[1]["evalue"], "")

    def test_read_error_empty_line(self):
        blastn_input_6 = [""]
        with self.assertRaises(Exception):
            next(BlastnOutput6NTReader(blastn_input_6))

    def test_read_error_too_many_columns(self):
        blastn_input_6 = [
            "# a comment",  # a comment
            "1	MK468611.1	100.0	126	0	0	1	126	8433	8308	1.1e-74	290.5	1	2	3",
        ]
        with self.assertRaises(Exception):
            next(BlastnOutput6NTReader(blastn_input_6))
    
    def test_read_error_wrong_data_type(self):
        blastn_input_6 = [
            "# a comment",  # a comment
            "1	MK468611.1	not_num	126	0	0	1	126	8433	8308	1.1e-74	290.5	11	22",
        ]
        with self.assertRaises(Exception):
            next(BlastnOutput6NTReader(blastn_input_6))

    def test_filtration(self):
        blastn_input_6 = [
            "# a comment",  # a comment
            "1	MK468611.1	100.0	126	0	0	1	126	8433	8308	1.1e-74	290.5	11	22",
            "2	MK468611.1	135.0	126	0	0	1	126	8433	8308	1.1e-74	290.5	11	22",  # pident too high
            "3 	MK468611.1	-0.25	126	0	0	1	126	8433	8308	1.1e-74	290.5	11	22",  # pident too low
            "4	MK468611.1	-0.25	126	0	0	1	126	8433	8308	NaN	290.5	11	22",  # NaN error
            "5	MK468611.1	-0.25	126	0	0	1	126	8433	8308	1	290.5	11	22",  # error too high
        ]

        rows = list(BlastnOutput6NTReader(blastn_input_6, filter_invalid=True))
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["qseqid"], "1")


class TestBlastnOutput6NTWriter(unittest.TestCase):
    def test_write(self):
        with TemporaryFile("w") as f:
            writer = BlastnOutput6NTWriter(f)
            writer.writerow({
                "qseqid": "1",
                "sseqid": "2",
                "pident": 1.0,
                "length": 5,
                "mismatch": 10,
                "gapopen": 15,
                "qstart": 20,
                "qend": 25,
                "sstart": 30,
                "send": 35,
                "evalue": 0.1,
                "bitscore": 1.0,
                "qlen": 1,
                "slen": 2,
            })
            writer.writerow({
                "qseqid": "2",
            })

    def test_write_error_bad_key(self):
        with TemporaryFile("w") as f:
            writer = BlastnOutput6Writer(f)
            with self.assertRaises(Exception):
                writer.writerow({"bad key": 1})


class TestBlastnOutput6NTRerankedReader(unittest.TestCase):
    def test_read(self):
        blastn_input_6 = [
            "# a comment",  # a comment
            "1	MK468611.1	100.0	126	0	0	1	126	8433	8308	1.1e-74	290.5	11	22	0.1	33",
            "2	MK468611.1	90.0	126	0	0	1	126	8433	8308		290.5	11	22	0.1	33", # missing evalue
        ]
        rows = list(BlastnOutput6NTRerankedReader(blastn_input_6))
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0]["qcov"], 0.1)
        self.assertEqual(rows[1]["evalue"], "")

    def test_read_error_empty_line(self):
        blastn_input_6 = [""]
        with self.assertRaises(Exception):
            next(BlastnOutput6NTRerankedReader(blastn_input_6))

    def test_read_error_too_many_columns(self):
        blastn_input_6 = [
            "# a comment",  # a comment
            "1	MK468611.1	100.0	126	0	0	1	126	8433	8308	1.1e-74	290.5	1	2   0.1	33	3",
        ]
        with self.assertRaises(Exception):
            next(BlastnOutput6NTRerankedReader(blastn_input_6))
    
    def test_read_error_wrong_data_type(self):
        blastn_input_6 = [
            "# a comment",  # a comment
            "1	MK468611.1	not_num	126	0	0	1	126	8433	8308	1.1e-74	290.5	11	22	0.1	33",
        ]
        with self.assertRaises(Exception):
            next(BlastnOutput6NTRerankedReader(blastn_input_6))

    def test_filtration(self):
        blastn_input_6 = [
            "# a comment",  # a comment
            "1	MK468611.1	100.0	126	0	0	1	126	8433	8308	1.1e-74	290.5	11	22	0.1	33",
            "2	MK468611.1	135.0	126	0	0	1	126	8433	8308	1.1e-74	290.5	11	22	0.1	33",  # pident too high
            "3 	MK468611.1	-0.25	126	0	0	1	126	8433	8308	1.1e-74	290.5	11	22	0.1	33",  # pident too low
            "4	MK468611.1	-0.25	126	0	0	1	126	8433	8308	NaN	290.5	11	22	0.1	33",  # NaN error
            "5	MK468611.1	-0.25	126	0	0	1	126	8433	8308	1	290.5	11	22	0.1	33",  # error too high
        ]

        rows = list(BlastnOutput6NTRerankedReader(blastn_input_6, filter_invalid=True))
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["qseqid"], "1")


class TestBlastnOutput6NTRerankedWriter(unittest.TestCase):
    def test_write(self):
        with TemporaryFile("w") as f:
            writer = BlastnOutput6NTRerankedWriter(f)
            writer.writerow({
                "qseqid": "1",
                "sseqid": "2",
                "pident": 1.0,
                "length": 5,
                "mismatch": 10,
                "gapopen": 15,
                "qstart": 20,
                "qend": 25,
                "sstart": 30,
                "send": 35,
                "evalue": 0.1,
                "bitscore": 1.0,
                "qlen": 1,
                "slen": 2,
                "qcov": 0.1,
                "hsp_count": 3,
            })
            writer.writerow({
                "qseqid": "2",
            })

    def test_write_error_bad_key(self):
        with TemporaryFile("w") as f:
            writer = BlastnOutput6Writer(f)
            with self.assertRaises(Exception):
                writer.writerow({"bad key": 1})


class TestHitSummaryReader(unittest.TestCase):
    def test_read(self):
        hit_summary_input = [
            "read_id_1	10	20	accession_id_1	30	40	50",
            "read_id_2	11	21		31	41	51", # missing accession_id
        ]
        rows = list(HitSummaryReader(hit_summary_input))
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0]["level"], 10)
        self.assertEqual(rows[1]["accession_id"], "")

    def test_read_error_empty_line(self):
        hit_summary_input = [""]
        with self.assertRaises(Exception):
            next(HitSummaryReader(hit_summary_input))

    def test_read_error_too_many_columns(self):
        hit_summary_input = [
            "read_id_1	10	20	accession_id_1	30	40	50	60",
        ]
        with self.assertRaises(Exception):
            next(HitSummaryReader(hit_summary_input))
    
    def test_read_error_wrong_data_type(self):
        hit_summary_input = [
            "read_id_1	10	not_num	accession_id_1	30	40	50",
        ]
        with self.assertRaises(Exception):
            next(HitSummaryReader(hit_summary_input))


class TestHitSummaryWriter(unittest.TestCase):
    def test_write(self):
        with TemporaryFile("w") as f:
            writer = HitSummaryWriter(f)
            writer.writerow({
                "read_id": "read_id_1",
                "level": 10,
                "taxid": 20,
                "accession_id": "accession_id_1",
                "species_taxid": 30,
                "genus_taxid": 40,
                "family_taxid": 60,
            })
            writer.writerow({
                "read_id": "2",
            })

    def test_write_error_bad_key(self):
        with TemporaryFile("w") as f:
            writer = HitSummaryWriter(f)
            with self.assertRaises(Exception):
                writer.writerow({"bad key": 1})


class TestHitSummaryMergedReader(unittest.TestCase):
    def test_read(self):
        hit_summary_input = [
            "read_id_1	10	20	accession_id_1	30	40	50	contig_id_1	contig_accession_id_1	60	70	80	from_assembly	souce_count_type",
            "read_id_2	11	21		31	41	51	contig_id_2	contig_accession_id_2	61	71	81	from_assembly	souce_count_type", # missing accession_id
        ]
        rows = list(HitSummaryMergedReader(hit_summary_input))
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0]["contig_species_taxid"], 60)
        self.assertEqual(rows[1]["accession_id"], "")

    def test_read_error_empty_line(self):
        hit_summary_input = [""]
        with self.assertRaises(Exception):
            next(HitSummaryMergedReader(hit_summary_input))

    def test_read_error_too_many_columns(self):
        hit_summary_input = [
            "read_id_1	10	20	accession_id_1	30	40	50	contig_id_1	contig_accession_id_1	60	70	80	from_assembly	souce_count_type	90",
        ]
        with self.assertRaises(Exception):
            next(HitSummaryMergedReader(hit_summary_input))
    
    def test_read_error_wrong_data_type(self):
        hit_summary_input = [
            "read_id_1	10	not_num	accession_id_1	30	40	50	contig_id_1	contig_accession_id_1	60	70	80	from_assembly	souce_count_type",
        ]
        with self.assertRaises(Exception):
            next(HitSummaryMergedReader(hit_summary_input))


class TestHitSummaryMergedWriter(unittest.TestCase):
    def test_write(self):
        with TemporaryFile("w") as f:
            writer = HitSummaryMergedWriter(f)
            writer.writerow({
                "read_id": "read_id_1",
                "level": 10,
                "taxid": 20,
                "accession_id": "accession_id_1",
                "species_taxid": 30,
                "genus_taxid": 40,
                "family_taxid": 60,
                "contig_id": "contig_id_1",
                "contig_accession_id": "contig_accession_id_1",
                "contig_species_taxid": 70,
                "contig_genus_taxid": 80,
                "contig_family_taxid": 90,
                "from_assembly": "from_assembly",
                "source_count_type": "source_count_type",
            })
            writer.writerow({
                "read_id": "2",
            })

    def test_write_error_bad_key(self):
        with TemporaryFile("w") as f:
            writer = HitSummaryMergedWriter(f)
            with self.assertRaises(Exception):
                writer.writerow({"bad key": 1})

