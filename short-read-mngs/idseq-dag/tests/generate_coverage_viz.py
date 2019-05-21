import unittest

from idseq_dag.steps.generate_coverage_viz import PipelineStepGenerateCoverageViz

# Test various static functions from PipelineStepGenerateCoverageViz.
# All subfunctions in generate_coverage_viz_data and generate_coverage_summary_data are tested.
class GenerateCoverageViz(unittest.TestCase):
   # Basic test of generate_coverage_viz_data
    def test_generate_coverage_viz_data_basic(self):
        accession_data = {
            "ACCESSION_1": {
                "total_length": 10,
                "name": "Test Name",
                "contigs": ["CONTIG_1", "CONTIG_2"],
                "reads": ["READS_1", "READS_2"]
            }
        }

        contig_data = {
            "CONTIG_1": {
                "total_length": 6,
                "accession": "ACCESSION_1",
                "query_start": 3,
                "query_end": 4,
                "subject_start": 6,
                "subject_end": 7,
                "coverage": [1, 2, 3, 4, 5, 6],
                "prop_mismatch": 0.1,
                "num_reads": 4,
                "alignment_length": 2,
                "percent_id": 99,
                "num_mismatches": 2,
                "num_gaps": 1,
                "byterange": [1, 100],
            },
             "CONTIG_2": {
                "total_length": 6,
                "accession": "ACCESSION_1",
                "query_start": 4,
                "query_end": 4,
                "subject_start": 8,
                "subject_end": 8,
                "coverage": [2, 2, 2, 2, 2, 2],
                "prop_mismatch": 0.1,
                "num_reads": 5,
                "alignment_length": 1,
                "percent_id": 98,
                "num_mismatches": 4,
                "num_gaps": 6,
                "byterange": [100, 200],
            }
        }

        read_data = {
            "READS_1": {
                "accession": "ACCESSION_1",
                "subject_start": 2,
                "subject_end": 3,
                "prop_mismatch": 0.3,
                "alignment_length": 2,
                "percent_id": 98,
                "num_mismatches": 1,
                "num_gaps": 0,
            },
            "READS_2": {
                "accession": "ACCESSION_1",
                "subject_start": 10,
                "subject_end": 10,
                "prop_mismatch": 0.3,
                "alignment_length": 1,
                "percent_id": 97,
                "num_mismatches": 0,
                "num_gaps": 1,
            }
        }

        coverage_viz_json = PipelineStepGenerateCoverageViz.generate_coverage_viz_data(
          accession_data, contig_data, read_data, 5
        )

        self.assertTrue("ACCESSION_1" in coverage_viz_json)

        accession_obj = coverage_viz_json["ACCESSION_1"]

        self.assertEqual(accession_obj["total_length"], 10)
        self.assertEqual(accession_obj["name"], "Test Name")

        # generate_hit_group_json is tested further elsewhere.
        self.assertEqual(accession_obj["hit_groups"], [
            [0, 1, 0, 2, 3, 2, .98, 1, 0, 1, []],
            [1, 0, 4, 6, 7, 2, .99, 2, 1, 3, [[1, 100]]],
            [1, 0, 5, 8, 8, 1, .98, 4, 6, 3, [[100, 200]]],
            [0, 1, 0, 10, 10, 1, .97, 0, 1, 4, []]
        ])

        # calculate_accession_coverage is tested further elsewhere.
        self.assertEqual(accession_obj["coverage"], [
            [0, 0.5, 0.5, 0, 1],
            [1, 0.5, 0.5, 0, 1],
            [2, 1.5, 0.5, 1, 0],
            [3, 3.0, 1, 2, 0],
            [4, 0.5, 0.5, 0, 1],
        ])
        self.assertEqual(accession_obj["coverage_bin_size"], 2)

        self.assertEqual(accession_obj["max_aligned_length"], 2)
        self.assertEqual(accession_obj["coverage_depth"], 1.2)
        self.assertEqual(accession_obj["coverage_breadth"], 0.6)
        self.assertEqual(accession_obj["avg_prop_mismatch"], 0.2)


   # Basic test of generate_coverage_viz_summary_data
    def test_generate_coverage_viz_summary_data_basic(self):
        accession_data = {
            "ACCESSION_1": {
                "reads": ["READS_1"],
                "contigs": ["CONTIGS_1"],
                "name": "Test Accession 1",
                "score": 100,
            },
            "ACCESSION_2": {
                "reads": ["READS_2", "READS_3"],
                "contigs": [],
                "name": "Test Accession 2",
                "score": 200,
            }
        }

        taxon_data = {
            "TAXON_1": {
              "accessions": ["ACCESSION_1", "ACCESSION_2"],
              "num_total_accessions": 100
            },
        }

        coverage_viz_obj = {
            "ACCESSION_1": {
                "coverage_depth": 10,
            },
            "ACCESSION_2": {
                "coverage_depth": 20,
            }
        }

        taxon_data_json = PipelineStepGenerateCoverageViz.generate_coverage_viz_summary_data(
            taxon_data, accession_data, coverage_viz_obj
        )

        self.assertEqual(len(taxon_data_json.keys()), 1)
        self.assertEqual(taxon_data_json["TAXON_1"], {
            "best_accessions": [
                {
                    "id": "ACCESSION_1",
                    "name": "Test Accession 1",
                    "num_contigs": 1,
                    "num_reads": 1,
                    "score": 100,
                    "coverage_depth": 10
                },
                {
                    "id": "ACCESSION_2",
                    "name": "Test Accession 2",
                    "num_contigs": 0,
                    "num_reads": 2,
                    "score": 200,
                    "coverage_depth": 20
                }
            ],
            "num_accessions": 100
        })


    # Basic test for calculate_accession_coverage where everything aligns nicely.
    def test_calculate_accession_coverage_basic(self):
        accession_data = {
            "total_length": 4,
            "contigs": ["CONTIG_1"],
            "reads": []
        }

        contig_data = {
            "CONTIG_1": {
                "total_length": 4,
                "accession": "ACCESSION_1",
                "query_start": 1,
                "query_end": 4,
                "subject_start": 1,
                "subject_end": 4,
                "coverage": [2] * 4,
            }
        }

        read_data = {}

        (coverage, bin_size) = PipelineStepGenerateCoverageViz.calculate_accession_coverage(
            "ACCESSION_1", accession_data, contig_data, read_data, 4
        )

        self.assertEqual(coverage, [
            [0, 2.0, 1, 1, 0],
            [1, 2.0, 1, 1, 0],
            [2, 2.0, 1, 1, 0],
            [3, 2.0, 1, 1, 0],
        ])
        self.assertEqual(bin_size, 1.0)


    # Test where contig is slightly smaller than accession.
    def test_calculate_accession_coverage_contig_small(self):
        accession_data = {
            "total_length": 10,
            "contigs": ["CONTIG_1"],
            "reads": []
        }

        contig_data = {
            "CONTIG_1": {
                "total_length": 9,
                "accession": "ACCESSION_1",
                "query_start": 1,
                "query_end": 9,
                "subject_start": 1,
                "subject_end": 10,
                "coverage": [2] * 5,
            }
        }

        read_data = {}

        (coverage, bin_size) = PipelineStepGenerateCoverageViz.calculate_accession_coverage(
            "ACCESSION_1", accession_data, contig_data, read_data, 10
        )

        self.assertEqual(coverage, [
            [0, 2.0, 1, 1, 0],
            [1, 2.0, 1, 1, 0],
            [2, 2.0, 1, 1, 0],
            [3, 2.0, 1, 1, 0],
            [4, 2.0, 1, 1, 0],
            [5, 2.0, 1, 1, 0],
            [6, 2.0, 1, 1, 0],
            [7, 2.0, 1, 1, 0],
            [8, 2.0, 1, 1, 0],
            [9, 2.0, 1, 1, 0],
        ])
        self.assertEqual(bin_size, 1.0)


    # Test where accession is slightly smaller than contig.
    def test_calculate_accession_coverage_accession_small(self):
        accession_data = {
            "total_length": 9,
            "contigs": ["CONTIG_1"],
            "reads": []
        }

        contig_data = {
            "CONTIG_1": {
                "total_length": 10,
                "accession": "ACCESSION_1",
                "query_start": 1,
                "query_end": 10,
                "subject_start": 1,
                "subject_end": 9,
                "coverage": [2] * 10,
            }
        }

        read_data = {}

        (coverage, bin_size) = PipelineStepGenerateCoverageViz.calculate_accession_coverage(
            "ACCESSION_1", accession_data, contig_data, read_data, 9
        )

        self.assertEqual(coverage, [
            [0, 2.0, 1, 1, 0],
            [1, 2.0, 1, 1, 0],
            [2, 2.0, 1, 1, 0],
            [3, 2.0, 1, 1, 0],
            [4, 2.0, 1, 1, 0],
            [5, 2.0, 1, 1, 0],
            [6, 2.0, 1, 1, 0],
            [7, 2.0, 1, 1, 0],
            [8, 2.0, 1, 1, 0],
        ])
        self.assertEqual(bin_size, 1.0)


    # Test where contig covers only part of accession, and coverage array is slightly small.
    def test_calculate_accession_coverage_partial_contig(self):
        accession_data = {
            "total_length": 10,
            "contigs": ["CONTIG_1"],
            "reads": []
        }

        contig_data = {
            "CONTIG_1": {
                "total_length": 7,
                "accession": "ACCESSION_1",
                "query_start": 1,
                "query_end": 7,
                "subject_start": 3,
                "subject_end": 9,
                "coverage": [2] * 6, # Coverage is 6 instead of 7.
            }
        }

        read_data = {}

        (coverage, bin_size) = PipelineStepGenerateCoverageViz.calculate_accession_coverage(
            "ACCESSION_1", accession_data, contig_data, read_data, 10
        )

        self.assertEqual(coverage, [
            [2, 2.0, 1, 1, 0],
            [3, 2.0, 1, 1, 0],
            [4, 2.0, 1, 1, 0],
            [5, 2.0, 1, 1, 0],
            [6, 2.0, 1, 1, 0],
            [7, 2.0, 1, 1, 0],
            [8, 2.0, 1, 1, 0],
        ])
        self.assertEqual(bin_size, 1.0)


    # Test where contig covers only part of accession and contig has varying coverage.
    def test_calculate_accession_coverage_offset_varying(self):
        accession_data = {
            "total_length": 10,
            "contigs": ["CONTIG_1"],
            "reads": []
        }

        contig_data = {
            "CONTIG_1": {
                "total_length": 5,
                "accession": "ACCESSION_1",
                "query_start": 1,
                "query_end": 5,
                "subject_start": 3,
                "subject_end": 9,
                "coverage": [1, 2, 3, 4, 5]
            }
        }

        read_data = {}

        (coverage, bin_size) = PipelineStepGenerateCoverageViz.calculate_accession_coverage(
            "ACCESSION_1", accession_data, contig_data, read_data, 10
        )

        self.assertEqual(coverage, [
            [2, 1.0, 1, 1, 0],
            [3, 1.5, 1, 1, 0],
            [4, 2.5, 1, 1, 0],
            [5, 3.0, 1, 1, 0],
            [6, 3.5, 1, 1, 0],
            [7, 4.5, 1, 1, 0],
            [8, 5.0, 1, 1, 0],
        ])
        self.assertEqual(bin_size, 1.0)


    # Previous test, but the alignment is reversed.
    def test_calculate_accession_coverage_offset_varying_reverse(self):
        accession_data = {
            "total_length": 10,
            "contigs": ["CONTIG_1"],
            "reads": []
        }

        contig_data = {
            "CONTIG_1": {
                "total_length": 5,
                "accession": "ACCESSION_1",
                "query_start": 1,
                "query_end": 5,
                "subject_start": 9,
                "subject_end": 3,
                "coverage": [1, 2, 3, 4, 5]
            }
        }

        read_data = {}

        (coverage, bin_size) = PipelineStepGenerateCoverageViz.calculate_accession_coverage(
            "ACCESSION_1", accession_data, contig_data, read_data, 10
        )

        self.assertEqual(coverage, [
            [2, 5.0, 1, 1, 0],
            [3, 4.5, 1, 1, 0],
            [4, 3.5, 1, 1, 0],
            [5, 3.0, 1, 1, 0],
            [6, 2.5, 1, 1, 0],
            [7, 1.5, 1, 1, 0],
            [8, 1.0, 1, 1, 0],
        ])
        self.assertEqual(bin_size, 1.0)


    # Test where accession size is slightly larger than max_bins.
    def test_calculate_accession_coverage_max_bins(self):
        accession_data = {
            "total_length": 10,
            "contigs": ["CONTIG_1"],
            "reads": []
        }

        contig_data = {
            "CONTIG_1": {
                "total_length": 5,
                "accession": "ACCESSION_1",
                "query_start": 1,
                "query_end": 5,
                "subject_start": 1,
                "subject_end": 10,
                "coverage": [2.0] * 5
            }
        }

        read_data = {}

        (coverage, bin_size) = PipelineStepGenerateCoverageViz.calculate_accession_coverage(
            "ACCESSION_1", accession_data, contig_data, read_data, 9
        )

        self.assertEqual(coverage, [
            [0, 2.0, 1, 1, 0],
            [1, 2.0, 1, 1, 0],
            [2, 2.0, 1, 1, 0],
            [3, 2.0, 1, 1, 0],
            [4, 2.0, 1, 1, 0],
            [5, 2.0, 1, 1, 0],
            [6, 2.0, 1, 1, 0],
            [7, 2.0, 1, 1, 0],
            [8, 2.0, 1, 1, 0],
        ])

        self.assertTrue((bin_size - 10 / 9) < 0.01)


    # Test where accession size is slightly larger than max_bins and contig has varying coverage.
    def test_calculate_accession_coverage_max_bins_varying(self):
        accession_data = {
            "total_length": 10,
            "contigs": ["CONTIG_1"],
            "reads": []
        }

        contig_data = {
            "CONTIG_1": {
                "total_length": 5,
                "accession": "ACCESSION_1",
                "query_start": 1,
                "query_end": 5,
                "subject_start": 1,
                "subject_end": 10,
                "coverage": [1, 2, 3, 4, 5]
            }
        }

        read_data = {}

        (coverage, bin_size) = PipelineStepGenerateCoverageViz.calculate_accession_coverage(
            "ACCESSION_1", accession_data, contig_data, read_data, 9
        )

        self.assertEqual(coverage, [
            [0, 1.0, 1, 1, 0],
            [1, 1.5, 1, 1, 0],
            [2, 2.0, 1, 1, 0],
            [3, 2.5, 1, 1, 0],
            [4, 3.0, 1, 1, 0],
            [5, 3.5, 1, 1, 0],
            [6, 4.0, 1, 1, 0],
            [7, 4.5, 1, 1, 0],
            [8, 5.0, 1, 1, 0],
        ])
        self.assertTrue((bin_size - 10 / 9) < 0.01)


    # Test where max bins is much smaller than accession.
    # This allows us to verify that the coverage is being reduced when the contig doesn't cover the whole bin.
    def test_calculate_accession_coverage_max_bins_small(self):
        accession_data = {
            "total_length": 10,
            "contigs": ["CONTIG_1"],
            "reads": []
        }

        contig_data = {
            "CONTIG_1": {
                "total_length": 6,
                "accession": "ACCESSION_1",
                "query_start": 1,
                "query_end": 4,
                "subject_start": 4,
                "subject_end": 7,
                "coverage": [1, 2, 3, 4, 5, 6]
            }
        }

        read_data = {}

        (coverage, bin_size) = PipelineStepGenerateCoverageViz.calculate_accession_coverage(
            "ACCESSION_1", accession_data, contig_data, read_data, 5
        )

        # The coverage is reduced (e.g. divided by 2) for bins with 0.5 and 2.0.
        # Since the contig only covers half of that bin.
        self.assertEqual(coverage, [
            [1, 0.5, 0.5, 1, 0],
            [2, 2.5, 1, 1, 0],
            [3, 2.0, 0.5, 1, 0],
        ])
        self.assertEqual(bin_size, 2)


    # Test multiple reads.
    def test_calculate_accession_coverage_multiple_reads(self):
        accession_data = {
            "total_length": 10,
            "contigs": [],
            "reads": ["READS_1", "READS_2"]
        }

        contig_data = {}

        read_data = {
            "READS_1": {
                "accession": "ACCESSION_1",
                "subject_start": 4,
                "subject_end": 7,
            },
            "READS_2": {
                "accession": "ACCESSION_1",
                "subject_start": 6,
                "subject_end": 9,
            }
        }

        (coverage, bin_size) = PipelineStepGenerateCoverageViz.calculate_accession_coverage(
            "ACCESSION_1", accession_data, contig_data, read_data, 10
        )

        self.assertEqual(coverage, [
            [3, 1.0, 1, 0, 1],
            [4, 1.0, 1, 0, 1],
            [5, 2.0, 1, 0, 2],
            [6, 2.0, 1, 0, 2],
            [7, 1.0, 1, 0, 1],
            [8, 1.0, 1, 0, 1],
        ])
        self.assertEqual(bin_size, 1)


    # Test where max bins is much smaller than contig.
    # This allows us to verify that the coverage is being reduced when the reads and contigs don't cover the whole bin.
    def test_calculate_accession_coverage_max_bins_small_reads(self):
        accession_data = {
            "total_length": 10,
            "contigs": ["CONTIG_1"],
            "reads": ["READS_1", "READS_2"]
        }

        contig_data = {
            "CONTIG_1": {
                "total_length": 6,
                "accession": "ACCESSION_1",
                "query_start": 3,
                "query_end": 4,
                "subject_start": 6,
                "subject_end": 7,
                "coverage": [1, 2, 3, 4, 5, 6]
            }
        }

        read_data = {
            "READS_1": {
                "accession": "ACCESSION_1",
                "subject_start": 2,
                "subject_end": 3,
            },
            "READS_2": {
                "accession": "ACCESSION_1",
                "subject_start": 10,
                "subject_end": 10,
            }
        }

        (coverage, bin_size) = PipelineStepGenerateCoverageViz.calculate_accession_coverage(
            "ACCESSION_1", accession_data, contig_data, read_data, 5
        )

        # [0.5, 0.5] from READS_1.
        # [1.5, 2.0] from CONTIG_1
        # last [0.5] from READS_2.
        self.assertEqual(coverage, [
            [0, 0.5, 0.5, 0, 1],
            [1, 0.5, 0.5, 0, 1],
            [2, 1.5, 0.5, 1, 0],
            [3, 2.0, 0.5, 1, 0],
            [4, 0.5, 0.5, 0, 1],
        ])
        self.assertEqual(bin_size, 2)


    # This test was created in response to a division by zero bug.
    # This bug is caused by a rounding error when calculating the bins in the accession that a contig overlaps.
    # This test verifies that this bug is no longer occurring.
    def test_calculate_accession_coverage_rounding_error(self):
        accession_data = {
            "total_length": 545,
            "contigs": ["CONTIG_1"],
            "reads": [],
        }

        contig_data = {
            "CONTIG_1": {
                "total_length": 400,
                "accession": "ACCESSION_1",
                "query_start": 320,
                "query_end": 322,
                "subject_start": 221,
                "subject_end": 219,
                "coverage": [1] * 400
            }
        }

        read_data = {}

        (coverage, bin_size) = PipelineStepGenerateCoverageViz.calculate_accession_coverage(
            "ACCESSION_1", accession_data, contig_data, read_data, 500
        )

        self.assertEqual(coverage, [
            [200, 1.0, 1.0, 1, 0],
            [201, 1.0, 1.0, 1, 0],
            [202, 0.75, 0.752, 1, 0]
        ])
        self.assertEqual(bin_size, 1.09)


    # Tests for calculate_accession_stats.
    def test_calculate_accession_stats_basic(self):
        accession_data = {
            "total_length": 10,
            "contigs": ["CONTIG_1"],
            "reads": []
        }

        contig_data = {
            "CONTIG_1": {
                "total_length": 10,
                "accession": "ACCESSION_1",
                "query_start": 1,
                "query_end": 10,
                "subject_start": 1,
                "subject_end": 10,
                "coverage": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                "prop_mismatch": 0.1
            }
        }

        read_data = {}

        stats = PipelineStepGenerateCoverageViz.calculate_accession_stats(
          accession_data, contig_data, {}
        )

        self.assertEqual(stats["max_aligned_length"], 10)
        self.assertEqual(stats["coverage_depth"], 5.5)
        self.assertEqual(stats["coverage_breadth"], 1)
        self.assertEqual(stats["avg_prop_mismatch"], 0.1)


    # Add multiple contigs and reads
    def test_calculate_accession_stats_multiple(self):
        accession_data = {
            "total_length": 100,
            "contigs": ["CONTIG_1", "CONTIG_2"],
            "reads": ["READ_1", "READ_2"]
        }

        contig_data = {
            "CONTIG_1": {
                "total_length": 100,
                "accession": "ACCESSION_1",
                "query_start": 1,
                "query_end": 40,
                "subject_start": 1,
                "subject_end": 40,
                "coverage": [2] * 40,
                "prop_mismatch": 0.01
            },
            "CONTIG_2": {
                "total_length": 100,
                "accession": "ACCESSION_1",
                "query_start": 1,
                "query_end": 20,
                "subject_start": 100,
                "subject_end": 81,
                "coverage": [4] * 20,
                "prop_mismatch": 0.02
            }
        }

        read_data = {
            "READ_1": {
                "total_length": 20,
                "accession": "ACCESSION_1",
                "query_start": 1,
                "query_end": 20,
                "subject_start": 1,
                "subject_end": 20,
                "prop_mismatch": 0.03
            },
            "READ_2": {
                "total_length": 20,
                "accession": "ACCESSION_1",
                "query_start": 1,
                "query_end": 20,
                "subject_start": 100,
                "subject_end": 81,
                "prop_mismatch": 0.04
            }
        }

        stats = PipelineStepGenerateCoverageViz.calculate_accession_stats(
          accession_data, contig_data, read_data
        )

        self.assertEqual(stats["max_aligned_length"], 40)
        self.assertEqual(stats["coverage_depth"], 2)
        self.assertEqual(stats["coverage_breadth"], 0.6)
        self.assertEqual(stats["avg_prop_mismatch"], 0.025)


    # Basic test of generate_hit_group_json
    def test_generate_hit_group_json_basic(self):
        accession_id = "ACCESSION_1"

        accession_data = {
            "total_length": 100,
            "name": "Test Name",
            "contigs": ["CONTIG_1", "CONTIG_2"],
            "reads": ["READS_1", "READS_2", "READS_3", "READS_4"]
        }

        contig_data = {
            "CONTIG_1": {
                "total_length": 25,
                "accession": "ACCESSION_1",
                "query_start": 1,
                "query_end": 25,
                "subject_start": 40,
                "subject_end": 64,
                "coverage": [2] * 25,
                "prop_mismatch": 0.1,
                "num_reads": 22,
                "alignment_length": 25,
                "percent_id": 100,
                "num_mismatches": 6,
                "num_gaps": 3,
                "byterange": [1, 100],
            },
             "CONTIG_2": {
                "total_length": 9,
                "accession": "ACCESSION_1",
                "query_start": 1,
                "query_end": 9,
                "subject_start": 85,
                "subject_end": 93,
                "coverage": [2, 2, 2, 2, 2, 2, 2, 2, 2],
                "num_reads": 15,
                "alignment_length": 9,
                "percent_id": 88,
                "num_mismatches": 4,
                "num_gaps": 5,
                "byterange": [100, 200],
            }
        }

        read_data = {
            "READS_1": {
                "accession": "ACCESSION_1",
                "subject_start": 1,
                "subject_end": 30,
                "alignment_length": 30,
                "percent_id": 98,
                "num_mismatches": 1,
                "num_gaps": 0,
            },
            "READS_2": {
                "accession": "ACCESSION_1",
                "subject_start": 69,
                "subject_end": 75,
                "alignment_length": 7,
                "percent_id": 97,
                "num_mismatches": 0,
                "num_gaps": 1,
            },
            "READS_3": {
                "accession": "ACCESSION_1",
                "subject_start": 83,
                "subject_end": 89,
                "alignment_length": 7,
                "percent_id": 90,
                "num_mismatches": 3,
                "num_gaps": 1,
            },
            "READS_4": {
                "accession": "ACCESSION_1",
                "subject_start": 90,
                "subject_end": 86,
                "alignment_length": 5,
                "percent_id": 95,
                "num_mismatches": 2,
                "num_gaps": 0,
            }
        }

        hit_group_json = PipelineStepGenerateCoverageViz.generate_hit_group_json(
          accession_data, accession_id, contig_data, read_data, 10
        )

        self.assertEqual(len(hit_group_json), 4)

        self.assertEqual(hit_group_json, [
            [0, 1, 0, 1, 30, 30, .98, 1, 0, 1, []],
            [1, 0, 22, 40, 64, 25, 1, 6, 3, 5, [[1, 100]]],
            [0, 1, 0, 69, 75, 7, .97, 0, 1, 7, []],
            # CONTIG_2, READ_3, and READ_4 are grouped together and their stats are averaged.
            [1, 2, 15, 83, 93, 7, .91, 3, 2, 8, [[100, 200]]],
        ])


    # Test that generate_hit_group_json aggregates contigs
    def test_generate_hit_group_json_contigs(self):
        accession_id = "ACCESSION_1"

        accession_data = {
           "total_length": 10,
            "name": "Test Name",
            "contigs": ["CONTIG_1", "CONTIG_2"],
            "reads": [],
        }


        contig_data = {
            "CONTIG_1": {
                "total_length": 6,
                "accession": "ACCESSION_1",
                "query_start": 3,
                "query_end": 4,
                "subject_start": 7,
                "subject_end": 7,
                "coverage": [1, 2, 3, 4, 5, 6],
                "prop_mismatch": 0.1,
                "num_reads": 4,
                "alignment_length": 2,
                "percent_id": 99,
                "num_mismatches": 2,
                "num_gaps": 1,
                "byterange": [1, 100],
            },
             "CONTIG_2": {
                "total_length": 6,
                "accession": "ACCESSION_1",
                "query_start": 4,
                "query_end": 4,
                "subject_start": 8,
                "subject_end": 8,
                "coverage": [2, 2, 2, 2, 2, 2],
                "prop_mismatch": 0.1,
                "num_reads": 5,
                "alignment_length": 1,
                "percent_id": 98,
                "num_mismatches": 4,
                "num_gaps": 6,
                "byterange": [100, 200],
            }
        }

        hit_group_json = PipelineStepGenerateCoverageViz.generate_hit_group_json(
          accession_data, accession_id, contig_data, [], 5
        )

        self.assertEqual(len(hit_group_json), 1)

        self.assertEqual(hit_group_json, [
            [2, 0, 9, 7, 8, 1.5, .985, 3, 3.5, 3, [[1, 100], [100, 200]]],
        ])


    # Basic test of remove_taxons_with_no_contigs
    def test_remove_taxons_with_no_contigs_basic(self):
        accession_data = {
            "ACCESSION_1": {
                "reads": ["READS_1"],
                "contigs": ["CONTIGS_1"]
            },
            "ACCESSION_2": {
                "reads": ["READS_2"],
                "contigs": []
            },
            "ACCESSION_3": {
                "reads": ["READS_3"],
                "contigs": []
            },
            "ACCESSION_4": {
                "reads": ["READS_4"],
                "contigs": []
            }
        }

        taxon_data = {
            "TAXON_1": {
              "accessions": ["ACCESSION_1", "ACCESSION_2"],
              "num_total_accessions": 2
            },
            "TAXON_2": {
              "accessions": ["ACCESSION_3", "ACCESSION_4"],
              "num_total_accessions": 2
            }
        }

        PipelineStepGenerateCoverageViz.remove_taxons_with_no_contigs(
          accession_data, taxon_data
        )

        self.assertEqual(len(accession_data.keys()), 2)
        self.assertTrue("ACCESSION_1" in accession_data)
        self.assertTrue("ACCESSION_2" in accession_data)

        self.assertEqual(len(taxon_data.keys()), 1)
        self.assertTrue("TAXON_1" in taxon_data)


    # Basic test of get_unassigned_reads_set
    def test_get_unassigned_reads_set_basic(self):
        accession_data = {
            "ACCESSION_1": {
                "reads": ["READS_1", "READS_3", "READS_4"],
                "contigs": ["CONTIGS_1"]
            },
            "ACCESSION_2": {
                "reads": ["READS_2"],
                "contigs": []
            }
        }

        unassigned_reads_set = PipelineStepGenerateCoverageViz.get_unassigned_reads_set(
          accession_data
        )

        self.assertEqual(len(unassigned_reads_set), 4)
        self.assertTrue("READS_1" in unassigned_reads_set)
        self.assertTrue("READS_2" in unassigned_reads_set)
        self.assertTrue("READS_3" in unassigned_reads_set)
        self.assertTrue("READS_4" in unassigned_reads_set)



    # Test that select best accessions correctly filters to top accessions.
    def test_select_best_accessions_per_taxon_basic(self):
        taxon_data = {
            "TAXON_1": {
              "accessions": ["ACCESSION_1", "ACCESSION_2", "ACCESSION_3"],
              "num_total_accessions": 3
            },
            "TAXON_2": {
              "accessions": ["ACCESSION_4", "ACCESSION_5", "ACCESSION_6"],
              "num_total_accessions": 3
            },
            "TAXON_3": {
              "accessions": ["ACCESSION_7"],
              "num_total_accessions": 1
            },
        }

        contig_data = {
            "CONTIGS_1": { "alignment_length": 100 },
            "CONTIGS_2": { "alignment_length": 100 },
            "CONTIGS_3": { "alignment_length": 100 },
            "CONTIGS_4": { "alignment_length": 120 },
            "CONTIGS_5": { "alignment_length": 120 },
            "CONTIGS_6": { "alignment_length": 300 },
            "CONTIGS_7": { "alignment_length": 200 },
        }

        read_data = {
            "READS_1": { "alignment_length": 150 },
            "READS_2": { "alignment_length": 150 },
            "READS_3": { "alignment_length": 150 },
            "READS_4": { "alignment_length": 150 },
        }

        accession_data = {
            "ACCESSION_1": {
                "reads": [],
                "contigs": ["CONTIGS_1", "CONTIGS_2", "CONTIGS_3"],
            },
            "ACCESSION_2": {
                "reads": [],
                "contigs": ["CONTIGS_4", "CONTIGS_5"],
            },
            "ACCESSION_3": {
                "reads": [],
                "contigs": ["CONTIGS_6"],
            },
            "ACCESSION_4": {
                "reads": ["READS_1"],
                "contigs": [],
            },
            "ACCESSION_5": {
                "reads": [],
                "contigs": ["CONTIGS_7"],
            },
            "ACCESSION_6": {
                "reads": ["READS_2", "READS_3"],
                "contigs": [],
            },
            "ACCESSION_7": {
                "reads": ["READS_4"],
                "contigs": [],
            },
        }

        (new_taxon_data, new_accession_data) = PipelineStepGenerateCoverageViz.select_best_accessions_per_taxon(
            taxon_data, accession_data, contig_data, read_data, 2
        )

        self.assertEqual(new_taxon_data, {
            "TAXON_1": {
              "accessions": ["ACCESSION_3", "ACCESSION_1", "ACCESSION_2"],
              "num_total_accessions": 3
            },
            "TAXON_2": {
              "accessions": ["ACCESSION_5", "ACCESSION_6"],
              "num_total_accessions": 3
            },
            "TAXON_3": {
              "accessions": ["ACCESSION_7"],
              "num_total_accessions": 1
            },
        })

        self.assertFalse("ACCESSION_4" in new_accession_data)

        self.assertTrue(new_accession_data["ACCESSION_3"]["score"] > new_accession_data["ACCESSION_1"]["score"])
        self.assertTrue(new_accession_data["ACCESSION_1"]["score"] > new_accession_data["ACCESSION_2"]["score"])
        self.assertTrue(new_accession_data["ACCESSION_5"]["score"] > new_accession_data["ACCESSION_6"]["score"])


    # Test that select best accessions correctly filters to top accessions.
    def test_select_best_accessions_per_taxon_high_read_count(self):
        taxon_data = {
            "TAXON_1": {
              "accessions": ["ACCESSION_1", "ACCESSION_2", "ACCESSION_3"],
              "num_total_accessions": 3
            },
            "TAXON_2": {
              "accessions": ["ACCESSION_4", "ACCESSION_5"],
              "num_total_accessions": 2
            },
        }

        contig_data = {
            "CONTIGS_1": { "alignment_length": 1 },
            "CONTIGS_2": { "alignment_length": 1 },
            "CONTIGS_3": { "alignment_length": 1 },
            "CONTIGS_4": { "alignment_length": 1 },
        }

        read_data = {
            "READS_1": { "alignment_length": 1 },
            "READS_2": { "alignment_length": 1 },
            "READS_3": { "alignment_length": 1 },
            "READS_4": { "alignment_length": 1 },
            "READS_5": { "alignment_length": 1 },
            "READS_6": { "alignment_length": 1 },
        }

        accession_data = {
            "ACCESSION_1": {
                "reads": ["READS_1", "READS_2", "READS_3"],
                "contigs": [],
            },
            "ACCESSION_2": {
                "reads": [],
                "contigs": ["CONTIGS_1"],
            },
            "ACCESSION_3": {
                "reads": [],
                "contigs": ["CONTIGS_2"],
            },
            "ACCESSION_4": {
                "reads": ["READS_4", "READS_5", "READS_6"],
                "contigs": ["CONTIGS_4"],
            },
            "ACCESSION_5": {
                "reads": [],
                "contigs": ["CONTIGS_3"],
            },
        }

        (new_taxon_data, new_accession_data) = PipelineStepGenerateCoverageViz.select_best_accessions_per_taxon(
            taxon_data, accession_data, contig_data, read_data, 2
        )

        self.assertEqual(new_taxon_data, {
            "TAXON_1": {
              "accessions": ["ACCESSION_2", "ACCESSION_3"],
              "num_total_accessions": 3
            },
            "TAXON_2": {
              "accessions": ["ACCESSION_4", "ACCESSION_5"],
              "num_total_accessions": 2
            },
        })

        self.assertFalse("ACCESSION_1" in new_accession_data)

        self.assertTrue(new_accession_data["ACCESSION_2"]["score"] == new_accession_data["ACCESSION_3"]["score"])
        self.assertTrue(new_accession_data["ACCESSION_4"]["score"] > new_accession_data["ACCESSION_5"]["score"])
