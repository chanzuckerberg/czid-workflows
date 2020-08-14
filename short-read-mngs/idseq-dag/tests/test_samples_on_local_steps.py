import os
import unittest

from idseq_dag.steps.combine_taxon_counts import PipelineStepCombineTaxonCounts
from idseq_dag.steps.generate_alignment_viz import \
    PipelineStepGenerateAlignmentViz
from idseq_dag.steps.generate_annotated_fasta import \
    PipelineStepGenerateAnnotatedFasta
from idseq_dag.steps.generate_taxid_fasta import PipelineStepGenerateTaxidFasta
from idseq_dag.steps.generate_taxid_locator import \
    PipelineStepGenerateTaxidLocator
from idseq_dag.steps.run_bowtie2 import PipelineStepRunBowtie2
from idseq_dag.steps.run_cdhitdup import PipelineStepRunCDHitDup
from idseq_dag.steps.run_gsnap_filter import PipelineStepRunGsnapFilter
from idseq_dag.steps.run_lzw import PipelineStepRunLZW
from idseq_dag.steps.run_priceseq import PipelineStepRunPriceSeq
from idseq_dag.steps.run_star import PipelineStepRunStar
from idseq_dag.steps.run_subsample import PipelineStepRunSubsample
from tests import test_utils
import idseq_dag.util.command as command


class TestSamplesOnLocalSteps(unittest.TestCase):
    def test_many_samples(self):
        dag_file = "examples/generic_test_dag.json"
        test_samples_base = "s3://idseq-samples-test/test-sets"
        samples_to_test = [
            "RR004_water_2_S23",
            # TODO: Set up the rest of these test sets
            # "OPS_006_VAP_RNA_TA_4153_P1037810_S177",
            # "X26_H2O-B_no-frag_rAPID_Nav-no-zymo_AMR-SA-532"
        ]
        output_dir_s3 = "s3://idseq-samples-test/test-run-outputs"

        for bundle in samples_to_test:
            bundle = os.path.join(test_samples_base, bundle)
            self.test_all_local_steps(dag_file, bundle, output_dir_s3)

    def test_all_local_steps(self, dag_file, test_bundle, output_dir_s3):
        step_pairs = [
            (PipelineStepRunStar, "star_out"),
            (PipelineStepRunPriceSeq, "priceseq_out"),
            (PipelineStepRunCDHitDup, "cdhitdup_out"),
            (PipelineStepRunLZW, "lzw_out"),
            (PipelineStepRunBowtie2, "bowtie2_out"),
            (PipelineStepRunSubsample, "subsampled_out"),
            (PipelineStepRunGsnapFilter, "gsnap_filter_out"),
            (PipelineStepCombineTaxonCounts, "taxon_count_out"),
            (PipelineStepGenerateAnnotatedFasta, "annotated_out"),
            (PipelineStepGenerateTaxidFasta, "taxid_fasta_out"),
            (PipelineStepGenerateTaxidLocator, "taxid_locator_out")]

        for step_class, step_name in step_pairs:
            test_utils.run_step_and_match_outputs(
                step_class, step_name, dag_file, test_bundle, output_dir_s3)

        test_utils.run_step_and_match_outputs(
            PipelineStepGenerateAlignmentViz, "alignment_viz_out", dag_file,
            test_bundle, output_dir_s3, "align_viz")

        command.execute("rm tmp*")
