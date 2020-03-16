import json
import os
from idseq_dag.engine.pipeline_step import PipelineStep
from idseq_dag.util.count import DAG_SURGERY_HACKS_FOR_READ_COUNTING
import idseq_dag.util.log as log

class PipelineStepCombineTaxonCounts(PipelineStep):
    '''
    Combine counts from gsnap and rapsearch
    '''
    def run(self):
        input_files = []
        for target in self.input_files_local:
            input_files.append(target[3])
        output_file = self.output_files_local()[0]
        self.combine_counts(input_files, output_file)

        # TODO:  Remove this hack as soon as the webapp has been updated so that
        # the regular inputs and outputs of this step are "_with_dcr.json".
        with_dcr = all(input_f.endswith("_with_dcr.json") for input_f in input_files)
        if not with_dcr:
            assert DAG_SURGERY_HACKS_FOR_READ_COUNTING
            input_files = [input_f.replace(".json", "_with_dcr.json") for input_f in input_files]
            output_file = output_file.replace(".json", "_with_dcr.json")
            if all(os.path.isfile(f) for f in input_files):
                self.combine_counts(input_files, output_file)
                self.additional_output_files_visible.append(output_file)
            else:
                # This is fine -- it means the assembler crashed, so we are not going to emit
                # these files in pipeline v3;  in pipeline v4 this code path will be gone
                log.write("Post-assembly combined taxon counts with DCR unavailable because assembler crashed.  See pre-assembly combined count.  This problem will diseappear after idseq-web fix to emit correct DAG json, long before v4.")

    def count_reads(self):
        pass

    @staticmethod
    def combine_counts(input_json_files, output_json_path):
        taxon_counts = []
        for input_file_name in input_json_files:
            with open(input_file_name, 'r') as input_file:
                data = json.load(input_file)
                taxon_counts += data["pipeline_output"]["taxon_counts_attributes"]
        output_dict = {"pipeline_output": {"taxon_counts_attributes": taxon_counts}}
        with open(output_json_path, 'w') as output_file:
            json.dump(output_dict, output_file)
