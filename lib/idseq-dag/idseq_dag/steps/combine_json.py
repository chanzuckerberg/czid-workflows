import json
from idseq_dag.engine.pipeline_step import PipelineStep
import idseq_dag.util.command as command
import idseq_dag.util.count as count

class PipelineStepCombineJson(PipelineStep):
    '''
    Combine json from multiple places
    '''
    def run(self):
        output_file = self.output_files_local()[0]
        self.combine_json(self.input_files_local, output_file)

    def count_reads(self):
        pass

    @staticmethod
    def combine_json(input_json_files, output_json_path):
        output_struct = None
        data_type = None

        for input_file_name in input_json_files:
            with open(input_file_name, 'r') as input_file:
                data = json.load(input_file)
                if not output_struct:
                    output_struct = data
                    data_type = type(data)
                else:
                    if data_type == list:
                        output_struct += data
                    elif data_type == dict:
                        output_struct.update(data)
                    else:
                        raise RuntimeError(f"combine json can't combine data_type {data_type}")
        with open(output_json_path, 'w') as output_file:
            json.dump(output_struct, output_file)
