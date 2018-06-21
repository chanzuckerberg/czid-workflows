import json

import os
import time

from idseq_dag.engine.pipeline_flow import PipelineFlow
import idseq_dag.util.command as command


class IdseqStepSetup(object):
    @staticmethod
    def get_step_object(step_class, step_name, paired=True):
        ''' return the PipelineStep with the default parameters ready for test. '''
        if paired:
            dag = IdseqStepSetup.paired_dag()
        else:
            dag = IdseqStepSetup.single_dag()
        step_info = {}
        for step in dag["steps"]:
            if step["out"] == step_name:
                step_info = step
                break
        if not step_info:
            raise ValueError("no steps correspond to %s" % step_name)

        # Download input data to local
        output_dir_s3 = os.path.join(dag["output_dir_s3"], "testrun_%s_%d_%d" % (step_name, int(paired),int(time.time())))
        result_dir_local = "/mnt/idseq/results/%s_%d/%d" % (step_name, int(paired), os.getpid())
        ref_dir_local = '/mnt/idseq/ref'
        command.execute("mkdir -p %s %s" % (result_dir_local, ref_dir_local))

        input_files = []
        for target in step_info["in"]:
            if target in dag["given_targets"]:
                input_dir_s3 = dag["given_targets"][target]["s3_dir"]
            else:
                input_dir_s3 = dag["output_dir_s3"]
            input_files.append(dag["targets"][target])
            PipelineFlow.fetch_input_files_from_s3(input_files[-1],
                                                   input_dir_s3,
                                                   result_dir_local)

        return step_class(
            name=step_name,
            input_files=input_files,
            output_files=dag["targets"][step_info["out"]],
            output_dir_local=result_dir_local,
            output_dir_s3=output_dir_s3,
            ref_dir_local=ref_dir_local,
            additional_files=step_info["additional_files"],
            additional_attributes=step_info["additional_attributes"]
        )

    @staticmethod
    def paired_dag():
        """Return a test dag based on paired fastqs"""
        with open("examples/paired_dag.json") as f:
            return json.load(f)

    @staticmethod
    def single_dag():
        """Return a test dag based on single fastq"""
        with open("examples/single_dag.json") as f:
            return json.load(f)
