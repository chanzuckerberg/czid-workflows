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
        step_info =  {}
        for step in dag["steps"]:
            if step["out"] == step_name:
                step_info = step
                break
        if not step_info:
            raise ValueError("no steps correspond to %s" % step_name)

        # Download input data to local
        output_dir_s3 = os.path.join(dag["output_dir_s3"], "testrun_bowtie2_%d_%d" % (int(paired),int(time.time())))
        result_dir_local = "/mnt/idseq/results/bowtie2_%d/%d" % (int(paired), os.getpid())
        command.execute("mkdir -p %s" % result_dir_local)

        input_files = []
        for node in step_info["in"]:
            if node in dag["head_nodes"]:
                input_dir_s3 = dag["head_nodes"]["s3_dir"]
            else:
                input_dir_s3 = dag["output_dir_s3"]
            input_files.append(dag["nodes"][node])
            PipelineFlow.fetch_input_files_from_s3(input_files[-1],
                                                   input_dir_s3,
                                                   result_dir_local)

        return step_class(
            name=step_name,
            input_files=input_files,
            output_files=dag["nodes"][step_info["out"]],
            output_dir_local=result_dir_local,
            output_dir_s3=output_dir_s3,
            ref_dir_local='/mnt/idseq/ref',
            additional_files=step_info["additional_files"],
            additional_attributes=step_info["additional_attributes"]
        )

    @staticmethod
    def paired_dag():
        ''' Return a test dag based on paired fastqs '''
        return json.loads('''
 {
  "output_dir_s3": "s3://idseq-samples-prod/test_samples/1/results",
  "nodes": {
    "fastqs": ["RR004_water_2_S23_R1_001.fastq.gz", "RR004_water_2_S23_R2_001.fastq.gz"],
    "star_out": ["unmapped.star.1.fq", "unmapped.star.2.fq"],
    "priceseq_out": ["priceseqfilter.unmapped.star.1.fasta", "priceseqfilter.unmapped.star.2.fasta"],
    "cdhitdup_out": ["cdhitdup.priceseqfilter.unmapped.star.1.fasta", "cdhitdup.priceseqfilter.unmapped.star.2.fasta"],
    "lzw_out": ["lzw.cdhitdup.priceseqfilter.unmapped.star.1.fasta", "lzw.cdhitdup.priceseqfilter.unmapped.star.2.fasta"],
    "bowtie_out": ["unmapped.bowtie2.lzw.cdhitdup.priceseqfilter.unmapped.star.1.fasta",
                   "unmapped.bowtie2.lzw.cdhitdup.priceseqfilter.unmapped.star.2.fasta",
                   "unmapped.bowtie2.lzw.cdhitdup.priceseqfilter.unmapped.star.merged.fasta"],
    "gsnap_filter_out": ["unmapped.gsnap_filter.bowtie2.lzw.cdhitdup.priceseqfilter.unmapped.star.1.fasta",
                         "unmapped.gsnap_filter.bowtie2.lzw.cdhitdup.priceseqfilter.unmapped.star.2.fasta",
                         "unmapped.gsnap_filter.bowtie2.lzw.cdhitdup.priceseqfilter.merged.star.1.fasta"]
  },

  "steps": [
    {
      "in" : ["fastqs"], "out": "star_out", "class": "PipelineStepRunStar", "module": "idseq_dag.steps.run_star",
      "additional_files": {"star_genome": "s3://idseq-database/host_filter/human/2018-02-15-utc-1518652800-unixtime__2018-02-15-utc-1518652800-unixtime/STAR_genome.tar"},
      "additional_attributes": {"truncate_reads_to": 75000000}
    },
    {
      "in" : ["star_out"], "out": "priceseq_out", "class": "PipelineStepRunPriceSeq", "module": "idseq_dag.steps.run_priceseq",
      "additional_files": {},
      "additional_attributes": {}
    },
    {
      "in" : ["priceseq_out"], "out": "cdhitdup_out", "class": "PipelineStepRunCDHitDup", "module": "idseq_dag.steps.run_cdhitdup",
      "additional_files": {},
      "additional_attributes": {}
    },
    {
      "in" : ["cdhitdup_out"], "out": "lzw_out", "class": "PipelineStepRunLZW", "module": "idseq_dag.steps.run_lzw",
      "additional_files": {},
      "additional_attributes": {"threshold": 0.45}
    },
    {
      "in" : ["lzw_out"], "out": "bowtie_out", "class": "PipelineStepRunBowtie2", "module": "idseq_dag.steps.run_bowtie2",
      "additional_files": {"bowtie2_genome": "s3://idseq-database/host_filter/human/2018-02-15-utc-1518652800-unixtime__2018-02-15-utc-1518652800-unixtime/bowtie2_genome.tar"},
      "additional_attributes": {}
    },
    {
      "in" : ["bowtie_out"], "out": "gsnap_filter_out", "class": "PipelineStepRunGsnapFilter", "module": "idseq_dag.steps.run_gsnap_filter",
      "additional_files": {"gsnap_genome": "s3://idseq-database/host_filter/human/2018-02-15-utc-1518652800-unixtime__2018-02-15-utc-1518652800-unixtime/hg38_pantro5_k16.tar"},
      "additional_attributes": {"truncate_reads_to": 10000000}
    }
  ],
  "head_nodes": {"fastqs": {"s3_dir":  "s3://idseq-samples-prod/test_samples/1/fastqs", "max_reads":75000000 } }
}
        ''')

    @staticmethod
    def single_dag():
        ''' Return a test dag based on single fastq '''
        return json.loads('''
 {
  "output_dir_s3": "s3://idseq-samples-prod/test_samples/1/results",
  "nodes": {
    "fastqs": ["RR004_water_2_S23_R1_001.fastq.gz"],
    "star_out": ["unmapped.star.1.fq"],
    "priceseq_out": ["priceseqfilter.unmapped.star.1.fasta"],
    "cdhitdup_out": ["cdhitdup.priceseqfilter.unmapped.star.1.fasta"],
    "lzw_out": ["lzw.cdhitdup.priceseqfilter.unmapped.star.1.fasta"],
    "bowtie_out": ["unmapped.bowtie2.lzw.cdhitdup.priceseqfilter.unmapped.star.1.fasta",
                   "unmapped.bowtie2.lzw.cdhitdup.priceseqfilter.unmapped.star.merged.fasta"],
    "gsnap_filter_out": ["unmapped.gsnap_filter.bowtie2.lzw.cdhitdup.priceseqfilter.unmapped.star.1.fasta",
                         "unmapped.gsnap_filter.bowtie2.lzw.cdhitdup.priceseqfilter.merged.star.1.fasta"]
  },

  "steps": [
    {
      "in" : ["fastqs"], "out": "star_out", "class": "PipelineStepRunStar", "module": "idseq_dag.steps.run_star",
      "additional_files": {"star_genome": "s3://idseq-database/host_filter/human/2018-02-15-utc-1518652800-unixtime__2018-02-15-utc-1518652800-unixtime/STAR_genome.tar"},
      "additional_attributes": {"truncate_reads_to": 75000000}
    },
    {
      "in" : ["star_out"], "out": "priceseq_out", "class": "PipelineStepRunPriceSeq", "module": "idseq_dag.steps.run_priceseq",
      "additional_files": {},
      "additional_attributes": {}
    },
    {
      "in" : ["priceseq_out"], "out": "cdhitdup_out", "class": "PipelineStepRunCDHitDup", "module": "idseq_dag.steps.run_cdhitdup",
      "additional_files": {},
      "additional_attributes": {}
    },
    {
      "in" : ["cdhitdup_out"], "out": "lzw_out", "class": "PipelineStepRunLZW", "module": "idseq_dag.steps.run_lzw",
      "additional_files": {},
      "additional_attributes": {"threshold": 0.45}
    },
    {
      "in" : ["lzw_out"], "out": "bowtie_out", "class": "PipelineStepRunBowtie2", "module": "idseq_dag.steps.run_bowtie2",
      "additional_files": {"bowtie2_genome": "s3://idseq-database/host_filter/human/2018-02-15-utc-1518652800-unixtime__2018-02-15-utc-1518652800-unixtime/bowtie2_genome.tar"},
      "additional_attributes": {}
    },
    {
      "in" : ["bowtie_out"], "out": "gsnap_filter_out", "class": "PipelineStepRunGsnapFilter", "module": "idseq_dag.steps.run_gsnap_filter",
      "additional_files": {"gsnap_genome": "s3://idseq-database/host_filter/human/2018-02-15-utc-1518652800-unixtime__2018-02-15-utc-1518652800-unixtime/hg38_pantro5_k16.tar"},
      "additional_attributes": {"truncate_reads_to": 10000000}
    }
  ],
  "head_nodes": {"fastqs": {"s3_dir":  "s3://idseq-samples-prod/test_samples/1/fastqs", "max_reads":75000000 } }

  }
        ''')

