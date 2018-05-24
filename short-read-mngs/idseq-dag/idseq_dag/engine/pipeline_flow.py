import json
import sys
import os
import multiprocess

class PipelineFlow:
  def __init__(self, lazy_run,
               nodes, steps, head_nodes,
               output_dir_s3,
               output_dir_local='/mnt/idseq/results',
               ref_dir_local='/mnt/idseq/ref'):
      '''
        # Here is an example of nodes/steps representation
        # Simplified for illustration
        # nodes: specify the list of files
        # steps: the PipelineStep s
        # head_nodes: list of input files that are given  (not required but nth)
        INPUT_DIR_S3 = 's3://idseq-samples-production/1/10/fastqs'
        FASTQS = ['input1.fastq', 'input2.fastq']
        STAR_OUT = ['unmapped1.fq', 'unmapped2.fq']
        DEDUP_OUT = ['dedup1.fa', 'dedup2.fa']
        BOWTIE_OUT = ['bowtie1.fa', 'bowtie2.fa']
        GSNAP_OUT = ['gsnap.m8']
        RAPSEARCH_OUT = ['rapsearch.m8']
        TAXON_COUNT_OUT = ['web_sample.json']
        nodes = [FASTQS, STAR_OUT, DEDUP_OUT, BOWTIE_OUT, GSNAP_OUT, RAPSEARCH_OUT, TAXON_COUNT_OUT]
        head_nodes = [(FASTQS, INPUT_DIR_S3)]
        steps = [ {in: [FASTQS], out: STAR_OUT, step: 'PipelineStepRunStar', additional_files: ['s3://idseq-database/wawawa'], additional_attributes: {} },
                  {in: [STAR_OUT], out: DEDUP_OUT, step: 'PipelineStepRunPriceSeq', additional_files: ['s3://idseq-database/wawawa'], additional_attributes: {} },
                  {in: [DEDUP_OUT], out: BOWTIE_OUT, step: 'PipelineStepRunBowtie2', additional_files: ['s3://idseq-database/wawawa'], additional_attributes: {} },
                  {in: [BOWTIE_OUT], out:  GSNAP_OUT, step: 'PipelineStepRunBowtie2', additional_files: ['s3://idseq-database/wawawa'], additional_attributes: {} },
                  {in: [BOWTIE_OUT], out:  RAPSEARCH_OUT, step: 'PipelineStepRunBowtie2', additional_files: ['s3://idseq-database/wawawa'], additional_attributes: {} },
                  {in: [GSNAP_OUT, RAPSEARCH_OUT], out: [TAXON_COUNT_OUT], step: 'PipelineStepRunBowtie2', additional_files: ['s3://idseq-database/wawawa'], additional_attributes: {} },
                 ]
    '''
    self.nodes = nodes
    self.edges = eges
    self.head_nodes = head_nodes

 def start():
   '''
          1. Traverse through the from the head nodes. Find the components that actually need to be run
          2. Come up with a large file download plan based on the components that need
              to be run
          3. start running based on the file dependency.
             At the end of the run all the nodes should be generated
   '''

