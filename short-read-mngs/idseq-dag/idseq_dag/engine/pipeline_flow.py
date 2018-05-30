import json
import sys
import os
import multiprocess
import idseq_dag

class PipelineFlow:
  def __init__(self, lazy_run,
               nodes, steps, head_nodes,
               output_dir_s3,
               output_dir_local='/mnt/idseq/results',
               ref_dir_local='/mnt/idseq/ref'):
    '''
        See examples/example_dag.json and
            idseq_dag.main.validate_dag_json for more details.
    '''
    self.lazy_run = lazy_run
    self.nodes = nodes
    self.steps = steps
    self.head_nodes = head_nodes
    self.output_dir_s3 = os.path.join(output_dir_s3, idseqdag.__version__)
    self.output_dir_local = output_dir_local
    self.ref_dir_local = ref_dir_local


 def plan(self):
     lazy_run = self.lazy_run

 def start(self):
   '''
          1. Traverse through the from the head nodes. Find the components that actually need to be run
          2. Come up with a large file download plan based on the components that need
              to be run
          3. start running based on the file dependency.
             At the end of the run all the nodes should be generated
   '''
   execution_plan = self.plan()

