import json
import sys
import os
import multiprocess
import threading
import idseq_dag
import importlib

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
        '''
            Traverse through the nodes and steps and calculate
            1. the large file download priority based on the how deep the step is
            2. if a step needs to be run based on the existence of output file and lazy run parameter
        '''
        covered_nodes = {}
        large_file_download_list = []
        step_list = []
        for hn in self.head_nodes:
            covered_nodes[hn[0]] = { 'round': 0, 'lazy_run': self.lazy_run, 's3_downlodable': True }
        steps_complete = set()
        while len(steps_complete) < len(self.steps):
            # run until all the steps can be run
            current_nodes = {}
            for step in self.steps:
                if step["step_id"] not in steps_complete:
                    step_can_be_run = True
                    round_max = 0
                    lazy_run = True
                    for node in step["in"]:
                        if node not in covered_nodes:
                            step_can_be_run = False
                            break
                        else:
                            round_max = max(covered_nodes[node]['round'], round_max)
                            if covered_nodes[node]['lazy_run'] == False:
                                    lazy_run = False
                    if step_can_be_run: # All the input is satisfied
                        steps_complete.add(step["step_id"])
                        file_list= self.nodes[step["out"]]
                        if lazy_run and idseq_dag.util.s3.check_s3_presnce_for_file_list(self.output_dir_s3, file_list):
                            # output can be lazily generated. touch the output
                            #idseq_dag.util.s3.touch_s3_file_list(self.output_dir_s3, file_list)
                            s3_downlodable = True
                        else:
                            # steps need to be run
                            lazy_run = False
                            s3_downlodable = False
                            step_list.append(step)
                            # The following can be changed to append if we want to get the round information
                            large_file_download_list += step["additional_files"].values()
                        # update nodes available for the next round
                        current_nodes[step["out"]] = { 'round': (round_max + 1), 'lazy_run': lazy_run, 's3_downlodable': s3_downlodable}
            covered_nodes.update(current_nodes)
        return (step_list, large_file_download_list, covered_nodes)

    def fetch_node_from_s3(node):
        ''' To Be Implemented. a .done file should be written to the result dir when the download is complete '''

    def start(self):
        # Come up with the plan
        (step_list, large_file_download_list, covered_nodes) = self.plan()

        for step in step_list: # download the files from s3 when necessary
            for node in step["in"]:
                node_info = covered_nodes[node]
                if node_info['s3_downloadable']:
                    threading.Thread(target=self.fetch_node_from_s3, args=(node,)).start()
        # use  fetch_from_s3 plus threading for the large_file_donwload for necessary steps
        for f in large_file_download_list:
            threading.Thread(target=fetch_from_s3, args=(f, self.output_dir_local,)).start()

        # Start initializing all the steps and start running them and wait until all of them are done
        step_instances = []
        for step in step_list:
            StepClass = getattr(importlib.import_module(step["module"]), step["step"])
            step_instance = StepClass(step["in"], step["out"],
                                      self.output_dir_local, self.output_dir_s3, self.ref_dir_local
                                      step["additional_files"], step["additional_attributes"])
            step_instance.start()
            step_instances.append(step_instance)
