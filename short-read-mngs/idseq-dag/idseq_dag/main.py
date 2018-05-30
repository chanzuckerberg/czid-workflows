import argparse
import json
import sys
import os
import traceback
import idseq_dag.util.s3
from idseq_dag.engine.pipeline_flow import PipelineFlow

def validate_dag_json(json_file):
    '''
    Validate the json format. see examples/*.json.
    Required fields are:
      "output_dir_s3": base results folder. a pipeline version number will be appended for real output folder.
      "nodes": lists of files that are given or would be generated
      "steps": steps that species actions to generate input and output
      "head_nodes": input files that are given

    '''
    try:
        dag = json.loads(open(json_file).read())
        output_dir = dag["output_dir_s3"]
        nodes = dag["nodes"]
        steps = dag["steps"]
        head_nodes = dag["head_nodes"]
        covered_nodes = set()
        for s in steps:
            # validate each step in/out are valid nodes
            output = nodes[s["out"]]
            covered_nodes.add(s["out"])
            for inode in s["in"]:
                input_n = nodes[inode]
        for hn in head_nodes:
            # validate the head nodes exist in s3
            node_name = hn[0]
            s3_path = hn[1]
            covered_nodes.add(node_name)
            for file_name in nodes[node_name]:
                s3_file = os.path.join(s3_path, file_name)
                if not idseq_dag.util.s3.check_s3_presence(s3_file):
                    print("%s file doesn't exist" % s3_file)
                    return False
        # Check that all nodes are covered
        # ALL Inputs Outputs VALIDATED
        for node_name in nodes.keys():
            if node_name not in covered_nodes:
                print("%s couldn't be generated from the steps" % node_name)
                return False
        return dag
    except:
        traceback.print_exc()
        return False

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('dag_json', help='pipeline dag in json file format')
    args = parser.parse_args()
    dag_json_file = args.dag_json
    dag = validate_dag_json(dag_json_file)
    if not dag:
        parser.print_help()
        sys.exit(1)
    print("everything is awesome. idseq dag is awesome")
    flow = PipelineFlow(lazy_run=True, nodes = dag["nodes"], steps = dag["steps"],
                        head_nodes = dag["head_nodes"], output_dir_s3 = dag["output_dir_s3"])
    flow.start()

