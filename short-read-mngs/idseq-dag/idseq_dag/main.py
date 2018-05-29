import argparse
import json
import sys
import os
import traceback
import idseq_dag.util.s3

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
        for s in steps:
            # validate each step in/out are valid nodes
            output = nodes[s["out"]]
            for inode in s["in"]:
                input_n = nodes[inode]
        for hn in head_nodes:
            # validate the head nodes exist in s3
            node_name = hn[0]
            s3_path = hn[1]
            for file_name in nodes[node_name]:
                s3_file = os.path.join(s3_path, file_name)
                if not idseq_dag.util.s3.check_s3_presence(s3_file):
                    print("%s file doesn't exist" % s3_file)
                    return False
        # ALL Inputs Outputs VALIDATED
        return True
    except:
        traceback.print_exc()
        return False

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('dag_json', help='pipeline dag in json file format')
    args = parser.parse_args()
    dag_json_file = args.dag_json
    if not validate_dag_json(dag_json_file):
        parser.print_help()
        sys.exit(1)
    print("everything is awesome. idseq dag is awesome")

