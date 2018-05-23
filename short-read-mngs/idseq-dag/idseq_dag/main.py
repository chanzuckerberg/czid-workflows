import argparse
import json
import sys

def validate_dag_json(json_file):
    return True
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('dag_json', help='pipeline dag in json file format')
    args = parser.parse_args()
    dag_json_file = args.dag_json
    if not validate_dag_json(dag_json_file):
        parser.print_help()
        sys.exit(1)
    print("idseq dag is awesome")

