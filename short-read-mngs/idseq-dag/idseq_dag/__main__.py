#!/usr/bin/env python3

import argparse
import os
import sys
import json
import contextlib
import importlib
import threading
import logging
import subprocess
import traceback

import idseq_dag
import idseq_dag.util.s3
import idseq_dag.util.count
from idseq_dag import __version__


def main():
    from idseq_dag.engine.pipeline_flow import PipelineFlow
    import idseq_dag.util.log as log
    log.configure_logger()

    parser = argparse.ArgumentParser()
    parser.add_argument('--version', action='version',
                        version='%(prog)s {version}'.format(version=__version__))
    parser.add_argument('dag_json', help='pipeline dag in json file format')
    parser.add_argument('--no-lazy-run', dest='lazy_run', action='store_false')
    parser.add_argument('--key-path-s3', dest='key_path_s3', help='ssh key')
    parser.add_argument('--no-versioned-output', dest='versioned_output', action='store_false')

    parser.set_defaults(lazy_run=True)
    parser.set_defaults(versioned_output=True)
    args = parser.parse_args()
    if args.key_path_s3:
        os.environ["KEY_PATH_S3"] = args.key_path_s3
    try:
        flow = PipelineFlow(lazy_run=args.lazy_run,
                            dag_json=args.dag_json,
                            versioned_output=args.versioned_output)
        log.write("everything is awesome. idseq dag is valid~")
    except:
        parser.print_help()
        raise
    log.write("start executing the dag")
    flow.start()

def count_input_reads(input_files, max_fragments):
    input_read_count = idseq_dag.util.count.reads_in_group(input_files[0],
                                                           max_fragments=max_fragments)
    counts_dict = dict(fastqs=input_read_count)
    if input_read_count == len(input_files) * max_fragments:
        counts_dict["truncated"] = input_read_count
    with open("fastqs.count", "w") as count_file:
        json.dump(counts_dict, count_file)

def run_step():
    root_logger = logging.getLogger()
    root_logger.setLevel(level=logging.INFO)
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(idseq_dag.util.log.JsonFormatter())
    root_logger.addHandler(stream_handler)

    parser = argparse.ArgumentParser()
    parser.add_argument('--workflow-name', required=True)
    parser.add_argument('--step-module', required=True)
    parser.add_argument('--step-class', required=True)
    parser.add_argument('--step-name', required=True)
    parser.add_argument('--input-files', type=json.loads, required=True)
    parser.add_argument('--output-files', type=json.loads, required=True)
    parser.add_argument('--output-dir-s3', required=True)
    parser.add_argument('--additional-files', type=json.loads, default={})
    parser.add_argument('--additional-attributes', type=json.loads, default={})
    args = parser.parse_args()

    # If the fasta/fastq input is unpaired, we get an empty string placeholder for it. Remove it.
    for target in args.input_files:
        if target[-1] == "":
            del target[-1]

    logging.info("idseq-dag %s running %s.%s", __version__, args.step_module, args.step_class)
    idseq_dag.util.s3.config["REF_DIR"] = os.getcwd()
    step = importlib.import_module(args.step_module)
    step_class = getattr(step, args.step_class)
    step_instance = step_class(
        name=args.step_name,
        input_files=[[None] * len(target) for target in args.input_files],
        output_files=args.output_files,
        output_dir_local=os.getcwd(),
        output_dir_s3=args.output_dir_s3,
        ref_dir_local=idseq_dag.util.s3.config["REF_DIR"],
        additional_files=args.additional_files,
        additional_attributes=args.additional_attributes,
        step_status_local=args.workflow_name + "_status.json",
        step_status_lock=contextlib.suppress()
    )
    step_instance.input_files_local = args.input_files

    with open(step_instance.step_status_local, "w") as status_file:
        json.dump(dict(), status_file)

    try:
        if args.step_class == "PipelineStepRunValidateInput":
            count_input_reads(input_files=step_instance.input_files_local,
                              max_fragments=step_instance.additional_attributes["truncate_fragments_to"])

        step_instance.update_status_json_file("running")
        step_instance.validate_input_files()
        with open(f"{args.step_name}.description.md", "wb") as outfile:
            # write step_description (which subclasses may generate dynamically) to local file
            outfile.write(step_instance.step_description().encode("utf-8"))
        step_instance.run()
        step_instance.count_reads()
        step_instance.save_counts()
        step_instance.update_status_json_file("uploaded")
    except Exception as e:
        try:
            # process exception for status reporting
            s = "user_errored" if isinstance(e, idseq_dag.engine.pipeline_step.InvalidInputFileError) else "pipeline_errored"
            step_instance.update_status_json_file(s)
        except Exception:
            logging.error("Failed to update status to '%s'", s)
        traceback.print_exc()
        exit(json.dumps(dict(wdl_error_message=True, error=type(e).__name__, cause=str(e), step_description_md=step_instance.step_description())))


if __name__ == "__main__":
    main()
