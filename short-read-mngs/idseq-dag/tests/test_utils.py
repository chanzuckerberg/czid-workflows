import os
import re
import filecmp

import idseq_dag.util.command as command
import idseq_dag.util.log as log
import idseq_dag.util.s3 as s3
from tests.idseq_step_setup import IdseqStepSetup


def should_match_exactly(expected, actual):
    to_compare = [expected, actual]
    download_to_compare(to_compare)
    compare_local_files(to_compare)


def should_match_sorted_fastq(expected, actual):
    # Download the files locally
    to_compare = [expected, actual]
    download_to_compare(to_compare)

    # Sort the fastqs
    for i, path in enumerate(to_compare):
        new_name = f"tmp-sorted-{i}-" + os.path.basename(path)
        command.execute(f"cat {path} | paste - - - - | sort -k1,1 -S 3G | tr '\t' '\n' > {new_name}")
        to_compare[i] = new_name

    compare_local_files(to_compare)


def should_match_sam(expected, actual):
    contents = [None, None]
    for i, path in enumerate([expected, actual]):
        if path.startswith("s3://"):
            contents[i] = s3_file_contents(path)
        else:
            contents[i] = cat_file_contents(path)
        # Remove lines starting with @PG from SAM file comparison because these
        # will differ each time. Ex: @PG ID:bowtie2 PN:bowtie2 VN:2.3.2
        # CL:"/usr/local/bin...
        contents[i] = re.sub("@PG.*\n", "", contents[i])

    if contents[0] == contents[1]:
        log.write(f"File {expected} is semantically the same as {actual}")
    else:
        raise ValueError(f"{expected} does not match {actual}")


def s3_file_contents(path):
    return command.execute_with_output(f"aws s3 cp {path} -")


def cat_file_contents(path):
    return command.execute_with_output("cat " + path)


def download_to_compare(to_compare):
    for i, path in enumerate(to_compare):
        if path.startswith("s3://"):
            local_name = f"tmp-{i}-" + os.path.basename(path)
            if os.path.isfile(local_name):
                os.remove(local_name)
            path = s3.fetch_from_s3(path, local_name, allow_s3mi=True)
            if path is None:
                raise RuntimeError(f"Fetch from S3 failed for {path}")
            to_compare[i] = path


def compare_local_files(to_compare):
    if filecmp.cmp(to_compare[0], to_compare[1]):
        log.write(f"File {to_compare[0]} matches {to_compare[1]}")
    else:
        raise ValueError(f"{to_compare[0]} does not match {to_compare[1]}")


def run_step_and_match_outputs(step_class,
                               step_name,
                               dag_file,
                               test_bundle,
                               output_dir_s3,
                               subfolder=None):
    # Run the step
    step = IdseqStepSetup.get_test_step_object(step_class, step_name, dag_file,
                                               test_bundle, output_dir_s3)
    step.start()
    step.wait_until_finished()

    # Get lists of expected and actual files
    if step.should_count_reads:
        step.output_files.append(f"{step_name}.count")

    expected_files = [os.path.join(test_bundle, f) for f in step.output_files]
    actual_files = step.output_files_local()

    if step.additional_output_files_hidden:
        actual_files += step.additional_output_files_hidden

        for f in step.additional_output_files_hidden:
            basename = os.path.basename(f)
            # Handle sub-folders like "align_viz"
            if subfolder and subfolder in f:
                to_append = os.path.join(test_bundle, subfolder, basename)
            else:
                to_append = os.path.join(test_bundle, basename)
            expected_files.append(to_append)

    # Check that results match
    for expected, actual in zip(expected_files, actual_files):
        log.write(f"Expected test file: {expected}")
        log.write(f"Actual test file: {actual}")
        if expected.endswith(".sam"):
            fn = should_match_sam
        elif expected.endswith((".fastq", ".fq")):
            fn = should_match_sorted_fastq
        else:
            fn = should_match_exactly

        fn(expected, actual)
