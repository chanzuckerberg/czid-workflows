import argparse
import json
import os


def main(spades_log: str, warnings_log: str):
    failure_info = {
        "log_present": False,
        "warnings_present": False,
        "warnings": None,
        "stack_trace": None,
        "errors": None,
    }
    failure_info = collect_spades_log(failure_info)
    failure_info = collect_warnings_log(failure_info)
    write_output(failure_info)


def collect_spades_log(failure_info: dict, spades_log: str):
    if os.path.isfile(spades_log):
        failure_info["log_present"] = True
        with open(spades_log) as log_file:
            log = log_file.readlines()
            error_lines = [line.replace("== Error ==", "*") for line in log if line.startswith("== Error ==")]
            if len(error_lines) > 0:
                failure_info["errors"] = "".join(error_lines)
            try:
                stack_trace_line_no = log.index("=== Stack Trace ===\n")
                stack_trace_end = log.index("\n", stack_trace_line_no)
                failure_info["stack_trace"] = "".join(log[stack_trace_line_no:stack_trace_end])
            except ValueError:
                pass
    return failure_info


def collect_warnings_log(failure_info: dict, warnings_log: str):
    if os.path.isfile(warnings_log):
        failure_info["warnings_present"] = True
        with open(warnings_log) as warnings_file:
            failure_info["warnings"] = "".join(warnings_file.readlines())
    return failure_info


def write_output(failure_info: dict):
    with open("spades_failure.json", "w") as output_file:
        output_file.write(json.dumps(failure_info))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Write SPAdes errors and warnings to a json file")
    parser.add_argument("spades_log", type=str, help="Path to spades.log")
    parser.add_argument("warnings_log", type=str, help="Path to warnings.log")
    args = parser.parse_args()
    main(args.spades_log, args.warnings_log)
