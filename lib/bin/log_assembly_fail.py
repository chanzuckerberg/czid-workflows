import json
import os

WARNINGS_LOG = "assembly/spades/warnings.log"
SPADES_LOG = "assembly/spades/spades.log"

failure_info = {
    "log_present": False,
    "warnings_present": False,
    "warnings": None,
    "stack_trace": None,
    "errors": None,
}

if os.path.isfile(WARNINGS_LOG):
    failure_info["warnings_present"] = True
    with open(WARNINGS_LOG) as warnings_file:
        failure_info["warnings"] = "".join(warnings_file.readlines())

if os.path.isfile(SPADES_LOG):
    failure_info["log_present"] = True
    with open(SPADES_LOG) as log_file:
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

with open("spades_failure.json", "w") as output_file:
    output_file.write(json.dumps(failure_info))
