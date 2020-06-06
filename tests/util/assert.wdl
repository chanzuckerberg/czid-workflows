version 1.0

task IsTrue {
    input {
        Boolean cond
        String msg
    }

    command <<<
        if [[ '~{cond}' != 'true' ]]; then
            # see miniwdl_run_with_custom_error.py for error reporting convention
            >&2 jq -nc --arg XXX "`<~{write_lines([msg])}`" '{"wdl_error_message": $XXX}'
            exit 1
        fi
    >>>

    runtime {
        docker: "stedolan/jq"
    }
}
