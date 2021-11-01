import json
import traceback

def print_exceptions(f):
    try:
        f()
    except Exception as e:
        traceback.print_exc()
        exit(json.dumps(dict(wdl_error_message=True, error=type(e).__name__, cause=str(e))))
