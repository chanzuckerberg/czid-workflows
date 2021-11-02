import json
import traceback


class IDseqDagError(Exception):
    def __init__(self, json):
        super().__init__()
        self.json = json

    def __str__(self):
        return str(self.json)


class InvalidInputFileError(IDseqDagError):
    pass


class InvalidOutputFileError(IDseqDagError):
    pass


class InsufficientReadsError(InvalidInputFileError):
    pass


class BrokenReadPairError(InvalidInputFileError):
    pass


class InvalidFileFormatError(InvalidInputFileError):
    pass


def print_exceptions(f):
    try:
        f()
    except Exception as e:
        traceback.print_exc()
        exit(json.dumps(dict(wdl_error_message=True, error=type(e).__name__, cause=str(e))))
