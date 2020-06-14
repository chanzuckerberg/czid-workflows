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
