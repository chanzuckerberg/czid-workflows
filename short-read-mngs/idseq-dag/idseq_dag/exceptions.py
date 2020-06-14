class InvalidInputFileError(Exception):
    def __init__(self, json):
        super().__init__()
        self.json = json


class InvalidOutputFileError(Exception):
    def __init__(self, json):
        super().__init__()
        self.json = json
