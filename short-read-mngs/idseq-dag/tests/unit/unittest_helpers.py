import re
import os

class MATCH_RE(object):
    """
    An assertion helper object that compares value to a regex.
    Usage:
      mock.assert_has_calls([
          call(MATCH_RE('ssh')),
          call(MATCH_RE(regex=re.compile('ssh'))
      ])
    """
    def __init__(self, pattern=None, regex=None):
        self.pattern = pattern
        if regex is not None:
            self.regex = regex
        else:
            self.regex = re.compile(pattern)

    def __eq__(self, other):
        return self.regex.match(other)

    def __ne__(self, other):
        return not self.regex.match(other)

    def __repr__(self):
        if self.pattern is not None:
            return f'<MATCH_RE({self.pattern})>'
        return f'<MATCH_RE(regex={self.regex})>'


def file_contents(file_path):
    with open(file_path, 'r') as f:
        return f.read()


def relative_file_path(file_path, relative_path):
    p1 = os.path.abspath(file_path)
    if os.path.isfile(file_path):
        p1 = os.path.dirname(p1)
    return os.path.abspath(os.path.join(p1, relative_path))
