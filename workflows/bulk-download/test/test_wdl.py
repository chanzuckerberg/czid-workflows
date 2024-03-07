import os
import zipfile
from test_util import WDLTestCase


def relpath(*args):
    """ helper to get a filepath relative to the current file """
    return os.path.join(os.path.dirname(__file__), *args)


class TestBulkDownloads(WDLTestCase):
    """Tests for bulk downloads"""
    wdl = relpath("..", "run.wdl")
    files = [relpath("host_filter_1.fastq"), relpath("host_filter_2.fastq")]

    def testConcatenate(self):
        """
        Test file concatenation
        """
        # Calculate expected concatenation
        concat_expected = ""
        for file in self.files:
            with open(file) as fp:
                concat_expected += fp.read()

        # Run task
        result = self.run_miniwdl(task="concatenate", task_input={"files": self.files})

        # Validate
        with open(result["outputs"]["concatenate.file"]) as fp_concat_observed:
            concat_observed = fp_concat_observed.read()
        self.assertEqual(concat_observed, concat_expected)

    def testZip(self):
        """
        Test zipping files into 1 .zip file
        """
        # Run task
        result = self.run_miniwdl(task="zip", task_input={"files": self.files})

        # Validate
        archive = zipfile.ZipFile(result["outputs"]["zip.file"], "r")
        for file in self.files:
            contents = archive.read(os.path.basename(file))
            with open(file) as fp:
                self.assertEqual(contents.decode("utf-8"), fp.read())
