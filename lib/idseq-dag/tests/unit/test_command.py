import unittest
from unittest.mock import patch, ANY
import os
import subprocess
import idseq_dag.util.command as command
import idseq_dag.util.command_patterns as command_patterns
from tests.unit.unittest_helpers import file_contents, relative_file_path, MATCH_RE

TMP_FILE_NAME = 'command_test_tmp_file.tmp'
TMP_FILE = f'/tmp/{TMP_FILE_NAME}'
TMP_FOLDER = '/tmp/command_test_tmp_folder'
TMP_PARENT_FOLDER = TMP_FOLDER
TMP_SRC_FOLDER = f'{TMP_PARENT_FOLDER}/command_test_tmp_src'
TMP_DEST_FOLDER = f'{TMP_PARENT_FOLDER}/command_test_tmp_dst'
TMP_SOURCE_FILE_PATH = os.path.join(TMP_SRC_FOLDER, TMP_FILE_NAME)
TESTFILE_ABC_TXT = relative_file_path(__file__, "util/dummy testfile abc.txt")
TESTFILE_BCD_TXT = relative_file_path(__file__, "util/dummy testfile bcd.txt")
TESTSCRIPT_HAPPY_PY = relative_file_path(__file__, "util/dummy testscript happy.py")


class CommandModuleMethodsExecuteCommands(unittest.TestCase):
    '''Tests for idseq_dag/util/command.py module methods to execute commands'''
    @patch('idseq_dag.util.command.log.write')
    def test_execute_legacy_format(self, _mock_log):
        '''WHEN command is a string, THEN execute and log a warning'''
        result = command.execute_with_output("echo 123")

        self.assertEqual(result, "123\n")
        _mock_log.assert_any_call(warning=True, message=MATCH_RE(".*legacy.*Use.*command_patterns"), obj_data=ANY)


    def test_execute_single_command_1(self):
        '''WHEN using SingleCommand with args that have special shell characters, THEN it doesn't execute subcommands'''
        cp1 = command_patterns.SingleCommand(
            cmd="echo",
            args=[
                "1",
                2,
                "$(pwd)",
                ";ls",
                "&& ls",
                "`pwd`",
                ">",
                "test.txt",
                "> test.txt",
                ">> test.txt",
                "&&",
                "ls",
                "$",
                "pwd"
            ]
        )

        result = command.execute_with_output(cp1)

        self.assertEqual(result, "1 2 $(pwd) ;ls && ls `pwd` > test.txt > test.txt >> test.txt && ls $ pwd\n")


    def test_execute_single_command_2(self):
        '''WHEN using SingleCommand with args that have special shell characters, THEN it doesn't execute subcommands'''
        assert " " in TESTFILE_ABC_TXT

        cp1 = command_patterns.SingleCommand(
            cmd="cat",
            args=[
                TESTFILE_ABC_TXT
            ]
        )

        result = command.execute_with_output(cp1)

        self.assertEqual(result, "abc")


    def test_execute_shell_script_command_1(self):
        '''WHEN using ShellScriptCommand with args that contain spaces or special characters, THEN it doesn't split them into separate arguments'''
        cp1 = command_patterns.ShellScriptCommand(
            script=r'''
                echo May the force be with you >> "$1"
                echo The truth is out there > "$1"
                echo Live longer and prosper >> "$1"
                grep "${@:3}" "$1" > "$2"
                cat "$2"
                rm "$1" "$2"
            ''',
            args=[
                r'''/tmp/tmp file with spaces, 'quotes', "double-quotes" and other bizarre characters `~&>.txt''',
                r'''/tmp/another tmp output file.txt''',
                "-e",
                "is out",
                "-e",
                "longer and prosper"
            ]
        )

        result = command.execute_with_output(cp1)

        self.assertEqual(result, "The truth is out there\nLive longer and prosper\n")


    def test_execute_shell_script_command_2(self):
        '''WHEN using SingleCommand with args that have special shell characters, THEN it doesn't execute subcommands'''
        cp1 = command_patterns.ShellScriptCommand(
            script=r'echo ${@: 1:$#-1} | sed ${@: $#}',
            args=[
                "1",
                2,
                "$(pwd)",
                ";ls",
                ";",
                "ls",
                "\n",
                "ls",
                "&& ls",
                "`pwd`",
                ">",
                "test.txt",
                "> test.txt",
                ">> test.txt",
                "&&",
                "ls",
                "$",
                "abc\nls"
                "pwd",
                "s/w/a/g"
            ]
        )

        result = command.execute_with_output(cp1)

        self.assertEqual(result, "1 2 $(pad) ;ls ; ls ls && ls `pad` > test.txt > test.txt >> test.txt && ls $ abc lspad\n")


    def test_execute_shell_script_command_3(self):
        '''WHEN using ShellScriptCommand with multiline script, THEN it works as expected'''
        cp1 = command_patterns.ShellScriptCommand(
            script=r'''
                a=123;
                echo May $a the force be with you \
                | sed "s/a/Z/g"
            ''',
            args=[
            ]
        )

        result = command.execute_with_output(cp1)

        self.assertEqual(result, "MZy 123 the force be with you\n")


    def test_specific_pattern_1(self):
        '''WHEN using ShellScriptCommand with multiline script, THEN it works as expected'''
        cp1 = command_patterns.ShellScriptCommand(
            script=r'''
                a=123;
                echo May $a the force be with you \
                | sed "s/a/Z/g"
            ''',
            args=[
            ]
        )

        result = command.execute_with_output(cp1)

        self.assertEqual(result, "MZy 123 the force be with you\n")


    def test_execute_python_cmd(self):
        '''WHEN using SingleCommand to invoke a .py file, THEN it works as expected'''
        cp1 = command_patterns.SingleCommand(
            cmd=TESTSCRIPT_HAPPY_PY,
            args=[
                "Hello!"
            ]
        )

        result = command.execute_with_output(cp1)

        self.assertEqual(result, "Python is happy to say: Hello!\n")


    def test_shellscript_with_param_array(self):
        '''WHEN using ShellScriptCommand to invoke a command with an array of parameters, THEN it works as expected'''
        cp1 = command_patterns.ShellScriptCommand(
            script=r'''paste "${slice_outputs[@]}"''',
            named_args={
                'slice_outputs': [
                    "-d",
                    r"\n",
                    TESTFILE_ABC_TXT,
                    TESTFILE_BCD_TXT
                ]
            }
        )

        result = command.execute_with_output(cp1)

        self.assertEqual(result, "abc\nbcd\n")


class CommandModuleMethodsWriteFileTestCase(unittest.TestCase):
    '''Tests for idseq_dag/util/command.py module methods that write content to files or create folders'''

    def tearDown(self):
        '''remove test temporary file used in test cases'''
        if os.path.exists(TMP_FILE):
            os.remove(TMP_FILE)
        subprocess.call(["rm", "-rf", TMP_FOLDER])


    def test_write_text_to_file(self):
        '''WHEN write_text_to_file is invoked, THEN write to the file the provided string'''
        command.write_text_to_file("done", TMP_FILE)

        self.assertEqual(file_contents(TMP_FILE), "done\n")


    def test_make_dirs(self):
        '''WHEN make_dirs is invoked, THEN it creates all required subdirs'''
        path = os.path.join(TMP_FOLDER, "1", "2")

        command.make_dirs(path)

        self.assertTrue(os.path.exists(path), f"folder {path} doesn't exist")


class CommandModuleMethodsFileOperationsTestCase(unittest.TestCase):
    '''Tests for idseq_dag/util/command.py module methods that manipulate files or create folders'''

    def setUp(self):
        '''creates move src and dest temp folders and a file to be copied/moved around'''
        subprocess.call(["mkdir", "-p", TMP_SRC_FOLDER])
        subprocess.call(["mkdir", "-p", TMP_DEST_FOLDER])
        subprocess.call(["touch", TMP_SOURCE_FILE_PATH])


    def tearDown(self):
        '''removes test temporary file used in test cases'''
        for p in [TMP_SRC_FOLDER, TMP_DEST_FOLDER]:
            subprocess.call(["rm", "-rf", p])


    def test_copy_file_to_folder(self):
        '''WHEN copy_file is invoked with a file path and a folder name, THEN it copies the file to the destination folder'''
        dest_file_path = os.path.join(TMP_DEST_FOLDER, TMP_FILE_NAME)

        command.copy_file(TMP_SOURCE_FILE_PATH, TMP_DEST_FOLDER)

        self.assertTrue(os.path.exists(dest_file_path), f"file {dest_file_path} doesn't exist")


    def test_copy_file_with_different_name(self):
        '''WHEN copy_file is invoked with a file path and a full path with a different file name, THEN it copies the file using the new file name'''
        new_file_name = TMP_FILE_NAME + ".new"
        dest_file_path = os.path.join(TMP_DEST_FOLDER, new_file_name)

        command.copy_file(TMP_SOURCE_FILE_PATH, dest_file_path)

        self.assertTrue(os.path.exists(dest_file_path), f"file {dest_file_path} doesn't exist")


    def test_move_file_to_folder(self):
        '''WHEN move_file is invoked with a file path and a folder name, THEN it copies the file to the destination folder'''
        dest_file_path = os.path.join(TMP_DEST_FOLDER, TMP_FILE_NAME)

        command.move_file(TMP_SOURCE_FILE_PATH, TMP_DEST_FOLDER)

        self.assertTrue(os.path.exists(dest_file_path), f"file {dest_file_path} doesn't exist")
        self.assertFalse(os.path.exists(TMP_SOURCE_FILE_PATH), f"file {TMP_SOURCE_FILE_PATH} shouldn't exist")


    def test_move_file_with_different_name(self):
        '''WHEN move_file is invoked with a file path and a full path with a different file name, THEN it copies the file using the new file name'''
        new_file_name = TMP_FILE_NAME + ".new"
        dest_file_path = os.path.join(TMP_DEST_FOLDER, new_file_name)

        command.move_file(TMP_SOURCE_FILE_PATH, dest_file_path)

        self.assertTrue(os.path.exists(dest_file_path), f"file {dest_file_path} doesn't exist")
        self.assertFalse(os.path.exists(TMP_SOURCE_FILE_PATH), f"file {TMP_SOURCE_FILE_PATH} shouldn't exist")


    def test_remove_file(self):
        '''WHEN remove_file is invoked, THEN it removes the file'''
        command.remove_file(TMP_SOURCE_FILE_PATH)

        self.assertFalse(os.path.exists(TMP_SOURCE_FILE_PATH), f"file {TMP_SOURCE_FILE_PATH} shouldn't exist")


    def test_remove_rf(self):
        '''WHEN remove_rf is invoked in a folder that contains files or other folders, THEN it removes the entire folder structure'''
        command.remove_rf(TMP_SRC_FOLDER)

        self.assertFalse(os.path.exists(TMP_SRC_FOLDER), f"file {TMP_SRC_FOLDER} shouldn't exist")


    def test_chmod(self):
        '''WHEN chmod is invoked, THEN it changes file permissions'''
        # Mark as write only
        command.chmod(TMP_SOURCE_FILE_PATH, 0o400)

        self.assertTrue(os.stat(TMP_SOURCE_FILE_PATH).st_mode & 0o777 == 0o400, f"file {TMP_SOURCE_FILE_PATH} should have 0o400 permissions")

        # Mark as read-write
        command.chmod(TMP_SOURCE_FILE_PATH, 0o600)

        self.assertTrue(os.stat(TMP_SOURCE_FILE_PATH).st_mode & 0o777 == 0o600, f"file {TMP_SOURCE_FILE_PATH} should have 0o600 permissions")


    def test_glob(self):
        '''WHEN glob is invoked, THEN it returns folders matching the pattern'''
        result = command.glob(f"{TMP_PARENT_FOLDER}/*")

        # assertCountEqual has a very misleading name, it actually compares two arrays disregarding item order - https://bugs.python.org/issue27071
        self.assertCountEqual(
            result,
            [
                TMP_SRC_FOLDER,
                TMP_DEST_FOLDER
            ]
        )


    def test_glob_with_limit(self):
        '''WHEN glob is invoked with a limit, THEN it returns folders matching the pattern'''
        result = command.glob(f"{TMP_PARENT_FOLDER}/*", max_results=1)

        either_one_of_those = {
            TMP_SRC_FOLDER,
            TMP_DEST_FOLDER
        }
        self.assertTrue(
            len(either_one_of_those.intersection(result))==1,
            f"result should be an array containing either {TMP_SRC_FOLDER} or {TMP_DEST_FOLDER}, but result contains {{ {result} }}"
        )


    def test_get_resource_filename(self):
        '''WHEN get_resource_filename is invoked with a relative file path to the root folder, THEN it returns the absolute path for that file.'''
        result = command.get_resource_filename("scripts/fastq-fasta-line-validation.awk")

        self.assertRegex(result, r"^.+/scripts/fastq-fasta-line-validation.awk")
