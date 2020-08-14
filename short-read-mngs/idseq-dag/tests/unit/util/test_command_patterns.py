import unittest
import subprocess
import os
import idseq_dag.util.command_patterns as command_patterns

TMP_DONE_FILE = '/tmp/123.done'
TMP_TEXT_FILE = '/tmp/123.txt'
TMP_FOLDER = '/tmp/single_command_temporary_folder/' # trailing slash is required

class SingleCommandTestCase(unittest.TestCase):
    '''Tests for idseq_dag/util/command_patterns.py SingleCommand class'''

    def test_open_1(self):
        '''WHEN args have special shell characters, THEN it doesn't execute subcommands'''
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

        p = cp1.open(stdout=subprocess.PIPE)
        stdout, stderr = p.communicate()

        self.assertFalse(stderr)
        self.assertEqual(stdout.decode(), "1 2 $(pwd) ;ls && ls `pwd` > test.txt > test.txt >> test.txt && ls $ pwd\n")


    def test_cd(self):
        '''WHEN using cd parameter, THEN it executes the command in the selected directory and resets it to previous dir before executing next command'''
        _shared_test_cd_parameter(
            test_context=self,
            pwd_command_pattern_with_cd=command_patterns.SingleCommand(
                cd=TMP_FOLDER,
                cmd="pwd",
                args=[]
            ),
            pwd_command_pattern_without_cd=command_patterns.SingleCommand(
                cmd="pwd",
                args=[]
            )
        )


class ScriptCommandTestCase(unittest.TestCase):
    def test_open_1(self):
        '''WHEN args have special shell characters, THEN it doesn't execute subcommands'''
        cp1 = command_patterns.ShellScriptCommand(
            script=r'''echo "${@:1:$#-1}" ${@:1:$#-1} | sed "${@:$#}";''',
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
                "abc\nls",
                "pwd",
                '"quotes"',
                "s/w/a/g"
            ]
        )

        p = cp1.open(stdout=subprocess.PIPE)
        stdout, stderr = p.communicate()

        self.assertFalse(stderr)
        self.assertEqual(
            stdout.decode(),
            "1 2 $(pad) ;ls ; ls \n ls && ls `pad` > test.txt > test.txt >> test.txt && ls $ abc\nls pad \"quotes\" "
            "1 2 $(pad) ;ls ; ls ls && ls `pad` > test.txt > test.txt >> test.txt && ls $ abc ls pad \"quotes\"\n"
        )

    def test_open_2(self):
        '''WHEN script uses shell variables within the script THEN it can get access to those variables but can use them from parameters'''
        cp1 = command_patterns.ShellScriptCommand(
            script=r'abc=123; echo abc=$abc, \$1=$1',
            args=[
                "$abc"
            ]
        )

        p = cp1.open(stdout=subprocess.PIPE)
        stdout, stderr = p.communicate()

        self.assertFalse(stderr)
        self.assertEqual(stdout.decode(), "abc=123, $1=$abc\n")

    def test_named_args(self):
        '''WHEN parameter named_args is used THEN variables are automatically expanded'''
        cp1 = command_patterns.ShellScriptCommand(
            script=r'''
                set -e;
                echo "original_string = \"${original_string}\"";
                echo "sed_patterns = \"${sed_patterns[@]}\"";
                echo "end_sed_options = \"${end_sed_options[@]}\"";
                echo "sed_patterns[0] = \"${sed_patterns[0]}\"";
                echo "sed_patterns[1] = \"${sed_patterns[1]}\"";
                echo "sed_patterns[2] = \"${sed_patterns[2]}\"";
                echo "sed_patterns[3] = \"${sed_patterns[3]}\"";
                echo "sed_patterns[4] = \"${sed_patterns[4]}\"";
                echo "sed_patterns[5] = \"${sed_patterns[5]}\"";
                echo "empty_array = \"${empty_array[@]}\"";
                echo "empty_str = \"${empty_str}\"";
                echo "${original_string}" | $sed_command "${sed_patterns[@]}" | sed "${end_sed_options[@]}" 
            ''',
            named_args={
                'original_string': "ABCDEF",
                'sed_command': 'sed',
                'end_sed_options': [
                    '-e',
                    's/X/Y/'
                ],
                'sed_patterns': [
                    '-e',
                    's/A/Z /',
                    '-e',
                    's/B/X/',
                    '-e',
                    's/;&`/^/'
                ],
                'empty_str': '',
                'empty_array': []
            }
        )

        p = cp1.open(stdout=subprocess.PIPE)
        stdout, stderr = p.communicate()

        self.assertFalse(stderr)
        self.assertEqual(
            stdout.decode(),
            'original_string = "ABCDEF"\n'
            'sed_patterns = "-e s/A/Z / -e s/B/X/ -e s/;&`/^/"\n'
            'end_sed_options = "-e s/X/Y/"\n'
            'sed_patterns[0] = "-e"\n'
            'sed_patterns[1] = "s/A/Z /"\n'
            'sed_patterns[2] = "-e"\n'
            'sed_patterns[3] = "s/B/X/"\n'
            'sed_patterns[4] = "-e"\n'
            'sed_patterns[5] = "s/;&`/^/"\n'
            'empty_array = ""\n'
            'empty_str = ""\n'
            'Z YCDEF\n'
        )

    def test_cd(self):
        '''WHEN using cd parameter, THEN it executes the command in the selected directory and resets it to previous dir before executing next command'''
        _shared_test_cd_parameter(
            test_context=self,
            pwd_command_pattern_with_cd=command_patterns.ShellScriptCommand(
                cd=TMP_FOLDER,
                script="pwd",
                args=[]
            ),
            pwd_command_pattern_without_cd=command_patterns.ShellScriptCommand(
                script="pwd",
                args=[]
            )
        )

# --- Shared test cases ---
def _shared_test_cd_parameter(test_context, pwd_command_pattern_with_cd, pwd_command_pattern_without_cd):
    '''WHEN using cd parameter, THEN it executes the command in the selected directory and resets it to previous dir before executing next command'''
    # Get current dir
    os.system(f"rm -rf {TMP_FOLDER}; mkdir -p {TMP_FOLDER};")
    current_dir = os.path.join(subprocess.check_output("pwd").decode().rstrip(), "")
    test_context.assertNotEqual(current_dir, TMP_FOLDER, "Current directory cannot be the same as TMP_FOLDER")

    # Test with CD
    p = pwd_command_pattern_with_cd.open(stdout=subprocess.PIPE)
    stdout, stderr = p.communicate()

    test_context.assertFalse(stderr)
    test_context.assertEqual(
        os.path.join(stdout.decode().rstrip(), ""),
        TMP_FOLDER,
        "Selected folder should match cd parameter"
    )

    # Test without CD
    p = pwd_command_pattern_without_cd.open(stdout=subprocess.PIPE)
    stdout, stderr = p.communicate()

    test_context.assertFalse(stderr)
    test_context.assertEqual(
        os.path.join(stdout.decode().rstrip(), ""),
        current_dir,
        "Folder has not been properly reseted between calls"
    )

    # teardown
    os.system(f"rm -rf {TMP_FOLDER};")
