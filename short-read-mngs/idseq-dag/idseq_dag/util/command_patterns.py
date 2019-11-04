from typing import List, Union, Dict, Any
import os
import json
import abc
from collections.abc import Iterable
import subprocess
import idseq_dag.util.log as log

class CommandPattern(abc.ABC):
    r'''
        Abstract class with common functionality for all command patterns.
        A command pattern should be used whenever you are invoking
        command.py#execute method instead of using the legacy string format.
    '''
    def __init__(self):
        # Saves extra info about which class created this command pattern. Very useful for debuging purposes.
        self.caller_info = log.get_caller_info(3)

    def __repr__(self):
        r'''
            Returns a more descriptive message when printing this object.
            for example:
                cp = command_patterns.SingleCommand(cmd="echo", args=["test"])
                print(cp)
            will return something like:
                <idseq_dag.util.command_patterns.SingleCommand object {"cmd": "echo", "args": ["test"], "cd": ".", "type": "SingleCommand", "caller_info": {"filename": "test_command_patterns.py", "method": "test_open_1", "f_lineno": 34}}>
            This is specially useful when CalledProcessError is raised, so you can get extra context about the command that has failed.
        '''
        json_data = json.dumps(self.as_dict(), default=log.default_json_serializer)
        class_name = ".".join([type(self).__module__, type(self).__qualname__])
        return f"<{class_name} object {json_data}>"

    def as_dict(self):
        r'''
            This method is used to retrieve a dictionary with all attributes from this class
            that will be sent to logs. Since this is an abstract class, descendants need
            to extend this method to add extra fields.
        '''
        return {'type': type(self).__name__, 'caller_info': self.caller_info}


class SingleCommand(CommandPattern):
    r'''
        This pattern is used to execute a single command.
        This is the preferred pattern and you should favor it over ShellScriptCommand,
        since you don't need to quote parameters when using SingleCommand.

        command.execute(
            command_patterns.SingleCommand(
                cmd="grep"
                args=[
                    "-e",
                    r'some input parameter that can include spaces, and other characters like ;&()" etc.',
                    "my_file.txt"
                ],
                cd="/work_dir"
            )
        )
    '''
    def __init__(self, cmd: str, args: List[Union[int, str]], cd: str = None):
        super().__init__()
        self.cmd = cmd
        self.args = args
        self.cd = cd

    def _command_args(self):
        return iter([self.cmd, *(str(a) for a in self.args)])

    def open(self, stdin: int = None, stdout: int = None, stderr: int = None) -> subprocess.Popen:
        popen_handler = subprocess.Popen(
            self._command_args(),
            shell=False,
            cwd=self.cd,
            stdin=stdin,
            stdout=stdout,
            stderr=stderr
        )
        return popen_handler

    def as_dict(self) -> Dict[str, Any]:
        result = {'cmd': self.cmd, 'args': self.args, 'cd': self.cd or os.curdir}
        result.update(super().as_dict())
        return result


class ShellScriptCommand(CommandPattern):
    r"""
        This pattern is used to execute a shell command that can include pipes, redirections and other things.
        It provides a way of passing parameters without using string interpolation, allowing you to properly send
        parameters containing spaces and other special characters.

        You should avoid using this pattern and try using SingleCommand instead whenever it is possible, since
        with SingleCommand you don't need to worry about properly escaping variables.

        This patterns has two modes.

        - Mode 1 - `args` array:
           `args` is an array of string (or int) parameters that will be passed to your script following
           regular shellscript standards ($1, $2, $3 ...). Each entry will be assigned to these variables
           and they can be used like in the example below.

            command.execute(
                command_patterns.ShellScriptCommand(
                    cmd=r'cat "$1" | grep -e "$2" > "$3";', # IMPORTANT: put quotes around these parameters
                    args=[
                        "my input file.txt",
                        "some grep pattern",
                        "output file.txt"
                    ],
                    cd="/work_dir"
                )
            )

          It is mandatory to use "quotes" around these parameters:

            Correct:
                cmd=r'''cat "$1" | grep -e "$2" > "$3";'''
            WRONG:
                cmd=r'''cat $1 | grep -e $2 > $3;''' # DON'T DO THAT

          otherwise bash will expand values containing spaces and other special characters
          in an unpredictable way (depending on your input), spreading them to multiple parameters.
          Using quotes is a safe way to guarantee these parameters will work in a consistent way no matter
          what is your input value.

          If you need to spread something into multiple parameters, you can use bash arrays:

          Correct:
            cmd=r'''cat "$1" | grep "${@:2:$#-1}" > "$5";''', # note the array reference is also quoted
            args=[
                "my input file.txt",
                "-e", "pattern 1",
                "-e", "pattern 2",
                "my output file.txt"
            ]
          WRONG:
            cmd=r'''cat $1 | grep -e ${@:2:$#-1} > $3;''' # DON'T DO THAT


        - Mode 2 - `named_args` object:
            This is way a more convenient way to expand values into variables. It will automatically define shell variables,
            including proper array expansion:

            Ex:
                cmd=r'''cat "${input_file}" | grep "${grep_pattern_array[@]}" > "${output_file}";''', # note the array expansion is also quoted
                named_args={
                    "input_file": "my input file.txt",
                    "grep_pattern_array": [
                        "-e", "pattern 1",
                        "-e", "pattern 2"
                    ],
                    "output_file": "my output file.txt"
                ]

            Behind the scenes, this method is declaring an script that behaves similar to:
                script.sh:
                    input_file="$1"
                    grep_pattern_array=( ${@:2:4} )
                    output_file=$5
                    cat "${input_file}" | grep "${grep_pattern_array[@]}" > "${output_file}";
                ...
                script.sh <...all parameters here...>

    """
    def __init__(self, script: str, args: List[Union[int, str]] = None, named_args: Dict[str, Union[int, str, List[Union[int, str]]]] = None, cd: str = None):
        super().__init__()
        assert (args is None) or (named_args is None), "You need to use either args or named_args"
        self.script = script
        self.args = args
        self.named_args = named_args
        self.popen_handler = None
        self.cd = cd

    def _script_named_args(self):
        args = []
        script_vars = []
        pos = 1
        for key, value in self.named_args.items():
            if isinstance(value, str) or isinstance(value, int):
                script_vars.append(f'{key}="${{{pos}}}"')  # ex: abc="${3}"
                args.append(value)
                pos += 1
            elif isinstance(value, Iterable):
                items_cnt = len(value)
                if items_cnt == 0:
                    shell_array = ""
                else:
                    shell_array = '${' + f"@:{pos}:{items_cnt}" + '}'  # ex: ${@:2:4}
                script_vars.append(f'{key}=( "{shell_array}" )')  # ex: abc=( "${@:2:4}" )
                args.extend(value)
                pos += items_cnt
            else:
                raise Exception(f"Invalid type for parameter {key}. Type: {type(value)}. Value: {value}.")

        script_named_args = "; ".join(script_vars) + "; "

        return script_named_args, args

    def _command_args(self):
        if self.args is not None:
            script_named_args = ""
            args = self.args
        elif self.named_args is not None:
            script_named_args, args = self._script_named_args()
        else:
            script_named_args = ""
            args = []

        return iter([script_named_args + self.script, "", *(str(a) for a in args)])

    def open(self, stdin: int = None, stdout: int = None, stderr: int = None) -> subprocess.Popen:
        self.popen_handler = subprocess.Popen(
            self._command_args(),
            shell=True,
            cwd=self.cd,
            stdin=stdin,
            stdout=stdout,
            stderr=stderr,
            executable="/bin/bash"
        )
        return self.popen_handler

    def as_dict(self) -> Dict[str, Any]:
        if self.named_args is not None:
            result = {'script': self.script, 'named_args': self.named_args, 'cd': self.cd or os.curdir}
        else:
            result = {'script': self.script, 'args': self.args, 'cd': self.cd or os.curdir}
        result.update(super().as_dict())
        return result
