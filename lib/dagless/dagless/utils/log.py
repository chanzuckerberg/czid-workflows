import logging
import multiprocessing
import os
import sys
import json
import datetime
import math
import functools
import secrets
from enum import Enum
from collections import namedtuple

from contextlib import contextmanager

print_lock = multiprocessing.RLock()  # print_lock is needed outside this module for voodoo magic in run_in_subprocess

LogContext = namedtuple('LogContext', ['values', 'extra_fields'])

def configure_logger(log_file=None):
    # reduce log noise from aws boto3
    for s in ['boto3', 'botocore', 's3transfer', 'urlib3']:
        logging.getLogger(s).setLevel(logging.WARNING)

    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)

    if log_file:
        handler = logging.FileHandler(log_file)
        formatter = JsonFormatter()
        handler.setFormatter(formatter)
        logger.addHandler(handler)

    # Echo to stdout so they get to CloudWatch
    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(logging.INFO)
    formatter = JsonFormatter()
    handler.setFormatter(formatter)
    logger.addHandler(handler)


def write(message: str = None, warning: bool = False, flush: bool = True,
          debug: bool = False, obj_data: object = None):
    '''
    Write a new log entry.

    Parameters:
    message (str): String message to log. This field is optional, since you might want to send obj_data instead.
    obj_data (object): Object that will be serialized as json. dict is the preferred format,
                       although any json serializable type might work (ex: bool). This field is optional,
                       since you might want to send obj_data instead.
    warning(boolean): Optional. Log this event using warning level. This level has precedence over debug.
    debug(boolean): Optional. Log this event using debug level.
    flush(boolean): Optional (default true). Flush stdout after logging.
    '''
    with print_lock:
        logger = logging.getLogger(__name__)
        if warning:
            logger.warning(msg=message, extra={"obj_data": obj_data})
        elif debug:
            logger.debug(msg=message, extra={"obj_data": obj_data})
        else:
            logger.info(msg=message, extra={"obj_data": obj_data})
        if flush:
            sys.stdout.flush()


def _get_current_time():
    return datetime.datetime.now()


def log_event(event_name: str, values: object = None, start_time: datetime.datetime = None,
              warning: bool = False, debug: bool = False, flush: bool = True,
              extra_fields: dict = {}):
    '''
    Write a new log event.

    Parameters:
    event_name (str): name of the event
    values(dict): Optional. Values associated to that event. Will be logged in json format.
    start_time(datetime): Optional. If given, it will calc the elapsed time from this datetime to now and write it to the log line.
    extra_fields(dict): Optional. If given, keys will be merged to the main event hash.
    warning(boolean): Optional. Log this event using warning level. This level has precedence over debug. (see log.write)
    debug(boolean): Optional. Log this event using debug level. (see log.write)
    flush(boolean): Optional (default true). Flush stdout after logging. (see log.write)

    Returns:
    datetime: Now. It can be used to pass to parameter start_time in a future log_event call.

    Example:
    start = log_event("downloaded_started", {"file":"abc"})
    # ... download your file here ...
    log_event("downloaded_completed", {"file":"abc"}, start_time=start)
    '''
    obj_data = {"event": event_name, **extra_fields}
    if values is not None:
        obj_data["values"] = values
    if start_time is not None:
        duration = (_get_current_time() - start_time).total_seconds()
        obj_data["duration_ms"] = math.floor(duration * 1000)
    write(obj_data=obj_data, warning=warning, debug=debug, flush=flush)
    return _get_current_time()


def log_function_execution(values: dict = None):
    '''
    Decorates a function and write to logs fn_start and fn_end events.

    Parameters:
    values(dict): Optional. Values associated to that event. Will be logged in json format.

    Example:

    from log import log_function_execution
    @log_function_execution()
    def my_function
        # ...

    The code above will generate two log entries:
    Event fn_start {"n": "my_function"}
    Event fn_end {"n": "my_function"} (0.3 sec)
    '''
    def decorator_fn(func):
        @functools.wraps(func)
        def wrapper_fn(*args, **kwargs):
            if values is None:
                val = {"name": func.__qualname__}
            else:
                val = {"name": func.__qualname__, "v": values}
            start = log_event("fn_start", val)
            try:
                result = func(*args, **kwargs)
                log_event("fn_end", val, start_time=start)
                return result
            except Exception as e:
                val["error_type"] = type(e).__name__
                val["error_args"] = e.args
                log_event("fn_error", val, start_time=start)
                raise e
        return wrapper_fn
    return decorator_fn


def get_caller_info(frame_depth):
    try:
        frame = sys._getframe(frame_depth)
        f_code = frame.f_code
        return {"filename": os.path.basename(f_code.co_filename), "method": f_code.co_name, "f_lineno": frame.f_lineno}
    except Exception as e:
        return {"error": "cannot retrieve caller info", **parse_exception(e)}


def parse_exception(e):
    return {"error_type": type(e).__name__, "error_args": e.args}


class LogContextMode(Enum):
    START_END_LOG_EVENTS = 'START_END_LOG_EVENTS'
    EXEC_LOG_EVENT = 'EXEC_LOG_EVENT'

@contextmanager
def log_context(context_name: str, values: dict = None, log_caller_info: bool = False, log_caller_info_depth: int = 0, log_context_mode: LogContextMode = LogContextMode.START_END_LOG_EVENTS):
    '''
    Log context manager to track started and completed events for a block of code

    Parameters:
    context_name (str): Name for this context.
    values(dict): Optional. Values associated to that event. Will be logged in json format.
    log_caller_info(bool): Log the caller file_name and method name.
    log_caller_info_depth(int): A number that indicates how deep it goes into the stack trace to get data about the caller.
                                Default is 0, which represents the method that invoked log_context method.
                                For example, if you increase this number to 1, this will log
                                the parent method that invoked the method where log_context is being invoked.
                                Ex:
                                    def method_a():
                                        method_b()
                                        # caller info will show method_a, because depth is 0
                                        log.log_context(context_name: 'ctx', log_caller_info: true, log_caller_info_depth=0)

                                    def method_b():
                                        # caller info will also show method_a, because depth is 1 and method_b has been invoked by method_a
                                        log.log_context(context_name: 'ctx', log_caller_info: true, log_caller_info_depth=1)
    log_context_mode(bool): defines the behavior of log_context
        - START_END_LOG_EVENTS: (default) writes a `ctx_start` event before start executing the wrapped code block, and a ctx_end when the code block succeeds.
        - EXEC_LOG_EVENT: only writes a `ctx_exec` event after executing the wrapped code block with success.
        In both modes, if an error happens, this function will write a ctx_error event.


    Example:
    # some_file.py
    from log import log_context

    def some_method
        # ... some code here ...
        with log_context("sub_step", {"abc": "123"}):
            # ... your code goes here ...
        with log_context("sub_step_2", {"def": 1}, log_context_mode: LogContextMode.EXEC_LOG_EVENT):
            # ... your code goes here ...

    The code above using log_context_mode START_END_LOG_EVENTS (default) will write the following entries:
        - Before start executing code block
            {"event": "ctx_start", "context_name": "sub_step", "uid": "BA31CF", "caller": {"filename": "some_file.py", "method": "some_method"}, "values": {"abc": 123}}
        - When code block succeeded
            {"event": "ctx_end", "context_name": "sub_step", "uid": "BA31CF", "caller": {"filename": "some_file.py", "method": "some_method"}, "values": {"abc": 123}, "duration_ms": 5006}
    The code above using log_context_mode EXEC_LOG_EVENT will write the following entries:
        - When block completed execution
            {"event": "ctx_exec", "context_name": "sub_step_2", "uid": "C31134", "caller": {"filename": "some_file.py", "method": "some_method"}, "values": {"def": 1}}
    In both cases:
        - When code block terminated with an error
            {"event": "ctx_error", "context_name": "sub_step", "uid": "BA31CF", "extra_fields": {"error_type": "RuntimeError", "error_args": ["Error message"]} "caller": {"filename": "some_file.py", "method": "some_method"}, "values": {"abc": 123}}
    '''
    context = LogContext(values=values, extra_fields={"context_name": context_name, "uid": secrets.token_hex(6)})

    if log_caller_info:
        context.extra_fields["caller"] = get_caller_info(3 + log_caller_info_depth)

    start = _get_current_time()
    if log_context_mode == LogContextMode.START_END_LOG_EVENTS:
        log_event("ctx_start", context.values, extra_fields=context.extra_fields)

    try:
        yield context
        if log_context_mode == LogContextMode.START_END_LOG_EVENTS:
            log_event("ctx_end", context.values, start_time=start, extra_fields=context.extra_fields)
        elif log_context_mode == LogContextMode.EXEC_LOG_EVENT:
            log_event("ctx_exec", values, start_time=start, extra_fields=context.extra_fields)
    except Exception as e:
        context.extra_fields.update(parse_exception(e))
        log_event("ctx_error", context.values, start_time=start, extra_fields=context.extra_fields)
        raise e


def default_json_serializer(obj):
    '''Default serializer that returns string <<non-sertializable: TYPE_NAME>> for that type'''
    return f"<<non-serializable: {type(obj).__qualname__}>>"


class JsonFormatter(logging.Formatter):
    default_time_format = '%Y-%m-%dT%H:%M:%S'
    default_msec_format = '%s.%03d'

    def format(self, record):
        obj = {"time": self.formatTime(record)}
        if record.msg is not None:
            obj['msg'] = super().format(record)
        if hasattr(record, 'obj_data') and (record.obj_data is not None):
            obj['data'] = record.obj_data
        obj.update({"thread": record.threadName, "pid": record.process, "level": record.levelname})
        return json.dumps(obj, default=default_json_serializer)

