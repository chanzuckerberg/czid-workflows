import logging
import multiprocessing
import os
import sys
import json
import datetime
import math
import functools
import secrets

from contextlib import contextmanager

__print_lock = multiprocessing.RLock()

def configure_logger(log_file=None):
    logger = logging.getLogger()
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
    with __print_lock:
        logger = logging.getLogger()
        if warning:
            logger.warning(msg=message, extra={"obj_data": obj_data})
        elif debug:
            logger.debug(msg=message, extra={"obj_data": obj_data})
        else:
            logger.info(msg=message, extra={"obj_data": obj_data})
        if flush:
            sys.stdout.flush()


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
        duration = (datetime.datetime.now()-start_time).total_seconds()
        obj_data["duration_ms"] = math.floor(duration * 1000)
    write(obj_data=obj_data, warning=warning, debug=debug, flush=flush)
    return datetime.datetime.now()


def log_execution(values: dict = None):
    '''
    Decorates a function and write to logs fn_start and fn_end events.

    Parameters:
    values(dict): Optional. Values associated to that event. Will be logged in json format.

    Example:

    from log import log_execution
    @log_execution()
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


@contextmanager
def log_context(context_name: str, values: dict = None, log_caller_info: bool = False):
    '''
    Log context manager to track started and completed events for a block of code

    Parameters:
    context_name (str): Name for this context.
    values(dict): Optional. Values associated to that event. Will be logged in json format.
    log_caller_info(bool): Log the caller file_name and method name.

    Example:
    # some_file.py
    from log import log_context

    def some_method
        # ... some code here ...
        with log_context("sub_step", {"abc": "123"}):
            # ... your code goes here ...

    The code above will generate two log entries:
    {"event": "ctx_start", "context_name": "sub_step", "caller": {"filename": "some_file.py", "method": "some_method"}, "values": {abc": 123}}
    {"event": "ctx_end", "context_name": "sub_step", "caller": {"filename": "some_file.py", "method": "some_method"}, "values": {abc": 123}, "duration_ms": 5006} (5.0 seconds)
    '''
    extra_fields = {"context_name": context_name, "uid": secrets.token_hex(6)}
    if log_caller_info:
        frame = sys._getframe(2)
        f_code = frame.f_code
        extra_fields["caller"] = {"filename": os.path.basename(f_code.co_filename), "method": f_code.co_name, "f_lineno": frame.f_lineno}
    start = log_event("ctx_start", values, extra_fields=extra_fields)
    try:
        yield
        log_event("ctx_end", values, start_time=start, extra_fields=extra_fields)
    except Exception as e:
        extra_fields["error_type"] = type(e).__name__
        extra_fields["error_args"] = e.args
        log_event("ctx_error", values, start_time=start, extra_fields=extra_fields)
        raise e


def set_up_stdout():
    # Unbuffer stdout and redirect stderr into stdout. This helps observe logged
    # events in realtime.
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)
    os.dup2(sys.stdout.fileno(), sys.stderr.fileno())


class JsonFormatter(logging.Formatter):
    default_time_format = '%Y-%m-%dT%H:%M:%S'
    default_msec_format = '%s.%03d'

    @staticmethod
    def _default_json_serializer(obj):
        '''Default serializer that returns string <<non-sertializable: TYPE_NAME>> for that type'''
        return f"<<non-serializable: {type(obj).__qualname__}>>"

    def format(self, record):
        obj = {"timestamp": self.formatTime(record)}
        if record.msg is not None:
            obj['msg'] = super().format(record)
        if record.obj_data is not None:
            obj['data'] = record.obj_data
        obj.update({"thread": record.threadName, "pid": record.process, "level": record.levelname})
        return json.dumps(obj, default=JsonFormatter._default_json_serializer)