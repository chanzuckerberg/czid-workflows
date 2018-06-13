import logging
import multiprocessing
import os
import sys

print_lock = multiprocessing.RLock()


def configure_logger(log_file=None):
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    if log_file:
        handler = logging.FileHandler(log_file)
        formatter = logging.Formatter("%(asctime)s: %(message)s")
        handler.setFormatter(formatter)
        logger.addHandler(handler)

    # Echo to stdout so they get to CloudWatch
    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(logging.INFO)
    formatter = logging.Formatter("%(asctime)s: %(message)s")
    handler.setFormatter(formatter)
    logger.addHandler(handler)


def write(message, warning=False, flush=True):
    logger = logging.getLogger()
    with print_lock:
        if warning:
            logger.warning(message)
        else:
            logger.info(message)
        if flush:
            sys.stdout.flush()


def set_up_stdout():
    # Unbuffer stdout and redirect stderr into stdout. This helps observe logged
    # events in realtime.
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)
    os.dup2(sys.stdout.fileno(), sys.stderr.fileno())
