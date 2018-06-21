import time
import threading
import subprocess
import os
import multiprocessing
import logging

import idseq_dag.util.command as command
import idseq_dag.util.log as log

# Peak network and storage perf for a typical small instance is saturated by
# just a few concurrent streams.
MAX_CONCURRENT_COPY_OPERATIONS = 8
IOSTREAM = multiprocessing.Semaphore(MAX_CONCURRENT_COPY_OPERATIONS)
# Make a second semaphore for uploads to reserve some capacity for downloads.
MAX_CONCURRENT_UPLOAD_OPERATIONS = 4
IOSTREAM_UPLOADS = multiprocessing.Semaphore(MAX_CONCURRENT_UPLOAD_OPERATIONS)


def check_s3_presence(s3_path):
    """True if s3_path exists. False otherwise."""
    try:
        o = command.execute_with_output("aws s3 ls %s" % s3_path)
        if o:
            return True
    except:
        pass
    return False


def check_s3_presence_for_file_list(s3_dir, file_list):
    for f in file_list:
        if not check_s3_presence(os.path.join(s3_dir, f)):
            return False
    return True


def touch_s3_file(s3_file_path):
    try:
        command.execute("aws s3 cp --metadata '{\"touched\":\"now\"}' %s %s" % (s3_file_path, s3_file_path))
        return True
    except:
        return False


def touch_s3_file_list(s3_dir, file_list):
    for f in file_list:
        touch_s3_file(os.path.join(s3_dir, f))


def install_s3mi(installed={}, mutex=threading.RLock()):  #pylint: disable=dangerous-default-value
    with mutex:
        if installed:  # Mutable default value persists
            return
        try:
            # This is typically a no-op.
            command.execute(
                "which s3mi || pip install git+git://github.com/chanzuckerberg/s3mi.git"
            )
            command.execute(
                "s3mi tweak-vm || echo s3mi tweak-vm is impossible under docker. Continuing..."
            )
        finally:
            installed['time'] = time.time()


def fetch_from_s3(src,
                  dst,
                  auto_unzip=False,
                  auto_untar=False,
                  allow_s3mi=False,
                  mutex=threading.RLock(),
                  locks={}):  #pylint: disable=dangerous-default-value
    """Fetch a file from S3 if needed using either s3mi or aws cp."""
    with mutex:
        if os.path.exists(dst) and os.path.isdir(dst):
            dst = os.path.join(dst, os.path.basename(src))
        # Right now untar and unzip are mutually exclusive
        # TODO(boris): we might change this
        unzip = auto_unzip and dst.endswith(".gz")
        untar = auto_untar and dst.endswith(".tar")
        if unzip:
            dst = dst[:-3]  # Remove .gz
        if untar:
            dst = dst[:-4] # Remove .tar
        abspath = os.path.abspath(dst)
        if abspath not in locks:
            locks[abspath] = threading.RLock()
        destination_lock = locks[abspath]

    with destination_lock:
        if os.path.exists(dst):
            # No need to fetch this file from s3, it has been just produced
            # on this instance.
            return dst

        try:
            destdir = os.path.dirname(dst)
            if destdir:
                os.makedirs(destdir)
        except OSError as e:
            # It's okay if the parent directory already exists, but all other
            # errors are fatal.
            if e.errno != os.errno.EEXIST:
                log.write("Error in creating destination directory.")
                raise

        with IOSTREAM:
            try:
                if allow_s3mi:
                    try:
                        install_s3mi()
                    except:
                        log.write("s3mi failed to install.")
                        allow_s3mi = False

                if untar:
                    command_params = f" | tar xvf - -C {destdir}"
                else:
                    pipe_filter = ""
                    if unzip:
                        pipe_filter = "| gzip -dc "
                    command_params = f"{pipe_filter} > {dst}"
                log.write(f"command_params: {command_params}")

                try:
                    assert allow_s3mi
                    cmd = f"s3mi cat {src} {command_params}"
                    command.execute(cmd)
                except:
                    log.write(
                        "Failed to download with s3mi. Trying with aws s3 cp..."
                    )
                    cmd = f"aws s3 cp --quiet {src} - {command_params}"
                    command.execute(cmd)
                return dst
            except subprocess.CalledProcessError:
                # Most likely the file doesn't exist in S3.
                log.write(
                    "Failed to fetch file from S3. Most likely does not exist."
                )
                # Delete potential empty file remnant
                if os.path.isfile(dst) and os.stat(dst).st_size == 0:
                    os.remove(dst)
                return None


@command.retry
def upload_with_retries(from_f, to_f):
    command.execute(f"aws s3 cp --quiet {from_f} {to_f}")


def upload(from_f, to_f, status, status_lock=threading.RLock()):
    try:
        with IOSTREAM_UPLOADS:  # Limit concurrent uploads so as not to stall the pipeline.
            with IOSTREAM:  # Still counts toward the general semaphore.
                upload_with_retries(from_f, to_f)
            with status_lock:
                status[from_f] = "success"
    except:
        with status_lock:
            status[from_f] = "error"
        raise


def upload_log_file(sample_s3_output_path, lock=threading.RLock()):
    with lock:
        logh = logging.getLogger().handlers[0]
        logh.flush()
        command.execute("aws s3 cp --quiet %s %s/;" % (logh.baseFilename,
                                                       sample_s3_output_path))
