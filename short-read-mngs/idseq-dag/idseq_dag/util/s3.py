import time
import threading
import subprocess
import os
import multiprocessing
import logging
import errno
from idseq_dag.util.trace_lock import TraceLock

import idseq_dag.util.command as command
import idseq_dag.util.log as log

# Peak network and storage perf for a typical small instance is saturated by
# just a few concurrent streams.
MAX_CONCURRENT_COPY_OPERATIONS = 8
IOSTREAM = multiprocessing.Semaphore(MAX_CONCURRENT_COPY_OPERATIONS)
# Make a second semaphore for uploads to reserve some capacity for downloads.
MAX_CONCURRENT_UPLOAD_OPERATIONS = 4
MAX_CONCURRENT_S3MI_DOWNLOADS = 2
# If a s3mi slot does not free up within MAX_S3MI_WAIT seconds, we use plain old aws s3 cp instead of s3mi.
MAX_S3MI_WAIT = 15
S3MI_SEM = multiprocessing.Semaphore(MAX_CONCURRENT_S3MI_DOWNLOADS)
IOSTREAM_UPLOADS = multiprocessing.Semaphore(MAX_CONCURRENT_UPLOAD_OPERATIONS)

def split_identifiers(s3_path):
    return s3_path[5:].split("/", 1)

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


def install_s3mi(installed={}, mutex=TraceLock("install_s3mi", threading.RLock())):  #pylint: disable=dangerous-default-value
    with mutex:
        if installed:  # Mutable default value persists
            return
        try:
            # This is typically a no-op.
            command.execute(
                "which s3mi || pip install git+git://github.com/chanzuckerberg/s3mi.git"
            )
            command.execute(
                "s3mi tweak-vm || echo s3mi tweak-vm sometimes fails under docker. Continuing..."
            )
        finally:
            installed['time'] = time.time()


DEFAULT_AUTO_UNZIP = False
DEFAULT_AUTO_UNTAR = False
DEFAULT_ALLOW_S3MI = False
ZIP_EXTENSIONS = {
    ".lz4": "lz4 -dc",
    ".bz2": "lbzip2 -dc",
    ".gz": "gzip -dc", # please avoid .gz if you can (slow to decompress)
}

def fetch_from_s3(src,   #pylint: disable=dangerous-default-value
                  dst,
                  auto_unzip=DEFAULT_AUTO_UNZIP,
                  auto_untar=DEFAULT_AUTO_UNTAR,
                  allow_s3mi=DEFAULT_ALLOW_S3MI,
                  mutex=TraceLock("fetch_from_s3", threading.RLock()),
                  locks={}):
    """Fetch a file from S3 if needed, using either s3mi or aws cp."""
    assert "prod/samples/549/" not in src, "Sorry, project 549 has been disabled."
    with mutex:
        if os.path.exists(dst) and os.path.isdir(dst):
            dst = os.path.join(dst, os.path.basename(src))
        unzip = ""
        if auto_unzip:
            for zext, zdecomp in ZIP_EXTENSIONS.items():
                if dst.lower().endswith(zext.lower()):
                    unzip += f" | {zdecomp}"
                    dst = dst[:-len(zext)]
                    break
        untar = auto_untar and dst.lower().endswith(".tar")
        if untar:
            dst = dst[:-4]  # Remove .tar
        abspath = os.path.abspath(dst)
        if abspath not in locks:
            locks[abspath] = TraceLock(f"fetch_from_s3: {abspath}", threading.RLock())
        destination_lock = locks[abspath]

    with destination_lock:
        # This check is a bit imperfect when untarring... unless you follow the discipline that
        # all contents of file foo.tar are under directory foo/...
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
            if e.errno != errno.EEXIST:
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

                if allow_s3mi:
                    wait_start = time.time()
                    allow_s3mi = S3MI_SEM.acquire(timeout=MAX_S3MI_WAIT)
                    wait_duration = time.time() - wait_start
                    if not allow_s3mi:
                        log.write(f"Failed to acquire S3MI semaphore after waiting {wait_duration} seconds for {src}.")
                    elif wait_duration >= 5:
                        log.write(f"Waited {wait_duration} seconds to acquire S3MI semaphore for {src}.")

                write_dst = f" > {dst}"
                if untar:
                    write_dst = f" | tar xvf - -C {destdir}"
                command_params = f"{unzip} {write_dst}"
                log.write(f"command_params: {command_params}")

                try:
                    assert allow_s3mi
                    cmd = f"set -o pipefail; s3mi cat {src} {command_params}"
                    command.execute(cmd)
                except:
                    if allow_s3mi:
                        allow_s3mi = False
                        S3MI_SEM.release()
                        log.write(
                            "Failed to download with s3mi. Trying with aws s3 cp..."
                        )
                    cmd = f"set -o pipefail; aws s3 cp --only-show-errors {src} - {command_params}"
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
            finally:
                if allow_s3mi:
                    S3MI_SEM.release()



def fetch_reference(src,
                    dst,
                    auto_unzip=True,
                    auto_untar=True,
                    allow_s3mi=DEFAULT_ALLOW_S3MI):
    '''
        This function behaves like fetch_from_s3 in most cases, with one excetpion:

            *  When auto_unzip=True, src is NOT a compressed file and NOT a tar file, and a compressed
               version of src exists in s3, this function downloads and uncompresses it.

        This function behaves identically to fetch_from_s3 when:

            * src is a tar file, or

            * auto_unzip=False, or

            * auto_unzip=True and src is a compressed file, or

            * auto_unzip=True, src is NOT a compressed file, and no compressed version of src exists in s3

        We generally try to keep all our references in S3 either lz4-compressed or tar'ed, because this
        guards against corruption errors that sometimes occur during download.  If we accidentally
        download a corrupt file, it would fail to unlz4 or to untar, and we would retry the download.
    '''
    if auto_unzip and not src.endswith(".tar") and not any(src.endswith(zext) for zext in ZIP_EXTENSIONS):
        # Try to guess which compressed version of the file exists in S3;  then download and decompress it.
        for zext in ZIP_EXTENSIONS:
            try:
                return fetch_from_s3(src + zext,
                                     dst,
                                     auto_unzip=auto_unzip,
                                     auto_untar=auto_untar,
                                     allow_s3mi=allow_s3mi)
            except:
                # It's okay - if all these fail, we'll try to fetch src directly.
                pass
    return fetch_from_s3(src,
                         dst,
                         auto_unzip=auto_unzip,
                         auto_untar=auto_untar,
                         allow_s3mi=allow_s3mi)


def fetch_byterange(first_byte, last_byte, bucket, key, output_file):
    get_range = f"aws s3api get-object " \
                f"--range bytes={first_byte}-{last_byte} " \
                f"--bucket {bucket} " \
                f"--key {key} {output_file}"
    get_range_proc = subprocess.Popen(get_range, shell=True, stdout=subprocess.DEVNULL)
    get_range_proc.wait()


@command.retry
def upload_with_retries(from_f, to_f):
    command.execute(f"aws s3 cp --only-show-errors {from_f} {to_f}")

@command.retry
def upload_folder_with_retries(from_f, to_f):
    command.execute(f"aws s3 cp --only-show-errors --recursive {from_f.rstrip('/')}/ {to_f.rstrip('/')}/")

def upload(from_f, to_f, status, status_lock=TraceLock("upload", threading.RLock())):
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


def upload_log_file(sample_s3_output_path, lock=TraceLock("upload_log_file", threading.RLock())):
    with lock:
        logh = logging.getLogger().handlers[0]
        logh.flush()
        command.execute(f"aws s3 cp --only-show-errors {logh.baseFilename} {sample_s3_output_path}/")
