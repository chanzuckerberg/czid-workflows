import time
import subprocess
import os
import multiprocessing
import errno
import re
import json
from urllib.parse import urlparse
import traceback
import boto3
import botocore.exceptions
import botocore.session
from idseq_dag.util.trace_lock import TraceLock
import idseq_dag.util.command_patterns as command_patterns

import idseq_dag.util.command as command
import idseq_dag.util.log as log

# Peak network and storage perf for a typical small instance is saturated by
# just a few concurrent streams.
MAX_CONCURRENT_COPY_OPERATIONS = 12
IOSTREAM = multiprocessing.Semaphore(MAX_CONCURRENT_COPY_OPERATIONS)
# Make a second semaphore for uploads to reserve some capacity for downloads.
MAX_CONCURRENT_UPLOAD_OPERATIONS = 8
MAX_CONCURRENT_S3MI_DOWNLOADS = 2
# If a s3mi slot does not free up within MAX_S3MI_WAIT seconds, we use plain old aws s3 cp instead of s3mi.
MAX_S3MI_WAIT = 15
S3MI_SEM = multiprocessing.Semaphore(MAX_CONCURRENT_S3MI_DOWNLOADS)
IOSTREAM_UPLOADS = multiprocessing.Semaphore(MAX_CONCURRENT_UPLOAD_OPERATIONS)

config = {
    # Configured in idseq_dag.engine.pipeline_flow.PipelineFlow
    "REF_DIR": None,
    "PURGE_SENTINEL": None
}

def split_identifiers(s3_path):
    return s3_path[5:].split("/", 1)


# Boto's default session is global and shared across threads.   As long as this is the
# only place that we use boto and we never call functions here from multiple processes
# (threads okay) and all accesses here are protected by this lock, we should be fine.
# But still easy to defeat this safety by just importing boto into another file.
# TODO:  Create our own private boto session protected by this lock.
botolock = multiprocessing.RLock()  # Please don't trace me, I am very lightweight.


def rate_limit_boto(average_delay=0.33, last_call={}):  # pylint: disable=dangerous-default-value
    '''Sleep to ensure at least average_delay has elapsed since the last call.'''
    with botolock:
        t_since_last_call = time.time() - last_call.get("time", 0)
        if t_since_last_call < average_delay:
            time.sleep(average_delay - t_since_last_call)
        last_call["time"] = time.time()


def _check_s3_presence(s3_path, allow_zero_byte_files):
    """True if s3_path exists. False otherwise."""
    with log.log_context(context_name="s3.check_s3_presence", values={'s3_path': s3_path}, log_context_mode=log.LogContextMode.EXEC_LOG_EVENT) as lc:
        parsed_url = urlparse(s3_path, allow_fragments=False)
        bucket = parsed_url.netloc
        key = parsed_url.path.lstrip('/')
        try:
            o = boto3.resource('s3').Object(
                bucket,
                key
            )
            size = o.content_length
            lc.values['size'] = size
            exists = (allow_zero_byte_files and size >= 0) or (not allow_zero_byte_files and size > 0)
        except botocore.exceptions.ClientError as e:
            if e.response['Error']['Code'] == "404":
                exists = False
            else:
                # Something else has gone wrong.
                raise
        lc.values['exists'] = exists
        return exists


def check_s3_presence(s3_path, allow_zero_byte_files=True):
    with botolock:
        rate_limit_boto()
        return _check_s3_presence(s3_path, allow_zero_byte_files)


def get_s3_object_by_path(s3_path):
    parsed_url = urlparse(s3_path, allow_fragments=False)
    bucket = parsed_url.netloc
    key = parsed_url.path.lstrip('/')
    try:
        o = boto3.resource('s3').Object(
            bucket,
            key
        )
        return o.get()['Body'].read()
    except botocore.exceptions.ClientError as e:
        if e.response['Error']['Code'] == "404" or e.response['Error']['Code'] == 'NoSuchKey':
            return None
        else:
            # raise all others
            raise


def check_s3_presence_for_file_list(s3_dir, file_list, allow_zero_byte_files=True):
    with botolock:
        for f in file_list:
            rate_limit_boto()
            if not _check_s3_presence(os.path.join(s3_dir, f), allow_zero_byte_files):
                return False
    return True


@command.retry
def _get_credentials():
    log.write("Refreshing credentials.")
    session = botocore.session.Session()
    credentials = session.get_credentials()
    return {
        "AWS_ACCESS_KEY_ID": credentials.access_key,
        "AWS_SECRET_ACCESS_KEY": credentials.secret_key,
        "AWS_SESSION_TOKEN": credentials.token,
        "AWS_DEFAULT_REGION": session.create_client("s3").meta.region_name
    }


def refreshed_credentials(credentials_mutex=multiprocessing.RLock(), credentials_cache={}):  # pylint: disable=dangerous-default-value
    with credentials_mutex:
        if credentials_cache.get("expiration_time", 0) < time.time() + 5 * 60:
            try:
                credentials_cache["vars"] = _get_credentials()
            except:
                log.write("ERROR:  Failed to refresh credentials with boto, even after retries.  Subcommands will have to do it themselves.")
                log.write(traceback.format_exc())
                return {}
            credentials_cache["expiration_time"] = time.time() + 15 * 60  # this is the default for most accounts
    return credentials_cache["vars"]


def find_oldest_reference(refdir):
    try:
        # To understand this ls command, please see the path forming explained in fetch_from_s3 when is_reference=True.
        ls = subprocess.run(f"ls -td {refdir}/*/* | tail -1", shell=True, executable="/bin/bash", check=True, capture_output=True)
    except:
        # This happens when there are no reference files to delete.
        return None
    result = ls.stdout.decode('utf-8').strip()
    assert config["PURGE_SENTINEL"] != None
    if result == config["PURGE_SENTINEL"]:
        log.write("WARNING:  Encountered purge sentinel before freeing sufficient space.  Job may run out of space.")
        return None
    return result


def need_more_space(refdir):
    try:
        df = subprocess.run(f"df -m {refdir}" + " | awk '{print $5}' | tail -1 | sed 's=%=='", shell=True, executable="/bin/bash", check=True, capture_output=True)
        percent_used = int(df.stdout.decode('utf-8').strip())
        log.write(f"Disk used space: {percent_used} percent.")
        return percent_used > 60
    except:
        log.write("Error:  Failed to determine available space on instance.  Will assume too little space is available, and will try to delete least recently used reference downloads.")
        traceback.format_exc()
        return True  # Better safe than sorry.


def really_make_space():
    refdir = config['REF_DIR']
    while need_more_space(refdir):
        lru = find_oldest_reference(refdir)
        if not lru:
            log.write("WARNING:  Too little available space on instance, and could not find any reference downloads to delete.")
            break
        command.remove_rf(lru)


def make_space(done={}, mutex=TraceLock("make_space", multiprocessing.RLock())):  # pylint: disable=dangerous-default-value
    with mutex:
        if not done:
            try:
                really_make_space()
            except:
                log.write("Error making space.  Please attend to this before instance storage fills up.")
                log.write(traceback.format_exc())
            done['time'] = time.time()


DEFAULT_AUTO_UNZIP = False
DEFAULT_AUTO_UNTAR = False
DEFAULT_ALLOW_S3MI = False
ZIP_EXTENSIONS = {
    ".lz4": "lz4 -dc",
    ".bz2": "lbzip2 -dc",
    ".gz": "gzip -dc",  # please avoid .gz if you can (slow to decompress)
}


# this should be a subset of ZIP_EXTENSIONS, and a very small one, because
# there is a high cost to each new entry (one more s3 op per reference)
REFERENCE_AUTOGUESS_ZIP_EXTENSIONS = {
    ".lz4": "lz4 -dc",
}


# WARNING: This will bypass download if src is a local path, see comment below for details
def fetch_from_s3(src,  # pylint: disable=dangerous-default-value
                  dst,
                  auto_unzip=DEFAULT_AUTO_UNZIP,
                  auto_untar=DEFAULT_AUTO_UNTAR,
                  allow_s3mi=DEFAULT_ALLOW_S3MI,
                  okay_if_missing=False,
                  is_reference=False,
                  touch_only=False,
                  mutex=TraceLock("fetch_from_s3", multiprocessing.RLock()),
                  locks={}):
    """Fetch a file from S3 if needed, using either s3mi or aws cp.

    IT IS NOT SAFE TO CALL THIS FUNCTION FROM MULTIPLE PROCESSES.
    It is totally fine to call it from multiple threads (it is designed for that).

    When is_reference=True, "dst" must be an existing directory.

    If src does not exist or there is a failure fetching it, the function returns None,
    without raising an exception.  If the download is successful, it returns the path
    to the downloaded file or folder.  If the download already exists, it is touched
    to update its timestamp.

    When touch_only=True, if the destination does not already exist, the function
    simply returns None (as if the download failed).  If the destination does exist,
    it is touched as usual.  This is useful in implementing an LRU cache policy.

    An exception is raised only if there is a coding error or equivalent problem,
    not if src simply doesn't exist.
    """
    # FIXME: this is a compatibility hack so we can replace this function
    #   We are removing ad-hoc s3 downloads from within steps and converting
    #   additional_files to wdl inputs. These files will be transparently
    #   downloaded by miniwdl. miniwdl will also handle the caching that
    #   is currently done here. This hack bypasses the s3 download if the
    #   source is already a local file, and returns the source (which is
    #   a local file path). This way, when we change the additional_files
    #   to inputs we can provide the local file path to the step instead
    #   of the s3 path and seamlessly transition without a coordinated
    #   change between idseq-dag and the idseq monorepo.
    if not src.startswith("s3://"):
        log.write(f"fetch_from_s3 is skipping download because source: {src} does not start with s3://")
        if not os.path.isfile(src):
            return None
        if auto_untar and src.endswith(".tar"):
            dst = src[:-4]
            if not os.path.isdir(dst):
                command.make_dirs(dst + ".untarring")
                script = 'tar xvf "${src}" -C "${tmp_destdir}"'
                named_args = {
                    "src": src,
                    "tmp_destdir": dst + ".untarring"
                }
                command.execute(
                    command_patterns.ShellScriptCommand(
                        script=script,
                        named_args=named_args
                    )
                )
                command.rename(dst + ".untarring/" + os.path.basename(dst), dst)
            return dst
        return src

    # Do not be mislead by the multiprocessing.RLock() above -- that just means it won't deadlock
    # if called from multiple processes but does not mean the behaivior will be correct.  It will
    # be incorrect, because the locks dict (cointaining per-file locks) cannot be shared across
    # processes, the way it can be shared across threads.

    if is_reference:
        assert config["REF_DIR"], "The is_reference code path becomes available only after initializing gloabal config['REF_DIR']"

    if os.path.exists(dst) and os.path.isdir(dst):
        dirname, basename = os.path.split(src)
        if is_reference or os.path.abspath(dst).startswith(config["REF_DIR"]):
            # Downloads to the reference dir are persisted from job to job, so we must include
            # version information from the full s3 path.
            #
            # The final destination for s3://path/to/source.db will look like /mnt/ref/s3__path__to/source.db
            # The final destination for s3://path/to/myarchive.tar will look like /mnt/ref/s3__path__to/myarchive/...
            #
            # We considered some other alternatives, for example /mnt/ref/s3__path__to__source.db, but unfortunately,
            # some tools incorporate the base name of their database input into the output filenames, so any approach
            # that changes the basename causes problems downstream.  An example such tool is srst2.
            is_reference = True
            if dirname.startswith("s3://"):
                dirname = dirname.replace("s3://", "s3__", 1)
            # If dirname contains slashes, it has to be flattened to single level.
            dirname = dirname.replace("/", "__")
            dst = os.path.join(dst, dirname, basename)
        else:
            dst = os.path.join(dst, basename)
    else:
        assert not is_reference, f"When fetching references, dst must be an existing directory: {dst}"

    unzip = ""
    if auto_unzip:
        file_without_ext, ext = os.path.splitext(dst)
        if ext in ZIP_EXTENSIONS:
            unzip = " | " + ZIP_EXTENSIONS[ext]  # this command will be used to decompress stdin to stdout
            dst = file_without_ext  # remove file extension from dst
    untar = auto_untar and dst.lower().endswith(".tar")
    if untar:
        dst = dst[:-4]  # Remove .tar

    # Downloads are staged under tmp_destdir.  Only after a download completes successfully it is moved to dst.
    destdir = os.path.dirname(dst)
    tmp_destdir = os.path.join(destdir, "tmp_downloads")
    tmp_dst = os.path.join(tmp_destdir, os.path.basename(dst))

    abspath = os.path.abspath(dst)
    with mutex:
        if abspath not in locks:
            locks[abspath] = TraceLock(f"fetch_from_s3: {abspath}", multiprocessing.RLock())
        destination_lock = locks[abspath]

    # shouldn't happen and makes it impossible to ensure that any dst that exists is complete and correct.
    assert tmp_dst != dst, f"Problematic use of fetch_from_s3 with tmp_dst==dst=='{dst}'"

    with destination_lock:
        # This check is a bit imperfect when untarring... unless you follow the discipline that
        # all contents of file foo.tar are under directory foo/... (which we do follow in IDseq)
        if os.path.exists(dst):
            command.touch(dst)
            return dst

        if touch_only:
            return None

        for (kind, ddir) in [("destinaiton", destdir), ("temporary download", tmp_destdir)]:
            try:
                if ddir:
                    command.make_dirs(ddir)
            except OSError as e:
                # It's okay if the parent directory already exists, but all other
                # errors fail the download.
                if e.errno != errno.EEXIST:
                    log.write(f"Error in creating {kind} directory.")
                    return None

        with IOSTREAM:
            try:
                if allow_s3mi:
                    wait_start = time.time()
                    allow_s3mi = S3MI_SEM.acquire(timeout=MAX_S3MI_WAIT)
                    wait_duration = time.time() - wait_start
                    if not allow_s3mi:
                        log.write(f"Failed to acquire S3MI semaphore after waiting {wait_duration} seconds for {src}.")
                    elif wait_duration >= 5:
                        log.write(f"Waited {wait_duration} seconds to acquire S3MI semaphore for {src}.")

                if untar:
                    write_dst = r''' | tar xvf - -C "${tmp_destdir}";'''
                    named_args = {'tmp_destdir': tmp_destdir}
                else:
                    write_dst = r''' > "${tmp_dst}";'''
                    named_args = {'tmp_dst': tmp_dst}
                command_params = f"{unzip} {write_dst}"

                named_args.update({'src': src})

                try_cli = not allow_s3mi
                if allow_s3mi:
                    if os.path.exists(tmp_dst):
                        command.remove_rf(tmp_dst)
                    try:
                        command.execute(
                            command_patterns.ShellScriptCommand(
                                script=r'set -o pipefail; s3mi cat --quiet "${src}" ' + command_params,
                                named_args=named_args
                            )
                        )
                    except subprocess.CalledProcessError:
                        try_cli = not okay_if_missing
                        allow_s3mi = False
                        S3MI_SEM.release()
                        if try_cli:
                            log.write(
                                "Failed to download with s3mi. Trying with aws s3 cp..."
                            )
                        else:
                            raise
                if try_cli:
                    if os.path.exists(tmp_dst):
                        command.remove_rf(tmp_dst)
                    if okay_if_missing:
                        script = r'set -o pipefail; aws s3 cp --quiet "${src}" - ' + command_params
                    else:
                        script = r'set -o pipefail; aws s3 cp --only-show-errors "${src}" - ' + command_params
                    command.execute(
                        command_patterns.ShellScriptCommand(
                            script=script,
                            named_args=named_args,
                            env=dict(os.environ, **refreshed_credentials())
                        )
                    )
                # Move staged download into final location.  Leave this last, so it only happens if no exception has occurred.
                # By this point we have already asserted that tmp_dst != dst.
                command.rename(tmp_dst, dst)
                return dst
            except BaseException as e:  # Deliberately super broad to make doubly certain that dst will be removed if there has been any exception
                if os.path.exists(dst):
                    command.remove_rf(dst)
                if not isinstance(e, subprocess.CalledProcessError):
                    # Coding error of some sort.  Best not hide it.
                    raise
                if okay_if_missing:
                    # We presume.
                    log.write(
                        "File most likely does not exist in S3."
                    )
                else:
                    log.write(
                        "Failed to fetch file from S3."
                    )
                return None
            finally:
                if allow_s3mi:
                    S3MI_SEM.release()
                if os.path.exists(tmp_dst):  # by this point we have asserted that tmp_dst != dst (and that assert may have failed, but so be it)
                    command.remove_rf(tmp_dst)

# WARNING: This will bypass download if src is a local path, see comment in fetch_from_s3 for details
def fetch_reference(src,  # pylint: disable=dangerous-default-value
                    dst,
                    auto_unzip=True,
                    auto_untar=True,
                    allow_s3mi=DEFAULT_ALLOW_S3MI,
                    touch_only=False):
    '''
        This function behaves like fetch_from_s3 in most cases, with one excetpion:

            *  When auto_unzip=True, src is NOT a compressed file and NOT a tar file, and a compressed
               version of src exists in s3, and touch_only=False, this function downloads and uncompresses it.

        This function behaves identically to fetch_from_s3 when:

            * src is a tar file, or

            * auto_unzip=False, or

            * auto_unzip=True and src is a compressed file, or

            * auto_unzip=True, src is NOT a compressed file, and no compressed version of src exists in s3, or

            * touch_only=True

        We generally try to keep all our references in S3 either lz4-compressed or tar'ed, because this
        guards against corruption errors that sometimes occur during download.  If we accidentally
        download a corrupt file, it would fail to unlz4 or to untar, and we would retry the download.

        We trust "aws s3 cp" to do its own integrity checking, so we only prefer the lz4 version of a file
        if we are allowed to use s3mi.
    '''
    if not touch_only and auto_unzip and not src.endswith(".tar") and not any(src.endswith(zext) for zext in ZIP_EXTENSIONS):
        # Try to guess which compressed version of the file exists in S3;  then download and decompress it.
        for zext in REFERENCE_AUTOGUESS_ZIP_EXTENSIONS:
            result = fetch_from_s3(src + zext,
                                   dst,
                                   auto_unzip=auto_unzip,
                                   auto_untar=auto_untar,
                                   allow_s3mi=allow_s3mi,
                                   okay_if_missing=True,  # It's okay if missing a compressed version.
                                   is_reference=True,
                                   touch_only=touch_only)
            if result:
                return result

    return fetch_from_s3(src,
                         dst,
                         auto_unzip=auto_unzip,
                         auto_untar=auto_untar,
                         allow_s3mi=allow_s3mi,
                         okay_if_missing=False,  # It's NOT okay if missing on the last attempt.
                         is_reference=True,
                         touch_only=touch_only)


@command.retry
def upload_with_retries(from_f, to_f, checksum=False):
    with IOSTREAM_UPLOADS:
        with IOSTREAM:
            args = []
            if checksum:
                args.append("--checksum")
            args.append(from_f)
            args.append(to_f)
            command.execute(
                command_patterns.SingleCommand(
                    cmd="s3parcp",
                    args=args,
                    env=dict(os.environ, **refreshed_credentials())
                )
            )

@command.retry
def upload_folder_with_retries(from_f, to_f, checksum=False):
    with IOSTREAM_UPLOADS:
        with IOSTREAM:
            args = ["--recursive"]
            if checksum:
                args.append("--checksum")
            args.append(os.path.join(from_f, ""))
            args.append(os.path.join(to_f, ""))
            command.execute(
                command_patterns.SingleCommand(
                    cmd="s3parcp",
                    args=args,
                    env=dict(os.environ, **refreshed_credentials())
                )
            )
