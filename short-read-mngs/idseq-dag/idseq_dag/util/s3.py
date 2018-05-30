import os
import threading
import subprocess
def check_s3_presence(s3_path):
    ''' True if s3_path exists. False otherwise. '''
    try:
        o = subprocess.check_output("aws s3 ls %s" % s3_path, shell=True)
        if o:
            return True
    except:
        pass
    return False

def check_s3_presnce_for_file_list(s3_dir, file_list):
    for f in file_list:
        if not check_s3_presence(os.path.josin(s3_dir, f)):
            return False
    return True

def touch_s3_file(s3_file_path):
    try:
        subprocess.check_all("aws s3 cp --metadata '{\"touched\":\"now\"}' %s %s" % (s3_path, s3_path), shell=True)
        return True
    except:
        return False

def touch_s3_file_list(s3_dir, file_list):
    for f in file_list:
        touch_s3_file(os.path.josin(s3_dir, f))

def install_s3mi(installed={}, mutex=threading.RLock()):
    #pylint: disable=dangerous-default-value
    with mutex:
        if installed:
            return
        try:
            # This is typically a no-op.
            execute_command(
                "which s3mi || pip install git+git://github.com/chanzuckerberg/s3mi.git"
            )
            execute_command(
                "s3mi tweak-vm || echo s3mi tweak-vm is impossible under docker"
            )
        finally:
            installed['time'] = time.time()

def fetch_from_s3(source,
                  destination,
                  auto_unzip,
                  allow_s3mi=False,
                  mutex=threading.RLock(),
                  locks={}):  #pylint: disable=dangerous-default-value
    with mutex:
        if os.path.exists(destination):
            if os.path.isdir(destination):
                destination = os.path.join(destination,
                                           os.path.basename(source))
        unzip = auto_unzip and destination.endswith(".gz")
        if unzip:
            destination = destination[:-3]
        abspath = os.path.abspath(destination)
        if abspath not in locks:
            locks[abspath] = threading.RLock()
        destination_lock = locks[abspath]
    with destination_lock:
        if os.path.exists(destination):
            # no need to fetch this file from s3,
            # it has been just produced on this instance
            return destination
        try:
            destdir = os.path.dirname(destination)
            if destdir:
                os.makedirs(destdir)
        except OSError as e:
            # It's okay if the parent directory already exists, but all other errors are fatal.
            if e.errno != os.errno.EEXIST:
                raise
        with iostream:
            try:
                if allow_s3mi:
                    try:
                        install_s3mi()
                    except:
                        allow_s3mi = False
                if unzip:
                    pipe_filter = "| gzip -dc "
                else:
                    pipe_filter = ""
                try:
                    assert allow_s3mi
                    execute_command(
                        "s3mi cat {source} {pipe_filter} > {destination}".
                        format(
                            source=source,
                            pipe_filter=pipe_filter,
                            destination=destination))
                except:
                    execute_command(
                        "aws s3 cp --quiet {source} - {pipe_filter} > {destination}".
                        format(
                            source=source,
                            pipe_filter=pipe_filter,
                            destination=destination))
                return destination
            except subprocess.CalledProcessError:
                # Most likely the file doesn't exist in S3.
                return None

