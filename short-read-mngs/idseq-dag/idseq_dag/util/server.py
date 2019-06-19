import os
import json
import time
import random
import threading
import traceback
import multiprocessing

from contextlib import contextmanager
from collections import defaultdict
from copy import deepcopy

from idseq_dag.util.periodic_thread import PeriodicThread

import idseq_dag.util.command as command
import idseq_dag.util.log as log


MIN_INTERVAL_BETWEEN_DESCRIBE_INSTANCES = 180
MAX_POLLING_LATENCY = 10  # seconds
MAX_INSTANCES_TO_POLL = 8
MAX_DISPATCHES_PER_MINUTE = 10


def _duct_tape(service_name):
     # TODO: remove this duct tape
    if service_name == 'rapsearch2':
        return 'rapsearch'
    return service_name


class ServiceUnreliable(RuntimeError):
    ''' Raised when the probability and cost of failures become too high relative to success. '''

    def __init__(self, service, chunk_id):
        super().__init__(f"Service {_duct_tape(service)} is unreliable, giving up for chunk {chunk_id}.")


class ChunkStatus:
    SUCCESS = "success"
    CRASH = "crash"
    TIMEOUT = "timeout"
    CORRUPT_OUTPUT = "corrupt_output"


class ChunkStatusTracker:

    def __init__(self, service):
        self.lock = multiprocessing.RLock()
        self.service = service
        # tally[status][server_ip] contains the list of (chunk_id, elapsed, try_number) for that server and status.
        # TODO:  A flatter map or list might be even better here.
        self.outcome_data = {
            ChunkStatus.SUCCESS: defaultdict(list),
            ChunkStatus.CRASH: defaultdict(list),
            ChunkStatus.TIMEOUT: defaultdict(list),
            ChunkStatus.CORRUPT_OUTPUT: defaultdict(list)
        }
        # chunks_waiting[try_number] is the set of chunk_ids that are waiting for an ASGinstance for this try_number
        # (initially each chunk is at try_number=1, then 2, etc...)
        # This helps us ensure that most chunks try at least once before any retries.
        self.chunks_waiting_lock = multiprocessing.RLock()
        self.chunks_waiting = defaultdict(set)

    def chunk_has_try_priority(self, _chunk_id, try_number):
        # Note try_number starts counting from 1.
        with self.chunks_waiting_lock:
            for prior_try in range(1, try_number):
                if self.chunks_waiting[prior_try]:
                    return False
        return True

    def register_chunk_waiting(self, chunk_id, try_number):
        with self.chunks_waiting_lock:
            self.chunks_waiting[try_number].add(chunk_id)

    def register_chunk_dispatched(self, chunk_id, try_number):
        with self.chunks_waiting_lock:
            self.chunks_waiting[try_number].discard(chunk_id)
            # TODO:  notify if we removed the last element of try_number,
            # when we remove the busy wait.

    def note_outcome(self, instance_ip, chunk_id, elapsed, status, try_number):
        with self.lock:
            self.outcome_data[status][instance_ip].append((chunk_id, elapsed, try_number))

    def tally_costs_by_outcome(self):
        with self.lock:
            outcome_counts = defaultdict(int)
            outcome_elapsed = defaultdict(float)
            for status, outcomes in self.outcome_data.items():
                for _instance, instance_outcomes in outcomes.items():
                    for _chunk_id, chunk_elapsed, _try_number in instance_outcomes:  # list of pairs
                        outcome_counts[status] += 1
                        outcome_elapsed[status] += chunk_elapsed
            return outcome_counts, outcome_elapsed

    def service_unreliable(self):
        outcome_counts, outcome_elapsed = self.tally_costs_by_outcome()
        total_count = sum(outcome_counts.values())
        total_elapsed = sum(outcome_elapsed.values())
        success_count = outcome_counts[ChunkStatus.SUCCESS]
        fail_count = total_count - success_count
        avg_elapsed_fail = (total_elapsed - outcome_elapsed[ChunkStatus.SUCCESS]) / fail_count if fail_count else 0.0
        avg_elapsed_success = outcome_elapsed[ChunkStatus.SUCCESS] / success_count if success_count else 0.0
        # If after 12 tries the probability of success is < 1/3, we give up.
        if total_count >= 12 and success_count < (total_count / 3):
            return True
        # If after 18 tries the average time to fail is large compared to the average time to succeed,
        # and failures are common, then it's too expensive to keep trying, so we give up.
        if total_count >= 18 and success_count < (total_count * 2 / 3) and avg_elapsed_fail >= 0.5 * avg_elapsed_success:
            return True
        # Ok, keep trying.
        return False

    def servers_with_rep(self, chunk_id):
        rep = defaultdict(int)
        terrible = set()
        with self.lock:
            for status, status_outcomes in self.outcome_data.items():
                if status == ChunkStatus.SUCCESS:
                    # We are only interested in the failures here.
                    continue
                for instance_ip, instance_outcomes in status_outcomes.items():
                    num_successes_for_ip = len(self.outcome_data[ChunkStatus.SUCCESS].get(instance_ip, []))
                    for failed_chunk_id, _elapsed, _try_number in instance_outcomes:  # this is a list
                        rep[instance_ip] += 1
                        # If this instance has had 3 failures and not a single success, we declare it "terrible"
                        # so it gets avoided.  If failures outnumber successes 4:1, also terrible.
                        if rep[instance_ip] >= 3 and num_successes_for_ip == 0:
                            terrible.add(instance_ip)
                        if rep[instance_ip] >= 4 * num_successes_for_ip and num_successes_for_ip > 0:
                            terrible.add(instance_ip)
                        # If chunk_id previously failed on this ip, we declare the ip "terrible" so we don't retry
                        # the chunk on the same server;  unless the failure was corrupt output, for which we allow
                        # retries on the same server.
                        if chunk_id == failed_chunk_id and status != ChunkStatus.CORRUPT_OUTPUT:
                            terrible.add(instance_ip)
        slightly_suspect = set(rep.keys()) - terrible
        # NOTE:  When the failure is a timeout, that doesn't mean that anything is wrong with the instance!
        # Some samples just need more time in rapsearch, and the appropriate behavior is to time-out their chunks
        # then rerun with smaller chunks after more stringent QC has removed most reads causing the slowness.
        # So "terrible" here means terrible FOR THIS SAMPLE's data;  not necessarily terrible for the other samples
        # that may be running at this time.  On the other hand, it COULD mean terrible for all samples -- if
        # the machine has become unhealthy -- but we just can't leap to that conclusion.
        return slightly_suspect, terrible

    def status_report(self, total_chunks):
        # This is called only after the last call to self.note_outcome(),
        # so we don't have to bother with self.lock.
        succeeded_first_try = set()
        succeeded_retry = set()
        attempted = set()
        elapsed_succcess = 0.0
        elapsed_failure = 0.0
        server_successes = defaultdict(int)
        server_failures = defaultdict(int)
        num_failures = 0
        for status, outcomes in self.outcome_data.items():
            for instance_ip, instance_outcomes in outcomes.items():
                for chunk_id, elapsed, try_number in instance_outcomes:
                    attempted.add(chunk_id)
                    if status == ChunkStatus.SUCCESS:
                        elapsed_succcess += elapsed
                        server_successes[instance_ip] += 1
                        if try_number == 1:
                            succeeded_first_try.add(chunk_id)
                        else:
                            succeeded_retry.add(chunk_id)
                    else:
                        elapsed_failure += elapsed
                        server_failures[instance_ip] += 1
                        num_failures += 1
        failed = attempted - succeeded_first_try - succeeded_retry
        elapsed_avg_success = 0.0
        elapsed_avg_failure = 0.0
        if succeeded_first_try or succeeded_retry:
            elapsed_avg_success = elapsed_succcess / (len(succeeded_first_try) + len(succeeded_retry))
        if num_failures:
            elapsed_avg_failure = elapsed_failure / num_failures
        server_stats = {}
        for instance_ip in sorted(set(server_successes.keys()) | set(server_failures.keys())):
            server_stats[instance_ip] = {
                'SUCCESSES': server_successes.get(instance_ip, 0),
                'FAILURES': server_failures.get(instance_ip, 0)  # this way we don't affect len(server_failures)
            }
        report = {
            'TOTAL_CHUNKS': total_chunks,
            'CHUNKS_ATTEMPTED': len(attempted),
            'CHUNKS_SUCCEEDED_ON_FIRST_TRY': len(succeeded_first_try),
            'CHUNKS_SUCCEEDED_ON_RETRY': len(succeeded_retry),
            'CHUNKS_FAILED': len(failed),
            'CHUNKS_ELAPSED_AVG_SUCCESS': elapsed_avg_success,
            'CHUNKS_ELAPSED_AVG_FAILURE': elapsed_avg_failure,
            'TOTAL_SERVERS': len(server_stats),
            'SERVERS_WITH_FAILURES': len(server_failures),
            'SERVER_STATS': server_stats
        }
        return report

    def log_stats(self, total_chunks):
        # This happens after the last call to tracker.note_outcome(), making it
        # safe to access tracker.outcome_data directly without holding a lock.
        # But, to be safe...
        with self.lock:
            report = self.status_report(total_chunks)
            outcomes = deepcopy(self.outcome_data)
        log.log_event(f"STATUS REPORT FOR {self.service}", values=report)
        log.log_event(f"FULL OUTCOME DATA FOR {self.service}", values=outcomes)


def chunk_status_tracker(service, #pylint: disable=dangerous-default-value
                         trackers={
                             "gsnap": ChunkStatusTracker("gsnap"),
                             "rapsearch": ChunkStatusTracker("rapsearch")
                         }):
    return trackers[_duct_tape(service)]


@command.retry
def get_server_ips_work(service_name, environment, draining_tag):
    ''' return a dict of relevant instance IPs to instance IDs '''
    value = "%s-%s" % (service_name, environment)
    describe_json = json.loads(
        command.execute_with_output(
            f"aws ec2 describe-instances --filters 'Name=tag:service,Values={value}' 'Name=instance-state-name,Values=running'"))
    server_ips = {
        instance["NetworkInterfaces"][0]["PrivateIpAddress"]: instance["InstanceId"]
        for reservation in describe_json["Reservations"]
        for instance in reservation["Instances"]
        if draining_tag not in [tag["Key"] for tag in instance["Tags"]]
    }
    return server_ips


def get_server_ips(service_name, #pylint: disable=dangerous-default-value
                   environment,
                   max_interval_between_describe_instances,
                   draining_tag,
                   aggressive=False,
                   cache={},
                   mutex=threading.RLock()):
    try:
        with mutex:
            if aggressive:
                period = MIN_INTERVAL_BETWEEN_DESCRIBE_INSTANCES
            else:
                period = max_interval_between_describe_instances
            now = time.time()
            cache_key = (service_name, environment)
            if cache_key not in cache or now - cache[cache_key][0] >= period:
                # this may raise an exception when the AWS account rate limit is exceeded due to many concurrent jobs
                cache[cache_key] = (now,
                                    get_server_ips_work(
                                        service_name, environment, draining_tag))
            return cache[cache_key][1]
    except:
        # return [] causes a sleep of wait_seconds before retrying (see below)
        traceback.print_exc()
        return []


def wait_for_server_ip_work(service_name, #pylint: disable=dangerous-default-value
                            key_path,
                            remote_username,
                            environment,
                            max_concurrent,
                            chunk_id,
                            max_interval_between_describe_instances,
                            draining_tag,
                            had_to_wait=[False]):
    while True:
        log.write(f"Chunk {chunk_id} of {service_name} is at third gate")
        instance_ip_id_dict = get_server_ips(
            service_name, environment,
            max_interval_between_describe_instances, draining_tag,
            aggressive=had_to_wait[0])
        all_servers = set(instance_ip_id_dict.keys())
        slightly_suspect, terrible = chunk_status_tracker(service_name).servers_with_rep(chunk_id)
        slightly_suspect &= all_servers
        terrible &= all_servers
        presumed_good = all_servers - terrible
        instance_ips = random.sample(presumed_good,
                                     min(MAX_INSTANCES_TO_POLL,
                                         len(presumed_good)))
        ip_nproc_dict = {}
        dict_mutex = threading.RLock()
        dict_writable = True

        def poll_server(ip):
            # ServerAliveInterval to fix issue with containers keeping open
            # an SSH connection even after worker machines had finished
            # running.
            commands = "ps aux | grep %s | grep -v bash || echo error" % service_name
            output = command.execute_with_output(
                command.remote(commands, key_path, remote_username, ip),
                timeout=MAX_POLLING_LATENCY).rstrip().split("\n")
            if output != ["error"]:
                with dict_mutex:
                    if dict_writable:
                        ip_nproc_dict[ip] = len(output) - 1

        poller_threads = []
        for ip in instance_ips:
            pt = threading.Thread(target=poll_server, args=[ip])
            pt.start()
            poller_threads.append(pt)
        for pt in poller_threads:
            pt.join(MAX_POLLING_LATENCY)
        with dict_mutex:
            # Any zombie threads won't be allowed to modify the dict.
            dict_writable = False
        if len(ip_nproc_dict) < len(poller_threads):
            log.write(
                "Only {} out of {} instances responded to polling;  {} threads are timing out.".
                format(
                    len(ip_nproc_dict), len(poller_threads),
                    len([pt for pt in poller_threads if pt.isAlive()])))
        log.write(f"Chunk {chunk_id} of {service_name} is at fourth gate")
        if not ip_nproc_dict:
            have_capacity = False
        else:
            min_nproc_ip = min(ip_nproc_dict, key=ip_nproc_dict.get)
            min_nproc = ip_nproc_dict[min_nproc_ip]
            have_capacity = (min_nproc < max_concurrent)
        if have_capacity:
            had_to_wait[0] = False
            # Make an urn where each ip occurs more times if it has more free slots, and not at all if it lacks free slots
            urn = []
            for ip, nproc in ip_nproc_dict.items():
                free_slots = max_concurrent - nproc
                if free_slots > 0:
                    # TODO:  If we are concerned about cost, we could bias this distribution toward the middle,
                    # so machines that have no jobs running are not preferred so much (this way we can
                    # reap them in the autoscaler sooner).   The only reason for not doing that yet is
                    # that we have a race condition (by design) here and that would increase the overload rate.
                    # TODO:  For reliability, we can factor in the historic success rate and efficiency
                    # of the machine in addition to its current load.
                    weight = 2**free_slots - 1
                    if ip in slightly_suspect:
                        # minimize the probability of choosing a server with a bad reputation
                        weight = 1
                    urn.extend([ip] * weight)
            min_nproc_ip = random.choice(urn)
            free_slots = max_concurrent - ip_nproc_dict[min_nproc_ip]
            log.write(f"{service_name} server {min_nproc_ip} has capacity {free_slots}. Kicking off ")
            return min_nproc_ip, instance_ip_id_dict[min_nproc_ip]
        else:
            if len(terrible) > len(all_servers) * 2 / 3.0:
                raise ServiceUnreliable(service_name, chunk_id)
            had_to_wait[0] = True
            wait_seconds = random.randint(
                max(20, MAX_POLLING_LATENCY), max(60, MAX_POLLING_LATENCY))
            log.write(f"{service_name} servers busy. Wait for {wait_seconds} seconds")
            time.sleep(wait_seconds)


def wait_for_server_ip(service_name, #pylint: disable=dangerous-default-value
                       key_path,
                       remote_username,
                       environment,
                       max_concurrent,
                       chunk_id,
                       max_interval_between_describe_instances,
                       draining_tag,
                       mutex=threading.RLock(),
                       mutexes={},
                       last_checks={}):
    # We rate limit these to ensure fairness across jobs regardless of job size
    service_name = _duct_tape(service_name)
    with mutex:
        if service_name not in mutexes:
            mutexes[service_name] = threading.RLock()  # TraceLock not required because "gate" prints are sufficient.
            last_checks[service_name] = [None]
        lc = last_checks[service_name]
        mx = mutexes[service_name]
    log.write(f"Chunk {chunk_id} of {service_name} is at second gate")
    with mx:
        if lc[0] is not None:
            sleep_time = (60.0 / MAX_DISPATCHES_PER_MINUTE) - (
                time.time() - lc[0])
            # At least pi seconds of sleep to keep this chunk from racing with the previous one to the same server.
            sleep_time = max(sleep_time, 3.141519)
            if sleep_time > 0:
                log.write(f"Sleeping for {sleep_time:3.1f} seconds to rate-limit wait_for_server_ip.")
                time.sleep(sleep_time)
        # if we had to wait above, that counts toward the rate limit delay
        lc[0] = time.time()
        result = wait_for_server_ip_work(service_name, key_path,
                                         remote_username, environment,
                                         max_concurrent, chunk_id,
                                         max_interval_between_describe_instances,
                                         draining_tag)
        return result


def unixtime_now():
    return int(time.time())


def build_job_tag(job_tag_prefix, chunk_id):
    batch_job_id = os.environ.get('AWS_BATCH_JOB_ID', 'local')
    job_tag_key = f"{job_tag_prefix}{batch_job_id}_chunk{chunk_id}"
    job_tag_value_func = unixtime_now
    return job_tag_key, job_tag_value_func


def create_tag(instance_iD, tag_key, tag_value_func):
    command.execute(f"aws ec2 create-tags --resources {instance_iD} --tags Key={tag_key},Value={tag_value_func()}")


def delete_tag_with_retries(instance_iD, tag_key):
    command.execute_with_retries(f"aws ec2 delete-tags --resources {instance_iD} --tags Key={tag_key}")


@contextmanager
def ASGInstance(service, key_path, remote_username, environment, chunk_id, try_number, additional_attributes):
    max_concurrent = additional_attributes["max_concurrent"]
    max_interval_between_describe_instances = additional_attributes.get("max_interval_between_describe_instances", 900)
    job_tag_prefix = additional_attributes.get("job_tag_prefix", "RunningIDseqBatchJob_")
    job_tag_refresh_seconds = additional_attributes.get("job_tag_refresh_seconds", 600)
    draining_tag = additional_attributes.get("draining_tag", "draining")
    tracker = chunk_status_tracker(service)
    if tracker.service_unreliable():
        raise ServiceUnreliable(service, chunk_id)
    tracker.register_chunk_waiting(chunk_id, try_number)
    try:
        t_start = time.time()
        t_print = t_start
        # This delay helps prioritize chunks that still haven't had their first try.
        if try_number > 1:
            time.sleep(5.0)
        while not tracker.chunk_has_try_priority(chunk_id, try_number):
            if time.time() - t_print >= 60:
                t_print = time.time()
                log.write(f"try {try_number} for chunk {chunk_id} has been waiting for {(t_print - t_start):3.1f} seconds for other chunks' prior tries to complete")
            time.sleep(15.0)
        instance_ip, instance_iD = wait_for_server_ip(service, key_path, remote_username, environment, max_concurrent, chunk_id,
                                                      max_interval_between_describe_instances, draining_tag)
    finally:
        tracker.register_chunk_dispatched(chunk_id, try_number)
    log.write(f"starting alignment for chunk {chunk_id} on {service} server {instance_ip}")
    job_tag_key, job_tag_value_func = build_job_tag(job_tag_prefix, chunk_id)
    t = PeriodicThread(target=create_tag, wait_seconds=job_tag_refresh_seconds, stop_latency_seconds=60,
                       args=(instance_iD, job_tag_key, job_tag_value_func))
    t.start()
    try:
        yield instance_ip
    finally:
        t.stop()
        t.join()
        delete_tag_with_retries(instance_iD, job_tag_key)


def run_test():
    print("Hello. Testing idseq_dag.utils.server")
    chunk_status_tracker("rapsearch2").note_outcome("10.10.11.190", 101, 118.0, ChunkStatus.SUCCESS, 1)
    chunk_status_tracker("rapsearch2").note_outcome("10.10.11.190", 102, 119.0, ChunkStatus.SUCCESS, 1)
    chunk_status_tracker("rapsearch2").note_outcome("10.10.11.191", 103, 109.0, ChunkStatus.SUCCESS, 1)
    chunk_status_tracker("rapsearch2").note_outcome("10.10.11.191", 104, 9.0, ChunkStatus.CRASH, 1)
    chunk_status_tracker("rapsearch2").note_outcome("10.10.11.192", 104, 140.0, ChunkStatus.CORRUPT_OUTPUT, 2)
    chunk_status_tracker("rapsearch2").note_outcome("10.10.11.192", 104, 100.0, ChunkStatus.SUCCESS, 3)
    chunk_status_tracker("rapsearch2").note_outcome("10.10.11.192", 105, 98.0, ChunkStatus.SUCCESS, 1)
    chunk_status_tracker("rapsearch2").note_outcome("10.10.11.191", 106, 18.0, ChunkStatus.CRASH, 1)
    chunk_status_tracker("rapsearch2").note_outcome("10.10.11.191", 107, 180.0, ChunkStatus.TIMEOUT, 1)
    chunk_status_tracker("rapsearch2").register_chunk_waiting(100, 1)
    chunk_status_tracker("rapsearch2").register_chunk_waiting(100, 1)
    chunk_status_tracker("rapsearch2").register_chunk_waiting(99, 2)
    assert True == chunk_status_tracker("rapsearch2").chunk_has_try_priority(100, 1)
    assert False == chunk_status_tracker("rapsearch2").chunk_has_try_priority(99, 2)
    chunk_status_tracker("rapsearch2").register_chunk_dispatched(100, 1)
    assert True == chunk_status_tracker("rapsearch2").chunk_has_try_priority(99, 2)
    input_chunks = {
        'rapsearch2': set(range(110)),
        'gsnap': set(range(30)),
    }
    chunk_status_tracker("gsnap").note_outcome("10.10.12.121", 6, 18.0, ChunkStatus.CRASH, 1)
    chunk_status_tracker("gsnap").note_outcome("10.10.12.120", 7, 180.0, ChunkStatus.TIMEOUT, 1)
    chunk_status_tracker("gsnap").note_outcome("10.10.12.121", 7, 180.0, ChunkStatus.CRASH, 2)
    chunk_status_tracker("gsnap").note_outcome("10.10.12.120", 6, 18.0, ChunkStatus.SUCCESS, 2)
    for chunk_id in range(9, 18):
        chunk_status_tracker("gsnap").note_outcome("10.10.12.120", chunk_id, 88.8, ChunkStatus.CRASH, 1)
        assert chunk_status_tracker("gsnap").service_unreliable() == (chunk_id >= 16), chunk_id
    for service in ("rapsearch2", "gsnap"):
        tracker = chunk_status_tracker(service)
        tracker_report = tracker.status_report(len(input_chunks[service]))
        print(f"STATUS REPORT FOR {service}: " + json.dumps(tracker_report, indent=4))
        print(f"FULL OUTCOME DATA FOR {service}: " + json.dumps(tracker.outcome_data, indent=4))
        for chunk_id in [104, 105, 106]:
            print(f"Servers with rep for {chunk_id}: {tracker.servers_with_rep(chunk_id)}")
        assert tracker.service_unreliable() == (service == "gsnap")


if __name__ == "__main__":
    run_test()
