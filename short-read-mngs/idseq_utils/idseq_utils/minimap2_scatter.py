import os
import shutil
import sys
import errno
from argparse import ArgumentParser
from os.path import abspath, basename, join
from subprocess import run, PIPE
from tempfile import TemporaryDirectory


class Minimap2MergeException(Exception):
    pass
################################################################################################################
#
# Minimap2
#
################################################################################################################


# TODO: add minimap2 make database here


def minimap2_alignment(cwd, par_tmpdir, cpus, database, out, queries):
    cmd = [
        "/usr/local/bin/minimap2",
        "-cx",
        "sr",
        f"-t {cpus}",
        "--split-map",
        out,
        f"{database}",
    ]
    for q in queries:
        print(q)
        cmd += [q]
    res = run(cmd, cwd=cwd, stdout=PIPE, stderr=PIPE)
    if res.returncode != 0:
        for line in res.stderr.decode().split("\n"):
            print(line)
        raise Exception(f"Command failed: {' '.join(cmd)}")


def minimap2_merge_cmd(cwd, par_tmpdir, chunks, queries):
    cmd = ["minimap2", "-cx", "sr", "--split-merge", "-o", f"{par_tmpdir}/out.paf"]
    for query in queries:
        cmd += [query]

    cmd += [
        ".",
    ]
    for chunk in chunks:
        cmd += [f"{chunk}"]

    res = run(cmd, cwd=cwd, stdout=PIPE, stderr=PIPE)
    if res.returncode != 0:
        for line in res.stderr.decode().split("\n"):
            print(line)
        raise Minimap2MergeException(f"Command {' '.join(cmd)} failed with result: {res.stderr.decode()}")


################################################################################################################
#
# Main
#
################################################################################################################


def zero_pad(n: int, i: int):
    tagged = str(n)
    return ("0" * (i - len(tagged))) + tagged


def make_par_dir(cwd: str, par_tmpdir: str):
    os.mkdir(join(cwd, par_tmpdir))
    p_dir = join(cwd, par_tmpdir, "parallelizer")
    os.mkdir(p_dir)
    with open(join(p_dir, "command"), "w"):
        pass
    with open(join(p_dir, "register"), "w"):
        pass
    with open(join(p_dir, "workers"), "w"):
        pass


def minimap2_chunk(db_chunk: str, output_dir: str, *query: str):
    """
    Run a single chunk of the database using minimap2-scatter
    This is no longer used, we should consider removing
    """

    # make output directory
    try:
        os.mkdir(output_dir)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    # get chunk

    # for diamond: chunk, n_seqs, n_letters = basename(db_chunk)[:-5].split("-")
    chunk = basename(db_chunk).split("_")[-1]  # example: nt.part_001

    # Don't really understand this temp dir thing?
    with TemporaryDirectory() as tmp_dir:
        make_par_dir(tmp_dir, "par-tmp")
        with open(join(tmp_dir, "par-tmp", f"align_todo_{zero_pad(0, 6)}"), "w") as f:
            f.writelines([f"Aligning chunk f{chunk}"])

        minimap2_alignment(
            cwd=tmp_dir,
            par_tmpdir="par-tmp",
            cpus=7,
            database=abspath(db_chunk),
            out=f"intermediate{chunk}",
            queries=[abspath(q) for q in query],
        )
        shutil.copy(
            join(tmp_dir, f"intermediate{chunk}"),
            join(output_dir, f"intermediate{chunk}"),
        )


def minimap2_merge(chunk_dir, out, *query):
    chunk_dir = abspath(chunk_dir)
    with TemporaryDirectory() as tmp_dir:
        make_par_dir(tmp_dir, "par-tmp")
        with open(join(tmp_dir, "par-tmp", f"join_todo_{zero_pad(0, 6)}"), "w") as f:
            f.write("TOKEN\n")
        chunks = []
        query_tmp = []
        for q in query:
            shutil.copy(q, join(tmp_dir, basename(q)))
            query_tmp.append(join(tmp_dir, basename(q)))
        for f in os.listdir(chunk_dir):
            shutil.copy(join(chunk_dir, f), join(tmp_dir, "par-tmp", f))
            chunks.append(join(tmp_dir, "par-tmp", f))
        minimap2_merge_cmd(tmp_dir, "par-tmp", chunks, query_tmp)
        shutil.copy(join(tmp_dir, "par-tmp", "out.paf"), out)


if __name__ == "__main__":
    parser = ArgumentParser()
    subparsers = parser.add_subparsers(title="commands", dest="command")

    minimap2_chunk_parser = subparsers.add_parser("minimap2-chunk")
    minimap2_chunk_parser.add_argument("--db", required=True)
    minimap2_chunk_parser.add_argument("--out-dir", required=True)
    minimap2_chunk_parser.add_argument("--query", required=True, action="append")

    minimap2_join_parser = subparsers.add_parser("minimap2-merge")
    minimap2_join_parser.add_argument("--chunk-dir", required=True)
    minimap2_join_parser.add_argument("--out", required=True)
    minimap2_join_parser.add_argument("--query", required=True, action="append")

    args = parser.parse_args(sys.argv[1:])
    if args.command == "minimap2-chunk":
        minimap2_chunk(args.db, args.out_dir, *args.query)
    elif args.command == "minimap2-merge":
        minimap2_merge(args.chunk_dir, args.out, *args.query)
