import os
import shutil
import sys
import errno

from argparse import ArgumentParser
from glob import glob
from multiprocessing import Pool
from os.path import abspath, basename, join
from subprocess import run, PIPE
from tempfile import NamedTemporaryFile, TemporaryDirectory
from typing import Iterable

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class DiamondBlastXException(Exception):
    pass


class DiamondJoinException(Exception):
    pass

################################################################################################################
#
# Diamond
#
################################################################################################################


def diamond_makedb(cwd: str, reference_fasta: str, database: str):
    cmd = [
        "diamond",
        "makedb",
        "--in",
        reference_fasta,
        "--db",
        database,
    ]
    run(cmd, check=True, cwd=cwd, stdout=PIPE, stderr=PIPE)


def diamond_blastx(
    cwd: str,
    par_tmpdir: str,
    block_size: float,
    database: str,
    out: str,
    queries: Iterable[str],
    chunk=False,
    join_chunks=0,
):
    cmd = [
        "diamond",
        "blastx",
        "--multiprocessing",
        "--parallel-tmpdir",
        par_tmpdir,
        "--block-size",
        str(block_size),
        "--db",
        database,
        "--out",
        out,
    ]
    for query in queries:
        cmd += ["--query", query]
    if chunk:
        cmd += ["--single-chunk"]
    if join_chunks > 0:
        cmd += ["--join-chunks", str(join_chunks)]
    res = run(cmd, cwd=cwd, stdout=PIPE, stderr=PIPE)
    if res.returncode != 0:
        for line in res.stderr.decode().split("\n"):
            print(line)
        raise DiamondBlastXException(f"Command {' '.join(cmd)} failed with error: {res.stderr.decode()}")


################################################################################################################
#
# Main
#
################################################################################################################


def _consume_iter(iterable: Iterable, n: int):
    for i, e in enumerate(iterable):
        yield e
        if i == n - 1:
            break


def align_chunk(ref_chunk: int, start: int, size: int, query_chunk: int):
    return f"{ref_chunk} {start} {size} # query_chunk={query_chunk}\n"


def zero_pad(n: int, m: int):
    tagged = str(n)
    return ("0" * (m - len(tagged))) + tagged


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


def make_db(reference_fasta: str, output_dir: str, chunks: int):
    os.mkdir(output_dir)
    chunk_size = (sum(1 for _ in SeqIO.parse(reference_fasta, "fasta")) // chunks) + 1
    seqs = SeqIO.parse(reference_fasta, "fasta")
    for i in range(chunks):
        print(f"STARTING  CHUNK {i}")
        fasta_name = f"{i}.fasta"
        SeqIO.write(_consume_iter(seqs, chunk_size), fasta_name, "fasta")
        n_seqs = n_letters = 0
        for seq in SeqIO.parse(fasta_name, "fasta"):
            n_seqs += 1
            n_letters += len(seq.seq)
        db_name = f"{i}-{n_seqs}-{n_letters}"
        print(f"INDEXING  CHUNK {i}")
        diamond_makedb(".", fasta_name, join(output_dir, db_name))
        os.remove(fasta_name)
        print(f"COMPLETED CHUNK {i}")


def blastx_chunk(db_chunk: str, output_dir: str, *query: str):
    try:
        os.mkdir(output_dir)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise
    chunk, n_seqs, n_letters = basename(db_chunk)[:-5].split("-")
    with TemporaryDirectory() as tmp_dir:
        make_par_dir(tmp_dir, "par-tmp")
        with open(join(tmp_dir, "par-tmp", f"align_todo_{zero_pad(0, 6)}"), "w") as f:
            f.writelines([align_chunk(int(chunk), 0, int(n_seqs), 0)])
        diamond_blastx(
            cwd=tmp_dir,
            par_tmpdir="par-tmp",
            block_size=int(n_letters) / 1e9,
            database=abspath(db_chunk),
            out="out.tsv",
            chunk=True,
            queries=(abspath(q) for q in query),
        )
        ref_block_name = f"ref_block_{zero_pad(0, 6)}_{zero_pad(int(chunk), 6)}"
        ref_dict_name = f"ref_dict_{zero_pad(0, 6)}_{zero_pad(int(chunk), 6)}"
        shutil.copy(
            join(tmp_dir, "par-tmp", ref_block_name), join(output_dir, ref_block_name)
        )
        shutil.copy(
            join(tmp_dir, "par-tmp", ref_dict_name), join(output_dir, ref_dict_name)
        )


def mock_reference_fasta(chunks: int, chunk_size: int):
    letters = chunk = i = 0
    while chunk < chunks:
        n = 100
        letters += n
        if letters > (1 + chunk) * chunk_size:
            chunk += 1
        yield SeqRecord(Seq("M" * n), "A")
        i += 1


def blastx_join(chunk_dir: str, out: str, *query: str):
    with TemporaryDirectory() as tmp_dir:
        make_par_dir(tmp_dir, "par-tmp")
        with open(join(tmp_dir, "par-tmp", f"join_todo_{zero_pad(0, 6)}"), "w") as f:
            f.write("TOKEN\n")

        for f in os.listdir(chunk_dir):
            shutil.copy(join(chunk_dir, f), join(tmp_dir, "par-tmp", f))

        chunks = len(os.listdir(chunk_dir)) // 2
        with NamedTemporaryFile() as ref_fasta, NamedTemporaryFile() as db:
            # make fake db to appease diamond
            SeqIO.write(SeqRecord(Seq("M"), "ID"), ref_fasta.name, "fasta")
            diamond_makedb(tmp_dir, ref_fasta.name, db.name)
            diamond_blastx(
                cwd=tmp_dir,
                par_tmpdir="par-tmp",
                block_size=1,
                database=db.name,
                out=out,
                join_chunks=chunks,
                queries=(abspath(q) for q in query),
            )

        with open(out, "w") as out_f:
            for out_chunk in glob(join(tmp_dir, f"{out}_*")):
                with open(out_chunk) as out_chunk_f:
                    out_f.writelines(out_chunk_f)
                os.remove(out_chunk)


if __name__ == "__main__":
    parser = ArgumentParser()
    subparsers = parser.add_subparsers(title="commands", dest="command")

    make_db_parser = subparsers.add_parser("make-db")
    make_db_parser.add_argument("--db", required=True)
    make_db_parser.add_argument("--in", required=True)
    make_db_parser.add_argument("--chunks", type=int, required=True)

    blastx_chunk_parser = subparsers.add_parser("blastx-chunk")
    blastx_chunk_parser.add_argument("--db", required=True)
    blastx_chunk_parser.add_argument("--out-dir", required=True)
    blastx_chunk_parser.add_argument("--query", required=True, action="append")

    blastx_chunks_parser = subparsers.add_parser("blastx-chunks")
    blastx_chunks_parser.add_argument("--db-dir", required=True)
    blastx_chunks_parser.add_argument("--out-dir", required=True)
    blastx_chunks_parser.add_argument("--query", required=True, action="append")

    blastx_join_parser = subparsers.add_parser("blastx-join")
    blastx_join_parser.add_argument("--chunk-dir", required=True)
    blastx_join_parser.add_argument("--out", required=True)
    blastx_join_parser.add_argument("--query", required=True, action="append")

    args = parser.parse_args(sys.argv[1:])
    if args.command == "make-db":
        make_db(args.__getattribute__("in"), args.db, args.chunks)
    elif args.command == "blastx-chunk":
        blastx_chunk(args.db, args.out_dir, *args.query)
    elif args.command == "blastx-chunks":

        def _blastx_chunk(db):
            print(f"STARTING:  {db}")
            res = blastx_chunk(join(args.db_dir, db), args.out_dir, *args.query)
            print(f"FINISHING: {db}")
            return res

        with Pool(48) as p:
            p.map(_blastx_chunk, os.listdir(args.db_dir))
    elif args.command == "blastx-join":
        blastx_join(args.chunk_dir, args.out, *args.query)
