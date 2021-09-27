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

################################################################################################################
#
# Minimap2
#
################################################################################################################


# TODO: add minimap2 make database here

def minimap2_alignment(cwd, par_tmpdir, cpus, database, out, queries):
    cmd =  [
        "/usr/local/bin/minimap2", 
        "-cx", "sr",
        f"-t {cpus}",
        "--split-map", out,
        f"{database}",
    ]
    print(queries)
    for q in queries:
        print(q)
        cmd += [q]
    res = run(cmd, cwd=cwd, stdout=PIPE, stderr=PIPE)
    if res.returncode != 0:
        for line in res.stderr.decode().split("\n"):
            print(line)
        raise Exception(f"Command failed: {' '.join(cmd)}")

def minimap2_merge_cmd(cwd, par_tmpdir, chunks, queries):
    cmd = [
        "minimap2", 
        "-cx", 
        "sr", 
        "--split-merge",
        "-o",
        f"{par_tmpdir}/out.paf"
    ]
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
        raise Exception(f"Command failed: {' '.join(cmd)}") 

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


def zero_pad(n: int, l: int):
    tagged = str(n)
    return ("0" * (l - len(tagged))) + tagged


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


def minimap2_chunk(db_chunk: str, output_dir: str, *query: str):
    """
    Run a single chunk of the database using minimap2-scatter 
    """

    # make output directory
    try:
        os.mkdir(output_dir)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    # get chunk

    # for diamond: chunk, n_seqs, n_letters = basename(db_chunk)[:-5].split("-")
    chunk = basename(db_chunk).split("_")[-1] # example: nt.part_001  


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
            queries=[abspath(q) for q in query]
        )
        shutil.copy(join(tmp_dir, f"intermediate{chunk}"), join(output_dir, f"intermediate{chunk}"))


def mock_reference_fasta(chunks: int, chunk_size: int):
    letters = chunk = i = 0
    while chunk < chunks:
        n = 100
        letters += n
        if letters > (1 + chunk) * chunk_size:
            chunk += 1
        yield SeqRecord(Seq("M" * n), "A")
        i += 1

def minimap2_merge(chunk_dir, out, *query):
    chunk_dir = abspath(chunk_dir)
    with TemporaryDirectory() as tmp_dir:
        make_par_dir(tmp_dir, "par-tmp")
        with open(join(tmp_dir, "par-tmp", f"join_todo_{zero_pad(0, 6)}"), "w") as f:
            f.write("TOKEN\n") 
        chunks = []
        for q in query:
            shutil.copy(q, join(tmp_dir, q))
        for f in os.listdir(chunk_dir):
            print("list", f)
            shutil.copy(join(chunk_dir, f), join(tmp_dir, "par-tmp", f)) 
            chunks.append(join(tmp_dir, "par-tmp", f))
        minimap2_merge_cmd(tmp_dir, "par-tmp", chunks, query)
        shutil.copy(join(tmp_dir, "par-tmp", "out.paf"), out)


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
                out="out.tsv",
                join_chunks=chunks,
                queries=(abspath(q) for q in query),
            )

        with open(out, "w") as out_f:
            for out_chunk in glob(join(tmp_dir, f"out.tsv_*")):
                with open(out_chunk) as out_chunk_f:
                    out_f.writelines(out_chunk_f)
                os.remove(out_chunk)


if __name__ == "__main__":
    parser = ArgumentParser()
    subparsers = parser.add_subparsers(title="commands", dest="command")

    minimap2_chunk_parser = subparsers.add_parser("minimap2-chunk")
    minimap2_chunk_parser.add_argument("--db", required=True)
    minimap2_chunk_parser.add_argument("--out-dir", required=True)
    minimap2_chunk_parser.add_argument("--query", required=True, action='append')

    minimap2_join_parser = subparsers.add_parser("minimap2-merge")
    minimap2_join_parser.add_argument("--chunk-dir", required=True)
    minimap2_join_parser.add_argument("--out", required=True)
    minimap2_join_parser.add_argument("--query", required=True, action='append')

    args = parser.parse_args(sys.argv[1:])
    if args.command == "minimap2-chunk":
        minimap2_chunk(args.db, args.out_dir, *args.query)
    elif args.command == "minimap2-merge":
        minimap2_merge(args.chunk_dir, args.out, *args.query)

