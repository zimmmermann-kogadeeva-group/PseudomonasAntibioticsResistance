#!/usr/bin/env python3

import argparse
import shutil
import subprocess
from pathlib import Path

from Bio import SearchIO


def run_blast(seq_file, db_file, blast_output):

    output_dir = Path(blast_output).parent
    output_dir.mkdir(parents=True, exist_ok=True)
    new_db_file = output_dir / Path(db_file).name
    shutil.copy(db_file, new_db_file)

    # make blast database
    run1 = subprocess.run(
        f"makeblastdb -in {new_db_file} -parse_seqids -dbtype nucl",
        shell=True,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if run1.returncode != 0:
        raise OSError(f"Failed to run makeblastdb: {run1.stderr.decode()}")

    # align input sequences with query using blastn with default options
    run2 = subprocess.run(
        f"blastn -query {seq_file} -db {new_db_file} -out {blast_output} -outfmt 5",
        shell=True,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if run2.returncode != 0:
        raise OSError(f"Failed to run blastn: {run2.stderr.decode()}")

    # Return a dictionary correcting the hit IDs which get an unnecessary
    # prefix from BLAST
    return {x.id: x for x in SearchIO.parse(blast_output, "blast-xml")}


def extract_seqs(blast_results, as_str=False):
    res = {k: v.hsps[0].aln.alignment.query.seq for k, v in blast_results.items()}
    if as_str:
        res = "\n".join([f">{k}\n{v}" for k, v in res.items()])
    return res


if __name__ == "__main__":
    parser = argparse.ArgumentParser("get_seqs_for_st.py")
    parser.add_argument("filename")
    parser.add_argument("dbfile")
    parser.add_argument("workdir")
    parser.add_argument("-o", "--output")

    args = parser.parse_args()

    blast_res = run_blast(args.filename, args.dbfile, args.workdir + "/blast.xml")

    res = extract_seqs(blast_res, as_str=True)

    if args.output is not None:
        with open(args.output, "w") as fh:
            fh.write(res)
    else:
        print(res)
