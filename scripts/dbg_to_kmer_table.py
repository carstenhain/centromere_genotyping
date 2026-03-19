#!/usr/bin/env python3

import argparse
import os

import pandas as pd


KMERS_PATH = os.path.join(os.path.abspath(os.path.dirname(__file__)), "..", "data", "kmers.tsv")


def main():
    parser = argparse.ArgumentParser(
        description="Convert a de Bruijn graph (DBG) file to a k-mer count table."
    )
    parser.add_argument("--name", required=True, help="Sample name")
    parser.add_argument("--file_path", required=True, help="Path to DBG input file")
    parser.add_argument("--output", required=True, help="Path to output TSV file")
    args = parser.parse_args()

    print(f"Processing sample: {args.name}")
    print(f"Input DBG file: {args.file_path}")

    try:
        kmers = pd.read_csv(KMERS_PATH, sep="\t")
    except FileNotFoundError:
        print(f"Error: kmers.tsv not found at {KMERS_PATH}")
        raise

    dbg = pd.read_csv(args.file_path, sep="\t")

    if "kmer" not in dbg.columns:
        raise ValueError(
            f"Input DBG file '{args.file_path}' does not contain a 'kmer' column. "
            "Expected columns: kmer (and optionally additional count columns)."
        )
    kmer_table = kmers.merge(dbg, how="left", on="kmer")
    kmer_table.insert(0, "sample", args.name)

    kmer_table.to_csv(args.output, sep="\t", index=False)
    print(f"K-mer table written to: {args.output}")


if __name__ == "__main__":
    main()
