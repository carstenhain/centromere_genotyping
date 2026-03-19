#!/usr/bin/env python3

import argparse
import os

import pandas as pd


KMERS_PATH = os.path.join(os.path.abspath(os.path.dirname(__file__)), "..", "data", "kmers.tsv")


def main():
    parser = argparse.ArgumentParser(
        description="Process a sample for centromere genotyping."
    )
    parser.add_argument("--name", required=True, help="Sample name")
    parser.add_argument("--file_path", required=True, help="Path to input file")
    args = parser.parse_args()

    print(f"Processing sample: {args.name}")
    print(f"Input file: {args.file_path}")

    try:
        kmers = pd.read_csv(KMERS_PATH, sep="\t")
    except FileNotFoundError:
        print(f"Error: kmers.tsv not found at {KMERS_PATH}")
        raise
    print(f"Number of rows in kmers.tsv: {len(kmers)}")


if __name__ == "__main__":
    main()
