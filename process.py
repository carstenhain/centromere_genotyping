#!/usr/bin/env python3

import argparse


def main():
    parser = argparse.ArgumentParser(
        description="Process a sample for centromere genotyping."
    )
    parser.add_argument("--name", required=True, help="Sample name")
    parser.add_argument("--file_path", required=True, help="Path to input file")
    args = parser.parse_args()

    print(f"Processing sample: {args.name}")
    print(f"Input file: {args.file_path}")


if __name__ == "__main__":
    main()
