import argparse
import os
import pandas as pd # type: ignore
import re
from tqdm import tqdm # type: ignore

from collections import defaultdict


def reverse_complement(dna_sequence):
    """
        Computes the reverse complement of a base sequence given in upper case

        Parameters:
        dna_sequence:                     DNA sequence in upper case (only ATCG)

        Returns:
        reverse complement of dna_sequence
    """
    complement = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C'
    }

    complement_sequence = ''.join(complement[base] for base in dna_sequence)
    reverse_complement_sequence = complement_sequence[::-1]
    return reverse_complement_sequence


def main():
    parser = argparse.ArgumentParser(
        description="Convert a de Bruijn graph (DBG) file, subsetted to the k-mers of interest, to a k-mer count table."
    )
    parser.add_argument("--name", required=True, help="Sample name")
    parser.add_argument("--zgrep_file", required=True, help="Path to zgrep preprocessed DBG input file")
    parser.add_argument("--kmer_list", required=True, help="Path to file containing list of k-mers")
    parser.add_argument("--output", required=True, help="Path to output TSV file")
    args = parser.parse_args()

    print(f"Processing sample: {args.name}")
    print(f"Input DBG file: {args.zgrep_file}")
    print(f"Output TSV file: {args.output}")

    ### reads contigs from the zgrep and extract sequence and km value from the header, store in dataframe
    unitig_data = {"UNITIG": [], "KM": []}
    unitig_lines = open(args.zgrep_file).readlines()
    for j in tqdm(range(0, len(unitig_lines)), total=len(unitig_lines), desc=f"Reads contigs for {args.name}"):
        if unitig_lines[j].startswith(">"):
            ## get unitig sequence
            unitig_data["UNITIG"].append(unitig_lines[j + 1].strip())
            ## get km from dbg header
            unitig_data["KM"].append(float(re.search(r'km:f:(\d+\.?\d*)', unitig_lines[j].strip()).group(1))) # type: ignore
    unitig_df = pd.DataFrame(unitig_data)

    ### extract kmers from the unitigs, also add reverse complement kmers to the table
    k = 61
    kmer_km = defaultdict(float)
    for _, row in tqdm(unitig_df.iterrows(), total=unitig_df.shape[0], desc=f"Extracting kmers for {args.name}"):
        unitig = row['UNITIG']
        for i in range(len(unitig) - k + 1):
            kmer = unitig[i:i + k]
            r_kmer = reverse_complement(unitig[i:i + k])

            kmer_km[kmer] += row['KM']
            kmer_km[r_kmer] += row['KM']

    ### transform into dataframe, drop kmer duplicates
    result_data = []
    for kmer in open(args.kmer_list):
        kmer = kmer.strip()
        if kmer in kmer_km:
            result_data.append({"KMER": kmer, args.name: kmer_km[kmer]})
        else:
            result_data.append({"KMER": kmer, args.name: 0})    
    result_df = pd.DataFrame(result_data).drop_duplicates(subset="KMER")
    
    ### build copy and convert kmer into their reverse complement form
    rc_result_df = result_df.copy()
    rc_result_df["KMER"] = rc_result_df["KMER"].apply(reverse_complement)
    
    ### concat, sort and set KMER as index
    sorted_complete_result_df = pd.concat([result_df, rc_result_df], axis=0).drop_duplicates(subset="KMER").sort_values(by="KMER")

    ### save to file
    sorted_complete_result_df.to_csv(args.output, sep="\t", index=False)

if __name__ == "__main__":
    main()
