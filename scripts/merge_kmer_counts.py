import argparse
import os
import pandas as pd # type: ignore
from tqdm import tqdm # type: ignore

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
        description="Merges kmer counts of multiple shards for a single samples together."
    )
    parser.add_argument("--samplesheet", required=True, help="Path to the samplesheet CSV file")
    parser.add_argument("--jellyfish_output_list", required=True, help="Path to the list of Jellyfish output files")

    args = parser.parse_args()

    kmer_counts = []
    shards = {}

    ### Iterate over each shard, parse sample name and add kmer counts to final table, also keep track of which shards belong to which sample in a dictionary
    jellyfish_files = open(args.jellyfish_output_list).readlines()
    for line in tqdm(jellyfish_files, total=len(jellyfish_files), desc="Loading kmer counts"):
        
        ### parse sample and column name from the file path and store all shards per sample in a dictionary
        sample_name = os.path.basename(line.strip()).split("_")[0]
        name = os.path.basename(line.strip()).split("_")[0] + "_" + os.path.basename(line.strip()).split(".")[1]
        if sample_name in shards:
            shards[sample_name].append(name)
        else:
            shards[sample_name] = [name]

        ### append the kmer counts
        kmer_counts.append(pd.read_csv(line.strip(), sep=" ", header=None, names=["KMER", name]).sort_values(by="KMER").set_index("KMER"))
    
    ### concat shards into a single dataframe
    kmer_counts_df = pd.concat(kmer_counts, axis=1)

    ### sum up counts of all shards for each sample 
    for sample in shards:
        kmer_counts_df[sample] = kmer_counts_df[shards[sample]].sum(axis=1)

    ### subset to complete samples only
    kmer_df = kmer_counts_df[shards.keys()].copy()
    ### add kmer as columns
    kmer_df["KMER"] = kmer_df.index
    kmer_df.reset_index(drop=True, inplace=True)
    ### drop duplciates
    kmer_df = kmer_df.drop_duplicates(subset="KMER")
    
    ### build copy and convert kmer into their reverse complement form
    rc_kmer_df = kmer_df.copy()
    rc_kmer_df["KMER"] = rc_kmer_df["KMER"].apply(reverse_complement)
    
    ### concat, sort and set KMER as index
    sorted_complete_kmer_df = pd.concat([kmer_df, rc_kmer_df], axis=0).sort_values(by="KMER").set_index("KMER")

    ### save to file
    sorted_complete_kmer_df.to_csv("reads_kmer_counts.tsv", sep="\t", index=True)
        
if __name__ == "__main__":
    main()
