import argparse
import os
import pandas as pd
from tqdm import tqdm

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
        kmer_counts.append(pd.read_csv(line.strip(), sep=" ", header=None, names=["KMER", name]).set_index("KMER"))
    
    ### concat shards into a single dataframe
    kmer_counts_df = pd.concat(kmer_counts, axis=1)

    ### sum up counts of all shards for each sample 
    for sample in shards:
        kmer_counts_df[sample] = kmer_counts_df[shards[sample]].sum(axis=1)

    ### write merged columns to file, include index as it is the kmer column
    kmer_counts_df[shards.keys()].to_csv("reads_kmer_counts.tsv", sep="\t", index=True)
        
if __name__ == "__main__":
    main()
