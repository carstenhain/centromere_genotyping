import argparse
import os
import polars as pl # type: ignore
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
    kmer_index = None
    shards = {}

    ### Iterate over each shard, parse sample name and add kmer counts to final table, also keep track of which shards belong to which sample in a dictionary
    # read all shards
    with open(args.jellyfish_output_list, "r") as f:
        jellyfish_files = [x.strip() for x in f.readlines()]

    # query each shard
    for line in tqdm(jellyfish_files, total=len(jellyfish_files), desc="Loading kmer counts"):
        
        ### parse sample and column name from the file path and store all shards per sample in a dictionary
        sample_name = os.path.basename(line.strip()).split("_")[0]
        name = os.path.basename(line.strip()).split("_")[0] + "_" + os.path.basename(line.strip()).split(".")[1]
        if sample_name in shards:
            shards[sample_name].append(name)
        else:
            shards[sample_name] = [name]

        ### load and append the kmer counts
        shard_df = pl.read_csv(line.strip(), separator=" ", has_header=False, new_columns=["KMER", name]).sort("KMER")
        if kmer_index is None:
            kmer_index = shard_df.select("KMER")
        kmer_counts.append(shard_df.select(name))

    ### concat shards into a single dataframe
    kmer_counts_df = pl.concat([kmer_index, *kmer_counts], how="horizontal")

    ### sum up counts of all shards for each sample 
    for sample in shards:
        kmer_counts_df = kmer_counts_df.with_columns(
            pl.sum_horizontal([pl.col(column).fill_null(0) for column in shards[sample]]).alias(sample)
        )

    ### subset to complete samples only
    kmer_df = kmer_counts_df.select(["KMER", *list(shards.keys())])
    ### drop duplicates
    kmer_df = kmer_df.unique(subset=["KMER"], keep="first")
    
    ### build copy and convert kmer into their reverse complement form
    rc_kmer_df = kmer_df.with_columns(
        pl.col("KMER").map_elements(reverse_complement, return_dtype=pl.Utf8)
    )
    
    ### concat and sort
    sorted_complete_kmer_df = pl.concat([kmer_df, rc_kmer_df], how="vertical").sort("KMER")

    ### save to file
    sorted_complete_kmer_df.write_csv("reads_kmer_counts.tsv", separator="\t")
        
if __name__ == "__main__":
    main()
