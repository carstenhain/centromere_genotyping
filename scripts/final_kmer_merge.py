import argparse
import pandas as pd # type: ignore

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--kmer_subtable_list', required=True)
    args = parser.parse_args()

    kmer_tables = []

    for line in open(args.kmer_subtable_list):
        tmp_df = pd.read_csv(line.strip(), sep="\t")
        for col in tmp_df.columns:
            if col == "KMER":
                continue
            tmp_df[col] = tmp_df[col].astype(float)
        kmer_tables.append(tmp_df.set_index("KMER").copy())

    pd.concat(kmer_tables, axis=1).to_csv("final_kmer_merged.tsv", sep="\t", index=True)

if __name__ == '__main__':
    main()
