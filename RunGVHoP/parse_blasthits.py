#!/usr/bin/python3

# parse_blasthits.py

# Hsin-Ying Chang <hyhazelchang@gmail.com>


import os
import argparse
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict


def main():
    parser = argparse.ArgumentParser(
            description=("Parse blast results to gene content matrix."),
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--blast_in",
                        type=str,
                        default="../output_blastp/",
                        help="The directory of blast results.")
    parser.add_argument("--in_file_ext",
                        default="tsv",
                        type=str)
    parser.add_argument("--column_names",
                        type=str,
                        default="../source_data/features/GVHoP_GVEUKs_all.tsv",
						required=True)
    parser.add_argument("--outfile",
                        type=str,
                        default="../source_data/example_inputs/ex_GVEUKs.tsv",
                        help="The directory of output file.")

    # Defining variables from input
    args = parser.parse_args()
    blast_in = args.blast_in
    in_file_ext = args.in_file_ext
    column_names = args.column_names
    outfile = args.outfile

    # Read data from input files
    blastfiles = [filename for filename in os.listdir(blast_in) if filename.endswith(in_file_ext)]

    # Check presence/absence of hits
    score_dict = defaultdict(dict)
    for file in blastfiles:
        filename = os.path.basename(file).replace('.'+in_file_ext, '')
        try:
            result = pd.read_csv(os.path.join(blast_in, file), sep='\t')
        except pd.errors.EmptyDataError:
            result = pd.DataFrame()
        if not result.empty:
            hits = result.iloc[:, 1].to_list()
            for i in range(len(hits)):
                if filename in score_dict[hits[i]].keys():
                    if result.iloc[i, 11] > score_dict[hits[i]][filename]:
                        score_dict[hits[i]][filename] = result.iloc[i, 11]
                else:
                    score_dict[hits[i]][filename] = result.iloc[i, 11]
    
    # Combine all rows into a DataFrame
    score_df = pd.DataFrame(score_dict)
    score_df = score_df.fillna(0).astype(int)
    score_df.columns = score_df.columns.astype(str)
    score_df.index.name = 'seq_name'
    print(score_df)
    wanted_columns = []
    for line in open(column_names, 'r'):
        wanted_columns += line.strip().split('\t')
    score_df = score_df.reindex(columns=wanted_columns, fill_value=0)

    # Write to output
    score_df.to_csv(outfile, sep='\t', index=True)

if __name__ == '__main__':
	main()