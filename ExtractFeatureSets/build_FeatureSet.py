#!/usr/bin/env python

import os, re, operator, argparse
import pandas as pd
import numpy as np
from collections import defaultdict


"""
########################################
## Giant Virus-Host Predictor (GVHoP) ##
########################################
build_FeatureSet.py
v1 2026

Author: Hsin-Ying Chang
Email: hyhazelchang@gmail.com

Usage:
# for hmmer results
ExtractFeatureSets/build_FeatureSet.py --in_dir=temp_out/output_hmm/ --in_file_ext=domout --outfile=source_data/example_inputs/ex_GVOGs.tsv --column_names=source_data/features/GVHoP_GVOGs_all.tsv --opt=hmmerhits
# for blast results
ExtractFeatureSets/build_FeatureSet.py --in_dir=temp_out/output_blastp/ --in_file_ext=tsv --outfile=source_data/example_inputs/ex_GVEUKs.tsv --column_names=source_data/features/GVHoP_GVEUKs_all.tsv --opt=blasthits

# This script is used to build feature sets from raw results of BLAST or HMMER searches. It parses the input files, checks for the presence/absence of hits, and generates a feature set in a tab-separated format.
"""


def blastparser(in_dir, in_file_ext):
	# Read data from input files
    blastfiles = [filename for filename in os.listdir(in_dir) if filename.endswith(in_file_ext)]
    # Check presence/absence of hits
    score_dict = defaultdict(dict)
    for file in blastfiles:
        filename = os.path.basename(file).replace('.'+in_file_ext, '')
        try:
            result = pd.read_csv(os.path.join(in_dir, file), sep='\t')
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
    return score_dict

def hmmparser(in_dir, in_file_ext):
	hits = []
	bit_dict = {}
	content_dict = defaultdict(dict)
	for filenames in os.listdir(in_dir):
		if filenames.endswith(in_file_ext):
			acc = re.sub('.' + in_file_ext, "", filenames)
			f = open(in_dir + "/" + filenames, 'r')
			o = open(in_dir + "/" + filenames + ".parsed", 'w')
			o.write("protein_id\taccession\tbest_hit\taln_start\taln_end\tscore\ttype\n")
			hit_dict = {}
			start_dict = {}
			end_dict = {}
			bit_dict = defaultdict(int)
			position_dict = defaultdict(list)
			for line in f.readlines():
				if line.startswith("#"):
					pass
				else:
					newline = re.sub(r'\s+', '\t', line)
					list1 = newline.split('\t')
					ids = list1[0]
					hit = re.sub("_", "", list1[3])
					#print(hit)
					score = float(list1[7])
					domain_evalue = float(list1[11])
					if score > bit_dict[ids] and domain_evalue < 1e-3:
						ids_hit = ids +"."+ hit
						start = int(list1[15])
						end   = int(list1[16])
						position_dict[ids_hit].append(start)
						position_dict[ids_hit].append(end)
						#print(ids_hit, score, domain_evalue, start, end)
						hit_dict[ids] = hit
						start_dict[ids] = start
						end_dict[ids] = end
						bit_dict[ids] = score
			bit_sorted = sorted(bit_dict.items(), key=operator.itemgetter(1), reverse=True)
			output_list = []
			for item in bit_sorted:
				entry = item[0]
				score = item[1]
				if score > 0:
					#print entry, item, filenames
					ids_hit = entry + "." + hit_dict[entry]
					output_list.append(entry +"\t"+ str(hit_dict[entry]) +"\t"+ str(min(position_dict[ids_hit])) +"\t"+ str(max(position_dict[ids_hit])) +"\t"+ str(bit_dict[entry]) )
			done = []
			for line in output_list:
				line1 = line.rstrip()
				tabs = line1.split("\t")
				ids = tabs[0]
				hits.append(ids)
				cog = tabs[1]
				start = tabs[2]
				end = tabs[3]
				aln_length = str(abs(float(end) - float(start)))
				score = tabs[4]
				nr = acc + "_" + cog
				if nr in done:
					content_dict[acc][cog] += 1
					o.write(ids +"\t"+ acc +"\t"+ cog +"\t"+ start +"\t"+ end +"\t"+ aln_length +"\t"+ score +"\tNH\n")
				else:
					content_dict[acc][cog] = 1
					o.write(ids +"\t"+ acc +"\t"+ cog +"\t"+ start +"\t"+ end +"\t"+ aln_length +"\t"+ score +"\tBH\n")
					done.append(nr)
			o.close()
	return content_dict

def main():
    parser = argparse.ArgumentParser(
            description=("Build feature sets from raw results."),
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--in_dir",
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
    parser.add_argument("--opt",
                        type=str,
                        default="blasthits",
                        help="Option for parsing method. blasthits or hmmerhits.")

    # Defining variables from input
    args = parser.parse_args()
    in_dir = args.in_dir
    in_file_ext = args.in_file_ext
    column_names = args.column_names
    outfile = args.outfile
    opt = args.opt

	# Parse blast or hmm results
    if opt == "blasthits":
        out_dict = blastparser(in_dir, in_file_ext)
    elif opt == "hmmerhits":
        out_dict = hmmparser(in_dir, in_file_ext)

    # Combine all rows into a DataFrame
    out_df = pd.DataFrame(out_dict)
    out_df = out_df.fillna(0).astype(int)
    out_df.columns = out_df.columns.astype(str)
    out_df.index.name = 'seq_name'
    print(out_df)
    with open(column_names, 'r') as f:
        wanted_columns = f.readline().strip().split('\t')
    out_df = out_df.reindex(columns=wanted_columns, fill_value=0)

    # Write to output
    out_df.to_csv(outfile, sep='\t', index=True)

if __name__ == '__main__':
	main()