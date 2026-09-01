#!/usr/bin/python3

# parse_hmmsearch.py

# Hsin-Ying Chang <hyhazelchang@gmail.com>
# v1 2025/08/01 some functions from ncldv_markersearch

# ref: Moniruzzaman et al. 2020 Nat. Commun.


import os, sys, subprocess, re, shlex, glob, operator, argparse
import pandas as pd
import numpy as np
from natsort import natsorted, ns
from collections import defaultdict
from operator import itemgetter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


""" Loop through and parse the checkm HMM output """
def hmm_parser(folder, suffix, output):
	combined_output = open(output, "w")
	combined_output.write("protein\tacc\thit\tstart\tend\taln_length\tscore\tcategory\n")
	hits = []
	bit_dict = {}
	for filenames in os.listdir(folder):
		if filenames.endswith(suffix):
			acc = re.sub('.'+suffix, "", filenames)
			f = open(folder+"/"+filenames, 'r')
			o = open(folder+"/"+filenames+".parsed", 'w')
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
					ids_hit = entry +"."+ hit_dict[entry]
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
				nr = acc +"_"+ cog
				if nr in done:
					combined_output.write(ids +"\t"+ acc +"\t"+ cog +"\t"+ start +"\t"+ end +"\t"+ aln_length +"\t"+ score +"\tNH\n")
					o.write(ids +"\t"+ acc +"\t"+ cog +"\t"+ start +"\t"+ end +"\t"+ score +"\tNH\n")
				else:
					combined_output.write(ids +"\t"+ acc +"\t"+ cog +"\t"+ start +"\t"+ end +"\t"+ aln_length +"\t"+ score +"\tBH\n")
					o.write(ids +"\t"+ acc +"\t"+ cog +"\t"+ start +"\t"+ end +"\t"+ score +"\tBH\n")
					done.append(nr)
			o.close()
			
def getprot(item):
	items = item.split("~")
	protein = items[0]
	cog = items[1]
	return(protein)

def parse_domout(path_to_parsed_hmmfile, acc, protein_dict, cog_name, protein2dups):
	parsed = open(path_to_parsed_hmmfile, "r")
	done = {}
	protein2coords = defaultdict(list)
	protein2align_length = {}

	main_hit = "NAN"
	rnap_hits = []
	main_hits = []
	protein2cog = defaultdict(lambda:"NA")
	protein2acc = {}
	protein2score = {}
	protein2category = {}
	protein2length = {}

	for n in parsed.readlines():
		line = n.rstrip()
		tabs = line.split("\t")
		protein = tabs[0]
		annot = tabs[2].rstrip(".trim")
		if annot == cog_name:
			rnap_hits.append(protein)
			id_hit = protein +"~"+ annot
			hmm_score = float(tabs[5])
			category = tabs[6]
			start = int(tabs[3])
			end =   int(tabs[4])
			align_length = abs(end - start)

			if hmm_score > 20 and align_length > 20:

				record = protein_dict[protein]
				prot_length = len(record.seq)

				nr = acc +"_"+ annot

				protein2cog[protein]    = annot
				protein2acc[protein]    = acc
				protein2score[protein]  = hmm_score
				protein2length[protein] = prot_length

				protein2coords[id_hit].append(start)
				protein2coords[id_hit].append(end)
				#if annot == "PolB":
				#	print(start, end, id_hit)

				protein2align_length[id_hit] = align_length

				if category == "BH" and annot == cog_name:
				#if annot == cog_name:
					#main_hit = protein
					protein2dups[id_hit] = "single_besthit"
					main_hits.append(protein)

				else:
					#main_hit = protein
					protein2dups[id_hit] = "secondary_hit"
					main_hits.append(protein)

				protein2category[protein] = category

	parsed.close()
	#main_hits = set(main_hits)
	#main_hits2 = []
	#main_hits2 = [main_hits2.append(x) for x in main_hits if x not in main_hits2]
	return main_hits, protein2cog, protein2acc, protein2score, protein2length, protein2category, protein2coords, protein2align_length


def get_proteinsonreplicon(proteinid, record_list, prox):
	contig_name = re.sub(r"_\d*$", "", proteinid)
	final_list = []
	index = []
	indexzero=0
	ind = int(0)
	#record_list = natsorted(seqdict.keys())
	#print record_list
	for record in record_list:
		#print record
		if contig_name in record:
			final_list.append(record)
			index.append(ind)
			if proteinid == record:
				indexzero = ind
			ind +=1

	#prox = int(5) # number of genes to look in front and in back of gene
	start = indexzero - prox
	end = indexzero + prox + 1

	if start < 0:
		start = 0
	if end > len(final_list):
		end = len(final_list)

	protein_list = final_list[start:end]

	#print proteinid, indexzero, protein_list	
	return(protein_list)

# main function that runs the program
def run_program(in_dir, out_dir, out_preffix, prox, in_file_ext, allhits, cogset, seq_dir, seq_file_ext):
	merged = open(out_dir + out_preffix + "_full_output.txt", "w")
	merged.write("New_protein_name\tgenome\thit\tprotein_length\tbit_score\tnum_proteins_merged\thit_type\tprotein_ids\thmm_aln_coords\n")	
	raw_output = out_dir + out_preffix + "_rawout.txt"
	with open(cogset, "r") as f:
		cog_set = f.readline().rstrip('\n').split('\t')
	if allhits:
		hitset = ['single_besthit', 'main_hit', 'secondary_hit']
	else:
		hitset = ['single_besthit', 'main_hit']
	hmm_parser(in_dir, in_file_ext, raw_output)
	print("Compiling results...")
	protein_tally = []
	df = pd.DataFrame()	
	for i in os.listdir(seq_dir):
		if i.endswith(seq_file_ext):
			markercount = defaultdict(float)
			domout_file = os.path.join(in_dir, i)			
			parsed = re.sub(seq_file_ext, "domout.parsed", domout_file)
			acc = re.sub('.'+seq_file_ext, "", i)		
			merged_protein_list = []
			protein_file = os.path.join(seq_dir, i)
			seq_handle = open(protein_file, "r")
			seq_dict = SeqIO.to_dict(SeqIO.parse(seq_handle, "fasta"))
			record_list = natsorted(seq_dict.keys())
			orf_set = [record.id for record in seq_dict.values()]
			prot2protlist = defaultdict(list)
			num_proteins = defaultdict(lambda:int(1))
			prot2loc = defaultdict(list)
			prot2locrange = defaultdict(list)
			# parse domout file and get protein hits and coordinates
			for cog in cog_set:
				protein2dups = defaultdict(lambda:"hits")
				main_hits, protein2cog, protein2acc, protein2score, protein2length, _, protein2coords, _ = parse_domout(parsed, acc, seq_dict, cog, protein2dups)
				for main_hit in main_hits:	
					if main_hit in protein_tally:
						pass
					else:
						protein_tally.append(main_hit)
						orf_set = get_proteinsonreplicon(main_hit, record_list, prox)
						if main_hit == "NAN":
							pass
						else:
							#print(main_hit)
							prot2protlist[main_hit].append(main_hit)
							#print(len(prot2protlist[main_hit]))
							id_hit1 = main_hit +"~"+ cog
							range1 = protein2coords[id_hit1]
							r1 = range(range1[0], range1[1])
							meanloc1 = np.mean(range1)
							locrange1 = str(range1[0]) +"-"+ str(range1[1])
							prot2locrange[main_hit].append(locrange1)
							prot2loc[main_hit].append(meanloc1)
							orf_set.remove(main_hit)
							for m in protein_tally:
								if m in orf_set:
									orf_set.remove(m)
							orf_set = set(orf_set)
							for d in orf_set:
								if protein2cog[d] == cog:
									id_hit2 = d +"~"+ cog
									range2 = protein2coords[id_hit2]
									r2 = range(range2[0], range2[1])
									meanloc2 = np.mean(range2)
									locrange2 = str(range2[0]) +"-"+ str(range2[1])									
									set1 = set(r1)
									inter = set1.intersection(r2)
									if int(len(inter)) > 10:
										protein2dups[id_hit2] = "secondary_hit"
										#print(id_hit1, id_hit2, r1, r2, range1, range2, len(inter))
									else:
										protein_tally.append(d)
										protein2dups[id_hit1] = "main_hit"
										protein2dups[id_hit2] = "secondary_hit"
										#print(main_hit, id_hit1, id_hit2, d, protein2dups[id_hit1], protein2dups[id_hit2])
										protein_tally.append(id_hit1)
										protein_tally.append(id_hit2)
										minrange = min(range1 + range2)
										maxrange = max(range1 + range2)
										protein2coords[id_hit1] = [minrange, maxrange]
										#print(minrange, maxrange)
										prot2locrange[main_hit].append(locrange2)
										#protein2align_length[id_hit1] = abs(maxrange - minrange)
										#protein2align_length[main_hit] = int(protein2align_length[main_hit]) + int(protein2align_length[d])
										protein2length[main_hit] = int(protein2length[main_hit]) + int(protein2length[d])
										protein2score[main_hit] = float(protein2score[main_hit]) + float(protein2score[d])
										prot2protlist[main_hit].append(d)
										prot2loc[main_hit].append(meanloc2)
										num_proteins[id_hit1] +=1
										#print(id_hit1, id_hit2, prot2protlist, prot2protlist[main_hit])
										merged_protein_list.append(id_hit2)
				all_besthits = [p for p in protein2dups.keys() if protein2dups[p] in ["single_besthit", "main_hit"]] 
				all_scores = [protein2score[getprot(p)] for p in all_besthits]
				if len(all_scores) > 0:
					max_index, _ = max(enumerate(all_scores), key=operator.itemgetter(1))
					best_hit = all_besthits[max_index]
					other_hits = [j for j in all_besthits if j != best_hit]
					for o in other_hits:
						protein2dups[o] = "secondary_hit"
				hit_tally = defaultdict(int)		
				#best_hit_bit = defaultdict(float)
				for item in protein2dups:
					if item in merged_protein_list:
						pass
					elif protein2dups[item] in hitset:
						items = item.split("~")
						protein = items[0]
						hit = items[1]
						protlist = prot2protlist[protein]
						loc_list = [float(loc) for loc in prot2loc[protein]]
						index_list = [i[0] for i in sorted(enumerate(loc_list), key=lambda x:x[1])]
						#sorted_loc_list = [i[1] for i in sorted(enumerate(loc_list), key=lambda x:x[1])]
						sorted_prot_list = [protlist[index] for index in index_list]
						prot_str = ";".join(sorted_prot_list)
						range_list = prot2locrange[protein]
						sorted_range_list = [range_list[index] for index in index_list]
						range_str = ";".join(sorted_range_list) 
						#loc_str = ";".join([str(n) for n in sorted_loc_list])
						acc = protein2acc[protein]
						final_name_str = re.sub("_", ".", acc) +"_"+ hit
						hit_tally[final_name_str] +=1
						ptally = str(hit_tally[final_name_str])
						final_name = final_name_str +".copy"+ptally
						merged.write(final_name +"\t"+ acc +"\t"+ hit +"\t"+ str(protein2length[protein]) +"\t"+ str(round(protein2score[protein], 1)) +"\t"+ str(num_proteins[item]) +"\t"+ protein2dups[item] +"\t"+ prot_str +"\t"+ range_str +"\n")
						markercount[hit] +=1
			s1 = pd.DataFrame(pd.Series(markercount, name = acc))
			df = pd.concat([df, s1], axis=1, sort=True)
	merged.close()
	df2 = df.transpose().fillna(0).astype("Int64")
	df2.fillna(0, inplace=True)
	df2 = df2.reindex(columns=cog_set, fill_value=0)
	df2.to_csv(f"{out_dir}{out_preffix}.tsv", sep="\t", index_label="genome")


def main():
	parser = argparse.ArgumentParser(
        description=("Parse hmmsearch data"),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
	)
	parser.add_argument('--in_dir', required=True, help='Input folder of domout file')
	parser.add_argument('--in_file_ext', required=True, help='.fna, .fa, or .fasta')
	parser.add_argument('--seq_dir', required=True, help='Input folder of seq file')
	parser.add_argument('--seq_file_ext', required=True, help='.fna, .fa, or .fasta')
	parser.add_argument('--out_dir', required=True, help='out_dir for output files')
	parser.add_argument('--out_preffix', required=True, help='out_dir name prefix for output files')
	parser.add_argument('--proximity', required=False, type=int, default=5, help='number of genes to look up- and downstream of hits to join genes (default=5)')
	parser.add_argument('--cogset', required=True, default=None)
	parser.add_argument('--allhits', type=bool, default=False, const=True, nargs='?', help='Provide all hits (default is to provide only best hits to each marker gene)')
	
	args = parser.parse_args()
	in_dir = args.in_dir
	in_file_ext = args.in_file_ext
	seq_dir = args.seq_dir
	seq_file_ext = args.seq_file_ext
	out_dir = args.out_dir
	out_preffix = args.out_preffix
	prox = args.proximity
	allhits = args.allhits
	cogset = args.cogset
	
    # make output directory
	os.makedirs(out_dir, exist_ok=True)
	
    # run program
	run_program(in_dir, out_dir, out_preffix, prox, in_file_ext, allhits, cogset, seq_dir, seq_file_ext)

if __name__ == "__main__":
    main()