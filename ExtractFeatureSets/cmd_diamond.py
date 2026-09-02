#!/usr/bin/env python

import argparse
import os
import glob


"""
########################################
## Giant Virus-Host Predictor (GVHoP) ##
########################################
cmd_diamond.py
v1 2026

Author: Hsin-Ying Chang
Email: hyhazelchang@gmail.com

Usage:
# for execute diamond
ExtractFeatureSets/cmd_diamond.py --task=blastp --in_dir=example_MAGs/ --out_dir=temp_out/output_blastp/ --sh_dir=temp_out/sh/blastp/ --in_file_ext=faa --out_file_ext=txt --db_dir=GVHoP_database_v1.0/GVHoP_GVEUKs.dmnd --outfmt=6 --threshold=8 --opt="--evalue 1e-5" --n_job=8

# This script is used to generate shell scripts for executing diamond searches. It takes input files, constructs the appropriate diamond commands, and writes them to shell scripts for parallel execution.
"""


def main():
    parser = argparse.ArgumentParser(
        description=("Make the shscript for blast execution."),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--task",
                        default="blastp",
                        required=True,
                        type=str,
                        help="Specify your task.")
    parser.add_argument("--in_dir",
                        default="../example_MAGs/",
                        required=True,
                        type=str,
                        help="Directory containing input files (the file should have fasta extension).")
    parser.add_argument("--out_dir",
                        default="../output_blastp/",
                        type=str,
                        help="Output directory.")
    parser.add_argument("--sh_dir",
                        default="../sh/blastp/",
                        type=str,
                        help="Directory for shell scripts.")
    parser.add_argument("--in_file_ext",
                        default="fasta",
                        type=str)
    parser.add_argument("--out_file_ext",
                        default="xml",
                        type=str)
    parser.add_argument("--db_dir",
                        default="../GVHoP_database_v1.0/GVHoP_GVEUKs.dmnd",
                        required=True,
                        type=str)
    parser.add_argument("--outfmt",
                        default="6",
                        type=str)    
    parser.add_argument("--threshold",
                        default="8",
                        type=str)  
    parser.add_argument("--opt",
                        default="--evalue 1e-5",
                        type=str)    
    parser.add_argument("--n_job",
                        default=1,
                        type=int)

    args = parser.parse_args()
    task = args.task
    in_dir = args.in_dir
    out_dir = args.out_dir
    sh_dir = args.sh_dir
    in_file_ext = args.in_file_ext
    out_file_ext = args.out_file_ext
    db_dir = args.db_dir
    outfmt = args.outfmt
    threshold = args.threshold
    opt = args.opt
    n_job = args.n_job

    # Find the sequence files
    in_files = glob.glob(in_dir + "*." + in_file_ext)

    # Make output directory
    os.makedirs(out_dir, exist_ok=True)

    # Get file names
    count = 0
    diamond_cmd = []
    for seq in in_files:
        count += 1
        seq_name = seq.replace(in_dir, "")
        seq_name = seq_name.replace(f".{in_file_ext}", "")
        diamond_cmd.append(f"diamond {task} -d {db_dir} -q {seq} -o {out_dir}{seq_name}.{out_file_ext} -f {outfmt} -p {threshold} {opt}")

    # print out job scripts
    os.makedirs(sh_dir, exist_ok=True)
    quo = int(count / n_job)
    mod = int(count % n_job)
    cmd_num = 0
    for n in range(n_job):
        job = open(sh_dir + "job" + str(n+1) + ".sh" , "w")
        if n + 1 <= mod:
            for _ in range(quo + 1):
                job.write(diamond_cmd[cmd_num] + "\n")
                cmd_num += 1
            job.close()
        else:
            for _ in range(quo):
                job.write(diamond_cmd[cmd_num] + "\n")
                cmd_num += 1
            job.close()
    print(f"The shell scripts for diamond execution have been generated in {sh_dir}.")
    print(f"Please check the generated scripts and run the execute_sh.py script to execute diamond.")
    print(f"Example: python execute_sh.py --sh_dir={sh_dir} --log_dir=../execute/diamond/")

if __name__ == "__main__":
    main()