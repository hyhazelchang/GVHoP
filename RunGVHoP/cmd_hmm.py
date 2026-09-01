#!/usr/bin/python3

# cmd_hmm.py

# Hsin-Ying Chang <hyhazelchang@gmail.com>

import argparse
import os
import glob

def main():
    parser = argparse.ArgumentParser(
        description=("Make the shscript for hmm execution."),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--in_dir",
                        default="../example_MAGs/",
                        required=True,
                        type=str,
                        help="Directory containing input files (the file should have fasta extension).")
    parser.add_argument("--db_dir",
                        default="../GVHoP_database_v1.0/GVHoP_GVOGs.hmm",
                        required=True,
                        type=str,
                        help="Directory containing database files.")
    parser.add_argument("--out_dir",
                        default="../output_hmm/",
                        type=str,
                        help="Output directory.")
    parser.add_argument("--sh_dir",
                        default="../sh/hmmsearch/",
                        type=str,
                        help="Directory for shell scripts.")
    parser.add_argument("--in_file_ext",
                        required=True,
                        default="faa",
                        type=str)
    parser.add_argument("--out_file_ext",
                        default="domout",
                        type=str)
    parser.add_argument("--opt",
                        default="hmmsearch --cpu 10 -E 1e-5 --domtblout",
                        required=True,
                        type=str)
    parser.add_argument("--n_job",
                        default=2,
                        type=int)

    args = parser.parse_args()
    in_dir = args.in_dir
    db_dir = args.db_dir
    out_dir = args.out_dir
    sh_dir = args.sh_dir
    in_file_ext = args.in_file_ext
    out_file_ext = args.out_file_ext
    opt = args.opt
    n_job = args.n_job

    # Find the input files
    in_files = glob.glob(in_dir + "*." + in_file_ext)

    # Make output directory
    os.makedirs(out_dir, exist_ok=True)

    # Get file names
    count = 0
    hmm_cmd = []
    for file in in_files:
        count += 1
        file_name = file.replace(in_dir, "")
        file_name = file_name.replace(f".{in_file_ext}", "")
        hmm_cmd.append(opt + " " + out_dir + file_name + "." + out_file_ext + " " + db_dir + " " + file)

    # print out job scripts
    os.makedirs(sh_dir, exist_ok=True)
    quo = int(count / n_job)
    mod = int(count % n_job)
    cmd_num = 0
    for n in range(n_job):
        job = open(sh_dir + "job" + str(n+1) + ".sh" , "w")
        if n + 1 <= mod:
            for _ in range(quo + 1):
                job.write(hmm_cmd[cmd_num] + "\n")
                cmd_num += 1
            job.close()
        else:
            for _ in range(quo):
                job.write(hmm_cmd[cmd_num] + "\n")
                cmd_num += 1
            job.close()
    print(f"The shell scripts for hmm execution have been generated in {sh_dir}.")
    print(f"Please check the generated scripts and run the execute_sh.py script to execute hmmsearch.")
    print(f"Example: python execute_sh.py --sh_dir={sh_dir} --log_dir=../execute/hmmsearch/")

if __name__ == "__main__":
    main()