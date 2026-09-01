#!/usr/bin/python3

# cmd_hmm.py

# Hsin-Ying Chang <hyhazelchang@gmail.com>
# v1 2025/08/01

# Usage: python3 cmd_hmm.py --hmm_dir=../hmm/bin/ --in_dir=../test_MAGs/ --out_dir=../outdir/hmms/ --sh_dir=../sh/hmmbuild/ --in_file_ext=fasta --out_file_ext=hmm --opt=hmmbuild --n_job=20


import argparse
import os
import glob

def main():
    parser = argparse.ArgumentParser(
        description=("Make the shscript for hmm execution."),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--hmm_dir",
                        default=None,
                        type=str,
                        help="HMM file directory. Please provide absolute path.")
    parser.add_argument("--in_dir",
                        default=None,
                        type=str,
                        help="Directory containing input files (the file should have fasta extension). Please provide absolute path.")
    parser.add_argument("--db_dir",
                        default=None,
                        required=False,
                        type=str,
                        help="Directory containing database files (the file should have fasta extension). Please provide absolute path.")
    parser.add_argument("--out_dir",
                        default=None,
                        type=str,
                        help="Output directory. Please provide absolute path.")
    parser.add_argument("--sh_dir",
                        default=None,
                        type=str,
                        help="Directory for shell scripts. Please provide absolute path.")
    parser.add_argument("--in_file_ext",
                        default="fasta",
                        type=str)
    parser.add_argument("--out_file_ext",
                        default="hmm",
                        type=str)
    parser.add_argument("--opt",
                        default=None,
                        type=str)    
    parser.add_argument("--n_job",
                        default=1,
                        type=int)

    args = parser.parse_args()
    hmm_dir = args.hmm_dir
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
        if db_dir:
            hmm_cmd.append(hmm_dir + opt + " " + out_dir + file_name + "." + out_file_ext + " " + db_dir + " " + file)
        else:
            hmm_cmd.append(hmm_dir + opt + " " + out_dir + file_name + "." + out_file_ext + " " + file)

    # print out job scripts
    os.system("mkdir -p " + sh_dir)
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

if __name__ == "__main__":
    main()
