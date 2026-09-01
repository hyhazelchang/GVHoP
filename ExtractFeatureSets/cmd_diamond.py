#!/usr/bin/python3

# cmd_diamond.py

# Hsin-Ying Chang <hyhazelchang@gmail.com>
# v1 2023/06/10

# Usage: python3 /RAID1/working/R220/xinchang/pyscripts/cmd_diamond.py --diamond_dir=../softwares/diamond --task=blastp --in_dir=../test_MAGs/ --out_dir=../outdir/blastp_nr/ --sh_dir=../sh/blast/ --in_file_ext=faa --out_file_ext=txt --db_dir=../db/GVHoP_GVEUKs --outfmt=6 --threshold=1 --opt="-evalue 1e-5" --n_job=1


import argparse
import os
import glob

def main():
    parser = argparse.ArgumentParser(
        description=("Make the shscript for blast execution."),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--diamond_dir",
                        default=None,
                        type=str,
                        help="DIAMOND file directory. Please provide absolute path.")
    parser.add_argument("--task",
                        default=None,
                        type=str,
                        help="Specify your task.")
    parser.add_argument("--in_dir",
                        default=None,
                        type=str,
                        help="Directory containing input files (the file should have fasta extension). Please provide absolute path.")
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
                        default="xml",
                        type=str)
    parser.add_argument("--db_dir",
                        default=None,
                        type=str)
    parser.add_argument("--outfmt",
                        default=None,
                        type=str)    
    parser.add_argument("--threshold",
                        default=None,
                        type=str)  
    parser.add_argument("--opt",
                        default=None,
                        type=str)    
    parser.add_argument("--n_job",
                        default=1,
                        type=int)

    args = parser.parse_args()
    diamond_dir = args.diamond_dir
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
    seqfiles = glob.glob(in_dir + "*." + in_file_ext)

    # Make output directory
    os.makedirs(out_dir, exist_ok=True)

    # Get file names
    count = 0
    diamond_cmd = []
    for seq in seqfiles:
        count += 1
        seq_name = seq.replace(in_dir, "")
        seq_name = seq_name.replace(f".{in_file_ext}", "")
        diamond_cmd.append(diamond_dir + " " + task + 
                           " -d " + db_dir + 
                           " -q " + seq + 
                           " -o " + out_dir + seq_name + "." + out_file_ext + 
                           " -f " + outfmt +
                           " -p " + threshold +
                           " " + opt)

    # print out job scripts
    os.system("mkdir -p " + sh_dir)
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

if __name__ == "__main__":
    main()
