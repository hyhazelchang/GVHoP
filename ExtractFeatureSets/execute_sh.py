#!/usr/bin/python3

# execute_sh.py

# Hsin-Ying Chang <hyhazelchang@gmail.com>
# v1 2023/05/10

# python3 execute_sh.py --sh_dir=../sh/hmmsearch/ --log_dir=../execute/hmmsearch/


import argparse
import os
import glob

def main():
    parser = argparse.ArgumentParser(
        description=("Execute the shell scripts."),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--sh_dir",
                        default=None,
                        type=str,
                        help="Directory containing shell scripts. Please provide absolute path.")
    parser.add_argument("--log_dir",
                        default=None,
                        type=str,
                        help="Directory for log files. Please provide absolute path.")
    
    args = parser.parse_args()
    sh_dir = args.sh_dir
    log_dir = args.log_dir

    # Create the log directory
    os.system("mkdir -p " + log_dir)

    # Change to executable mode of shell scripts and execute
    sh_scripts = glob.glob(sh_dir + "*.sh")
    for sh in sh_scripts:
        os.system("chmod +x " + sh)
        job_name = os.path.basename(sh).split(".")[0]
        os.system(sh + " > " + log_dir + job_name + ".log 2>&1 &")
        print("command: " + sh + " > " + log_dir + job_name + ".log 2>&1 &")

if __name__ == "__main__":
    main()
