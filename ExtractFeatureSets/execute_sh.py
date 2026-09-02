#!/usr/bin/env python

import argparse
import os
import glob
import subprocess
from concurrent.futures import ThreadPoolExecutor


"""
########################################
## Giant Virus-Host Predictor (GVHoP) ##
########################################
execute_sh.py
v1 2026

Author: Hsin-Ying Chang
Email: hyhazelchang@gmail.com

Usage:
# for execute shell scripts
ExtractFeatureSets/execute_sh.py --sh_dir=temp_out/sh/blastp/ --log_dir=../logs/

# This script is used to execute shell scripts generateds. It takes a directory of shell scripts and executes them in parallel, logging the output to specified log files.
"""


def run_script(sh, log_dir):
    os.chmod(sh, 0o755)
    job_name = os.path.basename(sh).split(".")[0]
    log_file = f"{log_dir}{job_name}.log"
    with open(log_file, "w") as f:
        subprocess.run(
            ["bash", sh],
            stdout=f,
            stderr=subprocess.STDOUT,
            check=True
        )
    print(f"finished: {sh}")

def main():
    parser = argparse.ArgumentParser(
        description=("Execute the shell scripts."),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--sh_dir",
                        default=None,
                        required=True,
                        type=str,
                        help="Directory containing shell scripts.")
    parser.add_argument("--log_dir",
                        default="./logs/",
                        type=str,
                        help="Directory for log files.")
    
    args = parser.parse_args()
    sh_dir = args.sh_dir
    log_dir = args.log_dir

    # Create the log directory
    os.makedirs(log_dir, exist_ok=True)

    # Change to executable mode of shell scripts and execute
    sh_scripts = glob.glob(sh_dir + "*.sh")
    with ThreadPoolExecutor(max_workers=8) as executor:
        executor.map(lambda sh: run_script(sh, log_dir), sh_scripts)

if __name__ == "__main__":
    main()