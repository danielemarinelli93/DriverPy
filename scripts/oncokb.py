import pandas as pd
import os
import subprocess

# Set config path
config_path = 'configs.txt'

# Set directories
merged_results_dir = 'working_dir/merged_results'

# Read configurations
with open(config_path, 'r') as file:
    # Read the contents of the file
    contents = file.readlines()

    for line in contents:
        if line.startswith('oncokb_dir'):
            oncokb_dir = line.split('=')[1].strip()
        elif line.startswith('oncokb_token'):
            oncokb_token = line.split('=')[1].strip()
        elif line.startswith('oncotree_code'):
            oncotree_code = line.split('=')[1].strip()
        elif line.startswith('genome_ver'):
            genome_ver = line.split('=')[1].strip()


############### OncoKB annotator run ###############
# Write the OncoKB annotator command
cmd = f"python3 {os.path.join(oncokb_dir, 'MafAnnotator.py')} -i {os.path.join(merged_results_dir, 'merged.tsv')} -o {os.path.join(merged_results_dir, 'merged_oncokb.tsv')} -t {oncotree_code} -q HGVSp_Short -r {genome_ver} -b {oncokb_token}"

# Use subprocess to run the command
subprocess.run(cmd, shell=True)