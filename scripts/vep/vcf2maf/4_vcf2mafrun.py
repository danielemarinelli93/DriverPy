import pandas as pd
import os
import subprocess
import requests
import shutil
import time

# Set config path
config_path = 'configs.txt'

# Set directories
working_dir = 'working_dir'
vep_output_dir = 'working_dir/vep_output/'
vcf2maf_output_dir = 'working_dir/vcf2maf_output/'

# Create directories
if not os.path.exists(vcf2maf_output_dir):
    os.makedirs(vcf2maf_output_dir)

# Read configurations
with open(config_path, 'r') as file:
    # Read the contents of the file
    contents = file.readlines()

    for line in contents:
        if line.startswith('input'):
            first_input = line.split('=')[1].strip()
        elif line.startswith('genome_ver'):
            genome_ver = line.split('=')[1].strip()
        elif line.startswith('vep_dir'):
            vep_dir = line.split('=')[1].strip()
        elif line.startswith('vep_cache_dir'):
            vep_cache_dir = line.split('=')[1].strip()
        elif line.startswith('vep_fasta'):
            vep_fasta = line.split('=')[1].strip()
        elif line.startswith('vcf2maf_dir'):
            vcf2maf_dir = line.split('=')[1].strip()
        elif line.startswith('retain_ann'):
            retain_ann = line.split('=')[1].strip()


############### VCF2MAF run ###############

# Get sorted file list
file_list = sorted(file for file in os.listdir(vep_output_dir) if not file.endswith('warnings.txt'))

def run_vcf2maf(input_dir):

    # Define the function to run VCF2MAF
    for input_file in file_list:
        
        # Set input file path
        input_file_path = os.path.join(vep_output_dir, input_file)

        # Set output file path
        output_file = os.path.join(vcf2maf_output_dir, input_file.replace('.vcf', '.maf'))

        # Run VCF2MAF
        command = [
            'perl',
            vcf2maf_dir,
            '--input-vcf', input_file_path,
            '--output-maf', output_file,
            '--ref-fasta', vep_fasta,
            '--tumor-id', input_file.replace('.vcf', ''),
            '--normal-id', '.',
            '--inhibit-vep',
            '--retain-ann', retain_ann
            ]

        # Run VCF2MAF using subprocess
        result = subprocess.run(command, capture_output=True, text=True)

        # Check if VCF2MAF was executed correctly
        if result.returncode == 0:
            print(f'vcf2maf executed successfully for {input_file_path}')
        else:
            print(f'vcf2maf executing VEP for {input_file_path}')
            print('Output:', result.stdout)
            print('Error:', result.stderr)

# Execute
run_vcf2maf(vep_output_dir)

# Define the function to append files for vcf2maf
def append_df(directory):
   
    # List all files in the directory
    file_list = os.listdir(directory)

    # Initialize an empty DataFrame to store the combined data
    merged = pd.DataFrame()

    for file in file_list:

        # Read file
        x = pd.read_csv(os.path.join(directory, file), sep='\t', skiprows=1)

        # Append df
        merged = merged.append(x, ignore_index=True)

    # Write merged output .tsv
    merged.to_csv(os.path.join(directory, 'merged.tsv'), sep='\t', index=False)

# Execute
append_df(vcf2maf_output_dir)