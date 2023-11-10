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
vep_input_dir = 'working_dir/vep_input/'
vep_output_dir = 'working_dir/vep_output/'

# Create directories
if not os.path.exists(vep_input_dir):
    os.makedirs(vep_input_dir)
if not os.path.exists(vep_output_dir):
    os.makedirs(vep_output_dir)

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
        elif line.startswith('LoF'):
            loftee = line.strip()
        elif line.startswith('SpliceAI'):
            spliceai = line.strip()


############### VEP run ###############

# Get the sorted list of input files from the input directory
file_list = sorted(os.listdir(vep_input_dir))

# Define function to run VEP
def run_vep(input_dir):

    # Loop through the input files in the directory
    for input_file in file_list:
        
        # Set input file path
        input_file_path = os.path.join(vep_input_dir, input_file)

        # Set the output file path
        output_file = os.path.join(vep_output_dir, input_file.replace('.tsv', '.vcf'))

        # Run VEP
        command = [
            os.path.join(vep_dir, 'vep'),
            '-i', input_file_path,
            '-o', output_file,
            '--dir_cache', vep_cache_dir,
            '--fasta', vep_fasta,
            '--offline',
            '--assembly', genome_ver,
            '--pick',
            '--cache',
            #'--dir_plugins', os.path.join(vep_cache_dir, 'Plugins/'), 
            '--plugin', loftee, 
            '--plugin', spliceai,
            '--vcf',
            '--fork', '4',
            '--hgvs',
            '--symbol',
            '--numbers',
            '--canonical',
            '--variant_class',
            '--no_stats'
        ]

        # Run VEP using subprocess
        result = subprocess.run(command, capture_output=True, text=True)

        # Check if VEP was executed correctly
        if result.returncode == 0:
            print(f'VEP executed successfully for {input_file_path}')
        else:
            print(f'Error executing VEP for {input_file_path}')
            print('Output:', result.stdout)
            print('Error:', result.stderr)

# Execute 
run_vep(vep_input_dir)
