import pandas as pd
import os
import subprocess

# Set config path
config_path = 'configs.txt'

# Set directories
vep_output_dir = 'working_dir/vep_output/'
cravat_output_dir = 'working_dir/cravat_output'

# Create directories
if not os.path.exists(cravat_output_dir):
    os.makedirs(cravat_output_dir)

# Read configurations
with open(config_path, 'r') as file:
    # Read the contents of the file
    contents = file.readlines()

    for line in contents:
        if line.startswith('cravat_dir'):
            cravat_dir = line.split('=')[1].strip()
        elif line.startswith('cravat_genome'):
            cravat_genome = line.split('=')[1].strip()
        elif line.startswith('cravat_anno'):
            cravat_anno = line.split('=')[1].strip()


############### openCRAVAT ###############

# Loop through all files in the input directory that end with '.vcf'
for filename in os.listdir(vep_output_dir):
    if filename.endswith('.vcf'):
        
        # Set the input and output file paths
        input_file = os.path.join(vep_output_dir, filename)

        # Run openCRAVAT
        cmd = f'{cravat_dir} run {input_file} -d {cravat_output_dir} -t tsv -l {cravat_genome} -a {cravat_anno}'
        subprocess.run(cmd, shell=True)

# Define the function to append files for openCRAVAT
def append_df(directory):
    
    # List all files in the directory
    file_list = [file for file in os.listdir(directory) if file.endswith('variant.tsv')]

    # Initialize an empty DataFrame to store the combined data
    merged = pd.DataFrame()

    for file in file_list:
        # Read file
        x = pd.read_csv(os.path.join(directory, file), sep='\t', comment='#')

        # Append df
        merged = merged.append(x, ignore_index=True)

    merged.to_csv(os.path.join(directory, 'merged.tsv'), sep='\t', index=False)

append_df(cravat_output_dir)
