import pandas as pd
import os

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
        elif line.startswith('LoF'):
            loftee = line.strip()
        elif line.startswith('SpliceAI'):
            spliceai = line.strip()


############### VEP input preparation ###############

# Read input dataframe removing duplicated rows
mut = pd.read_csv(first_input, sep='\t').drop_duplicates()

# Convert start position so integer
mut['start'] = mut['start'].astype(int)

# Set end position based on VEP requirements
mut['end'] = mut['start']
mut['end'] = mut['end'].astype(int)

def delins_pos(row):
    if row['var'].startswith('-'):
        row['end'] = row['start'] + len(row['ref']) - 1
    elif row['ref'].startswith('-'):
        row['start'] = row['end'] + 1
    return row

def dnv_pos(row):
    if len(row['var']) == 2 and len(row['ref']) == 2:
        row['end'] = row['start'] + 1
    return row

mut = mut.apply(dnv_pos, axis=1)
mut = mut.apply(delins_pos, axis=1)

# VEP requirement for allele
mut['allele'] = mut['ref'] + '/' + mut['var']

# Set strand column
mut['strand'] = '+'

# Select columns
mut = mut[['chr', 'start', 'end', 'allele', 'strand', 'tumour_id']]

# Arrange the input file for a quicker VEP run
mut = mut.sort_values(['chr', 'start'])

# Group the rows based on tumour_id
tumours = mut.groupby('tumour_id')

# Loop through the groups and write each group to a separate tsv file
for group_name, group_vep in tumours:
        
    # Define the output file path for this group
    output_file_path = os.path.join('working_dir/vep_input', f'{group_name}.tsv')
        
    # Write the group dataframe to the output file as a tsv
    group_vep.to_csv(output_file_path, sep='\t', index=False, header =False)
        
print('VEP input .tsv files successfully created')