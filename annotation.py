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
vcf2maf_output_dir = 'working_dir/vcf2maf_output/'
cgi_input_dir = 'working_dir/cgi_input'
cgi_output_dir = 'working_dir/cgi_output'
cravat_output_dir = 'working_dir/cravat_output'
merged_results_dir = 'working_dir/merged_results'

# Create directories
if not os.path.exists(working_dir):
    os.makedirs(working_dir)
if not os.path.exists(vep_input_dir):
    os.makedirs(vep_input_dir)
if not os.path.exists(vep_output_dir):
    os.makedirs(vep_output_dir)
if not os.path.exists(vcf2maf_output_dir):
    os.makedirs(vcf2maf_output_dir)
if not os.path.exists(cgi_input_dir):
    os.makedirs(cgi_input_dir)
if not os.path.exists(cgi_output_dir):
    os.makedirs(cgi_output_dir)
if not os.path.exists(cravat_output_dir):
    os.makedirs(cravat_output_dir)
if not os.path.exists(merged_results_dir):
    os.makedirs(merged_results_dir)

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
        elif line.startswith('CADD'):
            cadd = line.strip()
        elif line.startswith('SpliceAI'):
            spliceai = line.strip()
        elif line.startswith('vcf2maf_dir'):
            vcf2maf_dir = line.split('=')[1].strip()
        elif line.startswith('retain_ann'):
            retain_ann = line.split('=')[1].strip()
        elif line.startswith('oncokb_dir'):
            oncokb_dir = line.split('=')[1].strip()
        elif line.startswith('oncokb_token'):
            oncokb_token = line.split('=')[1].strip()
        elif line.startswith('oncotree_code'):
            oncotree_code = line.split('=')[1].strip()
        elif line.startswith('cravat_dir'):
            cravat_dir = line.split('=')[1].strip()
        elif line.startswith('cravat_anno'):
            cravat_anno = line.split('=')[1].strip()
        elif line.startswith('cgi_id'):
            cgi_id = line.split('=')[1].strip()
        elif line.startswith('cgi_title'):
            cgi_title = line.split('=')[1].strip()
        elif line.startswith('cgi_reference'):
            cgi_reference = line.split('=')[1].strip()
        elif line.startswith('cravat_dir'):
            cravat_dir = line.split('=')[1].strip()
            


############### Cancer genome interpreter ###############

# Read input file removing duplicated rows
mut = pd.read_csv(first_input, sep='\t').drop_duplicates()

# Convert start position so integer
mut['start'] = mut['start'].astype(int)

# Set column names
mut.columns = ['chr', 'pos', 'ref', 'alt', 'sample']

# Write .tsv input file
mut.to_csv(os.path.join(cgi_input_dir, 'cgi_input.tsv'), sep = '\t', index=False)

# Send input file to API
headers = {'Authorization': cgi_id}
payload = {'cancer_type': oncotree_code, 'title': first_input, 'reference': cgi_reference}
r = requests.post('https://www.cancergenomeinterpreter.org/api/v1',
                headers=headers,
                files={
                        'mutations': open('working_dir/cgi_input/cgi_input.tsv', 'rb')
                        },
                data=payload)
r.json()

# Get job_id
cgi_req = requests.get('https://www.cancergenomeinterpreter.org/api/v1', headers=headers)
cgi_ids = cgi_req.json()
cgi_job_id = cgi_ids[-1]



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
            '--pick',
            '--cache',
            '--offline',
            '--assembly', 'GRCh37',
            '--dir_plugins', os.path.join(vep_cache_dir, 'Plugins/'), 
            '--plugin', loftee, 
            '--plugin', spliceai,
            '--dir_cache', vep_cache_dir,
            '--fasta', os.path.join(vep_cache_dir, 'ref_fasta/Homo_sapiens.GRCh37.dna.toplevel.fa'),
            '--vcf',
            '--fork', '4',
            '--hgvs',
            '--symbol',
            '--numbers',
            '--canonical',
            '--variant_class',
            '--no_stats',
            '--force_overwrite'
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
            os.path.join(vcf2maf_dir, 'vcf2maf.pl'),
            '--input-vcf', input_file_path,
            '--output-maf', output_file,
            '--ref-fasta', os.path.join(vep_cache_dir, 'ref_fasta/Homo_sapiens.GRCh37.dna.toplevel.fa'),
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



############### openCRAVAT for CHASMplus calls ###############

# Loop through all files in the input directory that end with '.vcf'
for filename in os.listdir(vep_output_dir):
    if filename.endswith('.vcf'):
        
        # Set the input and output file paths
        input_file = os.path.join(vep_output_dir, filename)

        # Construct the command to run
        command = [
            cravat_dir,
            'run',
            input_file,
            '-d', cravat_output_dir,
            '-t', 'tsv',
            '-l', 'hg19',
            '-a', 'chasmplus_LUAD', 'chasmplus_LUSC', 'chasmplus', 'hg19'
        ]

        # Execute the command using subprocess
        subprocess.run(command)

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



############### Cancer genome interpreter output download ###############

# Set payload
payload = {'action': 'download'}

# Initialize a variable to track whether the job is completed
job_completed = False

while not job_completed:
    # Get CGI run status
    response = requests.get(f'https://www.cancergenomeinterpreter.org/api/v1/{cgi_job_id}', headers=headers)
    
    if response.status_code == 200:
        # If the job status is 200, it's completed
        job_completed = True
        print("CGI run completed successfully")
        
        # Download results
        r = requests.get(f'https://www.cancergenomeinterpreter.org/api/v1/{cgi_job_id}', headers=headers, params=payload)
        with open('cgi_res.zip', 'wb') as fd:
            fd.write(r._content)
        print('CGI results downloaded successfully')
        
        # Unzip output
        command = ['unzip', 'cgi_res.zip']
        subprocess.run(command, capture_output=True, text=True)

        file_list = ('alterations.tsv', 'biomarkers.tsv', 'summary.txt', 'input01.tsv', 'cgi_res.zip')

        # Move results to cgi output directory
        for file in file_list:
            shutil.move(file, 'working_dir/cgi_output/')
    else:
        # If the CGI run is not completed, wait for 15 minutes before checking again
        print(f"Job status: {response.status_code}. Retrying in 15 minutes...")
        time.sleep(900)


# Read openCRAVAT merged output
cravat = pd.read_csv(os.path.join(cravat_output_dir, 'merged.tsv'), sep='\t')

# Remove the 'chr' prefix
cravat['original_input.chrom'] = cravat['original_input.chrom'].str.replace('chr', '')

# Create join column
cravat['join'] = cravat[['original_input.chrom', 'original_input.pos', 'original_input.ref_base', 'original_input.alt_base', 'tags']].apply(lambda row: ' '.join(str(x) for x in row), axis=1)

# Select columns
filtered_cravat = cravat[['join', 'chasmplus.pval', 'chasmplus.score', 'chasmplus.transcript', 'chasmplus.all', 'chasmplus_LUAD.pval', 'chasmplus_LUAD.score', 'chasmplus_LUAD.transcript', 'chasmplus_LUAD.all', 'chasmplus_LUSC.pval', 'chasmplus_LUSC.score', 'chasmplus_LUSC.transcript', 'chasmplus_LUSC.all']]

# Write .tsv file
filtered_cravat.to_csv(os.path.join(merged_results_dir, 'filtered_cravat.tsv'), sep='\t', index=False)

# Read CGI alterations.tsv output
cgi = pd.read_csv(os.path.join(cgi_output_dir, 'alterations.tsv'), sep='\t')

# Create the join column
cgi['join'] = cgi[['CHROMOSOME', 'POSITION', 'REF', 'ALT', 'CGI-Sample ID']].apply(lambda row: ' '.join(str(x) for x in row), axis=1)

# Remove duplicated rows (i.e. rows with multiple calls)
filtered_cgi_size = cgi[cgi['join'].map(cgi.groupby('join').size()) == 1]

# Select columns
filtered_cgi = filtered_cgi_size[['join', 'CGI-Oncogenic Summary', 'CGI-Oncogenic Prediction']]

# Write .tsv file
filtered_cgi.to_csv(os.path.join(merged_results_dir, 'filtered_cgi.tsv'), sep='\t', index=False)

# Read VCF2maf merged output
vcf2maf = pd.read_csv(os.path.join(vcf2maf_output_dir, 'merged.tsv'), sep='\t')

# Select columns
filtered_vcf2maf = vcf2maf[['Hugo_Symbol', 'NCBI_Build', 'Chromosome', 'Start_Position', 'End_Position', 'Strand', 'Variant_Classification', 'Variant_Type', 'Reference_Allele', 'Tumor_Seq_Allele2', 'Tumor_Sample_Barcode', 'HGVSc', 'HGVSp', 'HGVSp_Short', 'Transcript_ID', 'EXON', 'INTRON', 'all_effects', 'Consequence', 'BIOTYPE', 'CANONICAL', 'IMPACT', 'LoF', 'LoF_filter', 'LoF_flags', 'SpliceAI_pred_SYMBOL', 'SpliceAI_pred_DS_AG', 'SpliceAI_pred_DS_AL', 'SpliceAI_pred_DS_DG', 'SpliceAI_pred_DS_DL']]

# Create the join column
filtered_vcf2maf['join'] = filtered_vcf2maf[['Chromosome', 'Start_Position', 'Reference_Allele', 'Tumor_Seq_Allele2', 'Tumor_Sample_Barcode']].apply(lambda row: ' '.join(str(x) for x in row), axis=1)

# Write .tsv file
filtered_vcf2maf.to_csv(os.path.join(merged_results_dir, 'filtered_vcf2maf.tsv'), sep='\t')

# Merge dataframes
merged_tmp = pd.merge(filtered_vcf2maf, filtered_cgi, on='join', how='left')
merged_final = pd.merge(merged_tmp, filtered_cravat, on='join', how='left')

# Write merged output
merged_final.to_csv(os.path.join(merged_results_dir, 'merged.tsv'), sep='\t')



############### OncoKB annotator run ###############
# Write the OncoKB annotator command
command = [
    "python",
    os.path.join(oncokb_dir, 'MafAnnotator.py'),
    '-i',
    os.path.join(merged_results_dir, 'merged.tsv'),
    '-o',
    os.path.join(merged_results_dir, 'merged_oncokb.tsv'),
    '-t',
    'oncotree_code',
    '-q',
    'HGVSp_Short',
    '-r',
    genome_ver,
    '-b',
    oncokb_token
]

# Use subprocess to run the command
process = subprocess.run(command, capture_output=True, text=True)



############### Compress output files ###############
def gzip_func(input_dir):
    
    # Get sorted file list
    file_list = sorted(os.listdir(input_dir))

    # Loop through the input files in the directory
    for input_file in file_list:
    
        # Construct the input file path
        input_file_path = os.path.join(input_dir, input_file)

        # Define the command to run
        command = ['gzip', '-9', input_file_path]

        # Run the command using subprocess
        subprocess.run(command, capture_output=True, text=True)

# Compress VEP output
gzip_func(vep_output_dir)
gzip_func(vcf2maf_output_dir)
gzip_func(cravat_output_dir)
gzip_func(cgi_input_dir)
gzip_func(cgi_output_dir)
gzip_func(merged_results_dir)
