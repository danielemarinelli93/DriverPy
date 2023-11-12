import pandas as pd
import os, subprocess, requests, shutil, time, argparse, logging, sys

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

# Read configurations
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
        elif line.startswith('LoF_hg37'):
            loftee_hg37 = line.strip()
        elif line.startswith('LoF_hg38'):
            loftee_hg38 = line.strip()
        elif line.startswith('SpliceAI_hg37'):
            spliceai_hg37 = line.strip()
        elif line.startswith('SpliceAI_hg38'):
            spliceai_hg38 = line.strip()
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
        elif line.startswith('cgi_reference'):
            cgi_reference = line.split('=')[1].strip()
        elif line.startswith('cravat_dir'):
            cravat_dir = line.split('=')[1].strip()

############### Cancer genome interpreter ###############

### CGI run
def cgi_run():

    # Create directories
    if not os.path.exists(cgi_input_dir):
        os.makedirs(cgi_input_dir)
    if not os.path.exists(cgi_output_dir):
        os.makedirs(cgi_output_dir)

    # Check if the CGI run was already started
    if os.path.exists(os.path.join(cgi_output_dir, 'cgi_job_id.txt')):

        # Read the CGI job id from the .txt file    
        with open(os.path.join(cgi_output_dir, 'cgi_job_id.txt'), 'r') as cgi_input:
            cgi_job_id = cgi_input.read().strip()
        
        # Get the list of existing job ids from the CGI API
        headers = {'Authorization': cgi_id}
        r = requests.get('https://www.cancergenomeinterpreter.org/api/v1', headers=headers)
        
        # Run CGI only if the same job id is not already running
        if cgi_job_id in r.text:
            logging.error(
                '\n'
                'CGI job not submitted: a job with the same id was already submitted'
                )
        if not cgi_job_id in r.text:
            logging.error(
                '\n'
                'CGI job not submitted: this previously submitted job is not anymore on the system\n'
                'please delete working_dir/cgi_output/cgi_job_id.txt and re-run CGI'                
            )

    if not os.path.exists(os.path.join(cgi_output_dir, 'cgi_job_id.txt')):

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
        payload = {'cancer_type': oncotree_code, 'title': first_input, 'reference': 'hg19' if genome_ver == 'GRCh37' else 'hg38' if genome_ver == 'GRCh38' else None}
        r = requests.post('https://www.cancergenomeinterpreter.org/api/v1',
                        headers=headers,
                        files={'mutations': open('working_dir/cgi_input/cgi_input.tsv', 'rb')},
                        data=payload)
        r.json()

        # Get job_id
        cgi_req = requests.get('https://www.cancergenomeinterpreter.org/api/v1', headers=headers)
        cgi_ids = cgi_req.json()
        cgi_job_id = cgi_ids[-1]

        # Write CGI job id to a text file
        with open(os.path.join(cgi_output_dir, 'cgi_job_id.txt'), 'w') as file:

            # Write the string to the file
            file.write(cgi_job_id)
        logging.info(
            '\n'
            'CGI job submitted'
            )

### CGI output download
def cgi_download():

    if not os.path.exists(os.path.join(cgi_output_dir, 'alterations.tsv')):

        # Check if the CGI run was already started
        if os.path.exists(os.path.join(cgi_output_dir, 'cgi_job_id.txt')):

            # Read the CGI job id from the .txt file    
            with open(os.path.join(cgi_output_dir, 'cgi_job_id.txt'), 'r') as cgi_input:
                cgi_job_id = cgi_input.read().strip()

            # Set payload
            payload = {'action': 'download'}

            # Initialize a variable to track whether the job is completed
            job_completed = False

            while not job_completed:
                # Get CGI run status
                headers = {'Authorization': cgi_id}
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
        
        else:
            print('CGI job id not available')


############### ENSEMBL-VEP & VCF2MAF & OncoKB ###############

### VEP run
def vep_run():
     
     # Create directories
    if not os.path.exists(vep_input_dir):
        os.makedirs(vep_input_dir)

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
            
    print('VEP input files successfully created')

    # Create directories
    if not os.path.exists(vep_output_dir):
        os.makedirs(vep_output_dir)

    # Loop through the input files in the directory
    file_list = sorted(os.listdir(vep_input_dir))
    
    for input_file in file_list:
        
        # Set input file path
        input_file_path = os.path.join(vep_input_dir, input_file)

        # Set the output file path
        output_file = os.path.join(vep_output_dir, input_file.replace('.tsv', '.vcf'))

        loftee = loftee_hg37 if genome_ver == 'GRCh37' else loftee_hg38 if genome_ver == 'GRCh38' else None
        spliceai = spliceai_hg37 if genome_ver == 'GRCh37' else spliceai_hg38 if genome_ver == 'GRCh38' else None

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

    # Create directories
    if not os.path.exists(vcf2maf_output_dir):
        os.makedirs(vcf2maf_output_dir)

    # Define the function to run VCF2MAF
    file_list = sorted(file for file in os.listdir(vep_output_dir) if not file.endswith('warnings.txt'))
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
    
    if not os.path.exists(os.path.join(vcf2maf_output_dir, 'merged.tsv')):
                          
        # List all files in the directory
        file_list = os.listdir(vcf2maf_output_dir)

        if len(file_list) == 1:
            
            for file in file_list:
                
                x = pd.read_csv(os.path.join(vcf2maf_output_dir, file), sep='\t', comment='#')
                
                x.to_csv(os.path.join(vcf2maf_output_dir, 'merged.tsv'), sep='\t', index=False)

        elif len(file_list) > 1:

            # Initialize an empty DataFrame to store the combined data
            merged = pd.DataFrame()

            for file in file_list:

                # Read file
                x = pd.read_csv(os.path.join(vcf2maf_output_dir, file), sep='\t', comment='#')

                # Append df
                merged = merged.append(x, ignore_index=True)

            # Write merged output .tsv
            merged.to_csv(os.path.join(vcf2maf_output_dir, 'merged.tsv'), sep='\t', index=False)
         
        elif len(file_list) == 0:
            logging.error(
                '\n'
                'No input files for VCF2MAF\n'
                'please check VEP run for errors'
            )
            sys.exit()

    # Write the OncoKB annotator command
    cmd = f"python3 {os.path.join(oncokb_dir, 'MafAnnotator.py')} -i {os.path.join(vcf2maf_output_dir, 'merged.tsv')} -o {os.path.join(vcf2maf_output_dir, 'merged_oncokb.tsv')} -t {oncotree_code} -q HGVSp_Short -r {genome_ver} -b {oncokb_token}"

    if not os.path.exists(os.path.join(vcf2maf_output_dir, 'merged_oncokb.tsv')):
        # Use subprocess to run the command
        subprocess.run(cmd, shell=True)


############### openCRAVAT ###############

### openCRAVAT run
def cravat_run():

    # Create directories
    if not os.path.exists(cravat_output_dir):
        os.makedirs(cravat_output_dir)

    # Loop through all files in the input directory that end with '.vcf'
    for filename in os.listdir(vep_output_dir):
        if filename.endswith('.vcf'):
            
            # Set the input and output file paths
            input_file = os.path.join(vep_output_dir, filename)

            # Run openCRAVAT
            cmd = f"{cravat_dir} run {input_file} -d {cravat_output_dir} -t tsv -l {'hg19' if genome_ver == 'GRCh37' else 'hg38' if genome_ver == 'GRCh38' else None} -a {cravat_anno} {'hg19' if genome_ver == 'GRCh37' else None}"
            subprocess.run(cmd, shell=True)

    # List all files in the directory
    file_list = [file for file in os.listdir(cravat_output_dir) if file.endswith('variant.tsv')]

    # Initialize an empty DataFrame to store the combined data
    merged = pd.DataFrame()

    for file in file_list:
        # Read file
        x = pd.read_csv(os.path.join(cravat_output_dir, file), sep='\t', comment='#')

        # Append df
        merged = merged.append(x, ignore_index=True)

    merged.to_csv(os.path.join(cravat_output_dir, 'merged.tsv'), sep='\t', index=False)


############### Merging results ###############
def merging_results():

    if os.path.exists(os.path.join(cravat_output_dir, 'merged.tsv')):
    
        # Read openCRAVAT merged output
        cravat = pd.read_csv(os.path.join(cravat_output_dir, 'merged.tsv'), sep='\t')

        # Remove the 'chr' prefix
        cravat['original_input.chrom'] = cravat['original_input.chrom'].str.replace('chr', '')

        if genome_ver == 'GRCh37': 
            
            # Create join column
            cravat['join'] = cravat[['original_input.chrom', 'original_input.pos', 'original_input.ref_base', 'original_input.alt_base', 'tags']].apply(lambda row: ' '.join(str(x) for x in row), axis=1)

            # Select columns
            filtered_cravat = cravat[['join', 'chasmplus.pval', 'chasmplus.score', 'chasmplus.transcript', 'chasmplus.all', 'chasmplus_LUAD.pval', 'chasmplus_LUAD.score', 'chasmplus_LUAD.transcript', 'chasmplus_LUAD.all', 'chasmplus_LUSC.pval', 'chasmplus_LUSC.score', 'chasmplus_LUSC.transcript', 'chasmplus_LUSC.all']]

            # Write .tsv file
            filtered_cravat.to_csv(os.path.join(merged_results_dir, 'filtered_cravat.tsv'), sep='\t', index=False)

        elif genome_ver == 'GRCh38':
            
            # Create join column
            cravat['join'] = cravat[['chrom', 'pos', 'ref_base', 'alt_base', 'tags']].apply(lambda row: ' '.join(str(x) for x in row), axis=1)

            # Select columns
            filtered_cravat = cravat[['join', 'chasmplus.pval', 'chasmplus.score', 'chasmplus.transcript', 'chasmplus.all', 'chasmplus_LUAD.pval', 'chasmplus_LUAD.score', 'chasmplus_LUAD.transcript', 'chasmplus_LUAD.all', 'chasmplus_LUSC.pval', 'chasmplus_LUSC.score', 'chasmplus_LUSC.transcript', 'chasmplus_LUSC.all']]

            # Write .tsv file
            filtered_cravat.to_csv(os.path.join(merged_results_dir, 'filtered_cravat.tsv'), sep='\t', index=False)

    if os.path.exists(os.path.join(cgi_output_dir, 'alterations.tsv')):

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

    if os.path.exists(os.path.join(vcf2maf_output_dir, 'merged_oncokb.tsv')):

        # Read VCF2maf merged output
        vcf2maf = pd.read_csv(os.path.join(vcf2maf_output_dir, 'merged_oncokb.tsv'), sep='\t')

        # Select columns
        filtered_vcf2maf = vcf2maf[['Hugo_Symbol', 'NCBI_Build', 'Chromosome', 'Start_Position', 'End_Position', 'Strand', 'Variant_Classification', 'Variant_Type', 'Reference_Allele', 'Tumor_Seq_Allele2', 'Tumor_Sample_Barcode', 'HGVSc', 'HGVSp', 'HGVSp_Short', 'Transcript_ID', 'EXON', 'INTRON', 'all_effects', 'Consequence', 'BIOTYPE', 'CANONICAL', 'IMPACT', 'LoF', 'LoF_filter', 'LoF_flags', 'SpliceAI_pred_SYMBOL', 'SpliceAI_pred_DS_AG', 'SpliceAI_pred_DS_AL', 'SpliceAI_pred_DS_DG', 'SpliceAI_pred_DS_DL', 'ANNOTATED', 'GENE_IN_ONCOKB', 'VARIANT_IN_ONCOKB',	'MUTATION_EFFECT', 'MUTATION_EFFECT_CITATIONS', 'ONCOGENIC', 'LEVEL_1', 'LEVEL_2', 'LEVEL_3A', 'LEVEL_3B', 'LEVEL_4', 'LEVEL_R1', 'LEVEL_R2', 'HIGHEST_LEVEL', 'HIGHEST_SENSITIVE_LEVEL', 'HIGHEST_RESISTANCE_LEVEL', 'TX_CITATIONS', 'LEVEL_Dx1', 'LEVEL_Dx2', 'LEVEL_Dx3', 'HIGHEST_DX_LEVEL', 'DX_CITATIONS', 'LEVEL_Px1', 'LEVEL_Px2', 'LEVEL_Px3', 'HIGHEST_PX_LEVEL', 'PX_CITATIONS']]

        # Create the join column
        filtered_vcf2maf['join'] = filtered_vcf2maf[['Chromosome', 'Start_Position', 'Reference_Allele', 'Tumor_Seq_Allele2', 'Tumor_Sample_Barcode']].apply(lambda row: ' '.join(str(x) for x in row), axis=1)

        # Write .tsv file
        filtered_vcf2maf.to_csv(os.path.join(merged_results_dir, 'filtered_vcf2maf.tsv'), sep='\t')

    if os.path.exists(os.path.join(merged_results_dir, 'filtered_cravat.tsv')) & os.path.exists(os.path.join(merged_results_dir, 'filtered_cgi.tsv')) & os.path.exists(os.path.join(merged_results_dir, 'filtered_vcf2maf.tsv')):
        
        # Merge dataframes
        merged_tmp = pd.merge(filtered_vcf2maf, filtered_cgi, on='join', how='left')
        merged_final = pd.merge(merged_tmp, filtered_cravat, on='join', how='left')

        # Write merged output
        merged_final.to_csv(os.path.join(merged_results_dir, 'merged.tsv'), sep='\t')

    else:
        logging.error(
            '\n'
            'Unable to merge. Please manually check run'
        )
        sys.exit()