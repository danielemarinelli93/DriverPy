import pandas as pd
import os, subprocess, requests, shutil, time, argparse, logging, sys, re

config_path = 'configs.txt'
working_dir = 'working_dir'
vep_input_dir = 'working_dir/vep_input/'
vep_output_dir = 'working_dir/vep_output/'
vcf2maf_output_dir = 'working_dir/vcf2maf_output/'
cgi_input_dir = 'working_dir/cgi_input'
cgi_output_dir = 'working_dir/cgi_output'
cravat_output_dir = 'working_dir/cravat_output'
merged_results_dir = 'working_dir/merged_results'

with open(config_path, 'r') as file:
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
        elif line.startswith('vep_fasta_hg37'):
            vep_fasta_hg37 = line.split('=')[1].strip()
        elif line.startswith('vep_fasta_hg38'):
            vep_fasta_hg38 = line.split('=')[1].strip()
        elif line.startswith('loftee_hg37'):
            loftee_hg37 = line.split('=')[1].strip()
        elif line.startswith('loftee_hg38'):
            loftee_hg38 = line.split('=')[1].strip()
        elif line.startswith('SpliceAI_hg37'):
            spliceai_hg37 = line.split('<-')[1].strip()
        elif line.startswith('SpliceAI_hg38'):
            spliceai_hg38 = line.split('<-')[1].strip()
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
        elif line.startswith('cgi_id'):
            cgi_id = line.split('=')[1].strip()
        elif line.startswith('cgi_reference'):
            cgi_reference = line.split('=')[1].strip()
        elif line.startswith('cravat_dir'):
            cravat_dir = line.split('=')[1].strip()


############### ENSEMBL-VEP ###############

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

    # Run VEP
    if not os.path.exists(vep_output_dir):
        os.makedirs(vep_output_dir)
    
    file_list = sorted(os.listdir(vep_input_dir)) 
    for input_file in file_list:
        input_file_path = os.path.join(vep_input_dir, input_file)
        output_file = os.path.join(vep_output_dir, input_file.replace('.tsv', '.vcf'))

        if genome_ver == 'GRCh37':
            vep_fasta = vep_fasta_hg37    
            loftee = loftee_hg37
            spliceai = spliceai_hg37

        elif genome_ver == 'GRCh38':
            vep_fasta = vep_fasta_hg38
            loftee = loftee_hg38
            spliceai = spliceai_hg38

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
            '--hgvs',
            '--symbol',
            '--numbers',
            '--canonical',
            '--variant_class',
            '--no_stats'
        ]
        subprocess.run(command, capture_output=True, text=True)

############### openCRAVAT ###############

### openCRAVAT run
def cravat_run():

    if not os.path.exists(cravat_output_dir):
        os.makedirs(cravat_output_dir)

    for filename in os.listdir(vep_output_dir):
        if filename.endswith('.vcf'):           
            input_file = os.path.join(vep_output_dir, filename)

            if genome_ver == 'GRCh37':
                cravat_genome = 'hg19'
            elif genome_ver == 'GRCh38':
                cravat_genome = 'hg38'

            if oncotree_code == 'NSCLC':
                cravat_anno = 'chasmplus_LUAD chasmplus_LUSC'
            elif oncotree_code == 'BRCA':
                cravat_anno = 'chasmplus_BRCA'
            
            # Run openCRAVAT
            cmd = f"{cravat_dir} run {input_file} -d {cravat_output_dir} -t vcf -l {cravat_genome} -a chasmplus {cravat_anno} {'hg19' if cravat_genome == 'hg19' else ''}"
            subprocess.run(cmd, shell=True)


############### Run VCF2MAF ###############
def vcf2maf_run():    
    if not os.path.exists(vcf2maf_output_dir):
        os.makedirs(vcf2maf_output_dir)

    # Initialize the variables outside the loop
    csq_columns = []
    info_columns = []

    # Process only the first .vcf file
    for file_name in os.listdir(cravat_output_dir):
        if file_name.endswith('.vcf'):
            with open(os.path.join(cravat_output_dir, file_name), 'r') as file:
                for line in file:
                    if line.startswith('##INFO=<ID=CSQ'):
                        csq_columns = re.search(r'Format: (.+)"', line).group(1).split('|')
                    elif line.startswith('##INFO=<ID=OC'):
                        info_column = re.search(r'ID=([^,]+),', line).group(1)
                        info_columns.append(info_column)
            break

    file_list = sorted(file for file in os.listdir(cravat_output_dir) if file.endswith('.vcf.vcf'))

    if genome_ver == 'GRCh37':
        vep_fasta = vep_fasta_hg37   
    elif genome_ver == 'GRCh38':
        vep_fasta = vep_fasta_hg38

    for input_file in file_list:
        input_file_path = os.path.join(cravat_output_dir, input_file)
        output_file = os.path.join(vcf2maf_output_dir, input_file.replace('.vcf.vcf', '.maf'))

        # Prepare the command
        command = [
            'perl',
            vcf2maf_dir,
            '--input-vcf', input_file_path,
            '--output-maf', output_file,
            '--ref-fasta', vep_fasta,
            '--tumor-id', input_file.replace('.vcf', ''),
            '--normal-id', '.',
            '--ncbi-build', genome_ver,
            '--inhibit-vep',
            '--retain-ann', ','.join(csq_columns),
            '--retain-info', ','.join(info_columns)
        ]
        # Execute the command
        subprocess.run(command, capture_output=True, text=True)

    if not os.path.exists(os.path.join(vcf2maf_output_dir, 'merged.maf')):
        
        file_list = os.listdir(vcf2maf_output_dir)
        if len(file_list) == 1:
            for file in file_list:
                x = pd.read_csv(os.path.join(vcf2maf_output_dir, file), sep='\t', comment='#')
                x.to_csv(os.path.join(vcf2maf_output_dir, 'merged.maf'), sep='\t', index=False)

        elif len(file_list) > 1:
            merged = pd.DataFrame()
            for file in file_list:
                x = pd.read_csv(os.path.join(vcf2maf_output_dir, file), sep='\t', comment='#')
                merged = merged.append(x, ignore_index=True)
            merged.to_csv(os.path.join(vcf2maf_output_dir, 'merged.maf'), sep='\t', index=False)
         
        elif len(file_list) == 0:
            logging.error(
                '\n'
                'No input files for VCF2MAF\n'
                'please check VEP run for errors'
            )

############### Run oncokb-annotator ###############
def oncokb_run():
    cmd = f"python3 {os.path.join(oncokb_dir, 'MafAnnotator.py')} -i {os.path.join(vcf2maf_output_dir, 'merged.maf')} -o {os.path.join(vcf2maf_output_dir, 'merged-oncokb.maf')} -t {oncotree_code} -q HGVSp_Short -r {genome_ver} -b {oncokb_token}"
    subprocess.run(cmd, shell=True)

############### Cancer genome interpreter ###############

def server_status(url):
    try:
        response = requests.get(url)
        return response.status_code == 200
    except requests.RequestException:
        return False

### CGI run
def cgi_run():
    
    # Check CGI API status
    cgi_api_url = 'https://www.cancergenomeinterpreter.org'

    if not server_status(cgi_api_url):
        logging.error('CGI server is down. Aborting script.')
        return

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

    if not os.path.exists(cgi_input_dir):
        logging.error('CGI server is down. No CGI results to download.')
        return

    job_completed = False
    while not job_completed:

        # Get CGI job id and status
        cgi_job_id = open(os.path.join(cgi_output_dir, 'cgi_job_id.txt'), 'r').read().strip()
        payload = {'action':'logs'}
        headers = {'Authorization': cgi_id}
        response = requests.get(f'https://www.cancergenomeinterpreter.org/api/v1/{cgi_job_id}', headers=headers, params=payload)
        cgi_log = response.json()
            
        # Download results only if the run is finished
        if cgi_log['status'] == 'Done':
            logging.info(
                '\n'
                'CGI run completed successfully'
                )
            job_completed = True
            payload={'action':'download'}
            r = requests.get(f'https://www.cancergenomeinterpreter.org/api/v1/{cgi_job_id}', headers=headers, params=payload)
            with open('cgi_res.zip', 'wb') as fd:
                fd.write(r._content)
            logging.info(
                '\n'
                'CGI results downloaded successfully'
                )      
            command = ['unzip', 'cgi_res.zip']
            subprocess.run(command, capture_output=True, text=True)
                
            file_list = ('alterations.tsv', 'biomarkers.tsv', 'summary.txt', 'input01.tsv', 'cgi_res.zip')
            for file in file_list:
                shutil.move(file, 'working_dir/cgi_output/')
            
            # Write filtered CGI output
            cgi = pd.read_csv(os.path.join(cgi_output_dir, 'alterations.tsv'), sep='\t')
            cgi['join'] = cgi[['CHROMOSOME', 'POSITION', 'REF', 'ALT', 'CGI-Sample ID']].apply(lambda row: ' '.join(str(x) for x in row), axis=1)
            filtered_cgi_size = cgi[cgi['join'].map(cgi.groupby('join').size()) == 1]
            filtered_cgi = filtered_cgi_size[['join', 'CGI-Oncogenic Summary', 'CGI-Oncogenic Prediction']]
            filtered_cgi.to_csv(os.path.join(cgi_output_dir, 'filtered_cgi.tsv'), sep='\t', index=False)
            
        if cgi_log['status'] == 'Running':
            logging.info(
                '\n'
                'CGI job running Retrying in 10 minutes...'
                )
            time.sleep(600)
            
        if cgi_log['status'] == 'Error':
            logging.error(
                '\n'
                'Error in the CGI run, please check and re-run'
            )
          
    if not os.path.exists(merged_results_dir):
        os.makedirs(merged_results_dir)

    vcf2maf = pd.read_csv(os.path.join(vcf2maf_output_dir, 'merged.maf'), sep='\t', dtype={'Start_Position': int, 'End_Position': int, 'Protein_position': int})
    cgi = pd.read_csv(os.path.join(cgi_output_dir, 'filtered_cgi.tsv'), sep='\t')
    vcf2maf['join'] = vcf2maf[['Chromosome', 'Start_Position', 'Reference_Allele', 'Tumor_Seq_Allele2', 'Tumor_Sample_Barcode']].apply(lambda row: ' '.join(str(x) for x in row), axis=1)
    merged = pd.merge(vcf2maf, cgi, on='join', how='left')
    merged.to_csv(os.path.join(merged_results_dir, 'merged.tsv'), sep='\t')

    cmd = f"python3 {os.path.join(oncokb_dir, 'MafAnnotator.py')} -i {os.path.join(merged_results_dir, 'merged.tsv')} -o {os.path.join(merged_results_dir, 'merged-oncokb.maf')} -t {oncotree_code} -q HGVSp_Short -r {genome_ver} -b {oncokb_token}"
    subprocess.run(cmd, shell=True)