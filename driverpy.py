import pandas as pd
import os, subprocess, requests, shutil, time, argparse, logging, sys, re

logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger('main')

config_path = '/nemo/lab/swantonc/working/griecoc/Tx842/drivers/scripts/driver_analysis/driver_anno/daniele_pipeline/2025-09-26_PEACErun/configs.txt'
working_dir = '/nemo/lab/swantonc/working/griecoc/Tx842/drivers/outputs/driver_anno/daniele_pipeline/working_dir/2025-11-07_PEACE/'
vep_input_dir = os.path.join(working_dir, 'vep_input/')
vep_output_dir = os.path.join(working_dir, 'vep_output/')
vcf2maf_output_dir = os.path.join(working_dir, 'vcf2maf_output/')
cgi_input_dir = os.path.join(working_dir, 'cgi_input/')
cgi_output_dir = os.path.join(working_dir, 'cgi_output/')
cravat_output_dir = os.path.join(working_dir, 'cravat_output/')
merged_results_dir = os.path.join(working_dir, 'merged_results/')

with open(config_path, 'r')import pandas as pd
import os, subprocess, requests, shutil, time, argparse, logging, sys, re

logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger('main')

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
        elif line.startswith('AlphaMissense_hg37'):
            AlphaMissense_hg37 = line.split('<-')[1].strip()
        elif line.startswith('AlphaMissense_hg38'):
            AlphaMissense_hg38 = line.split('<-')[1].strip()         
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
 as file:
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
        elif line.startswith('AlphaMissense_hg37'):
            AlphaMissense_hg37 = line.split('<-')[1].strip()
        elif line.startswith('AlphaMissense_hg38'):
            AlphaMissense_hg38 = line.split('<-')[1].strip()         
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
    try:
        # Create directories
        if not os.path.exists(vep_input_dir):
            os.makedirs(vep_input_dir)
            log.info(f"Created VEP input directory: {vep_input_dir}")

        # Read input dataframe removing duplicated rows
        try:
            mut = pd.read_csv(first_input, sep='\t').drop_duplicates()
            log.info(f"Loaded input file: {first_input}")
        except Exception as e:
            log.error(f"Failed to read input file {first_input}: {e}")
            return

        # Convert start position so integer
        try:
            mut['start'] = mut['start'].astype(int)
        except Exception as e:
            log.error(f"Failed to convert 'start' column to int: {e}")
            return

        # Set end position based on VEP requirements
        mut['end'] = mut['start']
        mut['end'] = mut['end'].astype(int)

        def delins_pos(row):
            # deletions
            if row['var'].startswith('-'):
                row['end'] = row['start'] + len(row['ref']) - 1
            # insertions
            elif row['ref'].startswith('-'):
                row['start'] = row['end'] + 1
            return row

        def mnv_pos(row):
            if len(row['var']) == len(row['ref']) and len(row['ref']) > 1:
                row['end'] = row['start'] + len(row['var']) - 1
            return row

        mut = mut.apply(mnv_pos, axis=1)
        mut = mut.apply(delins_pos, axis=1)

        # VEP requirement for allele
        mut['allele'] = mut['ref'] + '/' + mut['var']

        # Set strand column
        mut['strand'] = '+'

        # Select columns
        try:
            mut = mut[['chr', 'start', 'end', 'allele', 'strand', 'patient_tumour']]
            mut.rename(columns={'patient_tumour': 'tumour_id'}, inplace=True)
        except Exception as e:
            log.error(f"Error selecting/renaming columns: {e}")
            return

        # Arrange the input file for a quicker VEP run
        mut = mut.sort_values(['chr', 'start'])

        # Group the rows based on tumour_id
        tumours = mut.groupby('tumour_id')

        # Loop through the groups and write each group to a separate tsv file
        for group_name, group_vep in tumours:
            output_file_path = os.path.join(vep_input_dir, f'{group_name}.tsv')
            try:
                group_vep.to_csv(output_file_path, sep='\t', index=False, header=False)
                log.info(f"Created VEP input file: {output_file_path}")
            except Exception as e:
                log.error(f"Failed to write VEP input file {output_file_path}: {e}")

        log.info('VEP input files successfully created')

        # Run VEP
        if not os.path.exists(vep_output_dir):
            os.makedirs(vep_output_dir)
            log.info(f"Created VEP output directory: {vep_output_dir}")

        file_list = sorted(os.listdir(vep_input_dir))
        if not file_list:
            log.error("No VEP input files found to process.")
            return

        for input_file in file_list:
            input_file_path = os.path.join(vep_input_dir, input_file)
            output_file = os.path.join(vep_output_dir, input_file.replace('.tsv', '.vcf'))

            if genome_ver == 'GRCh37':
                vep_fasta = vep_fasta_hg37
                loftee = loftee_hg37
                spliceai = spliceai_hg37
                AlphaMissense = AlphaMissense_hg37
            elif genome_ver == 'GRCh38':
                vep_fasta = vep_fasta_hg38
                loftee = loftee_hg38
                spliceai = spliceai_hg38
                AlphaMissense = AlphaMissense_hg38
            else:
                log.error(f"Unsupported genome version: {genome_ver}")
                continue

            command = [
                'singularity', 'exec', '-B', '/nemo:/nemo',
                '/nemo/lab/swantonc/working/griecoc/singularity/ensembl-vepv114.sif', '/opt/vep/src/ensembl-vep/vep',
                '-i', input_file_path,
                '-o', output_file,
                '--dir_cache', vep_cache_dir,
                '--fasta', vep_fasta,
                '--offline',
                '--assembly', genome_ver,
                '--pick',
                '--pick_order', 'mane_select,canonical,rank,tsl,appris,biotype', '--ccds',
                '--cache',
                '--plugin', loftee,
                '--dir_plugins', '/nemo/lab/swantonc/working/griecoc/singularity/vep_data/loftee/plugin/GRCh37/',
                '--plugin', spliceai,
                '--plugin', AlphaMissense,
                '--vcf',
                '--hgvs',
                '--symbol',
                '--numbers',
                '--canonical',
                '--variant_class',
                '--no_stats',
                '--buffer_size', '100000',
            ]
            log.info(f"Running VEP for {input_file_path} -> {output_file}")
            result = subprocess.run(command, capture_output=True, text=True)
            if result.returncode != 0:
                log.error(f"VEP failed for {input_file_path}: {result.stderr}")
            else:
                log.info(f"VEP completed for {input_file_path}")

    except Exception as e:
        log.error(f"Unexpected error in vep_run: {e}")

############### openCRAVAT ###############

### openCRAVAT run
def cravat_run():
    try:
        if not os.path.exists(cravat_output_dir):
            os.makedirs(cravat_output_dir)
            log.info(f"Created openCRAVAT output directory: {cravat_output_dir}")

        vep_files = [f for f in os.listdir(vep_output_dir) if f.endswith('.vcf')]
        if not vep_files:
            log.error("No VEP output files found for openCRAVAT.")
            return

        for filename in vep_files:
            input_file = os.path.join(vep_output_dir, filename)

            if genome_ver == 'GRCh37':
                cravat_genome = 'hg19'
            elif genome_ver == 'GRCh38':
                cravat_genome = 'hg38'
            else:
                log.error(f"Unsupported genome version for openCRAVAT: {genome_ver}")
                continue

            if oncotree_code == 'NSCLC':
                cravat_anno = 'chasmplus_LUAD chasmplus_LUSC'
            elif oncotree_code == 'BRCA':
                cravat_anno = 'chasmplus_BRCA'
            else:
                cravat_anno = ''
                log.warning(f"No specific openCRAVAT annotation for oncotree_code: {oncotree_code}")

            cmd = f"{cravat_dir} run {input_file} -d {cravat_output_dir} -t vcf -l {cravat_genome} -a chasmplus {cravat_anno} {'hg19' if cravat_genome == 'hg19' else ''} --mp 1"
            log.info(f"Running openCRAVAT: {cmd}")
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            if result.returncode != 0:
                log.error(f"openCRAVAT failed for {input_file}: {result.stderr}")
            else:
                log.info(f"openCRAVAT completed for {input_file}")

    except Exception as e:
        log.error(f"Unexpected error in cravat_run: {e}")

############### Run VCF2MAF ###############
def vcf2maf_run():
    try:
        if not os.path.exists(vcf2maf_output_dir):
            os.makedirs(vcf2maf_output_dir)
            log.info(f"Created VCF2MAF output directory: {vcf2maf_output_dir}")

        csq_columns = []
        info_columns = []

        # Identify VEP and openCRAVAT VCF files
        vep_vcf_files = [f for f in os.listdir(vep_output_dir) if f.endswith('.vcf')]
        cravat_vcf_files = [f for f in os.listdir(cravat_output_dir) if f.endswith('.vcf.vcf')]

        # For each openCRAVAT VCF, create a .filled.vcf with missing VEP rows added (by first 5 columns)
        filled_vcf_files = []
        for cravat_file in cravat_vcf_files:
            vep_file = cravat_file.replace('.vcf.vcf', '.vcf')
            vep_path = os.path.join(vep_output_dir, vep_file)
            cravat_path = os.path.join(cravat_output_dir, cravat_file)
            filled_path = os.path.join(cravat_output_dir, cravat_file.replace('.vcf.vcf', '.filled.vcf'))

            # Read header and data from openCRAVAT VCF
            with open(cravat_path, 'r') as cravat_f:
                cravat_lines = cravat_f.readlines()
            cravat_header = [line for line in cravat_lines if line.startswith('#')]
            cravat_data = [line for line in cravat_lines if not line.startswith('#')]

            # Read data from VEP VCF
            if os.path.exists(vep_path):
                with open(vep_path, 'r') as vep_f:
                    vep_lines = vep_f.readlines()
                vep_data = [line for line in vep_lines if not line.startswith('#')]
            else:
                vep_data = []

            # Build sets of keys based on first 5 columns (tab-separated)
            def key5(line):
                return tuple(line.strip().split('\t')[:5])

            cravat_keys = set(key5(line) for line in cravat_data)
            vep_keys = set(key5(line) for line in vep_data)

            # Find missing data rows in openCRAVAT (by first 5 columns)
            missing_keys = vep_keys - cravat_keys
            missing_rows = [line for line in vep_data if key5(line) in missing_keys]

            # Write new .filled.vcf file
            with open(filled_path, 'w') as filled_f:
                for line in cravat_header:
                    filled_f.write(line)
                for line in sorted(cravat_data):
                    filled_f.write(line)
                for line in sorted(missing_rows):
                    filled_f.write(line)
            log.info(f"Created filled openCRAVAT VCF: {filled_path} ({len(missing_rows)} rows added)")
            filled_vcf_files.append(os.path.basename(filled_path))

        # Use .filled.vcf files for downstream processing
        if not filled_vcf_files:
            log.error("No .filled.vcf files found for VCF2MAF.")
            return

        # Extract CSQ/info columns from the first .filled.vcf file
        for file_name in filled_vcf_files:
            with open(os.path.join(cravat_output_dir, file_name), 'r') as file:
                for line in file:
                    if line.startswith('##INFO=<ID=CSQ'):
                        match = re.search(r'Format: (.+)"', line)
                        if match:
                            csq_columns = match.group(1).split('|')
                        else:
                            log.warning(f"Could not parse CSQ columns in {file_name}")
                    elif line.startswith('##INFO=<ID=OC'):
                        info_column = re.search(r'ID=([^,]+),', line).group(1)
                        info_columns.append(info_column)
            break

        # Run VCF2MAF on each .filled.vcf file
        for input_file in filled_vcf_files:
            input_file_path = os.path.join(cravat_output_dir, input_file)
            output_file = os.path.join(vcf2maf_output_dir, input_file.replace('.filled.vcf', '.maf'))

            if genome_ver == 'GRCh37':
                vep_fasta = vep_fasta_hg37
            elif genome_ver == 'GRCh38':
                vep_fasta = vep_fasta_hg38
            else:
                log.error(f"Unsupported genome version for VCF2MAF: {genome_ver}")
                return

            command = [
                'perl',
                vcf2maf_dir,
                '--input-vcf', input_file_path,
                '--output-maf', output_file,
                '--ref-fasta', vep_fasta,
                '--tumor-id', input_file.replace('.filled.vcf', ''),
                '--normal-id', '.',
                '--ncbi-build', genome_ver,
                '--inhibit-vep',
                '--retain-ann', ','.join(csq_columns),
                '--retain-info', ','.join(info_columns)
            ]
            log.info(f"Running VCF2MAF: {' '.join(command)}")
            result = subprocess.run(command, capture_output=True, text=True)
            if result.returncode != 0:
                log.error(f"VCF2MAF failed for {input_file_path}: {result.stderr}")
            else:
                log.info(f"VCF2MAF completed for {input_file_path}")

        # Merge all MAF files into a single MAF (unchanged)
        merged_maf_path = os.path.join(vcf2maf_output_dir, 'merged.maf')
        if not os.path.exists(merged_maf_path):
            file_list = [f for f in os.listdir(vcf2maf_output_dir) if f.endswith('.maf')]
            if len(file_list) == 1:
                try:
                    for file in file_list:
                        x = pd.read_csv(os.path.join(vcf2maf_output_dir, file), sep='\t', comment='#')
                        x.to_csv(merged_maf_path, sep='\t', index=False)
                        log.info(f"Single MAF file merged: {merged_maf_path}")
                except Exception as e:
                    log.error(f"Error merging single MAF file: {e}")
            elif len(file_list) > 1:
                dfs = []
                for file in file_list:
                    try:
                        x = pd.read_csv(os.path.join(vcf2maf_output_dir, file), sep='\t', comment='#')
                        dfs.append(x)
                    except Exception as e:
                        log.error(f"Error reading MAF file {file}: {e}")
                try:
                    merged = pd.concat(dfs, ignore_index=True)
                    merged.to_csv(merged_maf_path, sep='\t', index=False)
                    log.info(f"Multiple MAF files merged: {merged_maf_path}")
                except Exception as e:
                    log.error(f"Error writing merged MAF file: {e}")
            elif len(file_list) == 0:
                log.error(
                    '\n'
                    'No input files for VCF2MAF\n'
                    'please check VEP run for errors'
                )
    except Exception as e:
        log.error(f"Unexpected error in vcf2maf_run: {e}")

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
                f'please delete (os.path.join(cgi_output_dir, "cgi_job_id.txt")) and re-run CGI'
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
                files={'mutations': open(os.path.join(cgi_input_dir, 'cgi_input.tsv'), 'rb')},
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
                shutil.move(file, cgi_output_dir) #quotes in 'cgi_output_dir' meant this step was not working

            # Write filtered CGI output
            cgi = pd.read_csv(os.path.join(cgi_output_dir, 'alterations.tsv'), sep='\t')
            cgi['join'] = cgi[['CHROMOSOME', 'POSITION', 'REF', 'ALT', 'CGI-Sample ID']].apply(lambda row: ' '.join(str(x) for x in row), axis=1)
            filtered_cgi = cgi[cgi['join'].map(cgi.groupby('join').size()) == 1]
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

    vcf2maf = pd.read_csv(os.path.join(vcf2maf_output_dir, 'merged-oncokb.maf'), sep='\t')
    cgi = pd.read_csv(os.path.join(cgi_output_dir, 'filtered_cgi.tsv'), sep='\t')
    vcf2maf['join'] = vcf2maf[['Chromosome', 'Start_Position', 'Reference_Allele', 'Tumor_Seq_Allele2', 'Tumor_Sample_Barcode']].apply(lambda row: ' '.join(str(x) for x in row), axis=1)
    merged = pd.merge(vcf2maf, cgi, on='join', how='left')
    merged.to_csv(os.path.join(merged_results_dir, 'merged-oncokb-cgi.maf'), sep='\t', index=False)


###

def main(args):
    if args.help:
        logging.info(
            '\n'
            'DriverPy: a scalable tool to call driver mutations with multiple variant annotators\n'
            '\n'
            'for more details see: https://github.com/danielemarinelli93/DriverPy'
        )
        sys.exit()
    elif args.vep_run:
        logging.info(
            '\n'
            'Running ENSEMBL-VEP\n'
            'for more details see: https://www.ensembl.org/info/docs/tools/vep/script/index.html'
        )
        vep_run()
    elif args.cravat_run:
        logging.info(
            '\n'
            'Running openCRAVAT'
            'for more details see: https://open-cravat.readthedocs.io/en/latest/quickstart.html'
        )
        cravat_run()
    elif args.vcf2maf_run:
        logging.info(
            '\n'
            'Running VCF2MAF'
        )
        vcf2maf_run()
    elif args.all:
        logging.info(
            '\n'
            'Running VEP, openCRAVAT, VCF2MAF'
        )
        vep_run()
        cravat_run()
        vcf2maf_run()
    elif args.cgi_run:
        logging.info(
            '\n'
            'Running Cancer Genome Interpreter\n'
            f'Account/token: {cgi_id}'
        )        
        cgi_run()
    elif args.cgi_download:
        cgi_download()
    elif args.oncokb_run:
        logging.info(
            '\n'
            'Running the oncokb-annotator on vcf2maf merged output'
            'For more details see: https://github.com/oncokb/oncokb-annotator'
        )
        oncokb_run()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--help', action="store_true")
    parser.add_argument('--cgi_run', action='store_true')
    parser.add_argument('--cgi_download', action='store_true')
    parser.add_argument('--vep_run', action='store_true')
    parser.add_argument('--cravat_run', action='store_true')
    parser.add_argument('--vcf2maf_run', action='store_true')
    parser.add_argument('--oncokb_run', action='store_true')
    parser.add_argument('--all', action='store_true')
    
    args = parser.parse_args()
    main(args)
