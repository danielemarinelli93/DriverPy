import pandas as pd
import os, subprocess, requests, shutil, time, argparse, logging, sys

logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger('main')

from core import cgi_run, cgi_download, vep_run, cravat_run, merging_results

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
        elif line.startswith('vep_fasta'):
            vep_fasta = line.split('=')[1].strip()
        elif line.startswith('LoF'):
            loftee = line.strip()
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
        elif line.startswith('cgi_reference'):
            cgi_reference = line.split('=')[1].strip()
        elif line.startswith('cravat_dir'):
            cravat_dir = line.split('=')[1].strip()

def main(args):
    if args.help:
        logging.info(
            '\n'
            'DriverPy: a scalable tool to call driver mutations with multiple variant annotators\n'
            '\n'
            'for more details see: https://github.com/danielemarinelli93/DriverPy'
            'use --run_all to run CGI, VEP and openCRAVAT modules'
        )
        sys.exit()
    elif args.cgi_run:
        logging.info(
            '\n'
            'Running Cancer Genome Interpreter\n'
            f'Account/token: {cgi_id}'
        )        
        cgi_run()
    elif args.cgi_download:
        cgi_download()
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
    elif args.merging_results:
        logging.info(
            '\n'
            'Merging VEP, openCRAVAT and CGI results'
        )
        merging_results()
    elif args.run:
        logging.info(
            '\n'
            'Running all modules'
        )
        cgi_run()
        vep_run()
        cravat_run()
        cgi_download()
        merging_results()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--help', action="store_true")
    parser.add_argument('--cgi_run', action='store_true')
    parser.add_argument('--cgi_download', action='store_true')
    parser.add_argument('--vep_run', action='store_true')
    parser.add_argument('--cravat_run', action='store_true')
    parser.add_argument('--merging_results', action='store_true')
    parser.add_argument('--run', action='store_true')

    args = parser.parse_args()
    main(args)