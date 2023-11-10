import pandas as pd
import os
import requests

# Set config path
config_path = 'configs.txt'

# Set directories
working_dir = 'working_dir'
cgi_input_dir = 'working_dir/cgi_input'
cgi_output_dir = 'working_dir/cgi_output'

# Create directories
if not os.path.exists(working_dir):
    os.makedirs(working_dir)
if not os.path.exists(cgi_input_dir):
    os.makedirs(cgi_input_dir)
if not os.path.exists(cgi_output_dir):
    os.makedirs(cgi_output_dir)

# Read configurations
with open(config_path, 'r') as file:
    # Read the contents of the file
    contents = file.readlines()

    for line in contents:
        if line.startswith('input'):
            first_input = line.split('=')[1].strip()
        elif line.startswith('oncotree_code'):
            oncotree_code = line.split('=')[1].strip()
        elif line.startswith('cgi_id'):
            cgi_id = line.split('=')[1].strip()
        elif line.startswith('cgi_title'):
            cgi_title = line.split('=')[1].strip()
        elif line.startswith('cgi_reference'):
            cgi_reference = line.split('=')[1].strip()


############### Cancer genome interpreter ###############

# Check if the CGI run was already started
if os.path.exists(os.path.join(cgi_output_dir, 'cgi_job_id.txt')):

    # Read the CGI job id from the .txt file    
    with open(os.path.join(cgi_output_dir, 'cgi_job_id.txt'), 'r') as cgi_input:
        cgi_job_id = cgi_input.read().strip()
    
    # Get the list of existing job ids from the CGI API
    headers = {'Authorization': cgi_id}
    r = requests.get('https://www.cancergenomeinterpreter.org/api/v1', headers=headers)
    
    # Run CGI only if the same job id is not already running
    if cgi_job_id not in r.text:

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

        # Write CGI job id to a text file
        with open(os.path.join(cgi_output_dir, 'cgi_job_id.txt'), 'w') as file:

            # Write the string to the file
            file.write(cgi_job_id)
        print('CGI job submitted')

    else: 
        print('CGI job not submitted: a job with the same id was already submitted')