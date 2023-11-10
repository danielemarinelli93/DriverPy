import pandas as pd
import os
import requests
import subprocess
import requests
import shutil

# Set config path
config_path = 'configs.txt'

# Set directories
working_dir = 'working_dir'
cgi_input_dir = 'working_dir/cgi_input'
cgi_output_dir = 'working_dir/cgi_output'

# Read configurations
with open(config_path, 'r') as file:
    
    # Read the contents of the file
        for line in file.readlines():
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

# Read CGI job id
with open(os.path.join(cgi_output_dir, 'cgi_job_id.txt'), 'r') as file:
          
          for line in file.readlines():
               cgi_job_id = line.strip()

############### Cancer genome interpreter output download ###############

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