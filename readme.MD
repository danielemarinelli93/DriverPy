This python command-line script was designed as a scalable tool based on some of the most common variant annotators (VEP, openCRAVAT) to annotate cancer driver mutations. The use of multiple variant annotators yields maximal flexibility, as results from one tool might be available only for a specific annotator and the scripts can be updated with new plugins as they're released. 

Variant annotators/pipelines:
- ENSEMBL-VEP: https://www.ensembl.org/info/docs/tools/vep/index.html, v110.1
- openCRAVAT: https://opencravat.org/, v2.4.1
- VCF2MAF: https://github.com/mskcc/vcf2maf, v1.6.21
- oncoKB: https://github.com/oncokb/oncokb-annotator, v.3.4.0
- Cancer Genome Interpreter: https://www.cancergenomeinterpreter.org/

The core of the tool runs VEP on the 5-column .tsv input file ('chr', 'start', 'ref', 'var', 'tumour_id', according to https://www.ensembl.org/info/docs/tools/vep/vep_formats.html). 
A separate VEP input file for each unique 'tumour_id' is created in the 'vep_input' folder. The VEP module outputs a .vcf file for each unique 'tumour_id' in the 'vep_output' folder which is then used as an input for openCRAVAT. The updated .vcf files resulting from the openCRAVAT run are then converted to .maf files by VFC2MAF, retaining all VEP and openCRAVAT annotations. 
All outputs are in the 'working_dir' directory which is created in the directory from which the script is run. 

Since both VEP and openCRAVAT can rely on cache files, the core of the DriverPy tool only includes them. This version of DriverPy uses LOFTEE and SpliceAI plugins from VEP and CHASMplus, CHASMplus_LUAD, CHASMplus_LUSC and CHASMplus_BRCA modules from openCRAVAT. 
Additional plugins/modules can be integrated with ease by slightly modifying the script and the configs file.

Separate modules run Cancer Genome Interpreter (https://www.cancergenomeinterpreter.org/) and OncoKB (https://www.oncokb.org/) annotators, which are both API-based and theferore need identification (token and/or email). 
The oncokb-annotator runs on the output .maf file from the VCF2MAF run. 
Results from the Cancer Genome Interpreter (actionable calls, boostDM annotations) need to be merged with the .maf file from the VCF2MAF run; the easiest way to do it is through the genomic position and 'tumour_id'. 

I suggest to run DriverPy in a dedicated conda environment with python 3.7, either locally or on a server. 
The core of the tool through command line with the 'python3 main.py --all' command. 
All modules are stored in the "core.py" file.
VEP plugins are to be installed according to: https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html
openCRAVAT modules are to be installed according to: https://open-cravat.readthedocs.io/en/latest/quickstart.html.
The "configs.txt" file needs to be modified according to the location of binary and reference files for VEP/openCRAVAT. Additional VEP plugins/openCRAVAT modules are to be specified in the configs file.
To avoid conflicts with VCF2MAF, I suggest to unzip the VEP fasta files.

