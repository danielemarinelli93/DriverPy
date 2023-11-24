This python script was designed as a scalable tool based on some of the most common variant annotators (VEP, openCRAVAT) to annotate cancer driver mutations. The use of multiple variant annotators yields maximal flexibility, as results from one tool might be available only for a specific annotator and the scripts can be updated with new plugins as they're released. 

Variant annotators/pipelines:
- ENSEMBL-VEP: https://www.ensembl.org/info/docs/tools/vep/index.html
- openCRAVAT: https://opencravat.org/
- Cancer Genome Interpreter: https://www.cancergenomeinterpreter.org/home
- oncoKB: https://github.com/oncokb/oncokb-annotator

Annotations:
- SpliceAI
- LOFTEE
- CHASMplus
- boostDM
- oncoKB

------------

- The input is a 5-column .tsv file: 'chr', 'start', 'ref', 'var', 'tumour_id'. 
- The tool is to be initiated through the "main.py" script. All modules are stored in the "core.py" files. 
- VEP plugins are to be installed according to: https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html
- openCRAVAT modules are to be installed according to: https://open-cravat.readthedocs.io/en/latest/quickstart.html.
- The "configs.txt" file needs to be modified according to the location of binary and reference files for VEP/openCRAVAT. Additional VEP plugins/openCRAVAT modules are to be specified in the configs file. 

------------
To dos:
- add stats
- add final driver + actionable calls
- run vcf2maf after opencravat
- set working directory in script