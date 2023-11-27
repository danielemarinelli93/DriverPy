This python script was designed as a scalable tool based on some of the most common variant annotators (VEP, openCRAVAT) to annotate cancer driver mutations. The use of multiple variant annotators yields maximal flexibility, as results from one tool might be available only for a specific annotator and the scripts can be updated with new plugins as they're released. 

Variant annotators/pipelines:
- ENSEMBL-VEP: https://www.ensembl.org/info/docs/tools/vep/index.html, v110.1
- openCRAVAT: https://opencravat.org/, v2.4.1
- Cancer Genome Interpreter: https://www.cancergenomeinterpreter.org/home
- VCF2MAF: https://github.com/mskcc/vcf2maf, v1.6.21
- oncoKB: https://github.com/oncokb/oncokb-annotator, v.3.4.0

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
- To avoid conflicts with vcf2maf, the VEP fasta files should be unzipped

------------
To dos:
- add stats
- add final driver + actionable calls
- set working directory in configs
- get the APIs (oncoKB + CGI) as an optional add-on
