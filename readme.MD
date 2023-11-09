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

The input is a 5-column .tsv file: 'chr', 'start', 'ref', 'var', 'tumour_id'. 
The script was built and tested in a conda environment with vep, vcf2maf, oncokbannotator and opencravat installed with all their dependencies.