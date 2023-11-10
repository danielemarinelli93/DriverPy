import pandas as pd
import os
import shutil

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
if not os.path.exists(merged_results_dir):
    os.makedirs(merged_results_dir)

# Read openCRAVAT merged output
cravat = pd.read_csv(os.path.join(cravat_output_dir, 'merged.tsv'), sep='\t')

# Remove the 'chr' prefix
cravat['original_input.chrom'] = cravat['original_input.chrom'].str.replace('chr', '')

# Create join column
cravat['join'] = cravat[['original_input.chrom', 'original_input.pos', 'original_input.ref_base', 'original_input.alt_base', 'tags']].apply(lambda row: ' '.join(str(x) for x in row), axis=1)

# Select columns
filtered_cravat = cravat[['join', 'chasmplus.pval', 'chasmplus.score', 'chasmplus.transcript', 'chasmplus.all', 'chasmplus_LUAD.pval', 'chasmplus_LUAD.score', 'chasmplus_LUAD.transcript', 'chasmplus_LUAD.all', 'chasmplus_LUSC.pval', 'chasmplus_LUSC.score', 'chasmplus_LUSC.transcript', 'chasmplus_LUSC.all']]

# Write .tsv file
filtered_cravat.to_csv(os.path.join(merged_results_dir, 'filtered_cravat.tsv'), sep='\t', index=False)

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

# Read VCF2maf merged output
vcf2maf = pd.read_csv(os.path.join(vcf2maf_output_dir, 'merged.tsv'), sep='\t')

# Select columns
filtered_vcf2maf = vcf2maf[['Hugo_Symbol', 'NCBI_Build', 'Chromosome', 'Start_Position', 'End_Position', 'Strand', 'Variant_Classification', 'Variant_Type', 'Reference_Allele', 'Tumor_Seq_Allele2', 'Tumor_Sample_Barcode', 'HGVSc', 'HGVSp', 'HGVSp_Short', 'Transcript_ID', 'EXON', 'INTRON', 'all_effects', 'Consequence', 'BIOTYPE', 'CANONICAL', 'IMPACT', 'LoF', 'LoF_filter', 'LoF_flags', 'SpliceAI_pred_SYMBOL', 'SpliceAI_pred_DS_AG', 'SpliceAI_pred_DS_AL', 'SpliceAI_pred_DS_DG', 'SpliceAI_pred_DS_DL']]

# Create the join column
filtered_vcf2maf['join'] = filtered_vcf2maf[['Chromosome', 'Start_Position', 'Reference_Allele', 'Tumor_Seq_Allele2', 'Tumor_Sample_Barcode']].apply(lambda row: ' '.join(str(x) for x in row), axis=1)

# Write .tsv file
filtered_vcf2maf.to_csv(os.path.join(merged_results_dir, 'filtered_vcf2maf.tsv'), sep='\t')

# Merge dataframes
merged_tmp = pd.merge(filtered_vcf2maf, filtered_cgi, on='join', how='left')
merged_final = pd.merge(merged_tmp, filtered_cravat, on='join', how='left')

# Write merged output
merged_final.to_csv(os.path.join(merged_results_dir, 'merged.tsv'), sep='\t')
