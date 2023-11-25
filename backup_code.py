    # Run VCF2MAF
    if not os.path.exists(vcf2maf_output_dir):
        os.makedirs(vcf2maf_output_dir)

    file_list = sorted(file for file in os.listdir(vep_output_dir) if not file.endswith('warnings.txt'))
    for input_file in file_list:
        input_file_path = os.path.join(vep_output_dir, input_file)
        output_file = os.path.join(vcf2maf_output_dir, input_file.replace('.vcf', '.maf'))

        # Run VCF2MAF
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
            '--retain-ann', retain_ann
            ]
        subprocess.run(command, capture_output=True, text=True)
    
    if not os.path.exists(os.path.join(vcf2maf_output_dir, 'merged.tsv')):
        
        file_list = os.listdir(vcf2maf_output_dir)
        if len(file_list) == 1:
            for file in file_list:
                x = pd.read_csv(os.path.join(vcf2maf_output_dir, file), sep='\t', comment='#')
                x.to_csv(os.path.join(vcf2maf_output_dir, 'merged.tsv'), sep='\t', index=False)

        elif len(file_list) > 1:
            merged = pd.DataFrame()
            for file in file_list:
                x = pd.read_csv(os.path.join(vcf2maf_output_dir, file), sep='\t', comment='#')
                merged = merged.append(x, ignore_index=True)
            merged.to_csv(os.path.join(vcf2maf_output_dir, 'merged.tsv'), sep='\t', index=False)
         
        elif len(file_list) == 0:
            logging.error(
                '\n'
                'No input files for VCF2MAF\n'
                'please check VEP run for errors'
            )

    # Write filtered output
    filtered_vcf2maf = merged[['Hugo_Symbol', 'NCBI_Build', 'Chromosome', 'Start_Position', 'End_Position', 'Strand', 'Variant_Classification', 'Variant_Type', 'Reference_Allele', 'Tumor_Seq_Allele2', 'Tumor_Sample_Barcode', 'HGVSc', 'HGVSp', 'HGVSp_Short', 'Transcript_ID', 'EXON', 'INTRON', 'all_effects', 'Consequence', 'BIOTYPE', 'CANONICAL', 'IMPACT', 'LoF', 'LoF_filter', 'LoF_flags', 'SpliceAI_pred_SYMBOL', 'SpliceAI_pred_DS_AG', 'SpliceAI_pred_DS_AL', 'SpliceAI_pred_DS_DG', 'SpliceAI_pred_DS_DL']]
    filtered_vcf2maf['join'] = filtered_vcf2maf[['Chromosome', 'Start_Position', 'Reference_Allele', 'Tumor_Seq_Allele2', 'Tumor_Sample_Barcode']].apply(lambda row: ' '.join(str(x) for x in row), axis=1)
    filtered_vcf2maf.to_csv(os.path.join(vcf2maf_output_dir, 'filtered_vcf2maf.tsv'), sep='\t')

    # Run oncokb-annotator
    cmd = f"python3 {os.path.join(oncokb_dir, 'MafAnnotator.py')} -i {os.path.join(vcf2maf_output_dir, 'filtered_vcf2maf.tsv')} -o {os.path.join(vcf2maf_output_dir, 'merged_filtered_oncokb.tsv')} -t {oncotree_code} -q HGVSp_Short -r {genome_ver} -b {oncokb_token}"
    subprocess.run(cmd, shell=True)




    # Write merged and filtered dataframes
    file_list = [file for file in os.listdir(cravat_output_dir) if file.endswith('variant.tsv')]
    merged = pd.DataFrame()
    for file in file_list:
        x = pd.read_csv(os.path.join(cravat_output_dir, file), sep='\t', comment='#')
        merged = merged.append(x, ignore_index=True)
    merged.to_csv(os.path.join(cravat_output_dir, 'merged.tsv'), sep='\t', index=False)

    if genome_ver == 'GRCh37':   
        merged['original_input.chrom'] = merged['original_input.chrom'].str.replace('chr', '')
        merged['join'] = merged[['original_input.chrom', 'original_input.pos', 'original_input.ref_base', 'original_input.alt_base', 'tags']].apply(lambda row: ' '.join(str(x) for x in row), axis=1)
    elif genome_ver == 'GRCh38':
        merged['chrom'] = merged['chrom'].str.replace('chr', '')
        merged['join'] = merged[['chrom', 'pos', 'ref_base', 'alt_base', 'tags']].apply(lambda row: ' '.join(str(x) for x in row), axis=1)
    
    if oncotree_code == 'NSCLC':
        cravat_columns = ['chasmplus_LUAD.pval', 'chasmplus_LUAD.score', 'chasmplus_LUAD.transcript', 'chasmplus_LUAD.all', 'chasmplus_LUSC.pval', 'chasmplus_LUSC.score', 'chasmplus_LUSC.transcript', 'chasmplus_LUSC.all']
    elif oncotree_code == 'BRCA':
        cravat_columns = ['chasmplus_BRCA.pval', 'chasmplus_BRCA.score', 'chasmplus_BRCA.transcript', 'chasmplus_BRCA.all']
    
    filtered_cravat = merged[['join', 'chasmplus.pval', 'chasmplus.score', 'chasmplus.transcript', 'chasmplus.all'] + cravat_columns]
    filtered_cravat.to_csv(os.path.join(cravat_output_dir, 'filtered_cravat.tsv'), sep='\t', index=False)



