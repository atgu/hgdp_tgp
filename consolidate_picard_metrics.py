#!/usr/bin/env python3

__author__ = 'Lindo Nkambule'

import pandas as pd

bams_paths = pd.read_csv('gs://african-seq-data/gambian-genomes/text_files/all_394_bams.txt', sep='\t', header=None)
bam_files = []

for index, row in bams_paths.iterrows():
    bam_files.append(row[0])

all_files_stats_list = []

path = 'gs://african-seq-data/gambian-genomes/variant-calling-metrics'

for bam in bam_files:
    file = pd.read_csv(f'{path}/{bam}/{bam}.variant_calling_detail_metrics',
                       skiprows=6, nrows=7, sep='\t')

    file['SAMPLE_ALIAS'] = bam

    bam_stats_list = file.values.flatten().tolist()

    all_files_stats_list.append(bam_stats_list)

header_cols = ['SAMPLE', 'HET_HOMVAR_RATIO', 'PCT_GQ0_VARIANTS', 'TOTAL_GQ0_VARIANTS', 'TOTAL_HET_DEPTH', 'TOTAL_SNPS',
               'NUM_IN_DB_SNP', 'NOVEL_SNPS', 'FILTERED_SNPS', 'PCT_DBSNP', 'DBSNP_TITV', 'NOVEL_TITV', 'TOTAL_INDELS',
               'NOVEL_INDELS', 'FILTERED_INDELS', 'PCT_DBSNP_INDELS', 'NUM_IN_DB_SNP_INDELS', 'DBSNP_INS_DEL_RATIO',
               'NOVEL_INS_DEL_RATIO', 'TOTAL_MULTIALLELIC_SNPS', 'NUM_IN_DB_SNP_MULTIALLELIC', 'TOTAL_COMPLEX_INDELS',
               'NUM_IN_DB_SNP_COMPLEX_INDELS', 'SNP_REFERENCE_BIAS', 'NUM_SINGLETONS']

all_files_stats_df = pd.DataFrame(all_files_stats_list, columns = header_cols)

all_files_stats_df.to_csv('gs://african-seq-data/gambian-genomes/consolidated.variant_calling_detail_metrics',
                        sep='\t', index=False)
