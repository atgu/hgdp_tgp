# Author: Zan Koenig

import hail as hl

# Function from gnomAD library to apply genotype filters
from gnomad.utils.filtering import filter_to_adj

hl.init(tmp_dir='gs://hgdp-1kg/znk-temporary-files/')


# Functions
# Set up function to:
# apply gnomAD's sample, variant and genotype QC filters
# remove two contaminated samples identified using CHARR - https://pubmed.ncbi.nlm.nih.gov/37425834/
# remove the gnomAD sample that's added for QC purposes
# add gnomAD's HGDP+1kGP metadata with the updated population labels as a column field


def run_qc(mt):
    ## Apply sample QC filters to dataset
    # This filters to only samples that passed gnomAD's sample QC hard filters
    mt = mt.filter_cols(~mt.gnomad_sample_filters.hard_filtered)  # removed 31 samples

    ## Apply variant QC filters to dataset
    # This subsets to only PASS variants - those which passed gnomAD's variant QC
    # PASS variants have an entry in the filters field
    mt = mt.filter_rows(hl.len(mt.filters) != 0, keep=False)

    # Remove the two contaminated samples identified by CHARR and 'CHMI_CHMI3_WGS2'
    contaminated_samples = {'HGDP01371', 'LP6005441-DNA_A09'}
    contaminated_samples_list = hl.literal(contaminated_samples)
    mt = mt.filter_cols(~contaminated_samples_list.contains(mt['s']))

    # CHMI_CHMI3_WGS2 is a sample added by gnomAD for QC purposes and has no metadata info
    mt = mt.filter_cols(mt.s == 'CHMI_CHMI3_WGS2', keep=False)

    # Only keep the variants which are found in the samples that are left
    mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))

    # Read in and add the metadata with the updated population labels as a column field
    metadata = hl.import_table(metadata_path, impute=True, key='s')
    mt = mt.annotate_cols(meta_updated=metadata[mt.s])

    ## Apply genotype QC filters to the dataset
    # This is done using a function imported from gnomAD and is the last step in the QC process
    mt = filter_to_adj(mt)

    return mt


def remove_pca_outliers(mt, outlier_path):
    # Remove PCA outliers from the HGDP+1kGP dataset
    # Use hl.hadoop_open to read in the PCA outliers file into Hail from Google Cloud Storage
    with hl.utils.hadoop_open(outlier_path) as file:
        outliers = [line.rstrip('\n') for line in file]

    # Use hl.literal to convert the outliers list from a python object to a Hail expression so that it can be used to
    # filter out samples
    outlier_list = hl.literal(outliers)

    # Keep the samples which are not contained in the pca outlier list
    mt = mt.filter_cols(~outlier_list.contains(mt['s']))

    # Only keep the variants which are found in the samples that are left
    mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))

    return mt


# Reading in data

# Path for HGDP+1kGP dataset prior to applying gnomAD QC filters
pre_qc_path = 'gs://gcp-public-data--gnomad/release/3.1.2/mt/genomes/gnomad.genomes.v3.1.2.hgdp_1kg_subset_dense.mt'

# PCA outliers file
outliers_path = 'gs://gcp-public-data--gnomad/release/3.1/secondary_analyses/hgdp_1kg/pca/pca_outliers.txt'

# Path for gnomAD's HGDP+1kGP metadata with updated population labels
metadata_path = 'gs://hgdp-1kg/tutorial_datasets/metadata_and_qc/gnomad_meta_updated.tsv'

# # Path for the number of structural variants (SV) counts per genome
# sv_counts_path = 'gs://hgdp-1kg/tutorial_datasets/metadata_and_qc/sv_counts_per_genome.tsv'

# Paths to related and unrelated Matrix Tables (without outliers) written out in Notebook 2: PCA and Ancestry Analyses
unrelateds_path = 'gs://hgdp-1kg/tutorial_datasets/pca_results/unrelateds_without_outliers.mt'
relateds_path = 'gs://hgdp-1kg/tutorial_datasets/pca_results/relateds_without_outliers.mt'

# Path to sample list for downsampled dataset
downsample_path = 'gs://hgdp-1kg/hgdp_tgp/qc_and_figure_generation/downsample_sample_list.txt'

# Path for final output table in tsv format
final_table_path = 'gs://hgdp-1kg/hgdp_tgp/qc_and_figure_generation/s4_unrel_heterozygosity.tsv'

# Read in the HGDP+1kGP pre-QC mt
pre_qc_mt = hl.read_matrix_table(pre_qc_path)

# Run QC
print("\nRunning QC\n")
mt = run_qc(pre_qc_mt)

# Remove PCA outliers from the dataset
mt_without_outliers = remove_pca_outliers(mt, outliers_path)

# should be 159339147 variants, 4094 samples (as of 12/06/2023)
# mt_without_outliers.count()

# Annotating on the SV counts per genome so it can be used downstream for analyses

# Read in table with sv counts
# sv_counts = hl.import_table(sv_counts_path, impute=True, key='sample', delimiter='\t')

# Annotating the sv counts onto the main mt
# mt_without_outliers = mt_without_outliers.annotate_cols(sv_counts=sv_counts[mt_without_outliers.s])
# print("\nAnnotating SV counts\n")

# Running sample QC on the dataset in order to get the latest metrics
mt_without_outliers = hl.sample_qc(mt_without_outliers)

# Grab the column fields of the Matrix Table
col_table = mt_without_outliers.cols()

# Write a col table with only the columns needed for table 1
col_table = col_table.select(col_table.meta_updated['hgdp_tgp_meta.Genetic.region'],
                             col_table.meta_updated.population,
                             col_table.sample_qc.n_deletion,
                             col_table.sample_qc.n_insertion)

print("\nSelected fields for col table\n")

# Validity check - there should be 4096 samples
# col_table.count()

# Annotating table with relatedness info

# Need to get number of unrelateds annotated to the table
# Reading in the metadata file with the updated population labels
metadata = hl.import_table(metadata_path, impute=True, key='s')

# Read in the unrelated and related matrix tables which were written out in Notebook 2: PCA and Ancestry Analyses
unrelateds = hl.read_matrix_table(unrelateds_path)
relateds = hl.read_matrix_table(relateds_path)
print("\nread in unrel/rel data\n")

# Annotating the mts with the metadata with updated population labels
unrelateds = unrelateds.annotate_cols(meta_updated=metadata[unrelateds.s])
relateds = relateds.annotate_cols(meta_updated=metadata[relateds.s])

# Annotate both the unrelated and the related tables with a flag named 'unrelated'
# Set the unrelated flag to 'True' for those in the unrelated dataset and 'False' for those in the related dataset
unrelateds = unrelateds.annotate_cols(unrelated=True)
relateds = relateds.annotate_cols(unrelated=False)

# Use hl.cols() to obtain two tables with only the columns from the unrelated and related mts
unrelateds_cols = unrelateds.cols()
relateds_cols = relateds.cols()
print("\nCreating col tables for unrel/rel\n")

# Validity check
# print(unrelateds_cols.count(), relateds_cols.count())  # 3378 unrelated and 718 related samples = 4096 total samples

unrelateds_count = unrelateds_cols.aggregate(hl.agg.counter(unrelateds_cols.meta_updated.population))
relateds_count = relateds_cols.aggregate(hl.agg.counter(relateds_cols.meta_updated.population))

# Validity check - print out the number of unrelated and related individuals per population
# print(f"\nNumber of unrelated individuals per population: \
# {unrelateds_count}\n\nNumber of related individuals per population: {relateds_count}\n")

# Join the columns of the unrelated and related datasets
mt_combined = unrelateds.union_cols(relateds)

# Validity check - count the number of unrelateds (True values) in the mt to make sure it is as expected
# 3378 True and 718 False
# print(f'Count for relateds/unrelateds pre-QC: {mt_combined.aggregate_cols(hl.agg.counter(mt_combined.unrelated))}')
# Create a table with only the columns from the mt containing related information
# This is done since the final output will be a tsv and thus must be in table format
# Being a table of columns allows it to be annotated onto the existing mt_col_table as shown below
mt_combined_col_table = mt_combined.cols()

# Annotate the relatedness information onto the column table
col_table = col_table.annotate(unrelated=mt_combined_col_table[col_table.s].unrelated)

# Validity check - count the number of unrelateds (True values) in the mt to make sure it is as expected (post-QC)
# 3376 True and 718 False
# print(f'\nCount for relateds/unrelated post-QC: {col_table.aggregate(hl.agg.counter(col_table.unrelated))}\n')

# filtering out related individuals
col_table = col_table.filter(col_table.unrelated, keep=True)

# Calculating statistics for manuscript text

# calculating the mean coverage for the whole project
# col_table.aggregate(hl.agg.stats(col_table.mean_coverage))

# Calculating coverage metrics for each project

mt_cols = mt_without_outliers.cols()

# The code below is not needed to run since the downsample dataset is already written out
# # Downsampling the dataset to 6 individuals per population
#
# # grouping the table by population and genetic region, and then taking 6 samples from each population
# downsample_ht = col_table.group_by(
#     col_table['hgdp_tgp_meta.Genetic.region'], col_table.population).aggregate(
#     s=hl.agg.take(col_table.s, 4))
# print("\ndownsampled the dataset\n")
#
# # checking count, should be equal to num pops (80)
# # downsample_ht.count()
#
# # Exploding the arrays of six individuals per population so there are 468 samples in total
# ht_explode = downsample_ht.explode('s')
#
# # creating a list of the samples from the downsampled individual list
# downsample_list = ht_explode.s.collect()

with hl.utils.hadoop_open(downsample_path) as f:
    for line in f:
        ds_samples = line.split(',')

# Using hl.literal here to convert the list from a python object to a hail expression
# this is so that it can be used to filter out samples
downsample_list = hl.literal(ds_samples)

# Filtering the original matrix table so that there are only 6 individuals per population
downsample_mt = mt_without_outliers.filter_cols(downsample_list.contains(mt_without_outliers['s']))
print("\nfiltering qc'd mt to only downsampled indivs\n")

# Only keep the variants which are found in the samples that are left
downsample_mt = downsample_mt.filter_rows(hl.agg.any(downsample_mt.GT.is_non_ref()))

# # Checking the number of samples in the mt after subsetting to 6/indiv. should be 80*6=480
# downsample_mt.count()

# # Running sample QC on the dataset to get newly calculate n_singleton
# downsample_mt = hl.sample_qc(downsample_mt)
#
# Running sample QC on the dataset to get newly calculate AC
downsample_mt = hl.variant_qc(downsample_mt)

print("\nfiltering to AC > 1\n")
# Removing variants with AC < 1
downsample_mt = downsample_mt.filter_rows((hl.min(downsample_mt.variant_qc.AC) < 1),
                                          keep=False)

mt_k_test = downsample_mt.annotate_rows(pop_callstats=hl.agg.group_by(
    downsample_mt.meta_updated.population, hl.agg.call_stats(downsample_mt.GT, downsample_mt.alleles)))

mt = downsample_mt.group_cols_by(
    downsample_mt.hgdp_tgp_meta.genetic_region, downsample_mt.meta_updated.population
).aggregate(
    af=hl.agg.call_stats(downsample_mt.GT, downsample_mt.alleles).AF
)
# we have a variants by pop matrix now
mt = mt.annotate_entries(
    p=mt.af[0],
    q=mt.af[1],
)
mt = mt.annotate_entries(heterozygosity_by_pop=1 - (mt.p ** 2) - (mt.q ** 2))

# mt.heterozygosity_by_pop.show()

mt = mt.annotate_cols(heterozygosity_stats_by_pop=hl.agg.stats(mt.heterozygosity_by_pop))

output_het_ht = mt.cols()
output_het_ht = output_het_ht.flatten()

print("\nwriting out final table\n")
output_het_ht.export(final_table_path, header=True)
