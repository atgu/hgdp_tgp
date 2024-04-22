# script for running PCA (a python script version of tutorial nb2) - running in Jupyter nb was taking too long/nb crashing   

# author = Mary T. Yohannes

import hail as hl

# Functions from gnomAD library to apply genotype filters and project related samples  
from gnomad.utils.filtering import filter_to_adj
from gnomad.sample_qc.ancestry import pc_project

# Path for HGDP+1kGP dataset prior to applying gnomAD QC filters
pre_qc_path = 'gs://gcp-public-data--gnomad/release/3.1.2/mt/genomes/gnomad.genomes.v3.1.2.hgdp_1kg_subset_dense.mt'

# Path for gnomAD's HGDP+1kGP metadata for plotting 
metadata_path = 'gs://hgdp-1kg/tutorial_datasets/metadata_and_qc/gnomad_meta_updated.tsv'

# Save the filtered and LD pruned mt as an intermediate file since LD pruning takes a while to rerun
ld_pruned_path = 'gs://hgdp-1kg/tutorial_datasets/pca_preprocessing/ld_pruned.mt'

# Hail table of related sample IDs for separating unrelateds and relateds for PCA run 
related_sample_ids_path = 'gs://hgdp-1kg/tutorial_datasets/pca_preprocessing/related_sample_ids.ht'

# Path for with-outliers PCA results - global & subcontinental PCA 
pc_scores_with_outliers_path = 'gs://hgdp-1kg/tutorial_datasets/pca/pc_scores_with_outliers/'

# PCA outliers file 
outliers_path = 'gs://hgdp-1kg/tutorial_datasets/pca/pca_outliers.txt'

# Path for without-outliers PCA results - global & subcontinental PCA 
pc_scores_without_outliers_path = 'gs://hgdp-1kg/tutorial_datasets/pca/pc_scores_without_outliers/'

# Paths for unrelated and related datasets without outliers   
unrelateds_mt_without_outliers_path = 'gs://hgdp-1kg/tutorial_datasets/pca_results/unrelateds_without_outliers.mt'
relateds_mt_without_outliers_path = 'gs://hgdp-1kg/tutorial_datasets/pca_results/relateds_without_outliers.mt' 


# these are intermediate files that helped with running the script faster - not included in the main list of files 
unrelateds_mt_preoutlier = hl.read_matrix_table('gs://hgdp-1kg/tutorial_datasets/pca/temp_files/unrel_preout.mt')
relateds_mt_preoutlier = hl.read_matrix_table('gs://hgdp-1kg/tutorial_datasets/pca/temp_files/rel_preout.mt')

# Read in gnomAD's HGDP+1kGP metadata for plot annotation
metadata = hl.import_table(metadata_path, impute = True, key = 's')

def run_pca(mt: hl.MatrixTable):
    """
    Runs PCA on a dataset
    :param mt: dataset to run PCA on
    :return: loadings and pc scores of unrelated samples 
    """
    pca_evals, pca_scores, pca_loadings = hl.hwe_normalized_pca(mt.GT, k=20, compute_loadings=True)
    pca_mt = mt.annotate_rows(pca_af=hl.agg.mean(mt.GT.n_alt_alleles()) / 2)
    pca_loadings = pca_loadings.annotate(pca_af=pca_mt.rows()[pca_loadings.key].pca_af)
    pca_scores = pca_scores.transmute(**{f'PC{i}': pca_scores.scores[i - 1] for i in range(1, 21)})
    
    return pca_loadings, pca_scores 

def project_individuals(project_mt, pca_loadings, unrel_scores, out_path: str, reg_name:str, outlier_status:str):
    """
    Project samples into predefined PCA space
    :param project_mt: matrix table of related samples to project 
    :param pca_loadings: existing PCA space of unrelated samples 
    :param unrel_scores: unrelated samples' PC scores
    :param out_path: path for where to save PCA projection outputs
    :param reg_name: region name for saving output purposes
    :param outlier_status: is the dataset with or without outliers? 
    """
    ht_projections = pc_project(project_mt, pca_loadings)  
    ht_projections = ht_projections.transmute(**{f'PC{i}': ht_projections.scores[i - 1] for i in range(1, 21)}) 
    scores = unrel_scores.union(ht_projections) # combine the pc scores from both the unrelateds and relateds 
    scores.export(out_path + reg_name + '_scores_' + outlier_status + '.txt.bgz') # write output for plotting    

# Calculate PC scores 

# Dictionaries to hold unrelateds' PCA loadings and scores
loadings_dict = {}
unrel_scores_dict = {}

# Run PCA on unrelated samples as a whole
loadings_dict['GLOBAL'], unrel_scores_dict['GLOBAL'] = run_pca(unrelateds_mt_preoutlier)  

# Project related samples onto unrelated-samples' PC space 
project_individuals(relateds_mt_preoutlier, loadings_dict['GLOBAL'], unrel_scores_dict['GLOBAL'], pc_scores_with_outliers_path, 'GLOBAL', 'with_outliers')


# Calculate PC scores 

# Dictionaries to hold unrelateds' PCA loadings and scores
loadings_dict = {}
unrel_scores_dict = {}
regions = unrelateds_mt_preoutlier['hgdp_tgp_meta']['genetic_region'].collect() 
regions = list(dict.fromkeys(regions)) # convert into a list
# There are 7 regions: EUR, AFR, AMR, EAS, CSA, OCE, and MID

# For each region, run PCA on the unrelated samples 
for i in regions:  
    if i is not None: # exclude a none value
        # Filter the unrelateds per region
        subcont_unrelateds = unrelateds_mt_preoutlier.filter_cols(unrelateds_mt_preoutlier['hgdp_tgp_meta']['genetic_region'] == i) 

        # Run PCA
        loadings_dict[i], unrel_scores_dict[i] = run_pca(subcont_unrelateds)

        # Filter the related mt per region 
        subcont_relateds = relateds_mt_preoutlier.filter_cols(relateds_mt_preoutlier['hgdp_tgp_meta']['genetic_region'] == i)  

        # Project related samples onto unrelated-samples' PC space 
        project_individuals(subcont_relateds, loadings_dict[i], unrel_scores_dict[i], pc_scores_with_outliers_path, i, 'with_outliers')

