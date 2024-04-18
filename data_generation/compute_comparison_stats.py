import hail as hl
# Function from gnomAD library to apply genotype filters
from gnomad.utils.filtering import filter_to_adj

hl.init(tmp_dir='gs://hgdp-1kg/znk-temporary-files/')


def group_hist(hgdp_1kg, comparison, gnomad_bool):
    if not gnomad_bool:
        # Removing missing values from the dataset
        hgdp_1kg = hgdp_1kg.filter(hl.is_missing(hgdp_1kg.maf_hgdp_1kg), keep=False)
        # Annotating a field with bool of whether each var is in comparison dataset
        hgdp_1kg = hgdp_1kg.annotate(in_comparison=hl.is_defined(comparison[hgdp_1kg.locus, hgdp_1kg.alleles]))
        # Summary will be returned, grouped table with maf bins
        summary = hgdp_1kg.group_by(
            cat=hl.case()
            .when(hgdp_1kg.maf_hgdp_1kg < 1e-3, "0.01%-0.1%")
            .when(hgdp_1kg.maf_hgdp_1kg < 1e-2, "0.1%-1%")
            .when(hgdp_1kg.maf_hgdp_1kg < 1e-1, "1.0%-10%")
            .default("10-50%")
        ).aggregate(
            # number of variants in comparison
            n_var=hl.agg.count(),
            # number of variants in hgdp+1kgp also in comparison
            n_var_in_both=hl.agg.count_where(hgdp_1kg.in_comparison == True))

        comparison = comparison.annotate(not_in_hgdp_1kg=hl.is_missing(hgdp_1kg[comparison.locus, comparison.alleles]))
        num = comparison.aggregate(hl.agg.count_where(comparison.not_in_hgdp_1kg == True))

        summary = summary.union(
            hl.utils.range_table(1).key_by(cat="0%").select(n_var=num, n_var_in_both=0), unify=True)

    if gnomad_bool:
        # Annotating a field with bool of whether each var is in comparison dataset
        # go back to check if maf >0 in hgdp+1kg
        hgdp_1kg = hgdp_1kg.annotate(in_comparison=hgdp_1kg.gnomad_only_maf > 0)
        summary = hgdp_1kg.group_by(
            cat=hl.case()
            .when(hl.is_missing(hgdp_1kg.hdgp_1kg_only_maf), "missing")
            .when(hgdp_1kg.hdgp_1kg_only_maf < 1e-3, "0.01%-0.1%")
            .when(hgdp_1kg.hdgp_1kg_only_maf < 1e-2, "0.1%-1%")
            .when(hgdp_1kg.hdgp_1kg_only_maf < 1e-1, "1.0%-10%")
            .default("10-50%")
        ).aggregate(
            # number of variants in comparison
            n_var=hl.agg.count(),
            # number of variants in hgdp+1kgp also in comparison
            n_var_in_both=hl.agg.count_where(hgdp_1kg.in_comparison == True))

        comparison = comparison.annotate(not_in_hgdp_1kg=hl.is_missing(hgdp_1kg[comparison.locus, comparison.alleles]))
        num = comparison.aggregate(hl.agg.count_where(comparison.not_in_hgdp_1kg))
        # print(res)

        # print(hgdp_1kg.aggregate(hl.agg.counter(hl.is_defined(hgdp_1kg.gnomad_only_maf))))
        # print(hgdp_1kg.aggregate(hl.agg.counter(hgdp_1kg.gnomad_only_maf > 0)))

        summary = summary.union(
            hl.utils.range_table(1).key_by(cat="0%").select(n_var=num, n_var_in_both=0), unify=True)

    return summary


# Given a hail matrix table, and fieldname which is the struct within the AF array is located
# creates a new field which contains the minor allele frequency
# will create a hail table which only contains the row keys and the the maf field created
def get_maf_ht(mt, fieldname):
    return mt.select_rows(maf=hl.min(mt[fieldname].AF)).rows()


def run_qc(mt, metadata_path):
    # Set up function to:
    # apply gnomAD's sample, variant and genotype QC filters
    # remove two contaminated samples identified using CHARR - https://pubmed.ncbi.nlm.nih.gov/37425834/
    # remove the gnomAD sample that's added for QC purposes
    # add gnomAD's HGDP+1kGP metadata with the updated population labels as a column field

    # Apply sample QC filters to dataset
    # This filters to only samples that passed gnomAD's sample QC hard filters
    mt = mt.filter_cols(~mt.gnomad_sample_filters.hard_filtered)  # removed 31 samples

    # Apply variant QC filters to dataset
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

    # Apply genotype QC filters to the dataset
    # This is done using a function imported from gnomAD and is the last step in the QC process
    mt = filter_to_adj(mt)

    # removing X and Y chromosome data
    contigs_to_keep = []
    for i in range(1, 23):
        contigs_to_keep.append(f'chr{i}')
    mt = mt.filter_rows(hl.set(contigs_to_keep).contains(mt.locus.contig))

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


def filter_comparison(data, data_type):
    # Filter the comparison datasets using the filter fields on the respective mts or hts
    if data_type == "mt":
        # removing variants that failed the QC based on imputed filters column
        data = data.filter_rows(hl.len(data.filters) > 0, keep=False)
        # removing X and Y chromosome data
        contigs_to_keep = []
        for i in range(1, 23):
            contigs_to_keep.append(f'chr{i}')
        data = data.filter_rows(hl.set(contigs_to_keep).contains(data.locus.contig))
        # data = data.filter_rows(data.locus.contig != 'X')
        # data = data.filter_rows(data.locus.contig != 'Y')

    elif data_type == "ht":
        # removing variants that failed the QC based on imputed filters column
        data = data.filter(hl.len(data.filters) > 0, keep=False)
        # removing X and Y chromosome data
        contigs_to_keep = []
        for i in range(1, 23):
            contigs_to_keep.append(f'chr{i}')
        data = data.filter(hl.set(contigs_to_keep).contains(data.locus.contig))
        # data = data.filter(data.locus.contig != 'X')
        # data = data.filter(data.locus.contig != 'Y')

    else:
        print("please input 'ht' or 'mt' for data_type")

    return data


# Reading in Datasets
# Pathways for comparison datasets
# Setting variables for paths makes it easier to update paths in the future
phase3_1kg_path = "gs://hgdp-1kg/hgdp_tgp/comparison_data/ALL.chr*.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz"
nygc_1kg_path = "gs://hgdp-1kg/hgdp_tgp/comparison_data/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr{[1-9],[1-9][0-9]}.recalibrated_variants.vcf.gz"
berg_path = "gs://hgdp-1kg/hgdp_tgp/comparison_data/hgdp_wgs.20190516.full.chr{[1-9],[1-9][0-9]}.vcf.gz"
# gnomadv3_path = "gs://gcp-public-data--gnomad/release/3.1.2/ht/genomes/gnomad.genomes.v3.1.2.sites.ht"
hgdp_1kg_preQC_path = 'gs://gcp-public-data--gnomad/release/3.1.2/mt/genomes/gnomad.genomes.v3.1.2.hgdp_1kg_subset_dense.mt'
# Path for PCA outliers list
outliers_path = 'gs://hgdp-1kg/tutorial_datasets/pca/pca_outliers.txt'
# Path for gnomAD's HGDP+1kGP metadata with updated population labels
metadata_path = 'gs://hgdp-1kg/tutorial_datasets/metadata_and_qc/gnomad_meta_updated.tsv'

print("\nreading in datasets\n")
# Reading in comparison datasets
# use this to read in the list of vcfs
# mt_list = [hl.import_vcf(mt_path,force_bgz = True) for mt_path in mt_paths]
print('phase3')
phase3_1kg_mt = hl.import_vcf(phase3_1kg_path, force_bgz=True, reference_genome='GRCh38')
print('nygc')
nygc_1kg_mt = hl.import_vcf(nygc_1kg_path, reference_genome='GRCh38', force_bgz=True)
print('berg')
berg_mt = hl.import_vcf(berg_path, force_bgz=True, reference_genome='GRCh38')
# gnomadv3_ht = hl.read_table(gnomadv3_path)
print('hgdp_1kg')
hgdp_1kg_preQC_mt = hl.read_matrix_table(hgdp_1kg_preQC_path)


print("running qc and removing pca outliers\n")
# Run QC on the HGDP+1kGP dataset
hgdp_1kg_mt = run_qc(hgdp_1kg_preQC_mt, metadata_path)

# Remove PCA outliers from the dataset
mt_without_outliers = remove_pca_outliers(hgdp_1kg_mt, outliers_path)


print("splitting multi allelics\n")
# Splitting Multi-allelics
# will want to use hl.split_multi_hts
# splitting multialleleics for nygc 1kgp
nygc_split_mt = hl.split_multi_hts(nygc_1kg_mt)
# getting sample/var count for nygc after split multi
# nygc_split_mt.count()
# splitting multiallelics for bergstrom
berg_split_mt = hl.split_multi_hts(berg_mt, permit_shuffle=True)
# getting sample/var count after bergstrom split multi
# berg_split_mt.count()
# resetting the names to the ones that are used in the downstream steps of this script
nygc_1kg_mt = nygc_split_mt
berg_mt = berg_split_mt


print('filtering comparison datasets and removing X&Y chr\n')
# filtering comparison dataset variants and removing X and Y chromosome
print('nygc\n')
nygc_filt = filter_comparison(nygc_1kg_mt, 'mt')
print('phase3\n')
phase3_filt = filter_comparison(phase3_1kg_mt, 'mt')
print('berg\n')
berg_filt = filter_comparison(berg_mt, 'mt')
# print('gnomad\n')
# gnomad_filt = filter_comparison(gnomadv3_ht, 'ht')


print('getting maf tables\n')
# Running hl.variant_qc() on hgdp+1kgp dataset to get the AF
hgdp_1kg_mt_var = hl.variant_qc(mt_without_outliers)
# Creating maf tables for all the datasets
# done using get_maf_ht which takes the min of the allele frequency array in a dataset
# given the name of the field which contains and array with the allele frequencies
hgdp_1kg_maf = get_maf_ht(hgdp_1kg_mt_var, 'variant_qc')
hgdp_1kg_maf = hgdp_1kg_maf.rename({'maf': 'maf_hgdp_1kg'})

nygc_maf = get_maf_ht(nygc_filt, 'info')
nygc_maf = nygc_maf.rename({'maf': 'nygc_maf'})

berg_maf = get_maf_ht(berg_filt, 'info')
berg_maf = berg_maf.rename({'maf': 'berg_maf'})

phase3_maf = get_maf_ht(phase3_filt, 'info')
phase3_maf = phase3_maf.rename({'maf': 'phase3_maf'})

# # Getting the gnomADv3 maf table Cannot use the get_maf_ht function since it is already a table and the format is a
# # bit different Get a separate matrix table with only the gnomAD AF metrics needed in order to write out the gnomAD
# # comparison histogram adding a col with a bool for if samples are ~gnomad_high_quality total of True in this col is
# # equal to the number of samples which are in HGDP+1kGP but not in the gnomAD dataset
# # Use hl.hadoop_open to read in the PCA outliers file into Hail from Google Cloud Storage
# with hl.utils.hadoop_open(outliers_path) as file:
#     outliers = [line.rstrip('\n') for line in file]
# # Use hl.literal to convert the outliers list from a python object to a Hail expression so that it can be used to
# # filter out samples
# outlier_list = hl.literal(outliers)
#
# mt_pca_anno = hgdp_1kg_mt.annotate_cols(is_pca_outlier=outlier_list.contains(hgdp_1kg_mt['s']))
# mt = mt_pca_anno.annotate_cols(not_in_gnomad=~mt_pca_anno.gnomad_release)
# # calculating call stats for the whole hgdp_1kg matrix table as well as only for samples which were not in gnomad
# mt = mt.annotate_rows(hgdp_tgp_stats=hl.agg.call_stats(mt.GT, mt.alleles),
#                       not_in_gnomad_stats=hl.agg.filter(mt.not_in_gnomad == True,
#                                                         hl.agg.call_stats(mt.GT, mt.alleles)),
#                       pca_outlier_stats=hl.agg.filter(mt.is_pca_outlier, hl.agg.call_stats(mt.GT, mt.alleles)),
#                       non_pca_outlier_stats=hl.agg.filter(~mt.is_pca_outlier, hl.agg.call_stats(mt.GT, mt.alleles))
#                       )
# # potential way to filter out the non-pca outlier stats - have not run this yet
# tmp_ht = mt.filter_rows(mt.non_pca_outlier_stats.AF[1] > 0).rows()
#
# # print('writing out gnomad agg stats table\n')
# # mt.rows().write('gs://hgdp-1kg/hgdp_tgp/qc_and_figure_generation/comparison_hists/gnomad_agg_stats_v3.ht')
#
# print('calculating gnomad specific metrics\n')
# ht = hl.read_table('gs://hgdp-1kg/hgdp_tgp/qc_and_figure_generation/comparison_hists/gnomad_agg_stats_v3.ht')
# # this will give the frequency just for the samples which are only in gnomad (not in 1kgp or HGDP)
# ht = ht.annotate(gnomad_only_AC=ht.gnomad_freq[0].AC - (ht.hgdp_tgp_stats.AC - ht.not_in_gnomad_stats.AC))
# # adding a field with the calculated AN for gnomAD only
# ht = ht.annotate(gnomad_only_AN=ht.gnomad_freq[0].AN - (ht.hgdp_tgp_stats.AN - ht.not_in_gnomad_stats.AN))
# # Calculating the maf for variants in gnomAD only
# ht = ht.annotate(gnomad_only_maf=0.5 - hl.abs(
#     0.5 - (ht.gnomad_only_AC[1] / ht.gnomad_only_AN)))
# # calculating the maf for variants in HGDP+1kGP only
# ht = ht.annotate(hdgp_1kg_only_maf=hl.min(
#     (ht.not_in_gnomad_stats.AC - ht.pca_outlier_stats.AC) / (ht.not_in_gnomad_stats.AN - ht.pca_outlier_stats.AN)
# ))
# ht = ht.select_globals()
# # trying to create a ht with only the metrics needed for gnomAD comparison
# hgdp_1kg_gnomad_ht = ht.select('gnomad_only_maf',
#                                'hdgp_1kg_only_maf',
#                                'gnomad_only_AC')


# Changing the number of partitions for each of the datasets to speed up downstream steps
nygc_maf = nygc_maf.naive_coalesce(5000)
berg_maf = berg_maf.naive_coalesce(5000)
phase3_maf = phase3_maf.naive_coalesce(5000)
hgdp_1kg_maf = hgdp_1kg_maf.naive_coalesce(5000)
# gnomad_filt = gnomad_filt.naive_coalesce(5000)
# hgdp_1kg_gnomad_ht = hgdp_1kg_gnomad_ht.naive_coalesce(5000)

print('creating grouped histograms for all comparison datasets\n')
# # Creating a table with the aggregated values for gnomAD
# gnomAD_hist = group_hist(hgdp_1kg_gnomad_ht, gnomad_filt, True).persist()
# # gnomAD_hist.show()
# print('writing gnomad\n')
# gnomAD_hist.export('gs://hgdp-1kg/hgdp_tgp/qc_and_figure_generation/comparison_hists/gnomAD_hist_v1.tsv')

# Creating a table with the aggregated values - testing with bergstrom data 1st
bergstrom_hist = group_hist(hgdp_1kg_maf, berg_maf, False).persist()
# bergstrom_hist.show()
print('writing bergstrom\n')
bergstrom_hist.export('gs://hgdp-1kg/hgdp_tgp/qc_and_figure_generation/comparison_hists/berg_hist.tsv')

# Creating a table with the aggregated values for NYGC
nygc_hist = group_hist(hgdp_1kg_maf, nygc_maf, False).persist()
# nygc_hist.show()
print('writing nygc\n')
nygc_hist.export('gs://hgdp-1kg/hgdp_tgp/qc_and_figure_generation/comparison_hists/nygc_hist.tsv')

# Creating a table with the aggregated values for phase3 1kGP
phase3_hist = group_hist(hgdp_1kg_maf, phase3_maf, False).persist()
# phase3_hist.show()
print('writing phase3\n')
phase3_hist.export('gs://hgdp-1kg/hgdp_tgp/qc_and_figure_generation/comparison_hists/phase3_hist.tsv')
