import hail as hl
import re
from tabulate import tabulate

hl.init()



# path for sample metadata file
sample_metadata_path = 'gs://african-seq-data/hgdp_tgp/gnomad_meta_v1.tsv'

# path for variant qc info
var_metadata_path = 'gs://gcp-public-data--gnomad/release/3.1/ht/genomes/gnomad.genomes.v3.1.sites.ht'

# path for Konrad's densified matrix table
dense_mt_path = 'gs://hgdp_tgp/output/tgp_hgdp.mt'

# reading in sample qc metadata file
sample_meta = hl.import_table(sample_metadata_path, impute=True)

# reading in variant qc information
var_meta = hl.read_table(var_metadata_path)

# reading in densified matrix table
dense_mt = hl.read_matrix_table(dense_mt_path)


# Unflattening the Sample Metadata File
d = {}

row = sample_meta.row_value

for name in row:
    def recur(dict_ref, split_name):
        if (len(split_name) == 1):
            dict_ref[split_name[0]] = row[name]
            return
        existing = dict_ref.get(split_name[0])
        if existing is not None:
            assert isinstance(existing, dict), existing  # fails on foo.bar and foo.bar.baz
            recur(existing, split_name[1:])
        else:
            existing = {}
            dict_ref[split_name[0]] = existing
            recur(existing, split_name[1:])
    recur(d, name.split('.'))
def dict_to_struct(d):
    fields = {}
    for k, v in d.items():
        if isinstance(v, dict):
            v = dict_to_struct(v)
        fields[k] = v
    return hl.struct(**fields)

sample_meta = sample_meta.select(**dict_to_struct(d))
sample_meta = sample_meta.key_by('s')


# Correcting the field types from the flattened hail table
sample_meta = sample_meta.annotate(project_meta=sample_meta.project_meta.annotate(
    v2_age=hl.float64(sample_meta.project_meta.v2_age),
    v2_related=hl.bool(sample_meta.project_meta.v2_related),
    v2_internal=hl.bool(sample_meta.project_meta.v2_internal),
    v2_neuro=hl.bool(sample_meta.project_meta.v2_neuro),
    v2_control=hl.bool(sample_meta.project_meta.v2_control),
    v2_topmed=hl.bool(sample_meta.project_meta.v2_topmed),
    v2_high_quality=hl.bool(sample_meta.project_meta.v2_high_quality),
    v2_pcr_free=hl.bool(sample_meta.project_meta.v2_pcr_free),
    v2_release_2_0_2=hl.bool(sample_meta.project_meta.v2_release_2_0_2),
    age=hl.int32(sample_meta.project_meta.age),
    age_alt=hl.int32(sample_meta.project_meta.age_alt)),
    bam_metrics=sample_meta.bam_metrics.annotate(
    median_coverage=hl.float64(sample_meta.bam_metrics.median_coverage),
    median_insert_size=hl.float64(sample_meta.bam_metrics.median_insert_size)),
    sex_imputation=sample_meta.sex_imputation.annotate(
    chr20_mean_dp=hl.float32(sample_meta.sex_imputation.chr20_mean_dp),
    chrX_mean_dp=hl.float32(sample_meta.sex_imputation.chrX_mean_dp),
    chrY_mean_dp=hl.float32(sample_meta.sex_imputation.chrY_mean_dp),
    chrX_ploidy=hl.float32(sample_meta.sex_imputation.chrX_ploidy),
    chrY_ploidy=hl.float32(sample_meta.sex_imputation.chrY_ploidy),
    impute_sex_stats = sample_meta.sex_imputation.impute_sex_stats.annotate(
        n_called=hl.int64(sample_meta.sex_imputation.impute_sex_stats.n_called),
        expected_homs=hl.float64(sample_meta.sex_imputation.impute_sex_stats.expected_homs),
        observed_homs=hl.int64(sample_meta.sex_imputation.impute_sex_stats.observed_homs))),
    sample_qc=sample_meta.sample_qc.annotate(
    n_hom_ref=hl.int64(sample_meta.sample_qc.n_hom_ref),
    n_het=hl.int64(sample_meta.sample_qc.n_het),
    n_hom_var=hl.int64(sample_meta.sample_qc.n_hom_var),
    n_non_ref=hl.int64(sample_meta.sample_qc.n_non_ref),
    n_singleton=hl.int64(sample_meta.sample_qc.n_singleton),
    n_snp=hl.int64(sample_meta.sample_qc.n_snp),
    n_insertion=hl.int64(sample_meta.sample_qc.n_insertion),
    n_deletion=hl.int64(sample_meta.sample_qc.n_deletion),
    n_transition=hl.int64(sample_meta.sample_qc.n_transition),
    n_transversion=hl.int64(sample_meta.sample_qc.n_transversion),
    n_star=hl.int64(sample_meta.sample_qc.n_star)),
    population_inference=sample_meta.population_inference.annotate(
    pca_scores=hl.array(hl.parse_json(sample_meta.population_inference.pca_scores, hl.tarray(hl.tfloat64)))),
    sample_filters=sample_meta.sample_filters.annotate(
    hard_filters=hl.set(hl.parse_json(sample_meta.sample_filters.hard_filters, hl.tset(hl.tstr))),
    qc_metrics_filters=hl.set(hl.parse_json(sample_meta.sample_filters.qc_metrics_filters, hl.tset(hl.tstr)))),
    relatedness_inference=sample_meta.relatedness_inference.annotate(
    relationships=hl.set(hl.parse_json(sample_meta.relatedness_inference.relationships, hl.tset(hl.tstr)))))



# Merging the metadata table and dense mt
ht = sample_meta

mt = dense_mt.annotate_cols(**ht[dense_mt.s])
all_sample_filters = set(ht['sample_filters'])
# bad_sample_filters were removing whole populations (mostly AFR and OCE) that passed all other QC, bad filters
bad_sample_filters = {re.sub('fail_', '', x) for x in all_sample_filters if x.startswith('fail_')}
# this filters to only variants that passed all gnomad QC or only failed filters in bad_sample_filters
mt_filt = mt.filter_cols((mt['sample_filters']['qc_metrics_filters'].difference(bad_sample_filters).length() == 0))


# Merging the variant QC info with the dataset
mt = mt.annotate_rows(**var_meta[mt.locus, mt.alleles])

mt.write('gs://african-seq-data/hgdp_tgp/hgdp_tgp_dense_meta_filt.mt')


# Creating list with table information
table = []
pop_tmp = set(pop)
for superpop in pop_tmp:
    subpop_tmp = mt_filt.filter_cols(mt_filt['project_meta'].project_pop == superpop).project_meta.project_subpop.collect()
    for subpop in set(subpop_tmp):
        indiv_count = mt_filt.aggregate_cols(hl.agg.counter(mt_filt['project_meta'].project_subpop == subpop))
        cov = mt_filt.filter_cols(mt_filt)
        table.append([superpop,subpop,indiv_count[1]])



print(tabulate(table, headers =["Super Population", "Population", "# Individuals"]))