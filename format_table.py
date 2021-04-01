import hail as hl
import re

hl.init()


# Prints out the total sample/variant count for an mt
def printCount(mt):
    n = mt.count()
    print('Number of Samples: {}\nNumber of Variants: {}'.format(n[1],n[0]))


# path for Alicia's sample metadata file
sample_metadata_path = 'gs://african-seq-data/hgdp_tgp/gnomad_meta_v1.tsv'

# path for Julia's sample metadata file
jul_metadata_path = ('gs://hgdp_tgp/output/gnomad_v3.1_sample_qc_metadata_hgdp_tgp_subset.ht')

# path for variant qc info
var_metadata_path = 'gs://gcp-public-data--gnomad/release/3.1/ht/genomes/gnomad.genomes.v3.1.sites.ht'

# path for Konrad's densified matrix table
dense_mt_path = 'gs://hgdp_tgp/output/tgp_hgdp.mt'


# reading in Alicia's sample metadata file
sample_meta = hl.import_table(sample_metadata_path, impute=True)

# reading in Julia's sample metadata file
jul_meta = hl.read_table(jul_metadata_path)

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


new_meta = sample_meta.select(sample_meta.hgdp_tgp_meta, sample_meta.bergstrom)


ht = jul_meta.annotate(**new_meta[jul_meta.s])


ht.write('gs://african-seq-data/hgdp_tgp/hgdp_tgp_sample_metadata.ht')

ht = hl.read_table('gs://african-seq-data/hgdp_tgp/hgdp_tgp_sample_metadata.ht')


mt = dense_mt.annotate_cols(**ht[dense_mt.s])
all_sample_filters = set(ht['sample_filters'])
# bad_sample_filters were removing whole populations (mostly AFR and OCE) that passed all other QC, bad filters
bad_sample_filters = {re.sub('fail_', '', x) for x in all_sample_filters if x.startswith('fail_')}
# this filters to only variants that passed all gnomad QC or only failed filters in bad_sample_filters
mt_filt = mt.filter_cols((mt['sample_filters']['qc_metrics_filters'].difference(bad_sample_filters).length() == 0))

# Merging the variant QC info with the dataset
mt = mt.annotate_rows(**var_meta[mt.locus, mt.alleles])

mt.write('gs://african-seq-data/hgdp_tgp/hgdp_tgp_dense_meta_filt.mt')


# Reading in annotated dataset
mt = hl.read_matrix_table('gs://african-seq-data/hgdp_tgp/hgdp_tgp_dense_meta_filt.mt')


# Removing Outlier Samples
outlier_set = {'NA20314','NA20299','HG01880','HG01881',
               'HGDP00130','HGDP00013','HGDP00150',
               'HGDP00029','HGDP01298','HGDP01303',
               'LP6005443-DNA_B02','HGDP01300','HG01628',
               'HG01629','HG01630','HG01694','HG01696',
               'HGDP00621','HGDP01270','HGDP01271'}
set_to_remove = hl.literal(outlier_set)
mt = mt.filter_cols(~set_to_remove.contains(mt['s']))


# Generating a VCF
# changing types to float64 so they can be written out to vcf
mt_vcf = mt.annotate_rows(info = mt.info.annotate(
    QUALapprox=hl.float64(mt.info.QUALapprox),
    SB= mt.info.SB.map(lambda x: hl.float64(x)),
    MQ=hl.float64(mt.info.MQ),
    MQRankSum=hl.float64(mt.info.MQRankSum),
    VarDP=hl.float64(mt.info.VarDP),
    AS_ReadPosRankSum=hl.float64(mt.info.AS_ReadPosRankSum), 
    AS_pab_max=hl.float64(mt.info.AS_pab_max), 
    AS_QD=hl.float64(mt.info.AS_QD), 
    AS_MQ=hl.float64(mt.info.AS_MQ), 
    QD=hl.float64(mt.info.QD), 
    AS_MQRankSum=hl.float64(mt.info.AS_MQRankSum), 
    FS=hl.float64(mt.info.FS), 
    AS_FS=hl.float64(mt.info.AS_FS), 
    ReadPosRankSum=hl.float64(mt.info.ReadPosRankSum), 
    AS_QUALapprox=hl.float64(mt.info.AS_QUALapprox), 
    AS_SB_TABLE=mt.info.AS_SB_TABLE.map(lambda x: hl.float64(x)), 
    AS_VarDP=hl.float64(mt.info.AS_VarDP), 
    AS_SOR=hl.float64(mt.info.AS_SOR), 
    SOR=hl.float64(mt.info.SOR), 
    singleton=hl.bool(mt.info.singleton), 
    transmitted_singleton=hl.bool(mt.info.transmitted_singleton), 
    omni=hl.bool(mt.info.omni), 
    mills=hl.bool(mt.info.omni), 
    monoallelic=hl.bool(mt.info.monoallelic), 
    AS_VQSLOD=hl.float64(mt.info.AS_VQSLOD), 
    InbreedingCoeff=hl.float64(mt.info.InbreedingCoeff)
))

mt_vcf = mt_vcf.drop('gvcf_info')

hl.export_vcf(mt_vcf, 'gs://african-seq-data/hgdp_tgp/hgdp_tgp_dense_filt.vcf.bgz', parallel='separate_header')


# Creating a separate dataset with only the columns needed for table X
mt = hl.sample_qc(mt, name="new_sample_qc")

mt_table = mt.select_cols(mt.new_sample_qc.n_snp, 
                          mt.bam_metrics.mean_coverage, 
                          mt.new_sample_qc.n_singleton, 
                          mt.hgdp_tgp_meta.Population,
                          mt.hgdp_tgp_meta.Genetic.region
                         )


# Grabbing only the columns from the matrix table (outputs table of just columns)
col_table = mt_table.cols()

col_table.checkpoint('gs://african-seq-data/hgdp_tgp/table_x_checkpoint.ht')


# Reading in checkpoint col table
col_table = hl.read_table('gs://african-seq-data/hgdp_tgp/table_x_checkpoint.ht')


# Creating three separate tables for each column, and calculating stats for that column grouped by genetic region (pop)
n_snp = col_table.group_by(
    col_table.region, col_table.Population).aggregate(
    n_snp_stats = hl.agg.stats(col_table.n_snp))

n_singleton = col_table.group_by(
    col_table.region, col_table.Population).aggregate(
    n_singleton_stats = hl.agg.stats(col_table.n_singleton))

mean_coverage = col_table.group_by(
    col_table.region, col_table.Population).aggregate(
    cov_stats = hl.agg.stats(col_table.mean_coverage))

mean_coverage = mean_coverage.flatten()

n_snp = n_snp.flatten()

n_singleton = n_singleton.flatten()


mean_coverage = mean_coverage.select('region',
                                     'Population',
                                     'cov_stats.n',
                                     'cov_stats.mean',
                                     'cov_stats.stdev')

n_snp = n_snp.select('region',
                     'Population',
                     'n_snp_stats.n', 
                     'n_snp_stats.mean',
                     'n_snp_stats.stdev')

n_singleton = n_singleton.select('region',
                                 'Population',
                                 'n_singleton_stats.n',
                                 'n_singleton_stats.mean',
                                 'n_singleton_stats.stdev')


mean_coverage = mean_coverage.key_by('region','Population')

n_snp = n_snp.key_by('region', 'Population')

n_singleton = n_singleton.key_by('region', 'Population')


table = mean_coverage.annotate(n_snp = n_snp[mean_coverage.region, mean_coverage.Population])


table = table.annotate(n_singleton = n_singleton[table.region, table.Population])


# Flattening out the structs created from annotating the tables
table = table.flatten()


table = table.key_by(table.region, table.Population)


table.export('gs://african-seq-data/hgdp_tgp/table_x.tsv', header=True)

