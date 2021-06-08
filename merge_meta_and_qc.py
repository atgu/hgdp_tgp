import hail as hl
import re

hl.init()

# Takes a list of dicts and converts it to a struct format (works with nested structs too)
def dict_to_struct(d):
    fields = {}
    for k, v in d.items():
        if isinstance(v, dict):
            v = dict_to_struct(v)
        fields[k] = v
    return hl.struct(**fields)


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


# These bits below were written by Tim Poterba to help troubleshoot unflattening a ht with nested structure
# dict to hold struct names as well as nested field names
d = {}

# Getting just the row field names 
row = sample_meta.row_value

# returns a dict with the struct names as keys and their inner field names as values
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


# using the dict created from flattened struct, creating new structs now unflattened
sample_meta = sample_meta.select(**dict_to_struct(d))
sample_meta = sample_meta.key_by('s')


# grabbing the columns needed from Alicia's metadata
new_meta = sample_meta.select(sample_meta.hgdp_tgp_meta, sample_meta.bergstrom)

# creating a table with Julia's metadata and Alicia's metadata
ht = jul_meta.annotate(**new_meta[jul_meta.s])

ht.write('gs://african-seq-data/hgdp_tgp/hgdp_tgp_sample_metadata.ht')

# reading in table annotated with Alicia and Julia's respective metadata
ht = hl.read_table('gs://african-seq-data/hgdp_tgp/hgdp_tgp_sample_metadata.ht')

# filtering samples to those who should pass QC
mt = dense_mt.annotate_cols(**ht[dense_mt.s])
all_sample_filters = set(ht['sample_filters'])
# bad_sample_filters were removing whole populations (mostly AFR and OCE) that passed all other QC, bad filters
bad_sample_filters = {re.sub('fail_', '', x) for x in all_sample_filters if x.startswith('fail_')}
# this filters to only variants that passed all gnomad QC or only failed filters in bad_sample_filters
mt_filt = mt.filter_cols((mt['sample_filters']['qc_metrics_filters'].difference(bad_sample_filters).length() == 0))


# annotating partially filtered dataset with variant metadata
mt_filt = mt_filt.annotate_rows(**var_meta[mt_filt.locus, mt_filt.alleles])


mt_filt.write('gs://african-seq-data/hgdp_tgp/hgdp_tgp_dense_meta_filt.mt')
# Reading in the annotated & partially filtered dataset
mt = hl.read_matrix_table('gs://african-seq-data/hgdp_tgp/hgdp_tgp_dense_meta_filt.mt')


# This is a set of outliers found by Mary during PCA analyses as well as one duplicate sample
outlier_set = {'NA20314','NA20299','HG01880','HG01881',
               'HGDP00130','HGDP00013','HGDP00150',
               'HGDP00029','HGDP01298','HGDP01303',
               'LP6005443-DNA_B02','HGDP01300','HG01628',
               'HG01629','HG01630','HG01694','HG01696',
               'HGDP00621','HGDP01270','HGDP01271'}
set_to_remove = hl.literal(outlier_set)
mt = mt.filter_cols(~set_to_remove.contains(mt['s']))
mt = mt.distinct_by_col()


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

# writing out a vcf version of the dataset for downstream analyses
mt_vcf = mt_vcf.drop('gvcf_info')

hl.export_vcf(mt_vcf, 'gs://african-seq-data/hgdp_tgp/hgdp_tgp_postqc.vcf.bgz', parallel='separate_header')


# Subsetting the variants in the dataset to only PASS variants (those which passed variant QC)
mt = mt.filter_rows(hl.len(mt.filters) !=0  ,keep=False)


# writing out the postQC dataset with PCA sample outliers removed and subset to PASS variants
mt.write('gs://african-seq-data/hgdp_tgp/hgdp_tgp_postQC.mt')