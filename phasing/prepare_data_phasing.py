__author__ = 'Lindo Nkambule'

import hail as hl
from gnomad.utils.filtering import filter_to_adj
hl.init(tmp_dir='gs://hgdp-1kg/phasing/tmp/')  # to avoid export_vcf error about multiple datanodes


def run_qc(mt):
    """
    Apply gnomAD's sample, variant and genotype QC filters
    """

    # 1. Sample QC
    # This filters to only samples that passed gnomAD's sample QC hard filters
    mt = mt.filter_cols(~mt.gnomad_sample_filters.hard_filtered)  # removed 31 samples

    # 2. Variant QC
    # This subsets to only PASS variants - those which passed gnomAD's variant QC
    # PASS variants have an entry in the filters field
    mt = mt.filter_rows(hl.len(mt.filters) != 0, keep=False)

    # 3. Genotype QC filters to the dataset
    # This is done using a function imported from gnomAD and is the last step in the QC process
    mt = filter_to_adj(mt)

    return mt


# Path for HGDP+1kGP dataset prior to applying gnomAD QC filters
pre_qc_path = 'gs://gcp-public-data--gnomad/release/3.1.2/mt/genomes/gnomad.genomes.v3.1.2.hgdp_1kg_subset_dense.mt'
pre_qc_mt = hl.read_matrix_table(pre_qc_path)
print(f'Num of SNVs and samples prior to any analysis = {pre_qc_mt.count()}')

# 1. QC MT
post_qc_mt = run_qc(pre_qc_mt)
print(f'Num of SNVs and samples after applying QC filters = {post_qc_mt.count()}')

# 2. QUALapprox and AS_QUALapprox are coded as an int64 and VCF specs only allow up to int32
post_qc_mt = post_qc_mt.annotate_rows(info=post_qc_mt.info.annotate(AS_QUALapprox=hl.float64(post_qc_mt.info.AS_QUALapprox),
                                                                    QUALapprox=hl.float64(post_qc_mt.info.QUALapprox)))

# 3. Export VCFs per chromosome
chroms = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
for chrom in chroms:
    print(f'<============================{chrom}============================>')
    mt_chrom = hl.filter_intervals(
        post_qc_mt,
        [hl.parse_locus_interval(x, reference_genome='GRCh38') for x in [chrom]])

    hl.export_vcf(mt_chrom, f'gs://hgdp-1kg/phasing/input_vcfs/hgdp1kg_{chrom}.vcf.bgz', tabix=True)
