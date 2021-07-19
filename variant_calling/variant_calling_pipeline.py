#!/usr/bin/env python3

__author__ = 'Lindo Nkambule'

import hailtop.batch as hb
import hail as hl
import pandas as pd
import os

from variant_calling.get_file_size import bytes_to_gb
from variant_calling.scatter_interval import scatter_interval_list
from variant_calling.haplotype_caller import haplotype_caller_gatk
from variant_calling.merge_gvcfs import merge_vcf
from variant_calling.validate_gvcf import validate_vcf
from variant_calling.gvcf_index import index_gvcf
from variant_calling.collect_calling_metrics import collect_variant_calling_metrics


def var_call_pipeline(input_bams: str, ref_fasta: str, ref_ind: str, ref_dict: str, calling_int_list: str,
                      evaluation_int_list: str, dbsnp_vcf: str, dbsnp_vcf_ind: str, contamination: float = 0.0,
                      scatter_count: int = 50, out_dir: str = None, backend=None):

    ref_fasta_size = bytes_to_gb(ref_fasta)
    ref_dict_size = bytes_to_gb(ref_dict)
    ref_ind_size = bytes_to_gb(ref_ind)

    print("Scattering interval file list into chunks")
    scatter = hb.Batch(backend=backend,
                       name='scatter-interval-list')
    calling_interval_list = scatter.read_input(calling_int_list)
    scatter_intervals = scatter_interval_list(b=scatter, interval_list_file=calling_interval_list, out_dir=out_dir,
                                              scatter_count=scatter_count)
    scatter.run()

    interval_files = hl.utils.hadoop_ls(f'{out_dir}/scatter-intervals/**')

    var_call = hb.Batch(backend=backend,
                        name='variant-calling')
    print("Running variant calling")

    fasta = var_call.read_input_group(**{'fasta': ref_fasta,
                                         'fasta.fai': ref_ind,
                                         'dict': ref_dict})

    bams_paths = pd.read_csv(input_bams, sep='\t', header=None)
    bam_files = []
    for index, row in bams_paths.iterrows():
        bam_files.append((row[0], row[1]))

    for bam, bai in bam_files:
        base_bam = os.path.basename(bam)  # get file name without path
        bam_no_ext = os.path.splitext(base_bam)[0]
        in_bam = var_call.read_input_group(**{'bam': bam,
                                              'bam.bai': bai})

        # job storage
        potential_hc_divisor: int = scatter_count - 20
        hc_divisor: int = potential_hc_divisor if potential_hc_divisor > 1 else 1
        vc_disk_size = round(
            ((bytes_to_gb(bam) + 30) / hc_divisor) + ref_fasta_size + ref_dict_size + ref_ind_size) + 20

        for file in interval_files:
            interval_file = var_call.read_input(file['path'])
            # we have to use the first part of the scatter interval file name + bam basename + .g.vcf.gz so that the
            # gvcf don't have the same name
            interval_file_name = file['path']
            base_interval_file_name = os.path.basename(interval_file_name)
            base_interval_file_name_no_ext = os.path.splitext(base_interval_file_name)[0]
            var_calling = haplotype_caller_gatk(b=var_call, interval_list_file=interval_file, input_bam=in_bam,
                                                ref_fasta=fasta, out_dir=out_dir, contamination=contamination,
                                                bam_filename_no_ext=bam_no_ext, storage=vc_disk_size,
                                                interval_list_name=base_interval_file_name_no_ext)

    var_call.run()

    merge_gvcfs = hb.Batch(backend=backend,
                           name='merge-gvcfs')
    print("Merging GVCFs")

    for bam, bai in bam_files:
        base_bam = os.path.basename(bam)  # get file name without path
        bam_no_ext = os.path.splitext(base_bam)[0]

        # get gvcf files to merge
        gvcfs_to_merge = hl.utils.hadoop_ls(f'{out_dir}/variant-calling/{bam_no_ext}/*.vcf.gz')
        gvcfs_list = []
        gvcfs_sizes_sum = 0
        for file in gvcfs_to_merge:
            gvcfs_list.append(file['path'])
            gvcfs_sizes_sum += bytes_to_gb(file['path'])

        merge_disk_size = round(gvcfs_sizes_sum * 2.5) + 10

        merge_gvcfs_job = merge_vcf(b=merge_gvcfs, gvcf_list=gvcfs_list, storage=merge_disk_size,
                                    output_vcf_name=bam_no_ext, out_dir=out_dir)

    merge_gvcfs.run()

    validate_index_gvcf = hb.Batch(backend=backend,
                                   name='validate-index-gvcf')
    print("Validating and indexing GVCFs")

    fasta_validate = validate_index_gvcf.read_input_group(**{'fasta': ref_fasta,
                                                             'fasta.fai': ref_ind,
                                                             'dict': ref_dict})
    dbsnp_vcf_validate = validate_index_gvcf.read_input_group(**{'vcf': dbsnp_vcf,
                                                                 'vcf.idx': dbsnp_vcf_ind})
    calling_interval_list_validate = validate_index_gvcf.read_input(calling_interval_list)

    for bam, bai in bam_files:
        base_bam = os.path.basename(bam)
        bam_no_ext = os.path.splitext(base_bam)[0]

        in_vcf_size = bytes_to_gb(f'{out_dir}/merged-gvcf/{bam_no_ext}/{bam_no_ext}.g.vcf.gz')
        validate_ref_size = ref_fasta_size + ref_dict_size + ref_ind_size
        validate_disk_size = round(in_vcf_size + bytes_to_gb(dbsnp_vcf) + validate_ref_size) + 20

        in_vcf = validate_index_gvcf.read_input(f'{out_dir}/merged-gvcf/{bam_no_ext}/{bam_no_ext}.g.vcf.gz')

        validate_gvcf_out = validate_vcf(b=validate_index_gvcf, input_vcf=in_vcf, ref_fasta=fasta_validate,
                                         dbsnp_vcf_file=dbsnp_vcf_validate, storage=validate_disk_size,
                                         calling_int_file=calling_interval_list_validate,
                                         output_vcf_ind_name=bam_no_ext)

        gvcf_index_file = index_gvcf(b=validate_index_gvcf, input_vcf=in_vcf, output_vcf_ind_name=bam_no_ext,
                                     out_dir=out_dir)

    validate_index_gvcf.run()

    collect_vr_metric = hb.Batch(backend=backend,
                                 name='collect-variant-calling-metrics')
    print("Collecting variant calling metrics")

    fasta_collect = collect_vr_metric.read_input_group(**{'fasta': ref_fasta,
                                                          'fasta.fai': ref_ind,
                                                          'dict': ref_dict})
    dbsnp_vcf_collect = collect_vr_metric.read_input_group(**{'vcf': dbsnp_vcf,
                                                              'vcf.idx': dbsnp_vcf_ind})
    evaluation_interval_list_collect = collect_vr_metric.read_input(evaluation_int_list)

    for bam, bai in bam_files:
        base_bam = os.path.basename(bam)
        bam_no_ext = os.path.splitext(base_bam)[0]

        vcf_group = collect_vr_metric.read_input_group(
            **{'vcf': f'{out_dir}/merged-gvcf/{bam_no_ext}/{bam_no_ext}.g.vcf.gz',
               'vcf.tbi': f'{out_dir}/merged-gvcf/{bam_no_ext}/{bam_no_ext}.g.vcf.gz.tbi'})

        metrics_disk_size = round(bytes_to_gb(f'{out_dir}/merged-gvcf/{bam_no_ext}/{bam_no_ext}.g.vcf.gz') +
                                  bytes_to_gb(dbsnp_vcf)) + 20
        vc_metrics = collect_variant_calling_metrics(b=collect_vr_metric, input_vcf=vcf_group,
                                                     dbsnp_vcf_file=dbsnp_vcf_collect, ref_dict=fasta_collect,
                                                     evaluation_int_list=evaluation_interval_list_collect,
                                                     metrics_basename=bam_no_ext, storage=metrics_disk_size,
                                                     out_dir=out_dir)

    collect_vr_metric.run()


