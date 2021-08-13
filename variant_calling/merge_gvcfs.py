#!/usr/bin/env python3

__author__ = 'Lindo Nkambule'

import hailtop.batch as hb
from typing import List


def merge_vcf(b: hb.batch.Batch, gvcf_list: List = None, output_vcf_name: str = None,
              merge_vcfs_img: str = None, memory: int = 3, out_dir: str = None, storage: int = None):
    """
    Combine multiple VCFs or GVCFs from scattered HaplotypeCaller runs
    :param b: batch
    :param gvcf_list: list of GVCF files to merge
    :param output_vcf_name: output GVCF name
    :param merge_vcfs_img: image to use for the job
    :param storage: Storage to use fo the job
    :param out_dir: output directory
    :param memory: job memory
    :return:
    """

    docker_image = merge_vcfs_img if merge_vcfs_img else\
        'us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.3-1564508330'

    outname = output_vcf_name + '.g.vcf.gz'

    merge_vcf_i = ''

    for line in gvcf_list:
        input_gvcf = b.read_input(line)
        merge_vcf_i += f'I={input_gvcf} \t'

    merge_vcfs = b.new_job(name=output_vcf_name)
    merge_vcfs.image(docker_image)
    merge_vcfs.memory(f'{memory}Gi')
    merge_vcfs.storage(f'{storage}Gi')
    merge_vcfs.command(f'java -Xms2000m -jar /usr/gitc/picard.jar \
      MergeVcfs \
      {merge_vcf_i} \
      O={outname}')
    merge_vcfs.command(f'mv {outname} {merge_vcfs.ofile}')
    b.write_output(merge_vcfs.ofile, f'{out_dir}/merged-gvcf/{output_vcf_name}/{outname}')

    return merge_vcfs
