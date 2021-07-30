#!/usr/bin/env python3

__author__ = 'Lindo Nkambule'

import hailtop.batch as hb


def index_gvcf(b: hb.batch.Batch, input_vcf: hb.resource.ResourceFile, output_vcf_ind_name: str = None, memory: int = 3,
               storage: int = 5, docker_img: str = None, out_dir: str = None):
    """
    Index a GVCF file
    :param b: batch
    :param input_vcf: GVCF file to index
    :param output_vcf_ind_name: output GVCF index name
    :param memory: job memory
    :param storage: storage to use fo the job
    :param docker_img: image to use for the job
    :param out_dir: output directory
    :return:
    """

    docker_image = docker_img if docker_img else 'us.gcr.io/broad-gatk/gatk:4.2.0.0'
    outname = output_vcf_ind_name + '.g.vcf.gz.tbi'

    index_gvcf_file = b.new_job(name=f'index-{output_vcf_ind_name}')
    index_gvcf_file.image(docker_image)
    index_gvcf_file.memory(f'{memory}Gi')
    index_gvcf_file.storage(f'{storage}Gi')
    index_gvcf_file.command(f'gatk IndexFeatureFile \
         -I {input_vcf} \
         -O {outname}')
    index_gvcf_file.command(f'mv {outname} {index_gvcf_file.ofile}')
    b.write_output(index_gvcf_file.ofile, f'{out_dir}/merged-gvcf/{output_vcf_ind_name}/{outname}')

    return index_gvcf_file
