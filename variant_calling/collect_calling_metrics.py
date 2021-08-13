#!/usr/bin/env python3

__author__ = 'Lindo Nkambule'

import hailtop.batch as hb


def collect_variant_calling_metrics(b: hb.batch.Batch, input_vcf: hb.resource.ResourceGroup,
                                    dbsnp_vcf_file: hb.resource.ResourceGroup, ref_dict: hb.resource.ResourceGroup,
                                    evaluation_int_list: hb.resource.ResourceFile, metrics_basename: str = None,
                                    memory: int = 6, docker_img: str = None, storage: int = None, out_dir: str = None):

    """
    Call germline SNPs and indels
    :param b: batch
    :param input_vcf: GVCF file to collect variant calling metrics for
    :param dbsnp_vcf_file: DBSNP VCF and its index to use in collecting metrics
    :param ref_dict: reference dictionary file from a reference Group (fasta, index, and dict) Resource
    :param evaluation_int_list: evaluation list file
    :param metrics_basename: name to be used for output
    :param memory: job memory
    :param docker_img: image to use for the job
    :param storage: storage to use fo the job
    :param out_dir: output directory
    :return:
    """

    docker_image = docker_img if docker_img else 'us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.3-1564508330'

    collect_vc_metrics = b.new_job(name=metrics_basename)
    collect_vc_metrics.image(docker_image)
    collect_vc_metrics.memory(f'{memory}Gi')
    collect_vc_metrics.storage(f'{storage}Gi')
    collect_vc_metrics.command(f'java -Xms2000m -jar /usr/gitc/picard.jar \
      CollectVariantCallingMetrics \
      INPUT={input_vcf.vcf} \
      OUTPUT={metrics_basename} \
      DBSNP={dbsnp_vcf_file.vcf} \
      SEQUENCE_DICTIONARY={ref_dict.dict} \
      TARGET_INTERVALS={evaluation_int_list} \
      GVCF_INPUT=true')
    collect_vc_metrics.command(f'mv {metrics_basename}.variant_calling_detail_metrics {collect_vc_metrics.detail}')
    collect_vc_metrics.command(f'mv {metrics_basename}.variant_calling_summary_metrics {collect_vc_metrics.summary}')
    b.write_output(collect_vc_metrics.detail,
                   f'{out_dir}/variant-calling-metrics/{metrics_basename}/{metrics_basename}.variant_calling_detail_metrics')
    b.write_output(collect_vc_metrics.summary,
                   f'{out_dir}/variant-calling-metrics/{metrics_basename}/{metrics_basename}.variant_calling_summary_metrics')

    return collect_vc_metrics
