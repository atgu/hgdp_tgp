#!/usr/bin/env python3

__author__ = 'Lindo Nkambule'

import hailtop.batch as hb


def validate_vcf(b: hb.batch.Batch, input_vcf: hb.resource.ResourceFile, ref_fasta: hb.resource.ResourceGroup,
                 dbsnp_vcf_file: hb.resource.ResourceGroup, calling_int_file: hb.resource.ResourceFile,
                 validate_vcf_img: str = None, memory: int = 7, storage: int = None,
                 output_vcf_ind_name: str = None):
    """
    Validate the GVCF output of HaplotypeCaller
    :param b: batch
    :param input_vcf: GVCF file to validate
    :param ref_fasta: reference files, including fasta and index
    :param dbsnp_vcf_file: DBSNP VCF and its index to use in the validation
    :param calling_int_file: calling interval file
    :param validate_vcf_img: image to use for the job
    :param memory: job memory
    :param storage: storage to use fo the job
    :param output_vcf_ind_name: name to be used for the job
    :return:
    """

    docker_image = validate_vcf_img if validate_vcf_img else 'us.gcr.io/broad-gatk/gatk:4.2.0.0'

    validate_gvcf = b.new_job(name=output_vcf_ind_name)
    validate_gvcf.image(docker_image)
    validate_gvcf.memory(f'{memory}Gi')
    validate_gvcf.storage(f'{storage}Gi')
    validate_gvcf.command(f'gatk IndexFeatureFile \
             -I {input_vcf}')
    validate_gvcf.command(f'gatk --java-options -Xms6000m \
      ValidateVariants \
      -V {input_vcf} \
      -R {ref_fasta.fasta} \
      -L {calling_int_file} \
      -gvcf \
      --validation-type-to-exclude ALLELES \
      --dbsnp {dbsnp_vcf_file.vcf}')

    return validate_gvcf
