#!/usr/bin/env python3

__author__ = 'Lindo Nkambule'

import hailtop.batch as hb


def haplotype_caller_gatk(b: hb.batch.Batch, input_bam: hb.resource.ResourceGroup, ref_fasta: hb.resource.ResourceGroup,
                          interval_list_file: hb.resource.ResourceFile, bam_filename_no_ext: str = None,
                          out_dir: str = None, interval_list_name: str = None, storage: int = None,
                          contamination: float = None, gatk_img: str = None, memory: float = 6.5, ncpu: int = 2):
    """
    Call germline SNPs and indels
    :param b: batch
    :param input_bam: BAM file
    :param ref_fasta: reference files, including fasta and index
    :param interval_list_file: interval list file with intervals to run variant calling on
    :param bam_filename_no_ext: BAM filename without extension
    :param interval_list_name: interval list name, used to name output GVCF
    :param storage: storage to use fo the job
    :param contamination: fraction of contamination in sequencing data to aggressively remove
    :param gatk_img: image to use for the job
    :param out_dir: output directory
    :param memory: job memory
    :param ncpu: number of CPUs
    :return:
    """

    docker_image = gatk_img if gatk_img else 'us.gcr.io/broad-gatk/gatk:4.2.0.0'

    output_file_name = bam_filename_no_ext + '_' + interval_list_name + '.g.vcf.gz'

    variant_calling = b.new_job(name=bam_filename_no_ext)

    variant_calling.image(docker_image)
    variant_calling.cpu(ncpu)
    variant_calling.memory(f'{memory}Gi')
    variant_calling.storage(f'{storage}Gi')
    variant_calling.command(f'gatk --java-options "-Xms6000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
      HaplotypeCaller \
      -R {ref_fasta.fasta} \
      -I {input_bam.bam} \
      -L {interval_list_file} \
      -O {variant_calling.ofile}\
      -contamination {contamination} \
      -G StandardAnnotation -G StandardHCAnnotation -G AS_StandardAnnotation \
      -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
      -ERC GVCF')
    # variant_calling.command(f'mv {output_file_name} {variant_calling.ofile}')
    b.write_output(variant_calling.ofile, f'{out_dir}/variant-calling/{bam_filename_no_ext}/{output_file_name}')

    return variant_calling
