#!/usr/bin/env python3

__author__ = 'Lindo Nkambule'

import argparse
import hailtop.batch as hb
import hail as hl
import pandas as pd
import os


def bytes_to_gb(in_file):
    """ Convert the size from bytes to GB"""

    file_info = hl.utils.hadoop_stat(in_file)
    size_bytes = file_info['size_bytes']
    size_gigs = size_bytes / (1024 * 1024 * 1024)

    return size_gigs


def scatter_interval_list(b: hb.batch.Batch, interval_list_file: hb.resource.ResourceFile, scatter_count: int = 50,
                          break_bands_at_multiples_of: int = 1000000, scatter_img: str = None, memory: int = 2,
                          out_dir: str = None):
    # break the calling interval list into sub-intervals
    docker_image = scatter_img if scatter_img else 'us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.3-1564508330'
    job_memory = str(memory) + 'Gi'

    scatter_list = b.new_job(name='scatter-interval-list')

    scatter_list.image(docker_image)
    scatter_list.cpu(4)  # this should be lower, check DSP pipeline
    scatter_list.memory(job_memory)
    scatter_list.command('mkdir /scatter_intervals')
    scatter_list.command(f'java -Xms1g -jar /usr/gitc/picard.jar \
      IntervalListTools \
      SCATTER_COUNT={scatter_count} \
      SUBDIVISION_MODE=BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
      UNIQUE=true \
      SORT=true \
      BREAK_BANDS_AT_MULTIPLES_OF={break_bands_at_multiples_of} \
      INPUT={interval_list_file} \
      OUTPUT=/scatter_intervals')
    scatter_list.command('''
    cat > my_script.py <<EOF
import sys
import os
import glob

intervals = sorted(glob.glob('/scatter_intervals/*/*.interval_list'))
for i, interval in enumerate(intervals):
      (directory, filename) = os.path.split(interval)
      newName = os.path.join(directory, str(i + 1) + filename)
      os.rename(interval, newName)
EOF
python3 my_script.py
    ''')
    scatter_list.command(f'mv /scatter_intervals {scatter_list.outfiles}')
    b.write_output(scatter_list.outfiles, f'{out_dir}/scatter_intervals')

    # We return the `scatter_list` Job object that can be used in downstream jobs.
    return scatter_list


def haplotype_caller_gatk(b: hb.batch.Batch, input_bam: hb.resource.ResourceGroup, ref_fasta: hb.resource.ResourceGroup,
                          interval_list_file: hb.resource.ResourceFile, bam_filename_no_ext: str = None,
                          out_dir: str = None, interval_list_name: str = None,
                          contamination: float = None, gatk_img: str = None, memory: float = 6.5, ncpu: int = 2):
    docker_image = gatk_img if gatk_img else 'us.gcr.io/broad-gatk/gatk:4.0.10.1'
    job_memory = str(memory) + 'Gi'

    output_file_name = bam_filename_no_ext + '_' + interval_list_name + '.g.vcf.gz'
    # print(output_file_name)

    variant_calling = b.new_job(name='variant-calling')

    variant_calling.image(docker_image)
    variant_calling.cpu(ncpu)
    variant_calling.memory(job_memory)
    variant_calling.command(f'gatk --java-options "-Xms6000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
      HaplotypeCaller \
      -R {ref_fasta.fasta} \
      -I {input_bam.bam} \
      -L {interval_list_file} \
      -O {output_file_name} \
      -contamination {contamination} \
      -G StandardAnnotation -G StandardHCAnnotation -G AS_StandardAnnotation \
      -new-qual \
      -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
      -ERC GVCF')
    variant_calling.command(f'mv {output_file_name} {variant_calling.ofile}')
    b.write_output(variant_calling.ofile, f'{out_dir}/variant_calling/{bam_filename_no_ext}/{output_file_name}')

    # We return the `variant_calling` Job object that can be used in downstream jobs.
    return variant_calling


def merge_vcf(b: hb.batch.Batch, inputs_vcfs_list_file: hb.resource.ResourceFile, output_vcf_name: str = None,
              merge_vcfs_img: str = None, memory: int = 3, out_dir: str = None):
    # Combine multiple VCFs or GVCFs from scattered HaplotypeCaller runs
    docker_image = merge_vcfs_img if merge_vcfs_img else 'us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.3-1564508330'
    job_memory = str(memory) + 'Gi'
    outname = output_vcf_name + '.g.vcf.gz'

    # disk_size = bytes_to_gb((inputs_vcfs_list * 2.5)) + 10

    merge_vcfs = b.new_job(name='merge-vcfs')
    merge_vcfs.image(docker_image)
    merge_vcfs.memory(job_memory)
    # merge_vcfs.memory(disk_size)
    merge_vcfs.command(f'java -Xms2000m -jar /usr/gitc/picard.jar \
      MergeVcfs \
      INPUT={inputs_vcfs_list_file} \
      OUTPUT={outname}')
    merge_vcfs.command(f'mv {outname} {merge_vcfs.ofile}')
    b.write_output(merge_vcfs.ofile, f'{out_dir}/merged_gvcf/{output_vcf_name}/{outname}')

    return merge_vcfs


def validate_vcf(b: hb.batch.Batch, input_vcf: hb.resource.ResourceFile, ref_fasta: hb.resource.ResourceGroup,
                 dbsnp_vcf_file: hb.resource.ResourceGroup, calling_int_file: hb.resource.ResourceFile,
                 validate_vcf_img: str = None, memory: int = 7):
    # Validate the (g)VCF output of HaplotypeCaller

    docker_image = validate_vcf_img if validate_vcf_img else 'us.gcr.io/broad-gatk/gatk:4.0.10.1'
    job_memory = str(memory) + 'Gi'

    # ref_size = bytes_to_gb(ref_fasta) + bytes_to_gb(ref_fasta_index) + bytes_to_gb(ref_dict)
    # disk_size = bytes_to_gb(input_vcf) + bytes_to_gb(dbsnp_vcf) + ref_size + 20

    validate_gvcf = b.new_job(name='validate-vcf')
    validate_gvcf.image(docker_image)
    validate_gvcf.memory(job_memory)
    # validate_gvcf.storage(disk_size)
    validate_gvcf.command(f'gatk --java-options -Xms6000m \
      ValidateVariants \
      -V {input_vcf} \
      -R {ref_fasta.fasta} \
      -L {calling_int_file} \
      -gvcf \
      --validation-type-to-exclude ALLELES \
      --dbsnp {dbsnp_vcf_file.vcf}')

    return validate_gvcf


def collect_variant_calling_metrics(b: hb.batch.Batch, input_vcf: hb.resource.ResourceGroup,
                                    dbsnp_vcf_file: hb.resource.ResourceGroup, ref_dict: hb.resource.ResourceGroup,
                                    evaluation_int_list: hb.resource.ResourceFile, metrics_basename: str = None,
                                    memory: int = 3, docker_img: str = None, out_dir: str = None):
    docker_image = docker_img if docker_img else 'us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.3-1564508330'
    job_memory = str(memory) + 'Gi'

    # disk_size = bytes_to_gb(input_vcf) + bytes_to_gb(dbsnp_vcf) + 20

    collect_vc_metrics = b.new_job(name='collect-variant-calling-metrics')
    collect_vc_metrics.image(docker_image)
    collect_vc_metrics.memory(job_memory)
    # collect_vc_metrics.storage(disk_size)
    collect_vc_metrics.command(f'java -Xms2000m -jar /usr/gitc/picard.jar \
      CollectVariantCallingMetrics \
      INPUT={input_vcf.vcf} \
      OUTPUT={metrics_basename} \
      DBSNP={dbsnp_vcf_file.vcf} \
      SEQUENCE_DICTIONARY={ref_dict.dict} \
      TARGET_INTERVALS={evaluation_int_list} \
      GVCF_INPUT=true')
    collect_vc_metrics.command(f'mv -R {metrics_basename}.variant_calling_* {collect_vc_metrics.ofile}')
    b.write_output(collect_vc_metrics.ofile,
                                    f'{out_dir}/variant-calling-metrics/{metrics_basename}.variant_calling_*')

    return collect_vc_metrics


if __name__ == '__main__':
    parser = argparse.ArgumentParser('Variant Calling Pipeline')
    parser.add_argument('--input-bams', required=True)
    parser.add_argument('--ref-fasta',
                        default='gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta')
    parser.add_argument('--ref-index',
                        default='gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai')
    parser.add_argument('--ref-dict',
                        default='gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict')
    parser.add_argument('--billing-project', default='diverse-pop-seq-ref')
    parser.add_argument('--bucket', default='african-seq-data')
    parser.add_argument('--local', action='store_true')
    parser.add_argument('--calling-interval-list',
                        default='gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.interval_list')
    parser.add_argument('--evaluation-interval-list',
                        default=
                        'gs://gcp-public-data--broad-references/hg38/v0/wgs_evaluation_regions.hg38.interval_list')
    parser.add_argument('--dbsnp-vcf',
                        default='gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf')
    parser.add_argument('--dbsnp-vcf-ind',
                        default=
                        'gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx')
    parser.add_argument('--contamination', type=float, default=0.0)
    parser.add_argument('--out-dir', default='gs://african-seq-data')

    args = parser.parse_args()

    if args.local:
        backend = hb.LocalBackend()
    else:
        backend = hb.ServiceBackend(billing_project=args.billing_project,
                                    bucket=args.bucket)

    scatter = hb.Batch(backend=backend,
                       name='scatter-interval-list')
    calling_interval_list = scatter.read_input(args.calling_interval_list)
    evaluation_interval_list = scatter.read_input(args.evaluation_interval_list)
    scatter_intervals = scatter_interval_list(b=scatter, interval_list_file=calling_interval_list, out_dir=args.out_dir)
    scatter.run()

    interval_files = hl.utils.hadoop_ls(f'{args.out_dir}/scatter_intervals/**')

    var_call = hb.Batch(backend=backend,
                        name='variant-calling')
    fasta = var_call.read_input_group(**{'fasta': args.ref_fasta,
                                         'fasta.fai': args.ref_index,
                                         'dict': args.ref_dict})
    dbsnp_vcf = var_call.read_input_group(**{'vcf': args.dbsnp_vcf,
                                             'vcf.idx': args.dbsnp_vcf_ind})

    bams_paths = pd.read_csv(args.input_bams, sep='\t', header=None)
    bam_files = []
    for index, row in bams_paths.iterrows():
        bam_files.append([row[0], row[1]])

    for bam in bam_files:
        base_bam = os.path.basename(bam[0])  # get file name without path
        bam_no_ext = os.path.splitext(base_bam)[0]
        in_bam = var_call.read_input_group(**{'bam': bam[0],
                                              'bam.bai': bam[1]})
        for file in interval_files:
            interval_file = var_call.read_input(file['path'])
            # we have to use the first part of the scatter interval file name + bam basename + .g.vcf.gz so that the
            # gvcf don't have the same name
            interval_file_name = file['path']
            base_interval_file_name = os.path.basename(interval_file_name)
            base_interval_file_name_no_ext = os.path.splitext(base_interval_file_name)[0]
            var_calling = haplotype_caller_gatk(b=var_call, interval_list_file=interval_file, input_bam=in_bam,
                                                ref_fasta=fasta, out_dir=args.out_dir, contamination=args.contamination,
                                                bam_filename_no_ext=bam_no_ext,
                                                interval_list_name=base_interval_file_name_no_ext)

    var_call.run()

    for bam in bam_files:
        base_bam = os.path.basename(bam[0])  # get file name without path
        bam_no_ext = os.path.splitext(base_bam)[0]
        gvcfs_to_merge = hl.utils.hadoop_ls(f'{args.out_dir}/variant_calling/{bam_no_ext}/*')
        with hl.hadoop_open(f'{args.out_dir}/variant_calling/{bam_no_ext}/gvcfs_to_merge.list', 'w') as f:
            for file in gvcfs_to_merge:
                f.write(file['path'] + '\n')
        # hl.hadoop_copy(f'{bam_no_ext}_gvcfs_to_merge.list', f'{args.out_dir}/variant_calling/{bam_no_ext}/')

    merge_gvcfs = hb.Batch(backend=backend,
                           name='merge-gvcfs')

    for bam in bam_files:
        base_bam = os.path.basename(bam[0])  # get file name without path
        bam_no_ext = os.path.splitext(base_bam)[0]
        gvcfs_to_merge_file = merge_gvcfs.read_input(f'{args.out_dir}/variant_calling/{bam_no_ext}/gvcfs_to_merge.list')
        merge_gvcfs_job = merge_vcf(b=merge_gvcfs, inputs_vcfs_list_file=gvcfs_to_merge_file,
                                    output_vcf_name=bam_no_ext, out_dir=args.out_dir)

    # This is where the pipeline currently fails and throws
    # Exception in thread "main" java.nio.file.ProviderNotFoundException: Provider "gs" not found
    # This happens because the file {args.out_dir}/variant_calling/{bam_no_ext}/gvcfs_to_merge.list has a list of gvcs,
    # one each line, to be merged with a Google Storage path i.e. path starts with a gs:// and picard is not able to
    # read files from the cloud
    # TODO: Ask Jackie how we can 'localize' the files
    merge_gvcfs.run()

        # merge_gvcfs = merge_vcf(b=batch, inputs_vcfs_list_file=gvcfs_to_merge_file, output_vcf_name=bam_no_ext)
        # vcf_validation = validate_vcf(b=batch, input_vcf=merge_gvcfs.ofile, ref_fasta=fasta,dbsnp_vcf_file=dbsnp_vcf,
        # calling_int_file=calling_interval_list)
        # collect_var_call_metrics = collect_variant_calling_metrics(b=batch, input_vcf=merge_gvcfs.ofile,
        # dbsnp_vcf_file=dbsnp_vcf, ref_dict=fasta,
        # metrics_basename=bam_no_ext,
        # evaluation_int_list=evaluation_interval_list)