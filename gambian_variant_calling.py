import argparse
import hailtop.batch as hb
import hail as hl


def bytes_to_gb(file):
    """ Convert the size from bytes to GB"""

    file_info = hl.utils.hadoop_stat(file)
    size_bytes = file_info['size_bytes']
    size_gigs = size_bytes / (1024 * 1024 * 1024)

    return size_gigs


def scatter_interval_list(b: hb.batch.Batch, interval_list: hb.resource.ResourceFile, scatter_count: int = 50,
                          break_bands_at_multiples_of: int = 1000000, scatter_img: str = None, memory: int = 2,
                          out_dir: str = None):
    # break the calling interval list into sub-intervals
    docker_image = scatter_img if scatter_img else 'us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.3-1564508330'
    job_memory = str(memory) + 'Gi'

    scatter_list = b.new_job(name='scatter-interval-list')

    scatter_list.image(docker_image)
    scatter_list.cpu(4)  # this should be lower, check DSP pipeline
    scatter_list.memory(job_memory)
    scatter_list.command(f'java -Xms1g -jar /usr/gitc/picard.jar \
      IntervalListTools \
      SCATTER_COUNT=~{scatter_count} \
      SUBDIVISION_MODE=BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
      UNIQUE=true \
      SORT=true \
      BREAK_BANDS_AT_MULTIPLES_OF={break_bands_at_multiples_of} \
      INPUT=~{interval_list} \
      OUTPUT={out_dir}')

    # We return the `scatter_list` Job object that can be used in downstream jobs.
    return scatter_list


def haplotype_caller_gatk(b: hb.batch.Batch, input_bam: str = None, interval_list: str = None,
                          vcf_basename: str = None, ref_fasta: str = None, ref_dict: str = None, ref_index: str = None,
                          contamination: float = None, gatk_img: str = None, memory: float = 6.5, ncpu: int = 2):

    docker_image = gatk_img if gatk_img else 'us.gcr.io/broad-gatk/gatk:4.0.10.1'
    job_memory = str(memory) + 'Gi'
    output_file_name = vcf_basename + '.g.vcf.gz'

    variant_calling = b.new_job(name='variant-calling')

    variant_calling.image(docker_image)
    variant_calling.cpu(ncpu)
    variant_calling.memory(job_memory)
    variant_calling.command(f'gatk --java-options "-Xms6000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
      HaplotypeCaller \
      -R ~{ref_fasta} \
      -I ~{input_bam} \
      -L ~{interval_list} \
      -O ~{output_file_name} \
      -contamination {contamination} \
      -G StandardAnnotation -G StandardHCAnnotation -G AS_StandardAnnotation \
      -new-qual \
      -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
      -ERC GVCF')

    # We return the `variant_calling` Job object that can be used in downstream jobs.
    return variant_calling


def merge_vcf(b: hb.batch.Batch, inputs_vcfs_list: list = None, output_vcf_name: str = None,
              merge_vcfs_img: str = None, memory: int = 3):
    # Combine multiple VCFs or GVCFs from scattered HaplotypeCaller runs
    docker_image = merge_vcfs_img if merge_vcfs_img else 'us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.3-1564508330'
    job_memory = str(memory) + 'Gi'

    # disk_size = bytes_to_gb((inputs_vcfs_list * 2.5)) + 10

    merge_vcfs = b.new_job(name='merge-vcfs')
    merge_vcfs.image(docker_image)
    merge_vcfs.memory(job_memory)
    # merge_vcfs.memory(disk_size)
    merge_vcfs.command(f'ava -Xms2000m -jar /usr/gitc/picard.jar \
      MergeVcfs \
      INPUT={inputs_vcfs_list} \
      OUTPUT={output_vcf_name}')

    return merge_vcfs


def validate_vcf(b: hb.batch.Batch, input_vcf: str = None, input_vcf_index: str = None, ref_fasta: str = None,
                 ref_fasta_index: str = None, ref_dict: str = None, dbsnp_vcf_file: str = None,
                 dbsnp_vcf_index: str = None, calling_int_list: str = None,
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
      -R {ref_fasta} \
      -L {calling_int_list} \
      -gvcf \
      --validation-type-to-exclude ALLELES \
      --dbsnp {dbsnp_vcf_file}')

    return validate_gvcf


def collect_variant_calling_metrics(b: hb.batch.Batch, input_vcf: str = None, input_vcf_index: str = None,
                                    metrics_basename: str = None, dbsnp_vcf_file: str = None, ref_dict: str = None,
                                    dbsnp_vcf_index: str = None, evaluation_int_list: str = None,
                                    memory: int = 3, docker_img: str = None):
    docker_image = docker_img if docker_img else 'us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.3-1564508330'
    job_memory = str(memory) + 'Gi'

    # disk_size = bytes_to_gb(input_vcf) + bytes_to_gb(dbsnp_vcf) + 20

    collect_vc_metrics = b.new_job(name='collect-variant-calling-metrics')
    collect_vc_metrics.image(docker_image)
    collect_vc_metrics.memory(job_memory)
    # collect_vc_metrics.storage(disk_size)
    collect_vc_metrics.command(f'java -Xms2000m -jar /usr/gitc/picard.jar \
      CollectVariantCallingMetrics \
      INPUT={input_vcf} \
      OUTPUT={metrics_basename} \
      DBSNP={dbsnp_vcf_file} \
      SEQUENCE_DICTIONARY=~{ref_dict} \
      TARGET_INTERVALS={evaluation_int_list} \
      GVCF_INPUT=true')

    return collect_vc_metrics


if __name__ == '__main__':
    parser = argparse.ArgumentParser('Variant Calling Pipeline')
    parser.add_argument('--ref-fasta', required=True,
                        default='gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta')
    parser.add_argument('--ref-index', required=True,
                        default='gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai')
    parser.add_argument('--ref-dict', required=True,
                        default='gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict')
    parser.add_argument('--billing-project', required=True)
    parser.add_argument('--bucket', required=True)
    parser.add_argument('--local', action='store_true')
    parser.add_argument('--calling-interval-list', required=True,
                        default='gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.interval_list')
    parser.add_argument('--dbsnp-vcf', required=True,
                        default='gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf')
    parser.add_argument('--dbsnp-vcf-ind', required=True,
                        default=
                        'gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx')

    args = parser.parse_args()

    if args.local:
        backend = hb.LocalBackend()
    else:
        backend = hb.ServiceBackend(billing_project=args.billing_project,
                                    bucket=args.bucket)

    batch = hb.Batch(backend=backend,
                     name='Variant Calling Pipeline')

    # Define inputs
    vcf = batch.read_input(args.vcf)
    calling_interval_list = batch.read_input(args.calling_interval_list)
    fasta = batch.read_input_group(**{'fasta': args.ref_fas,
                                      'fasta.fai': args.ref_index,
                                      'dict': args.ref_dict})
    dbsnp_vcf = batch.read_input_group(**{'vcf': args.dbsnp_vcf,
                                          'vcf.idx': args.dbsnp_vcf_ind})
