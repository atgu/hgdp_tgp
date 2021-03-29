#!/usr/bin/env python3

__author__ = 'Lindo Nkambule'

import argparse
import hailtop.batch as hb
import hail as hl
import os
import pandas as pd


def bytes_to_gb(file):
    """ Convert the size from bytes to GB"""

    file_info = hl.utils.hadoop_stat(file)
    size_bytes = file_info['size_bytes']
    size_gigs = size_bytes / (1024 * 1024 * 1024)

    return size_gigs


def cram_to_bam(b: hb.batch.Batch, input_cram_file: str = None, ref_fasta: str = None, ref_dict: str = None,
                ref_ind: str = None, bam_out_name: str = None, memory: int = 15, samtools_image: str = None,
                out_dir: str = None):
    docker_image = samtools_image if samtools_image else 'gcr.io/genomics-tools/samtools'

    out_bam_name = bam_out_name + '.bam'

    output_bam_size: float = bytes_to_gb(input_cram_file) / 0.40
    ref_size: float = bytes_to_gb(ref_fasta) + bytes_to_gb(ref_ind)
    disk_size: int = round(bytes_to_gb(input_cram_file) + output_bam_size + ref_size) + 25

    job_memory = str(memory) + 'Gi'
    job_storage = str(disk_size) + 'Gi'

    crams_to_bams = b.new_job(name=out_bam_name)
    in_cram = b.read_input(input_cram_file)
    fasta = b.read_input_group(**{'fasta': ref_fasta,
                                  'fasta.fai': ref_ind,
                                  'dict': ref_dict})

    crams_to_bams.memory(job_memory)
    crams_to_bams.image(docker_image)
    crams_to_bams.storage(job_storage)
    crams_to_bams.command(f'samtools view -b -T {fasta.fasta} -o {out_bam_name} {in_cram}')
    crams_to_bams.command(f'samtools index {out_bam_name}')
    crams_to_bams.command(f'mv {out_bam_name} {crams_to_bams.bamout}')
    crams_to_bams.command(f'mv {out_bam_name}.bai {crams_to_bams.bamind}')
    b.write_output(crams_to_bams.bamout, f'{out_dir}/BAMS/{out_bam_name}')
    b.write_output(crams_to_bams.bamind, f'{out_dir}/BAMS/{out_bam_name}.bai')

    return crams_to_bams


if __name__ == '__main__':
    parser = argparse.ArgumentParser('CRAMs to BAMs')
    parser.add_argument('--input-crams', required=True)
    parser.add_argument('--ref-fasta',
                        default='gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta')
    parser.add_argument('--ref-index',
                        default='gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai')
    parser.add_argument('--ref-dict',
                        default='gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict')
    parser.add_argument('--billing-project', default='diverse-pop-seq-ref')
    parser.add_argument('--bucket', default='african-seq-data')
    parser.add_argument('--local', action='store_true')
    parser.add_argument('--out-dir', default='gs://african-seq-data')

    args = parser.parse_args()

    if args.local:
        backend = hb.LocalBackend()
    else:
        backend = hb.ServiceBackend(billing_project=args.billing_project,
                                    bucket=args.bucket)

    crams_to_bams_job = hb.Batch(backend=backend, name='cram-to-bam')

    crams_paths = pd.read_csv(args.input_crams, sep='\t', header=None)
    crams_paths.columns = ['cram']
    cram_files = []
    for cram in crams_paths.index:
        cram_files.append(crams_paths['cram'][cram])

    for cram in cram_files:
        base_cram = os.path.basename(cram)  # get file name without path
        cram_no_ext = os.path.splitext(base_cram)[0]
        convert_cram_to_bam = cram_to_bam(b=crams_to_bams_job, input_cram_file=cram, ref_fasta=args.ref_fasta,
                                          ref_dict=args.ref_dict, ref_ind=args.ref_index,
                                          bam_out_name=cram_no_ext, out_dir=args.out_dir)

    crams_to_bams_job.run()


# input-crams should be a text file, each line with a full path to a CRAM file. Example below
# gs://african-seq-data/gambian-genomes/SC_GMFUL5306338.alt_bwamem_GRCh38DH.20151208.FULA.gambian_lowcov.cram
# gs://african-seq-data/gambian-genomes/SC_GMFUL5306339.alt_bwamem_GRCh38DH.20151208.FULA.gambian_lowcov.cram
# e.g. python3 cram_to_bam.py --input-crams gs://african-seq-data/crams_to_convert.txt
