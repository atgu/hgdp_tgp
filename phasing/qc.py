__author__ = 'Lindo Nkambule'

import argparse
import hail as hl
import hailtop.batch as hb
from hailtop.batch.job import Job


def get_file_size(file):
    """Get file size"""

    file_info = hl.utils.hadoop_stat(file)
    size_bytes = file_info['size_bytes']
    size_gigs = size_bytes / (1024 * 1024 * 1024)

    return size_gigs


def qc_vcf(
        b: hb.batch.Batch = None,
        input_vcf: hb.ResourceGroup = None,
        chrom: str = None,
        ncpu: int = 16,
        memory: str = 'standard',
        out_dir: str = None,
        storage: int = None,
        img: str = 'docker.io/lindonkambule/shapeit5_2023-03-23_a4a1818:latest',
) -> Job:
    """"""
    j = b.new_job(name=f'QC: {chrom}')

    j.image(img)
    j.cpu(ncpu)
    j.memory(memory)
    j.regions(['us-central1'])
    j.storage(f'{storage}Gi')

    j.declare_resource_group(
        qced_chunk={
            'qced.bcf': '{root}.bcf',
            'qced.bcf.csi': '{root}.bcf.csi'
        }
    )

    # check number of variants
    j.command(f"""
                echo Initial number of variants before QC
                bcftools query -f '%POS\n' {input_vcf['bcf']} | wc -l
                """
              )

    # HWE, GQ, and ExcessHet filters
    j.command(f"""
                echo 1. QC
                bcftools +fill-tags {input_vcf['bcf']} -Ou -- -t all | bcftools view -i'HWE>=1e-30 && F_MISSING<=0.1 && ExcHet >=0.5 && ExcHet <=1.5' --output {j.qced_chunk['qced.bcf']}
                echo 2. Indexing
                bcftools index {j.qced_chunk['qced.bcf']} --output {j.qced_chunk['qced.bcf.csi']} --threads {ncpu-2}
                echo 3. Number of variants after  QC
                bcftools query -f '%POS\n' {j.qced_chunk['qced.bcf']} | wc -l
                """)

    b.write_output(j.qced_chunk,
                   f'{out_dir}/qced_bcfs/hgdp1kgp_{chrom}')

    return j


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--billing-project', type=str, required=True)
    parser.add_argument('--work-dir', type=str, required=True)

    args = parser.parse_args()

    backend = hb.ServiceBackend(billing_project=args.billing_project,
                                remote_tmpdir=f'{args.work_dir}/tmp/')
    batch = hb.Batch(backend=backend,
                     name='hgdp1kg-qc')

    # we have already qc'ed chr2
    # chroms = [i for i in range(1, 23) if i != 2]
    # for i in chroms:
    for i in range(1, 23):
        vcf_path = f'{args.work_dir}/filtered_bcfs/hgdp1kg_chr{i}_filtered.bcf'
        vcf_size = round(get_file_size(vcf_path))
        chrom_vcf = batch.read_input_group(**{'bcf': vcf_path,
                                              'bcf.csi': f'{vcf_path}.csi'})

        qc_vcf(b=batch, input_vcf=chrom_vcf, chrom=f'chr{i}', out_dir=args.work_dir, storage=vcf_size*2+10)

    batch.run()


if __name__ == '__main__':
    main()
