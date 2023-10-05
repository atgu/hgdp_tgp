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


def remove_singletons(
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
    j = b.new_job(name=f'Filter: {chrom}')

    j.image(img)
    j.cpu(ncpu)
    j.memory(memory)
    j.regions(['us-central1'])
    j.storage(f'{storage}Gi')

    j.declare_resource_group(
        qced_chunk={
            'bcf': '{root}.bcf',
            'bcf.csi': '{root}.bcf.csi',
            'stats.txt': '{root}.stats.txt'
        }
    )

    # check number of variants before filtering out singletons (MAC<2)
    j.command(f"""
                echo Initial number of variants before QC
                bcftools query -f '%POS\n' {input_vcf['bcf']} | wc -l
                """
              )

    # HWE, GQ, and ExcessHet filters
    j.command(f"""
                echo 1. QC
                bcftools view -i'MAC>=2' --output {j.qced_chunk['bcf']} {input_vcf['bcf']} 
                echo 2. Indexing
                bcftools index {j.qced_chunk['bcf']} --output {j.qced_chunk['bcf.csi']} --threads {ncpu}
                echo 3. Number of variants after  QC
                bcftools query -f '%POS\n' {j.qced_chunk['bcf']} | wc -l
                echo 3. Writing BCF stats to a file
                bcftools stats {j.qced_chunk['bcf']} > {j.qced_chunk['stats.txt']}
                """)

    b.write_output(j.qced_chunk,
                   f'{out_dir}/shapeit5/filtered_phased_SNVs_INDELs/hgdp1kgp_{chrom}.filtered.SNV_INDEL.phased.shapeit5')

    return j


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--billing-project', type=str, required=True)
    parser.add_argument('--work-dir', type=str, required=True)

    args = parser.parse_args()

    backend = hb.ServiceBackend(billing_project=args.billing_project,
                                remote_tmpdir=f'{args.work_dir}/tmp/')
    batch = hb.Batch(backend=backend,
                     name='hgdp1kg-remove-singletons')

    for i in range(1, 23):
        vcf_path = f'{args.work_dir}/shapeit5/phase_rare/hgdp1kgp_chr{i}.full.shapeit5_rare.bcf'
        vcf_size = round(get_file_size(vcf_path))
        chrom_vcf = batch.read_input_group(**{'bcf': vcf_path,
                                              'bcf.csi': f'{vcf_path}.csi'})

        remove_singletons(b=batch, input_vcf=chrom_vcf, chrom=f'chr{i}',
                          out_dir=args.work_dir, storage=round(vcf_size*0.7+10))

    batch.run()


if __name__ == '__main__':
    main()
