__author__ = 'Lindo Nkambule'

import argparse
import hailtop.batch as hb
import hailtop.fs as hfs
from hailtop.batch.job import Job


def get_file_size(file: str):
    """
    Convert the size from bytes to GiB
    :param file: path to file, str
    :return: file size in GiB
    """

    file_info = hfs.stat(file)   # returns a named tuple
    size_gigs = file_info.size / (1024 * 1024 * 1024)

    return size_gigs


def qc_vcf(
        b: hb.batch.Batch = None,
        input_vcf: hb.ResourceGroup = None,
        females_list: hb.ResourceFile = None,
        chrom: str = None,
        ncpu: int = 16,
        memory: str = 'standard',
        out_dir: str = None,
        storage: int = None,
        img: str = 'docker.io/lindonkambule/gwaspy_phase_impute:latest',
) -> Job:
    """"""
    j = b.new_job(name=f'QC: {chrom}')

    j.image(img)
    j.cpu(ncpu)
    j.memory(memory)
    j._preemptible = False  # use non-preemptibles
    j.regions(['us-central1'])
    j.storage(f'{storage}Gi')

    j.declare_resource_group(
        qced_chunk={
            'bcf': '{root}.bcf',
            'bcf.csi': '{root}.bcf.csi'
        }
    )

    # check number of variants
    j.command(f"""
                echo Initial number of variants before QC
                bcftools query -f '%POS\n' {input_vcf['bcf']} | wc -l
                """
              )

    # HWE, GQ, and ExcessHet filters
    # bcftools HWE and ExcHet are tests, where 1=good and 0=bad
    # similar to UKBB phasing, we compute HWE and ExcHet using females only
    # https://github.com/odelaneau/shapeit5/blob/main/tasks/phasingUKB_200k_release/chrX/script_stats_qc.chrX.sh
    j.command(f"""
                echo 1. Annotating VCF
                bcftools +fill-tags {input_vcf['bcf']} -Ou -o out1.bcf -- -S {females_list} -t HWE,ExcHet
                bcftools +fill-tags out1.bcf -Ou -o out2.bcf -- -t F_MISSING
                rm out1.bcf
                echo 2. Filtering VCF
                bcftools view -Ou --output {j.qced_chunk['bcf']} -i'HWE>=1e-30 && F_MISSING<=0.1 && ExcHet >=0.5 && ExcHet <=1.5' out2.bcf
                rm out2.bcf
                echo 3. Indexing
                bcftools index {j.qced_chunk['bcf']} --output {j.qced_chunk['bcf.csi']} --threads {ncpu}
                echo 4. Number of variants after  QC
                bcftools query -f '%POS\n' {j.qced_chunk['bcf']} | wc -l
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
                     name='hgdp1kg-qc-chrX')

    females = batch.read_input(f'{args.work_dir}/hgdp1kg.females_with_pops.tsv')

    chroms = ['X']
    for i in chroms:
        # vcf_path = f'{args.work_dir}/hgdp1kg_chr{i}.vcf.bgz'
        vcf_path = f'{args.work_dir}/hgdp1kgp_chr{i}.bcf'
        vcf_size = round(get_file_size(vcf_path))
        chrom_vcf = batch.read_input_group(**{'bcf': vcf_path,
                                              'bcf.csi': f'{vcf_path}.csi'})

        qc_vcf(b=batch, input_vcf=chrom_vcf, females_list=females, chrom=f'chr{i}', out_dir=args.work_dir,
               storage=vcf_size*2+30)

    batch.run()


if __name__ == '__main__':
    main()
