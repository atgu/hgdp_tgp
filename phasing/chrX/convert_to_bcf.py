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


def convert_to_bcf(
        b: hb.batch.Batch = None,
        input_vcf: hb.ResourceGroup = None,
        chrom: str = None,
        ncpu: int = 16,
        memory: str = 'standard',
        out_dir: str = None,
        storage: int = None,
        img: str = 'docker.io/lindonkambule/gwaspy_phase_impute:latest'
) -> Job:
    j = b.new_job(name=f'Convert to BCF: {chrom}')

    j.image(img)
    j.cpu(ncpu)
    j.memory(memory)
    j._preemptible = False  # use non-preemptibles
    j.regions(['us-central1'])
    j.storage(f'{storage}Gi')

    j.declare_resource_group(
        output_bcf={
            'bcf': '{root}.bcf',
            'bcf.csi': '{root}.bcf.csi'
        }
    )

    j.command(f"""
                echo 1. Removing annotations
                bcftools annotate -x INFO,^FORMAT/GT {input_vcf['bcf']} -Ob --output {j.output_bcf['bcf']} --threads {ncpu}
                echo 2. Indexing
                bcftools index {j.output_bcf['bcf']} --output {j.output_bcf['bcf.csi']} --threads {ncpu} --force
                """)
    b.write_output(j.output_bcf,
                    f'{out_dir}/debugging/hgdp1kgp_{chrom}_pruned_annotations')

    return j


def fixploidy(
        b: hb.batch.Batch = None,
        input_vcf: hb.ResourceGroup = None,
        chrom: str = None,
        ncpu: int = 16,
        memory: str = 'standard',
        out_dir: str = None,
        storage: int = None,
        img: str = 'docker.io/lindonkambule/gwaspy_phase_impute:latest',
        fix_ploidy: str = 'updated',
        sample_info_file: hb.ResourceFile = None,
        ploidy_file: hb.ResourceFile = None
) -> Job:
    fixploidy_ = 'Updated fixploidy' if fix_ploidy=='updated' else 'Old fixploidy'
    j = b.new_job(name=f'{fixploidy_}: {chrom}')

    j.image(img)
    j.cpu(ncpu)
    j.memory(memory)
    j._preemptible = False  # use non-preemptibles
    j.regions(['us-central1'])
    j.storage(f'{storage}Gi')

    j.declare_resource_group(
        output_bcf={
            'bcf': '{root}.bcf',
            'bcf.csi': '{root}.bcf.csi'
        }
    )

    if fix_ploidy=='old':
        j.command(f"""
                    echo 1. Fixing Ploidy
                    bcftools +fixploidy {input_vcf['bcf']} > fixed.vcf
                    echo 2. Converting
                    bcftools +fill-tags fixed.vcf -Ob -o {j.output_bcf['bcf']} -- -t MAF
                    echo 3. Indexing
                    bcftools index {j.output_bcf['bcf']} --output {j.output_bcf['bcf.csi']} --threads {ncpu}
                    """)
        b.write_output(j.output_bcf,
                       f'{out_dir}/debugging/hgdp1kgp_{chrom}_fixploidy_old')
    else:
        j.command(f"""
                    echo 1. Fixing Ploidy
                    bcftools +fixploidy {input_vcf['bcf']} -- -d 2 -p {ploidy_file} -s {sample_info_file} > fixed.vcf
                    echo 2. Converting
                    bcftools +fill-tags fixed.vcf -Ob -o {j.output_bcf['bcf']} -- -t MAF
                    echo 3. Indexing
                    bcftools index {j.output_bcf['bcf']} --output {j.output_bcf['bcf.csi']} --threads {ncpu}
                    """)

        b.write_output(j.output_bcf,
                       f'{out_dir}/debugging/hgdp1kgp_{chrom}_fixploidy_updated')

    return j


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--billing-project', type=str, required=True)
    parser.add_argument('--work-dir', type=str, required=True)

    args = parser.parse_args()

    backend = hb.ServiceBackend(billing_project=args.billing_project,
                                remote_tmpdir=f'{args.work_dir}/tmp/')
    batch = hb.Batch(backend=backend,
                     name='hgdp1kg-chrX-convert2bcf-debugging')
    sample_sex_info_file = batch.read_input(f'{args.work_dir}/hgdp1kg.with.sex.txt')
    ploidy_file = batch.read_input(f'{args.work_dir}/ploidy.txt')

    # prune annotations we do not need from VCF to reduce file size
    vcf_path = f'{args.work_dir}/hgdp1kg_chrX.vcf.bgz'
    vcf_size = round(get_file_size(vcf_path))
    chrom_vcf = batch.read_input_group(**{'bcf': vcf_path,
                                          'bcf.csi': f'{vcf_path}.tbi'})
    converted_bcf = convert_to_bcf(b=batch, input_vcf=chrom_vcf, chrom=f'chrX',
                                   out_dir=args.work_dir, storage=int(vcf_size*2+20)).output_bcf

    ploidy = ['old', 'updated']
    for val in ploidy:
        fixploidy(b=batch, input_vcf=converted_bcf, chrom=f'chrX',
                       out_dir=args.work_dir, storage=int(vcf_size*2+20), fix_ploidy=val,
                       sample_info_file=sample_sex_info_file, ploidy_file=ploidy_file)

    batch.run()


if __name__ == '__main__':
    main()

