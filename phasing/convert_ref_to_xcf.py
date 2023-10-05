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


def convert_ref_to_xcf(
        b: hb.batch.Batch = None,
        ref: hb.ResourceGroup = None,
        chrom: str = None,
        ncpu: int = 8,
        memory: str = 'standard',
        out_dir: str = None,
        storage: int = None,
        img: str = 'docker.io/lindonkambule/shapeit5_2023-05-05_d6ce1e2:v5.1.1',
) -> Job:
    """Convert reference chromosome files to XCF"""
    j = b.new_job(name=f'convert_xcf: {chrom}')

    j.declare_resource_group(
        converted_chrom={
            'bcf': '{root}.bcf',
            'bcf.csi': '{root}.bcf.csi',
            'bin': '{root}.bin',
            'fam': '{root}.fam'
        }
    )

    j.image(img)
    j.cpu(ncpu)
    j.memory(memory)
    j.regions(['us-central1'])
    j.storage(f'{storage}Gi')

    j.command(f"""
                xcftools_static view -i {ref['bcf']} -o {j.converted_chrom['bcf']} -O sh -r {chrom} -T{ncpu} -m 0.03125
                """
              )

    b.write_output(j.converted_chrom,
                   f'{out_dir}/shapeit5/filtered_phased_SNVs_INDELs_XCF/hgdp1kgp_{chrom}.filtered.SNV_INDEL.phased.shapeit5_xcf')

    return j


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--billing-project', type=str, required=True)
    parser.add_argument('--work-dir', type=str, required=True)

    args = parser.parse_args()

    backend = hb.ServiceBackend(billing_project=args.billing_project,
                                remote_tmpdir=f'{args.work_dir}/tmp/')
    batch = hb.Batch(backend=backend,
                     name='hgdp1kg-convert-to-xcf')

    for i in range(1, 23):
        vcf_path = f'{args.work_dir}/shapeit5/filtered_phased_SNVs_INDELs/hgdp1kgp_chr{i}.filtered.SNV_INDEL.phased.shapeit5.bcf'
        vcf_size = round(get_file_size(vcf_path))
        chrom_vcf = batch.read_input_group(**{'bcf': vcf_path,
                                              'bcf.csi': f'{vcf_path}.csi'})

        convert_ref_to_xcf(b=batch, ref=chrom_vcf, chrom=f'chr{i}',
                           out_dir=args.work_dir, storage=round(vcf_size*2+2))

    batch.run()


if __name__ == '__main__':
    main()
