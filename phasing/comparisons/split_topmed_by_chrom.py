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


def subset_to_chrom(
        b: hb.batch.Batch = None,
        input_vcf: hb.ResourceGroup = None,
        chrom: str = None,
        ncpu: int = 8,
        memory: str = 'standard',
        outfile: str = None,
        out_dir: str = None,
        storage: int = None,
        img: str = 'docker.io/lindonkambule/shapeit5_2023-03-23_a4a1818:latest',
) -> Job:
    """Subset the TOPMed imputed files by chromosome and convert to BCF"""
    j = b.new_job(name=f'{outfile}: {chrom}')

    j.image(img)
    j.cpu(ncpu)
    j.memory(memory)
    j.regions(['us-central1'])
    j.storage(f'{storage}Gi')

    j.declare_resource_group(
        chrom_subset_outfile={
            'bcf': '{root}.bcf',
            'bcf.csi': '{root}.bcf.csi'
        }
    )

    # split input file by chromosome
    j.command(f"""
                bcftools view {input_vcf['vcf']} --regions {chrom} \
                --output-type b --output {j.chrom_subset_outfile['bcf']} \
                --threads {ncpu-1}
                """
              )

    # index chrom BCF
    j.command(f"""
                echo INDEXING FILTERED FILE
                bcftools index {j.chrom_subset_outfile['bcf']} \
                --output {j.chrom_subset_outfile['bcf.csi']} \
                --threads {ncpu-1}
                """
              )

    b.write_output(j.chrom_subset_outfile,
                   f'{out_dir}/imputation/{outfile}')

    return j


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--billing-project', required=True)
    parser.add_argument('--work-dir', type=str, required=True)

    args = parser.parse_args()

    backend = hb.ServiceBackend(billing_project=args.billing_project,
                                remote_tmpdir=f'{args.work_dir}/tmp/')
    batch = hb.Batch(backend=backend,
                     name='split-topmed-imputed-vcf')

    # compared 5 arrays: H3Africa, GSA, MEGA, PsychChip, Omni2.5. Only worth looking at H3Africa, GSA, and MEGA.
    arrays = ['GSA', 'H3Africa', 'MEGA']

    for arr in arrays:
        vcf_path = f'{args.work_dir}/NeuroGap_30x_{arr}_topmed.dose.vcf.gz'
        vcf_size = round(get_file_size(vcf_path))
        chrom_vcf = batch.read_input_group(**{'vcf': vcf_path,
                                              'vcf.ind': f'{vcf_path}.csi'})

        for i in range(1, 23):
            outfilename = f'NeuroGap_30x_{arr}_chr{i}_topmed.dose'
            subset_to_chrom(b=batch, input_vcf=chrom_vcf, chrom=f'chr{i}', storage=vcf_size+2, outfile=outfilename,
                            out_dir=args.work_dir)

    batch.run()


if __name__ == '__main__':
    main()
