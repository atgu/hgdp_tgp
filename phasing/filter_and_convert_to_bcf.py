__author__ = 'Lindo Nkambule'

import argparse
import hailtop.batch as hb
from hailtop.batch.job import Job


def filter_vcf_and_convert(
        b: hb.batch.Batch = None,
        input_vcf: hb.ResourceGroup = None,
        samples_to_filter: hb.ResourceFile = None,
        chrom: str = None,
        out_dir: str = None,
        ncpu: int = 16,
        memory: str = 'highmem',
        storage: int = 500,
        img: str = 'docker.io/lindonkambule/shapeit5_2023-03-23_a4a1818:latest',
) -> Job:
    """This function is for removing 29 (24 PCA outliers+5 duplicate) samples and convert to BCF"""
    j = b.new_job(name=f'filter: {chrom}')

    j.declare_resource_group(
        filtered_outfile={
            'bcf': '{root}.bcf',
            'bcf.csi': '{root}.bcf.csi'
        }
    )

    j.image(img)
    j.cpu(ncpu)
    j.memory(memory)
    j.regions(['us-central1'])
    j.storage(f'{storage}Gi')

    # remove 5 duplicates and 24 PCA outliers
    j.command(f"""
                echo REMOVING DUPLICATES AND PCA OUTLIERS
                bcftools view --samples-file ^{samples_to_filter} {input_vcf['vcf']} \
                --output-type b --output {j.filtered_outfile['bcf']} \
                --threads {ncpu-2}
                """
              )

    # index
    j.command(f"""
                echo INDEXING FILTERED FILE
                bcftools index {j.filtered_outfile['bcf']} \
                --output {j.filtered_outfile['bcf.csi']} \
                --threads {ncpu-2}
                """
              )

    # check number of samples (should be 4091)
    j.command(f"""
                echo CHECKING THE NUMBER OF SAMPLES IN FILTERED FILE
                bcftools query -l {j.filtered_outfile['bcf']} | wc -l
                """
              )

    # check number of variants
    j.command(f"""
                echo CHECKING THE NUMBER OF VARIANTS IN FILTERED FILE
                bcftools query -f '%POS\n' {j.filtered_outfile['bcf']} | wc -l
                """
              )

    b.write_output(j.filtered_outfile,
                   f'{out_dir}/hgdp1kg_{chrom}_filtered')

    return j


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--billing-project', required=True)
    parser.add_argument('--out-dir', type=str, default='gs://hgdp-1kg/phasing/filtered_bcfs')

    args = parser.parse_args()

    backend = hb.ServiceBackend(billing_project=args.billing_project,
                                remote_tmpdir=f'{args.out_dir}/tmp/')
    batch = hb.Batch(backend=backend,
                     name='filter-dups-pcaoutliers-hgdp1kg')

    dups_outliers_list = batch.read_input('gs://hgdp-1kg/phasing/5duplicate_and_24pcaoutliers_samples.tsv')

    for i in range(1, 5):
        vcf_path = f'gs://hgdp-1kg/phasing/input_vcfs/hgdp1kg_chr{i}.vcf.bgz'
        chrom_vcf = batch.read_input_group(**{'vcf': vcf_path,
                                              'vcf.tbi': f'{vcf_path}.tbi'})

        filter_vcf_and_convert(b=batch, input_vcf=chrom_vcf, samples_to_filter=dups_outliers_list, chrom=f'chr{i}',
                               out_dir=args.out_dir)

    batch.run()


if __name__ == '__main__':
    main()
