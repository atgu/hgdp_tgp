__author__ = 'Lindo Nkambule'

import argparse
import hailtop.fs as hfs
import hailtop.batch as hb
from hailtop.batch.job import Job
import pandas as pd
from typing import List


def get_file_size(file: str):
    """
    Convert the size from bytes to GiB
    :param file: path to file, str
    :return: file size in GiB
    """

    file_info = hfs.stat(file)   # returns a named tuple
    size_gigs = file_info.size / (1024 * 1024 * 1024)

    return size_gigs


def filter_bcf(
        b: hb.batch.Batch = None,
        input_vcf: hb.ResourceGroup = None,
        males: hb.ResourceFile = None,
        chrom: str = None,
        region: str = None,
        ncpu: int = 16,
        memory: str = 'standard',
        out_dir: str = None,
        storage: int = None,
        img: str = 'docker.io/lindonkambule/gwaspy_phase_impute:latest',
) -> Job:
    """"""
    j = b.new_job(name=f'Filter BCF: {region if region else chrom}')

    j.image(img)
    j.cpu(ncpu)
    j.memory(memory)
    j.regions(['us-central1'])
    j.storage(f'{storage}Gi')

    j.declare_resource_group(
        qced_chunk={
            'bcf': '{root}.bcf',
            'bcf.csi': '{root}.bcf.csi'
        }
    )

    # check number of variants before filtering out singletons
    j.command(f"""
                echo Initial number of variants before QC
                bcftools query -f '%POS\n' {input_vcf['bcf']} | wc -l
                """
              )

    # For males in the non-PAR, there were some variants (7667/4745481) that had missing genotypes but SHAPEIT5 imputed
    # them as hets (1|0). Here we removed these
    # chrX:2781480-155701382
    if chrom == 'chrX_non_par':
        j.command(f"""
                    echo 1. Getting regions where there are no Het males
                    bcftools view {input_vcf['bcf']} -r {region} -Ob -o chunk.bcf
                    bcftools index chunk.bcf
                    bcftools view chunk.bcf -S {males} -r chrX:2781480-155701382 -Ob -o nonpar.bcf
                    bcftools index nonpar.bcf
                    echo variants in non-PAR
                    bcftools query -f '%POS\n' nonpar.bcf | wc -l
                    """
                  )
        j.command(f"""
                    bcftools query -e'GT="het"' -f'%CHROM\t%POS\n' nonpar.bcf > regions.tsv
                    """
                  )
        j.command(f"""
                    echo variants in non-PAR excluding ones where males are imputed as Hets
                    wc -l regions.tsv
                    rm nonpar.bcf*
                    
                    echo 2. QC
                    bcftools view -i'MAC>=2' -R regions.tsv --output {j.qced_chunk['bcf']} chunk.bcf
                    echo 3. Indexing
                    bcftools index {j.qced_chunk['bcf']} --output {j.qced_chunk['bcf.csi']} --threads {ncpu}
                    echo 4. Number of variants after  QC
                    bcftools query -f '%POS\n' {j.qced_chunk['bcf']} | wc -l
                    """)
    else:
        j.command(f"""
                    echo 1. QC
                    bcftools view -i'MAC>=2' --output {j.qced_chunk['bcf']} {input_vcf['bcf']} 
                    echo 2. Indexing
                    bcftools index {j.qced_chunk['bcf']} --output {j.qced_chunk['bcf.csi']} --threads {ncpu}
                    echo 3. Number of variants after  QC
                    bcftools query -f '%POS\n' {j.qced_chunk['bcf']} | wc -l
                    """)

    if out_dir:
        b.write_output(j.qced_chunk,
                       f'{out_dir}/shapeit5/filtered_phased_SNVs_INDELs/hgdp1kgp_{chrom}.filtered.SNV_INDEL.phased.shapeit5')

    return j


def concatenate_filtered_chunks(
        b: hb.batch.Batch,
        filtered_variants_chunks_list: List[hb.ResourceGroup] = None,
        chrom: str = None,
        memory: str = 'standard',
        out_dir: str = None,
        storage: int = None,
        ncpu: int = 4,
        img: str = 'docker.io/lindonkambule/gwaspy_phase_impute:latest',
    ) -> Job:
    j = b.new_job(name=f'concatenate: chrX non-PAR')

    chunks = '\n'.join([f'{v["bcf"]}' for v in filtered_variants_chunks_list])

    j.declare_resource_group(
            concatenated_chrom={
                'bcf': '{root}.bcf',
                'bcf.csi': '{root}.bcf.csi'
            }
    )

    j.image(img)
    j.cpu(ncpu)
    j.memory(memory)
    j.regions(['us-central1'])
    j.storage(f'{storage}Gi')
    j.command(f'echo "{chunks}" > list_concatenate.txt')
    j.command(f"""
                bcftools concat -n -f list_concatenate.txt -o {j.concatenated_chrom['bcf']}
                """
                )

    j.command(f"""
                bcftools index {j.concatenated_chrom['bcf']} \
                --output {j.concatenated_chrom['bcf.csi']} \
                --threads {ncpu}
                """
            )

    b.write_output(j.concatenated_chrom,
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
                     name='hgdp1kg-chrX-filter')

    males = batch.read_input(f'{args.work_dir}/hgdp1kg.males.txt')

    chrx_file = [f'{args.work_dir}/shapeit5/phase_common/hgdp1kgp_chrX_par1.shapeit5_common.bcf',
                 f'{args.work_dir}/shapeit5/phase_common/hgdp1kgp_chrX_par2.shapeit5_common.bcf',
                 f'{args.work_dir}/shapeit5/phase_rare/hgdp1kgp_chrX_non_par.full.shapeit5_rare.bcf']

    names = ['chrX_par1', 'chrX_par2', 'chrX_non_par']

    for chrom, vcf_path in zip(names, chrx_file):
        vcf_size = round(get_file_size(vcf_path))
        chrom_vcf = batch.read_input_group(**{'bcf': vcf_path,
                                              'bcf.csi': f'{vcf_path}.csi'})

        if chrom == 'chrX_non_par':
            chunks_file = pd.read_csv(
                'https://raw.githubusercontent.com/odelaneau/shapeit5/main/resources/chunks/b38/UKB_WGS_200k/small_chunks_6cM/chunks_chrX.txt',
                sep='\t', header=None,
                names=['index', 'chrom', 'irg', 'org', 'col5', 'col6', 'col7', 'col8'])
            regions = chunks_file['org'].tolist()

            filtered_chunks = [
                filter_bcf(b=batch, input_vcf=chrom_vcf, males=males, chrom=chrom, region=regions[i],
                       storage=round(vcf_size*0.7+10)
                           ).qced_chunk
                for i in range(len(regions))
            ]

            concatenate_filtered_chunks(
                b=batch,
                filtered_variants_chunks_list=filtered_chunks,
                chrom=chrom,
                out_dir=args.work_dir,
                storage=round(vcf_size + 10)
            )

        else:
            filter_bcf(b=batch, input_vcf=chrom_vcf, males=males, chrom=chrom,
                       out_dir=args.work_dir, storage=round(vcf_size * 0.7 + 10))
    batch.run()


if __name__ == '__main__':
    main()
