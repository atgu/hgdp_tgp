__author__ = 'Lindo Nkambule'

import argparse
import os
import hail as hl
import hailtop.batch as hb
import pandas as pd
from hailtop.batch.job import Job
from typing import List, Union


def get_file_size(file):
    """Get file size"""

    file_info = hl.utils.hadoop_stat(file)
    size_bytes = file_info['size_bytes']
    size_gigs = size_bytes / (1024 * 1024 * 1024)

    return size_gigs


def run_phasing_batch(
        output_path: str = None,
        backend: Union[hb.ServiceBackend, hb.LocalBackend] = None,):

    def phase_common(
            b: hb.batch.Batch = None,
            input_vcf: hb.ResourceGroup = None,
            maf: float = 0.001,
            region: str = None,
            genetic_map_file: hb.ResourceFile = None,
            pedigree: hb.ResourceFile = None,
            ncpu: int = 16,
            memory: str = 'highmem',
            storage: int = None,
            img: str = 'docker.io/lindonkambule/shapeit5_2023-03-23_a4a1818:latest',
    ) -> Job:
        """Phase common variants by chunks"""
        j = b.new_job(name=f'phase_common: {region}')

        j.declare_resource_group(
            phased_common_chunk={
                'chunk.bcf': '{root}.bcf',
                'chunk.bcf.csi': '{root}.bcf.csi',
                'chunk.log': '{root}.log'
            }
        )

        j.image(img)
        j.cpu(ncpu)
        j.memory(memory)
        j.regions(['us-central1'])
        j.storage(f'{storage}Gi')

        # phase common variants
        j.command(f"""
                    phase_common_static --input {input_vcf['bcf']} \
                    --map {genetic_map_file} \
                    --output {j.phased_common_chunk['chunk.bcf']} \
                    --thread {ncpu-2} \
                    --log {j.phased_common_chunk['chunk.log']} \
                    --filter-maf {maf} \
                    --region {region} \
                    --pedigree {pedigree}
                    """
                  )

        return j

    def ligate_common_chunks(
            b: hb.batch.Batch,
            common_variants_chunks_list: List[hb.ResourceGroup] = None,
            pedigree: hb.ResourceFile = None,
            output_vcf_name: str = None,
            chrom: str = None,
            memory: str = 'standard',
            out_dir: str = None,
            storage: int = None,
            ncpu: int = 16,
            img: str = 'docker.io/lindonkambule/shapeit5_2023-03-23_a4a1818:latest',
    ) -> Job:
        # requires a VCF/BCF with its index
        j = b.new_job(name=f'ligate_common: {chrom}')

        phased_common_chunks = '\n'.join([f'{v["chunk.bcf"]}' for v in common_variants_chunks_list])

        j.declare_resource_group(
            ligated_chrom={
                'bcf': '{root}.bcf',
                'bcf.csi': '{root}.bcf.csi',
                'log': '{root}.log'
            }
        )

        j.image(img)
        j.cpu(ncpu)
        j.memory(memory)
        j.regions(['us-central1'])
        j.storage(f'{storage}Gi')
        j.command(f'echo "{phased_common_chunks}" > common_chunks_list_ligate.txt')
        j.command(f"""
                    ligate_static --input common_chunks_list_ligate.txt \
                    --pedigree {pedigree} \
                    --output {j.ligated_chrom['bcf']} \
                    --thread {ncpu-2} \
                    --log {j.ligated_chrom['log']} \
                    --index
                    """
                  )

        b.write_output(j.ligated_chrom,
                       f'{out_dir}/phase_common/{output_vcf_name}_{chrom}.shapeit5_common')

        return j

    def phase_rare(
            b: hb.batch.Batch = None,
            input_vcf: hb.ResourceGroup = None,
            scaffold_vcf: hb.ResourceGroup = None,
            scaffold_region: str = None,
            input_region: str = None,
            genetic_map_file: hb.ResourceFile = None,
            pedigree: hb.ResourceFile = None,
            ncpu: int = 16,
            memory: str = 'highmem',
            storage: int = None,
            img: str = 'docker.io/lindonkambule/shapeit5_2023-03-23_a4a1818:latest',
    ) -> Job:
        """Phase common variants by chunks"""
        j = b.new_job(name=f'phase_rare: {scaffold_region}')

        j.declare_resource_group(
            phased_rare_chunk={
                'chunk.bcf': '{root}.bcf',
                'chunk.bcf.csi': '{root}.bcf.csi',
                'chunk.log': '{root}.log'
            }
        )

        j.image(img)
        j.cpu(ncpu)
        j.memory(memory)
        j.regions(['us-central1'])
        j.storage(f'{storage}Gi')

        # phase rare variants
        # --input-region should be column 4 from GLIMPSE chunk output
        # --scaffold-region should be column 3 from GLIMPSE chunk output
        j.command(f"""
                    phase_rare_static \
                    --input {input_vcf['bcf']} --input-region {input_region} \
                    --scaffold {scaffold_vcf['bcf']} --scaffold-region {scaffold_region} \
                    --map {genetic_map_file} \
                    --pedigree {pedigree} \
                    --output {j.phased_rare_chunk['chunk.bcf']} \
                    --thread {ncpu-2} \
                    --log {j.phased_rare_chunk['chunk.log']}
                    """
                  )

        return j

    def concatenate_rare_chunks(
            b: hb.batch.Batch,
            rare_variants_chunks_list: List[hb.ResourceGroup] = None,
            output_vcf_name: str = None,
            chrom: str = None,
            memory: str = 'standard',
            out_dir: str = None,
            storage: int = None,
            ncpu: int = 16,
            img: str = 'docker.io/lindonkambule/shapeit5_2023-03-23_a4a1818:latest',
    ) -> Job:
        # requires a VCF/BCF with its index
        j = b.new_job(name=f'concatenate_rare: {chrom}')

        phased_rare_chunks = '\n'.join([f'{v["chunk.bcf"]}' for v in rare_variants_chunks_list])

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
        j.command(f'echo "{phased_rare_chunks}" > list_concatenate.txt')
        j.command(f"""
                    bcftools concat -n -f list_concatenate.txt -o {j.concatenated_chrom['bcf']}
                    """
                  )

        j.command(f"""
                    bcftools index {j.concatenated_chrom['bcf']} \
                    --output {j.concatenated_chrom['bcf.csi']} \
                    --threads {ncpu-2}
                    """
                  )

        b.write_output(j.concatenated_chrom,
                       f'{out_dir}/phase_rare/{output_vcf_name}_{chrom}.full.shapeit5_rare')

        return j

    batch = hb.Batch(backend=backend,
                     name='shapeit5-phase-hgdp1kg')

    # CHANGE PATH ONCE FAM FILE HAS BEEB FORMATTED
    ped_file = batch.read_input('gs://hgdp-1kg/phasing/hgdp1kg_pedigree.fam')

    for i in range(1, 23):
        # read chrom input files
        vcf_path = f'gs://hgdp-1kg/phasing/filtered_bcfs/hgdp1kg_chr{i}_filtered.bcf'
        chrom_vcf = batch.read_input_group(**{'bcf': vcf_path,
                                              'bcf.csi': f'{vcf_path}.csi'})
        vcf_size = round(get_file_size(vcf_path))

        map_file = batch.read_input(f'gs://hgdp-1kg/phasing/maps/b38/chr{i}.b38.gmap.gz')

        # 1A. Phase common chunks
        common_chunks_file = pd.read_csv(f'gs://hgdp-1kg/phasing/chunks/b38/20cM/chunks_chr{i}.txt',
                                         sep='\t', header=None,
                                         names=['index', 'chrom', 'irg', 'org'])
        common_regions = common_chunks_file['irg'].values.tolist()
        # common_regions = [(irg, org) for irg, org in zip(common_chunks_file['irg'], common_chunks_file['org'])]
        common_chunks_phased = [
            phase_common(
                b=batch,
                input_vcf=chrom_vcf,
                pedigree=ped_file,
                region=common_regions[i],
                genetic_map_file=map_file,
                storage=round(vcf_size*1.5)
            ).phased_common_chunk
            # for i in range(2)
            for i in range(len(common_regions))
        ]

        # 1B. Ligate phased common chunks into one chromosome
        common_phased_ligated_scaffold = ligate_common_chunks(
            b=batch,
            common_variants_chunks_list=common_chunks_phased,
            pedigree=ped_file,
            output_vcf_name='hgdp1kg',
            chrom=f'chr{i}',
            out_dir=output_path,
            storage=round(vcf_size*2)
        ).ligated_chrom

        # 2A. Phase rare chunks
        rare_chunks_file = pd.read_csv(f'gs://hgdp-1kg/phasing/chunks/b38/4cM/chunks_chr{i}.txt',
                                       sep='\t', header=None,
                                       names=['index', 'chrom', 'irg', 'org', 'col5', 'col6', 'col7', 'col8'])
        # irg (3rd col) is SCAFFOLD_REG and org (4th col) is INPUT_REG
        # https://github.com/odelaneau/shapeit5/tree/main/resources/chunks/b38
        rare_regions = [(irg, org) for irg, org in zip(rare_chunks_file['irg'], rare_chunks_file['org'])]
        rare_chunks_phased = [
            phase_rare(
                b=batch,
                input_vcf=chrom_vcf,
                scaffold_vcf=common_phased_ligated_scaffold,
                scaffold_region=rare_regions[i][0],  # irg (3rd col in chunks)
                input_region=rare_regions[i][1],  # org (4th col in chunks)
                genetic_map_file=map_file,
                pedigree=ped_file,
                storage=round(vcf_size*2*1.5)  # we have two input files (unphased VCF+scaffold) and one output
            ).phased_rare_chunk
            # for i in range(2)
            for i in range(len(rare_regions))
        ]

        # 2B. Concatenate phased rare chunks
        concatenated_chunks = concatenate_rare_chunks(
            b=batch,
            rare_variants_chunks_list=rare_chunks_phased,
            output_vcf_name='hgdp1kg',
            chrom=f'chr{i}',
            out_dir=output_path,
            storage=round(vcf_size*2)
        ).ligated_chunk

    batch.run()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--billing-project', required=True)
    parser.add_argument('--out-dir', type=str, default='gs://hgdp-1kg/phasing/shapeit5')

    args = parser.parse_args()

    backend = hb.ServiceBackend(billing_project=args.billing_project,
                                remote_tmpdir=f'{args.out_dir}/tmp/')

    run_phasing_batch(output_path=args.out_dir,
                      backend=backend)


if __name__ == '__main__':
    main()
