__author__ = 'Lindo Nkambule'

import argparse
import hailtop.fs as hfs
import hailtop.batch as hb
import pandas as pd
from hailtop.batch.job import Job
from typing import List, Union, Optional


def size(file: str):
    """
    Convert the size from bytes to GiB
    :param file: path to file, str
    :return: file size in GiB
    """

    file_info = hfs.stat(file)   # returns a named tuple
    size_gigs = file_info.size / (1024 * 1024 * 1024)

    return size_gigs


def run_phasing_batch(
        output_path: str = None,
        backend: Union[hb.ServiceBackend, hb.LocalBackend] = None):

    def annotate_vcf(
            b: hb.batch.Batch = None,
            vcf: hb.ResourceGroup = None,
            region: str = None,
            ncpu: int = 8,
            memory: str = 'standard',
            storage: int = None,
            img: str = 'docker.io/lindonkambule/gwaspy_phase_impute:latest',
    ) -> Job:
        j = b.new_job(name=f'Add AC, AN tags: {region}')

        j.image(img)
        j.cpu(ncpu)
        j.memory(memory)
        j.regions(['us-central1'])
        j.storage(f'{storage}Gi')

        j.declare_resource_group(
            annotated_vcf={
                'bcf': '{root}.bcf',
                'bcf.csi': '{root}.bcf.csi'
            }
        )

        j.command(f"""
                    bcftools +fill-tags {vcf['bcf']} -Ou -o {j.annotated_vcf['bcf']} -- -t AN,AC
                    bcftools index {j.annotated_vcf['bcf']} --output {j.annotated_vcf['bcf.csi']} --threads {ncpu}
                    """)

        return j

    def phase_common(
            b: hb.batch.Batch = None,
            input_vcf: hb.ResourceGroup = None,
            maf: float = 0.001,
            region: str = None,
            haploids_file: Optional[hb.ResourceFile] = None,
            pedigree: hb.ResourceFile = None,
            ncpu: int = 8,
            memory: str = 'standard',
            storage: int = None,
            img: str = 'docker.io/lindonkambule/gwaspy_phase_impute:latest',
            out_dir: str = None
    ) -> Job:
        """Phase common variants by chunks"""
        if region == 'par1':
            phase_region = 'chrX:10001-2781479'
            outfilename = 'hgdp1kgp_chrX_par1'
            genetic_map_file = '/root/gwaspy/resources/maps/b38/chrX_par1.b38.gmap.gz'
            region_type = 'PAR1'
        elif region == 'par2':
            phase_region = 'chrX:155701383-156030895'
            outfilename = 'hgdp1kgp_chrX_par2'
            genetic_map_file = '/root/gwaspy/resources/maps/b38/chrX_par2.b38.gmap.gz'
            region_type = 'PAR2'
        else:
            phase_region = region
            outfilename = None
            genetic_map_file = '/root/gwaspy/resources/maps/b38/chrX.b38.gmap.gz'
            region_type = 'Non-PAR'

        j = b.new_job(name=f'phase_common {region_type}: {phase_region}')

        j.declare_resource_group(
            phased_common_chunk={
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

        j.command(f"""
                    phase_common_static --input {input_vcf['bcf']} \
                    --map {genetic_map_file} \
                    --output {j.phased_common_chunk['bcf']} \
                    --thread {ncpu} \
                    --log {j.phased_common_chunk['log']} \
                    --filter-maf {maf} \
                    --region {phase_region} \
                    {f"--haploids {haploids_file}" if haploids_file else ''} \
                    --pedigree {pedigree}
                    """
                  )

        if out_dir:
            b.write_output(j.phased_common_chunk,
                           f'{out_dir}/phase_common/{outfilename}.shapeit5_common')

        return j

    def ligate_common_chunks(
            b: hb.batch.Batch,
            common_variants_chunks_list: List[hb.ResourceGroup] = None,
            pedigree: hb.ResourceFile = None,
            memory: str = 'standard',
            out_dir: str = None,
            storage: int = None,
            ncpu: int = 4,
            img: str = 'docker.io/lindonkambule/gwaspy_phase_impute:latest',
    ) -> Job:
        # requires a VCF/BCF with its index
        j = b.new_job(name=f'ligate_common: chrX non-PAR')

        phased_common_chunks = '\n'.join([f'{v["bcf"]}' for v in common_variants_chunks_list])

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
                    --thread {ncpu} \
                    --log {j.ligated_chrom['log']} \
                    --index
                    """
                  )

        b.write_output(j.ligated_chrom,
                       f'{out_dir}/phase_common/hgdp1kgp_chrX_non_par.shapeit5_common')

        return j

    def phase_rare(
            b: hb.batch.Batch = None,
            input_vcf: hb.ResourceGroup = None,
            scaffold_vcf: hb.ResourceGroup = None,
            scaffold_region: str = None,
            input_region: str = None,
            pedigree: hb.ResourceFile = None,
            haploids_file: hb.ResourceFile = None,
            ncpu: int = 4,
            memory: str = 'standard',
            storage: int = None,
            img: str = 'docker.io/lindonkambule/gwaspy_phase_impute:latest',
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

        j.command(f"""
                    phase_rare_static \
                    --input {input_vcf['bcf']} --input-region {input_region} \
                    --scaffold {scaffold_vcf['bcf']} --scaffold-region {scaffold_region} \
                    --map /root/gwaspy/resources/maps/b38/chrX.b38.gmap.gz \
                    --pedigree {pedigree} \
                    --haploids {haploids_file} \
                    --output {j.phased_rare_chunk['chunk.bcf']} \
                    --thread {ncpu} \
                    --log {j.phased_rare_chunk['chunk.log']}
                    """
                  )

        return j

    def concatenate_rare_chunks(
            b: hb.batch.Batch,
            rare_variants_chunks_list: List[hb.ResourceGroup] = None,
            memory: str = 'standard',
            out_dir: str = None,
            storage: int = None,
            ncpu: int = 4,
            img: str = 'docker.io/lindonkambule/gwaspy_phase_impute:latest',
    ) -> Job:
        # requires a VCF/BCF with its index
        j = b.new_job(name=f'concatenate_rare: chrX non-PAR')

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
                    --threads {ncpu-1}
                    """
                  )

        b.write_output(j.concatenated_chrom,
                       f'{out_dir}/phase_rare/hgdp1kgp_chrX_non_par.full.shapeit5_rare')

        return j

    batch = hb.Batch(backend=backend,
                     name='shapeit5-phase-hgdp1kg-chrX')

    ped_file = batch.read_input(f'{output_path}/hgdp1kg_pedigree.fam')
    males_file = batch.read_input(f'{output_path}/hgdp1kg.males.txt')

    regions = ['par1', 'par2', 'non_par']

    vcf_path = f'{output_path}/qced_bcfs/hgdp1kgp_chrX.bcf'
    chrom_vcf_in = batch.read_input_group(**{'bcf': vcf_path,
                                             'bcf.csi': f'{vcf_path}.csi'})
    vcf_size = round(size(vcf_path))

    for reg in regions:
        # vcf_path = f'{output_path}/qced_bcfs/hgdp1kgp_chrX_{reg}.bcf'
        # chrom_vcf_in = batch.read_input_group(**{'bcf': vcf_path,
        #                                          'bcf.csi': f'{vcf_path}.csi'})
        # vcf_size = round(size(vcf_path))

        chrom_vcf = annotate_vcf(
            b=batch,
            vcf=chrom_vcf_in,
            region=reg,
            storage=round(vcf_size*1.5 + 2)
        ).annotated_vcf

        # A. Phase PAR1 and PAR2 regions without chunking
        if reg == 'par1' or reg == 'par2':
            phase_common(
                b=batch,
                input_vcf=chrom_vcf,
                maf=0.0,
                pedigree=ped_file,
                region=reg,
                storage=round(vcf_size*1.5),
                out_dir=f'{output_path}/shapeit5'
            )

        # B. Phase NON-PAR region in chunks
        else:
            # B1. Phase common chunks
            common_chunks_file = pd.read_csv(
                'https://raw.githubusercontent.com/odelaneau/shapeit5/main/resources/chunks/b38/UKB_WGS_200k/large_chunks_25cM/chunks_chrX.txt',
                sep='\t', header=None,
                names=['index', 'chrom', 'irg', 'org'])
            common_regions = common_chunks_file['irg'].values.tolist()

            common_chunks_phased = [
                phase_common(
                    b=batch,
                    input_vcf=chrom_vcf,
                    maf=0.001,
                    pedigree=ped_file,
                    region=common_regions[i],
                    haploids_file=males_file,
                    storage=round(vcf_size*1.5)
                ).phased_common_chunk
                for i in range(len(common_regions))
            ]

            # B2. Ligate phased common chunks into one chromosome
            common_phased_ligated_scaffold = ligate_common_chunks(
                b=batch,
                common_variants_chunks_list=common_chunks_phased,
                pedigree=ped_file,
                out_dir=f'{output_path}/shapeit5',
                storage=round(vcf_size*0.1)
            ).ligated_chrom

            # B3. Phase rare chunks
            rare_chunks_file = pd.read_csv(
                'https://raw.githubusercontent.com/odelaneau/shapeit5/main/resources/chunks/b38/UKB_WGS_200k/small_chunks_6cM/chunks_chrX.txt',
                sep='\t', header=None,
                names=['index', 'chrom', 'irg', 'org', 'col5', 'col6', 'col7', 'col8'])

            rare_regions = [(irg, org) for irg, org in zip(rare_chunks_file['irg'], rare_chunks_file['org'])]
            rare_chunks_phased = [
                phase_rare(
                    b=batch,
                    input_vcf=chrom_vcf,
                    scaffold_vcf=common_phased_ligated_scaffold,
                    scaffold_region=rare_regions[i][0],  # irg (3rd col in chunks)
                    input_region=rare_regions[i][1],  # org (4th col in chunks)
                    pedigree=ped_file,
                    haploids_file=males_file,
                    storage=round(vcf_size*1.5*1.5)  # we have two input files (unphased VCF+scaffold) and one output
                ).phased_rare_chunk
                for i in range(len(rare_regions))
            ]

            # B4. Concatenate phased rare chunks
            concatenate_rare_chunks(
                b=batch,
                rare_variants_chunks_list=rare_chunks_phased,
                out_dir=f'{output_path}/shapeit5',
                storage=round(vcf_size*0.1)
            )

    batch.run()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--billing-project', type=str, required=True)
    parser.add_argument('--work-dir', type=str, required=True)

    args = parser.parse_args()

    backend = hb.ServiceBackend(billing_project=args.billing_project,
                                remote_tmpdir=f'{args.work_dir}/tmp/')

    run_phasing_batch(output_path=args.work_dir,
                      backend=backend)


if __name__ == '__main__':
    main()
