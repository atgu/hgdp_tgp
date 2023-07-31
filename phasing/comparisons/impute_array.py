__author__ = 'Lindo Nkambule'

import argparse
import hail as hl
import hailtop.batch as hb
import pandas as pd

from hailtop.batch.job import Job
from typing import List


def get_file_size(file):
    """Get file size"""

    file_info = hl.utils.hadoop_stat(file)
    size_bytes = file_info['size_bytes']
    size_gigs = size_bytes / (1024 * 1024 * 1024)

    return size_gigs


def phase_array_chrom(
        b: hb.batch.Batch = None,
        input_vcf: hb.ResourceGroup = None,
        ref: hb.ResourceGroup = None,
        genetic_map_file: hb.ResourceFile = None,
        header_file: hb.ResourceFile = None,
        chrom: str = None,
        array_name: str = None,
        ncpu: int = 16,
        memory: str = 'standard',
        storage: int = None,
        img: str = 'docker.io/lindonkambule/shapeit5_2023-05-05_d6ce1e2:v5.1.1',
) -> Job:
    """Phase array variants by chromosome"""
    j = b.new_job(name=f'phase {array_name}: {chrom}')

    j.declare_resource_group(
        phased_chrom={
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

    # SHAPEIT5 requires input to have AC and AN tags
    j.command(f"""
                echo FIXING HEADER
                bcftools reheader -h {header_file} {input_vcf['vcf']} > header_fixed.vcf
                bcftools +fill-tags header_fixed.vcf -Ob -o input_annotated.bcf -- -t AN,AC
                bcftools index input_annotated.bcf
                """
              )

    # phase array variants by chrom
    j.command(f"""
                echo PHASING
                phase_common_static --input input_annotated.bcf \
                --reference {ref['vcf']} \
                --map {genetic_map_file} \
                --region {chrom} \
                --output {j.phased_chrom['bcf']} \
                --log {j.phased_chrom['log']} \
                --thread {ncpu-1}
                """
              )

    return j


def imputation(
        b: hb.batch.Batch,
        phased_vcf: hb.ResourceGroup = None,
        ref_vcf: hb.ResourceGroup = None,
        region: str = None,
        genetic_map_file: hb.ResourceFile = None,
        array_name: str = None,
        ncpu: int = 8,
        memory: str = 'highmem',
        storage: int = None,
        img: str = 'docker.io/lindonkambule/gwaspy:v1',
) -> Job:

    j = b.new_job(name=f'impute {array_name}: {region}')

    j.declare_resource_group(
        imputed_chunk={
            'chunk.bcf': '{root}.bcf',
            'chunk.bcf.csi': '{root}.bcf.csi'
        }
    )

    j.cpu(ncpu)
    j.memory(memory)
    j.storage(f'{storage}Gi')
    j.image(img)

    j.command(f"""
                impute5_1.1.5_static \
                --h {ref_vcf['vcf']} \
                --m {genetic_map_file} \
                --g {phased_vcf['bcf']} \
                --r {region} \
                --out-gp-field \
                --o {j.imputed_chunk['chunk.bcf']} \
                --threads {ncpu-1}
                """
              )

    j.command(f"""
                bcftools index {j.imputed_chunk['chunk.bcf']} \
                --output {j.imputed_chunk['chunk.bcf.csi']} \
                --threads {ncpu-1}
                """
              )

    return j


def concatenate_chunks(
        b: hb.batch.Batch,
        chunks_list: List[hb.ResourceGroup] = None,
        output_vcf_name: str = None,
        memory: str = 'standard',
        out_dir: str = None,
        storage: int = None,
        ncpu: int = 4,
        img: str = 'docker.io/lindonkambule/shapeit5_2023-05-05_d6ce1e2:v5.1.1',
) -> Job:
    """Concatenate multiple phased or imputed chunks"""
    # requires a VCF/BCF with its index
    j = b.new_job(name=f'concatenate_imputed: {output_vcf_name}')

    phased_rare_chunks = '\n'.join([f'{v["chunk.bcf"]}' for v in chunks_list])

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
                   f'{out_dir}/comparisons/imputation_accuracy/array/imputation/{output_vcf_name}.imputed')

    return j


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--reference', type=str, required=True)
    parser.add_argument('--billing-project', type=str, required=True)
    parser.add_argument('--work-dir', type=str, required=True)

    args = parser.parse_args()

    backend = hb.ServiceBackend(billing_project=args.billing_project,
                                remote_tmpdir=f'{args.work_dir}/tmp/')
    batch = hb.Batch(backend=backend,
                     name=f'array-imputation-{args.reference}-all-chroms')

    # compared 5 arrays: H3Africa, GSA, MEGA, PsychChip, Omni2.5. Only worth looking at H3Africa, GSA, and MEGA.
    # arrays = ['GSA', 'H3Africa', 'MEGA', 'PsychChip', 'Omni2.5']
    arrays = ['GSA', 'H3Africa', 'MEGA']

    if args.reference == 'hgdp1kgp':
        print('<========== Using HGDP1KGP reference ==========>')
    else:
        print('<========== Using NYGC 1kGP reference ==========>')

    for i in range(1, 23):
        # read input files
        if args.reference == 'hgdp1kgp':
            ref_vcf_p = f'{args.work_dir}/shapeit5/phase_rare/hgdp1kgp_chr{i}.full.shapeit5_rare.bcf'
            ref_in_vcf = batch.read_input_group(**{'vcf': ref_vcf_p,
                                                   'vcf.ind': f'{ref_vcf_p}.csi'})
        else:
            ref_vcf_p = f'{args.work_dir}/comparisons/1kGP_ref_panel/1kGP_high_coverage_Illumina.chr{i}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz'
            ref_in_vcf = batch.read_input_group(**{'vcf': ref_vcf_p,
                                                   'vcf.ind': f'{ref_vcf_p}.tbi'})

        ref_vcf_size = round(get_file_size(ref_vcf_p))

        map_file = batch.read_input(f'{args.work_dir}/maps/b38/chr{i}.b38.gmap.gz')

        for array in arrays:
            array_vcf_p = f'{args.work_dir}/comparisons/array_data/NeuroGap_30x_{array}.vcf.gz'
            array_vcf = batch.read_input_group(**{'vcf': array_vcf_p,
                                                  'vcf.ind': f'{array_vcf_p}.csi'})

            array_vcf_size = round(get_file_size(array_vcf_p))
            array_hdr = batch.read_input(f'{args.work_dir}/comparisons/array_data/NeuroGap_30x_{array}_hdr.txt')

            # 1. Phase array data using a reference
            phased_chrom_array = phase_array_chrom(
                b=batch,
                input_vcf=array_vcf,
                ref=ref_in_vcf,
                genetic_map_file=map_file,
                header_file=array_hdr,
                chrom=f'chr{i}',
                array_name=array,
                storage=round(ref_vcf_size+array_vcf_size+2)
            ).phased_chrom

            imputation_chunks = pd.read_csv(f'{args.work_dir}/chunks/b38/4cM/chunks_chr{i}.txt',
                                            sep='\t', header=None,
                                            names=['index', 'chrom', 'irg', 'org', 'col5', 'col6', 'col7', 'col8'])
            # https://github.com/odelaneau/shapeit5/tree/main/resources/chunks/b38
            imp_chunks_no_buffer = imputation_chunks['org'].tolist()  # chunks with overlapping positions

            # 2. Impute genotypes
            imputed_chunks = [
                imputation(
                    b=batch,
                    phased_vcf=phased_chrom_array,
                    ref_vcf=ref_in_vcf,
                    region=imp_chunks_no_buffer[i],  # 4th column (with no buffer between chunks)
                    genetic_map_file=map_file,
                    array_name=array,
                    storage=round(ref_vcf_size+array_vcf_size+5)
                ).imputed_chunk
                for i in range(len(imp_chunks_no_buffer))
            ]

            # 3. Concatenate imputed chunks
            concatenated_imputed_chunks = concatenate_chunks(
                b=batch,
                chunks_list=imputed_chunks,
                output_vcf_name=f'NeuroGap_30x_{array}_chr{i}_{args.reference}',
                out_dir=args.work_dir,
                storage=round(ref_vcf_size+array_vcf_size+10)
            ).concatenated_chrom

    batch.run()


if __name__ == '__main__':
    main()
