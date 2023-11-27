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


def phase_impute(
        b: hb.batch.Batch = None,
        input_vcf: hb.ResourceGroup = None,
        reference: hb.ResourceGroup = None,
        ref_name: str = None,
        region: str = None,
        output_region: str = None,
        genetic_map_file: hb.ResourceFile = None,
        coverage: str = None,
        ncpu: int = 8,
        memory: str = 'standard',
        storage: int = None,
        img: str = 'docker.io/simrub/glimpse:v1.1.1-c27e90d_20210521',
) -> Job:
    # requires a VCF/BCF with its index
    j = b.new_job(name=f'{ref_name} impute {coverage}: {output_region}')

    j.declare_resource_group(
        imputed_chunk={
            'chunk.bcf': '{root}.bcf',
            'chunk.log': '{root}.log'
        }
    )

    j.image(img)
    j.cpu(ncpu)
    j.memory(memory)
    j.regions(['us-central1'])
    j.storage(f'{storage}Gi')

    j.command(f"""
                GLIMPSE_phase_v1.1.1 --input {input_vcf['vcf']} \
                --reference {reference['vcf']} \
                --thread {ncpu} \
                --input-region {region} \
                --map {genetic_map_file} \
                --output-region {output_region} \
                --output {j.imputed_chunk['chunk.bcf']} \
                --log {j.imputed_chunk['chunk.log']}
                """
              )

    return j


def ligate_chunks(
        b: hb.batch.Batch,
        gvcf_list: List[hb.ResourceGroup] = None,
        output_vcf_name: str = None,
        memory: str = 'lowmem',
        coverage: str = None,
        out_dir: str = None,
        storage: int = None,
        ncpu: int = 16,
        img: str = 'docker.io/simrub/glimpse:v1.1.1-c27e90d_20210521',
) -> Job:
    # requires a VCF/BCF with its index
    j = b.new_job(name=f'ligate: {output_vcf_name}')

    imp_chunks = '\n'.join([f'{v["chunk.bcf"]}' for v in gvcf_list])

    j.declare_resource_group(
        ligated_chunk={
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
    j.command(f'echo "{imp_chunks}" > list.imputed.txt')
    j.command(f"""
                GLIMPSE_ligate_v1.1.1 --input list.imputed.txt \
                --output {j.ligated_chunk['bcf']} \
                --log {j.ligated_chunk['log']}
                """
              )

    # index
    j.command(f"""
                bcftools index --force {j.ligated_chunk['bcf']} \
                --output {j.ligated_chunk['bcf.csi']} \
                --threads {ncpu-2}
                """
              )

    b.write_output(j.ligated_chunk,
                   f'{out_dir}/comparisons/imputation_accuracy/coverage/filtered_hgdp1kgp/imputation/{coverage}/{output_vcf_name}.imputed')

    return j


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--reference', type=str, required=True)
    parser.add_argument('--billing-project', required=True)
    parser.add_argument('--work-dir', type=str, required=True)

    args = parser.parse_args()

    backend = hb.ServiceBackend(billing_project=args.billing_project,
                                remote_tmpdir=f'{args.work_dir}/tmp/')
    batch = hb.Batch(backend=backend,
                     name=f'lowcov-imputation-hgdp_1kgp_filtered_ref-all-chroms')

    # compared N depths of coverage (0.5X, 1X, 2X, 4X, 6X, 10X, 20X).
    # Only worth looking at 0.5X, 1X, 2X, and 4X here since we won't be sequencing higher depths.
    depths = ['0.5X', '1.0X', '2.0X', '4.0X']

    if args.reference == 'hgdp1kgp':
        print('<========== Using HGDP1KGP reference ==========>')
    else:
        print('<========== Using NYGC 1kGP reference ==========>')

    for i in range(1, 23):
        # read input files
        if args.reference == 'hgdp1kgp':
            # ref_vcf_p = f'{args.work_dir}/shapeit5/phase_rare/hgdp1kgp_chr{i}.full.shapeit5_rare.bcf'
            ref_vcf_p = f'{args.work_dir}/shapeit5/filtered_phased_SNVs_INDELs/hgdp1kgp_chr{i}.filtered.SNV_INDEL.phased.shapeit5.bcf'
            ref_in_vcf = batch.read_input_group(**{'vcf': ref_vcf_p,
                                                   'vcf.ind': f'{ref_vcf_p}.csi'})
        else:
            ref_vcf_p = f'{args.work_dir}/comparisons/1kGP_ref_panel/1kGP_high_coverage_Illumina.chr{i}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz'
            ref_in_vcf = batch.read_input_group(**{'vcf': ref_vcf_p,
                                                   'vcf.ind': f'{ref_vcf_p}.tbi'})

        ref_vcf_size = round(get_file_size(ref_vcf_p))

        map_file = batch.read_input(f'{args.work_dir}/maps/b38/chr{i}.b38.gmap.gz')

        for depth in depths:
            cov_vcf_p = f'{args.work_dir}/comparisons/lowcov_data/neurogap_{depth}.vcf.gz'
            cov_vcf = batch.read_input_group(**{'vcf': cov_vcf_p,
                                                'vcf.ind': f'{cov_vcf_p}.tbi'})

            cov_vcf_size = round(get_file_size(cov_vcf_p))
            chunks_file = pd.read_csv(f'{args.work_dir}/chunks/b38/4cM/chunks_chr{i}.txt',
                                      sep='\t', header=None,
                                      names=['index', 'chrom', 'irg', 'org', 'col5', 'col6', 'col7', 'col8'])
            # irg (3rd col) is SCAFFOLD_REG and org (4th col) is INPUT_REG
            # https://github.com/odelaneau/shapeit5/tree/main/resources/chunks/b38
            regions = [(irg, org) for irg, org in zip(chunks_file['irg'], chunks_file['org'])]

            imputed_chunks = [
                phase_impute(
                    b=batch,
                    input_vcf=cov_vcf,
                    reference=ref_in_vcf,
                    ref_name=args.reference,
                    region=regions[i][0],
                    output_region=regions[i][1],
                    genetic_map_file=map_file,
                    coverage=depth,
                    storage=round(ref_vcf_size + cov_vcf_size + 2)
                ).imputed_chunk
                for i in range(len(regions))
            ]

            # merge phased/imputed chunks
            outfilename = f'neurogap_{depth}_chr{i}_{args.reference}'
            merged = ligate_chunks(
                b=batch,
                gvcf_list=imputed_chunks,
                output_vcf_name=outfilename,
                coverage=depth,
                out_dir=args.work_dir,
                storage=round(cov_vcf_size*0.1)
            ).ligated_chunk

    batch.run()


if __name__ == '__main__':
    main()

# For imputation using 1kGP, Hail Batch randomly failed at the ligation step for two chromosomes:
# 1.0X chr5
# 0.5X chr6

# had to re-run imputation for the 2 chromosomes (at the two coverages)