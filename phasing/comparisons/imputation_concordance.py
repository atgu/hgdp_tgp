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


def concordance(
        b: hb.batch.Batch = None,
        imputed_vcf: hb.ResourceGroup = None,
        validation_vcf: hb.ResourceGroup = None,
        gnomad_frequencies_vcf: hb.ResourceGroup = None,
        chrom: str = None,
        variant_type: str = None,
        out_file_name: str = None,
        out_dir: str = None,
        ncpu: int = 8,
        memory: str = 'standard',
        storage: int = None,
        data_type: str = None,
        img: str = 'docker.io/simrub/glimpse:v2.0.0-27-g0919952_20221207',
) -> Job:
    """Compute the genotyping error rate"""
    j = b.new_job(name=f'{variant_type} concordance: {out_file_name}')

    j.image(img)
    j.cpu(ncpu)
    j.memory(memory)
    j.regions(['us-central1'])
    j.storage(f'{storage}Gi')

    j.command(f"""
                echo splitting VCF to {variant_type} only
                bcftools view -v {variant_type} {imputed_vcf['vcf']} -Ob -o imputed_filtered_by_type.bcf --threads {ncpu}
                bcftools index imputed_filtered_by_type.bcf
                """
              )

    # https://github.com/odelaneau/GLIMPSE/blob/master/tutorial/concordance.lst
    j.command(f"""
                echo "{chrom} {gnomad_frequencies_vcf['vcf']} {validation_vcf['vcf']} imputed_filtered_by_type.bcf" > \
                concordance.txt
                """
              )

    if data_type == 'array':
        # cannot use --gt-val with (1) --min-val-dp or (2) --min-val-gl. Filter the validation first
        j.command(f"""
                    GLIMPSE2_concordance --gt-val --input concordance.txt --output concordance_output \
                    --bins 0.00000 0.0002 0.0005 0.001 0.002 0.005 0.01 0.02 0.05 0.10 0.15 0.20 0.3 0.4 0.50 0.60 \
                    --af-tag AF_afr \
                    --thread {ncpu}
                    """
                  )
    else:
        j.command(f"""
                    GLIMPSE2_concordance --input concordance.txt --min-val-dp 8 --output concordance_output \
                    --min-val-gl 0.9999 \
                    --bins 0.00000 0.0002 0.0005 0.001 0.002 0.005 0.01 0.02 0.05 0.10 0.15 0.20 0.3 0.4 0.50 0.60 \
                    --af-tag AF_afr \
                    --thread {ncpu}
                    """
                  )

    j.command(f'mv concordance_output.rsquare.grp.txt.gz {j.ofile}')
    b.write_output(j.ofile,
                   f'{out_dir}/{variant_type}/{out_file_name}.{variant_type}.rsquare.grp.txt.gz')

    return j


def array_concordance(
        b: hb.batch.Batch = None,
        ref: str = None,
        out_dir: str = None,
):
    """Compute the genotyping error rate for imputed array data"""
    # arrays = ['GSA', 'H3Africa', 'MEGA', 'PsychChip', 'Omni2.5']
    arrays = ['GSA', 'H3Africa', 'MEGA']

    for i in range(1, 23):
        # read input files
        # gnomad sites VCF
        gnomad_sites_p = f'gs://gcp-public-data--gnomad/release/3.1/vcf/genomes/gnomad.genomes.v3.1.sites.chr{i}.vcf.bgz'
        gnomad_sites_size = round(get_file_size(gnomad_sites_p))
        gnomad_sites_vcf = b.read_input_group(**{'vcf': gnomad_sites_p,
                                                 'vcf.ind': f'{gnomad_sites_p}.tbi'})

        # validation VCF
        validation_p = f'{out_dir}/NeuroGap_30x_Pilot_Callset.vcf.gz'
        validation_size = round(get_file_size(validation_p))
        validation_vcf = b.read_input_group(**{'vcf': validation_p,
                                               'vcf.ind': f'{validation_p}.tbi'})

        imputed_calls_path = f'{out_dir}/imputation_accuracy/array/filtered_hgdp1kgp/imputation'
        for arr in arrays:
            # imputed VCF
            if ref == 'hgdp1kgp' or ref == '1kgp':
                imputed_p = f'{imputed_calls_path}/NeuroGap_30x_{arr}_chr{i}_{ref}.imputed.bcf'
            elif ref == 'topmed':
                imputed_p = f'{imputed_calls_path}/NeuroGap_30x_{arr}_chr{i}_{ref}.dose.bcf'
            else:
                raise Exception(f'Unknown reference panel {ref} provided')

            imputed_vcf = b.read_input_group(**{'vcf': imputed_p,
                                                'vcf.ind': f'{imputed_p}.csi'})
            imputed_vcf_size = round(get_file_size(imputed_p))

            outfilename = f'NeuroGap_30x_{arr}_chr{i}_{ref}.imputed'

            variants = ['snps', 'indels']

            for var in variants:
                storage_buffer = 25 if var == "snps" else 20
                concordance(b=b,
                            imputed_vcf=imputed_vcf,
                            validation_vcf=validation_vcf,
                            gnomad_frequencies_vcf=gnomad_sites_vcf,
                            chrom=f'chr{i}',
                            variant_type=var,
                            out_file_name=outfilename,
                            out_dir=f'{out_dir}/imputation_accuracy/array/filtered_hgdp1kgp/concordance',
                            storage=round(gnomad_sites_size+validation_size+imputed_vcf_size+storage_buffer),
                            data_type='array')


def lowcov_concordance(
        b: hb.batch.Batch = None,
        ref: str = None,
        out_dir: str = None,
):
    """Compute the genotyping error rate for imputed array data"""
    depths = ['0.5X', '1.0X', '2.0X', '4.0X']

    for i in range(1, 23):
        # read input files
        # gnomad sites VCF
        gnomad_sites_p = f'gs://gcp-public-data--gnomad/release/3.1/vcf/genomes/gnomad.genomes.v3.1.sites.chr{i}.vcf.bgz'
        gnomad_sites_size = round(get_file_size(gnomad_sites_p))
        gnomad_sites_vcf = b.read_input_group(**{'vcf': gnomad_sites_p,
                                                 'vcf.ind': f'{gnomad_sites_p}.tbi'})

        # validation VCF
        validation_p = f'{out_dir}/NeuroGap_30x_Pilot_Callset.vcf.gz'
        validation_size = round(get_file_size(validation_p))
        validation_vcf = b.read_input_group(**{'vcf': validation_p,
                                               'vcf.ind': f'{validation_p}.tbi'})

        for depth in depths:
            # imputed VCF
            imputed_p = f'{out_dir}/imputation_accuracy/coverage/filtered_hgdp1kgp/imputation/{depth}/neurogap_{depth}_chr{i}_{ref}.imputed.bcf'

            imputed_vcf = b.read_input_group(**{'vcf': imputed_p,
                                                'vcf.ind': f'{imputed_p}.csi'})
            imputed_vcf_size = round(get_file_size(imputed_p))

            outfilename = f'neurogap_{depth}_chr{i}_{ref}.imputed'

            variants = ['snps', 'indels']

            for var in variants:
                storage_buffer = 25 if var == "snps" else 20
                concordance(b=b,
                            imputed_vcf=imputed_vcf,
                            validation_vcf=validation_vcf,
                            gnomad_frequencies_vcf=gnomad_sites_vcf,
                            chrom=f'chr{i}',
                            variant_type=var,
                            out_file_name=outfilename,
                            out_dir=f'{out_dir}/imputation_accuracy/coverage/filtered_hgdp1kgp/concordance',
                            storage=round(gnomad_sites_size+validation_size+imputed_vcf_size+storage_buffer),
                            data_type='lowcov')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--data-type', type=str, default='array', choices=['array', 'lowcov'])
    parser.add_argument('--billing-project', required=True)
    parser.add_argument('--work-dir', type=str, required=True)

    args = parser.parse_args()

    backend = hb.ServiceBackend(billing_project=args.billing_project,
                                remote_tmpdir=f'{args.work_dir}/tmp/')
    batch = hb.Batch(backend=backend,
                     name=f'{args.data_type}-imputation-accuracy-hgdp_1kgp_filtered_ref')

    references = ['hgdp1kgp', '1kgp', 'topmed'] if args.data_type == 'array' else ['hgdp1kgp', '1kgp']
    for reference in references:
        if args.data_type == 'array':
            array_concordance(
                b=batch,
                ref=reference,
                out_dir=args.work_dir)
        else:
            lowcov_concordance(
                b=batch,
                ref=reference,
                out_dir=args.work_dir)

    batch.run()


if __name__ == '__main__':
    main()

# 1. array data
# number of jobs = 22 (chroms) * 3 (arrays) * 3 (ref panels) * 2 (1snps & 1indels) = 396

# 2. low-coverage data
# number of jobs = 22 (chroms) * 4 (coverages) * 2 (ref panels) * 2 (1snps & 1indels) = 352
