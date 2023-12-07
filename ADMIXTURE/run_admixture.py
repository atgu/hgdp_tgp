__author__ = 'Zan Koenig'

import hailtop.batch as hb
import hail as hl
import argparse
import re
import subprocess
import os
import math


def run_admixture(#b: hb.batch.Batch,
        b,
        infile,
        k,
        run,
        out_dir,
        cpu: float = 1,
        memory: str = 'standard',
        img: str = 'docker.io/zkoenig/admixture:v3'):
    plink_files = b.read_input_group(bed=f'{infile}.bed',
                           bim=f'{infile}.bim',
                           fam=f'{infile}.fam')
                           
    admixture = b.new_job(name=f'admixture-{infile}-{k}')
    admixture.cpu(cpu)
    admixture.memory(memory)
    # calculating storage from the estimated data sizes, the *2 is to give some buffer room
    admixture.storage('10Gi')
    admixture.image(img)
    output_file_name = infile.split('/')[-1]

    # on first run, run cross validation
    if run == 1:
        cmd = f'admixture --cv=5 {plink_files.bed} {k}'
    else:
        cmd = f'admixture {plink_files.bed} {k}'
    admixture.command(cmd)
    admixture.command(f'mv {output_file_name}.{k}.P {admixture.pfile}')
    admixture.command(f'mv {output_file_name}.{k}.Q {admixture.qfile}')
    # writing each p and q file into a folder, will result in one folder/run
    b.write_output(admixture.pfile, f'{out_dir}/run{run}/{output_file_name}.{k}.P')
    b.write_output(admixture.qfile, f'{out_dir}/run{run}/{output_file_name}.{k}.Q')


def main(args):
    # setting the billing project and bucket where the output data will be stored
    backend = hb.ServiceBackend(billing_project='diverse-pop-seq-ref',
                                bucket='hgdp-1kg')
    # setting the name of the batch job
    b = hb.Batch(backend=backend, name='admixture')

    # specifying the number of k, (2,11) == up to 10 k
    for k in range(2, 11):
        # specifying how many runs will happen for each k (1,11) == 10 runs
        for run in range(1, 11):
            run_admixture(b, args.bed_file, k, run, args.output_path)

    b.run()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    # input file path, in plink format (bed, bim, fam)
    parser.add_argument('--bed_file',
                        default='gs://hgdp-1kg/hgdp_tgp/qc_and_figure_generation/admixture_plink_dataset/post_qc_ld_unrel_only')
    # output file path, main directory which will contain run folders
    parser.add_argument('--output_path', default='gs://hgdp-1kg/hgdp_tgp/qc_and_figure_generation/admixture_results')

    args = parser.parse_args()

    main(args)
