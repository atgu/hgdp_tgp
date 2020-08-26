#!/usr/bin/env python3

__author__ = 'konradk'

import argparse
import hail as hl
from gnomad.utils import slack
from hail.experimental.vcf_combiner.vcf_combiner import combine_gvcfs


def get_reference_mt_path(dataset: str = 'tgp_hgdp', sparse: bool = False, extension='mt'):
    return f'gs://hgdp_tgp/output/{dataset}{"_sparse" if sparse else ""}.{extension}'


def get_header_path(dataset: str = 'tgp'):
    return f'gs://hgdp_tgp/misc/{dataset}_header.vcf'


def get_samples_path(dataset: str = 'tgp'):
    return f'gs://hgdp_tgp/misc/sample_names_{dataset}.txt'


def get_sample_names_from_list_of_files(input_files, output_fname):
    sample_names_dict = hl.grep('#CHROM', input_files, max_count=100000, show=False)
    sample_names = []
    for fname, lines in sample_names_dict.items():
        sample_names.append('\t'.join([lines[0].strip().split('\t')[-1], fname]))
    with hl.hadoop_open(output_fname, 'w') as f:
        f.write('\n'.join(sample_names))


def get_sample_list_in_order(sample_path, inputs):
    all_samples = {}
    with hl.hadoop_open(sample_path) as f:
        for line in f:
            sample_id, path = line.strip().split('\t')
            all_samples[path] = sample_id
    sample_list = []
    for input_file in inputs:
        input_file = input_file.rsplit('/', 1)[1]
        if input_file not in all_samples:
            raise ValueError(f'{input_file} is not in {sample_path}')
        sample_list.append(all_samples[input_file])
    return sample_list


def main(args):
    hl.init(default_reference='GRCh38')

    hgdp_inputs = []
    tgp_inputs = []
    with hl.hadoop_open('gs://hgdp_tgp/misc/tgp_plus_hgdp_30x_reblocked_gvcfs.txt', 'r') as f:
        for line in f:
            line = line.strip()
            hgdp_inputs.append(line) if 'HGDP' in line else tgp_inputs.append(line)

    temp_bucket = 'gs://gnomad-tmp/tgp_hgdp'

    if args.get_sample_names:
        get_sample_names_from_list_of_files(tgp_inputs, get_samples_path('tgp'))
        get_sample_names_from_list_of_files(hgdp_inputs, get_samples_path('hgdp'))

    if args.create_sparse_mt:
        sample_names = get_sample_list_in_order(get_samples_path('tgp'), tgp_inputs)
        hl.experimental.run_combiner(tgp_inputs, out_file=get_reference_mt_path('tgp', sparse=True),
                                     tmp_path=temp_bucket, overwrite=args.overwrite,
                                     header=get_header_path('tgp'), sample_names=sample_names,
                                     use_genome_default_intervals=True)
        sample_names = get_sample_list_in_order(get_samples_path('hgdp'), hgdp_inputs)
        hl.experimental.run_combiner(hgdp_inputs, out_file=get_reference_mt_path('hgdp', sparse=True),
                                     tmp_path=temp_bucket, overwrite=args.overwrite,
                                     header=get_header_path('hgdp'), sample_names=sample_names,
                                     use_genome_default_intervals=True)
        tgp_mt = hl.read_matrix_table(get_reference_mt_path('tgp', sparse=True))
        tgp_mt = tgp_mt.annotate_entries(
            gvcf_info=tgp_mt.gvcf_info.drop('MQ0', 'VariantType')).drop('AB', 'MQ0')
        hgdp_mt = hl.read_matrix_table(get_reference_mt_path('hgdp', sparse=True))
        hgdp_mt = hgdp_mt.annotate_entries(
            gvcf_info=hgdp_mt.gvcf_info.select(*tgp_mt.gvcf_info))
        mt = combine_gvcfs([tgp_mt, hgdp_mt])
        mt.write(get_reference_mt_path(sparse=True), overwrite=args.overwrite)

    if args.densify_mt:
        mt = hl.read_matrix_table(get_reference_mt_path(sparse=True)).key_rows_by('locus', 'alleles')
        mt = hl.experimental.densify(hl.experimental.sparse_split_multi(mt))
        mt = mt.filter_rows(hl.len(mt.alleles) > 1)
        mt.naive_coalesce(5000).write(get_reference_mt_path(), args.overwrite)

    mt = hl.read_matrix_table(get_reference_mt_path()).drop('gvcf_info')
    hl.export_vcf(mt, get_reference_mt_path(extension='vcf.bgz'), parallel='header_per_shard')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    parser.add_argument('--get_sample_names', help='Get sample names from gVCFs', action='store_true')
    parser.add_argument('--create_sparse_mt', help='Create sparse MT callset', action='store_true')
    parser.add_argument('--densify_mt', help='Densify MT', action='store_true')
    parser.add_argument('--slack_channel', help='Send message to Slack channel/user', default='@konradjk')
    args = parser.parse_args()

    if args.slack_channel:
        from slack_token_pkg.slack_creds import slack_token
        with slack.slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)