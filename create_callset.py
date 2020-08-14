#!/usr/bin/env python3

__author__ = 'konradk'

import argparse
import hail as hl
from gnomad.utils import slack


def get_hgdp_tgp_path(sparse: bool = False):
    return f'gs://hgdp_tgp/output/hgdp_tgp{"_sparse" if sparse else ""}.mt'


def main(args):
    hl.init(default_reference='GRCh38')

    if args.create_sparse_mt:
        path_to_input_list = 'gs://hgdp_tgp/misc/tgp_plus_hgdp_30x_reblocked_gvcfs.txt'

        inputs = []
        with hl.hadoop_open(path_to_input_list, 'r') as f:
            for line in f:
                inputs.append(line.strip())

        temp_bucket = 'gs://gnomad-tmp/hgdp_tgp'
        hl.experimental.run_combiner(inputs, out_file=get_hgdp_tgp_path(sparse=True),
                                     tmp_path=temp_bucket, overwrite=args.overwrite)

    if args.densify_mt:
        mt = hl.read_matrix_table(get_hgdp_tgp_path(sparse=True))
        mt = hl.experimental.densify(mt)
        mt.write(get_hgdp_tgp_path(), args.overwrite)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
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