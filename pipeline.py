#!/usr/bin/env python3

__author__ = 'Lindo Nkambule'

import argparse
import hailtop.batch as hb
from variant_calling.variant_calling_pipeline import var_call_pipeline


if __name__ == '__main__':
    parser = argparse.ArgumentParser('Variant Calling Pipeline')
    parser.add_argument('--input-bams', required=True)
    parser.add_argument('--ref-fasta',
                        default='gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta')
    parser.add_argument('--ref-index',
                        default='gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai')
    parser.add_argument('--ref-dict',
                        default='gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict')
    parser.add_argument('--billing-project', default='diverse-pop-seq-ref')
    parser.add_argument('--bucket', default='african-seq-data')
    parser.add_argument('--local', action='store_true')
    parser.add_argument('--calling-interval-list',
                        default='gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.interval_list')
    parser.add_argument('--evaluation-interval-list',
                        default=
                        'gs://gcp-public-data--broad-references/hg38/v0/wgs_evaluation_regions.hg38.interval_list')
    parser.add_argument('--dbsnp-vcf',
                        default='gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf')
    parser.add_argument('--dbsnp-vcf-ind',
                        default=
                        'gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx')
    parser.add_argument('--contamination', type=float, default=0.0)
    parser.add_argument('--scatter-count', type=int, default=50)
    parser.add_argument('--run', type=str, default='vc', choices=['all', 'align', 'vc'])
    parser.add_argument('--out-dir', default='gs://african-seq-data')

    args = parser.parse_args()

    if args.local:
        backend = hb.LocalBackend()
    else:
        backend = hb.ServiceBackend(billing_project=args.billing_project,
                                    bucket=args.bucket)

    if args.run == "vc":
        print('Running Variant Calling Only')
        var_call_pipeline(input_bams=args.input_bams, ref_fasta=args.ref_fasta, ref_ind=args.ref_index,
                                  ref_dict=args.ref_dict, calling_int_list=args.calling_interval_list,
                                  evaluation_int_list=args.evaluation_interval_list, dbsnp_vcf=args.dbsnp_vcf,
                                  dbsnp_vcf_ind=args.dbsnp_vcf_ind, contamination=args.contamination,
                                  scatter_count=args.scatter_count, out_dir=args.out_dir, backend=backend)
    elif args.run == 'align':
        print('Running Alignment Only')
    else:
        print('Running both Alignment & Variant Calling')

