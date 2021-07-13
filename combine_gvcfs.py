__author__ = 'Lindo Nkambule'

import argparse
import hail as hl


def combine_gvcfs(gvcf_input_list: str, out_path: str = 'gs://african-seq-data/gambian-genomes/',
                  out_mt_name: str = 'gambian_genomes_merged_gvcfs', temp_bucket: str = 'gs://african-seq-data',
                  reference: str = 'GRCh38', use_genome_default_intervals: bool = True, overwrite: bool = True,
                  key_by_locus_and_alleles: bool = True):
    """
    Combine single-sample GVCFs into a multi-sample matrix table

    :param gvcf_input_list: MT to filter
    :param out_path: path to where multi-sample MT will be written
    :param out_mt_name: name to use for output MT
    :param temp_bucket: bucket for storing intermediate files
    :param reference: reference genome to use
    :param use_genome_default_intervals: import GVCFs with uniform partition intervals of default size for whole-genome
    data
    :param overwrite: overwrite MT if it exists
    :param key_by_locus_and_alleles: key by both locus and alleles in the final output
    """

    inputs = []
    with hl.hadoop_open(gvcf_input_list, 'r') as f:
        for line in f:
            inputs.append(line.strip())

    output_file = f'{out_path}{out_mt_name}.mt'  # output destination

    hl.experimental.run_combiner(
        inputs, out_file=output_file,
        tmp_path=temp_bucket,
        reference_genome=reference,
        use_genome_default_intervals=use_genome_default_intervals,
        overwrite=overwrite,
        key_by_locus_and_alleles=key_by_locus_and_alleles)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--gvcfs-list', type=str, required=True)
    parser.add_argument('--out-dir', type=str, default='gs://african-seq-data/gambian-genomes/')
    parser.add_argument('--out-name', default='gambian_genomes_merged_gvcfs')
    parser.add_argument('--tmp-bucket', default='gs://african-seq-data')
    parser.add_argument('--reference', type=str, default='GRCh38')
    parser.add_argument('--overwrite', action='store_true')

    args = parser.parse_args()

    print("Combining GVCFs into a multi-sample matrix table")
    combine_gvcfs(gvcf_input_list=args.gvcfs_list, out_path=args.out_dir, out_mt_name=args.out_name,
                  temp_bucket=args.tmp_bucket, reference=args.reference, overwrite=args.overwite)

    print("Done!")


if __name__ == '__main__':
    main()
