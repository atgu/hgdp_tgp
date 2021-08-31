#!/usr/bin/env python3

__author__ = 'Lindo Nkambule'

import hailtop.batch as hb


def scatter_interval_list(b: hb.batch.Batch, interval_list_file: hb.resource.ResourceFile, scatter_count: int = 50,
                          break_bands_at_multiples_of: int = 1000000, scatter_img: str = None, memory: int = 2,
                          out_dir: str = None):
    """
    break the calling interval list into sub-intervals
    :param b: batch
    :param interval_list_file: one or more interval lists
    :param scatter_count: the number of files into which to scatter the resulting list by locus
    :param break_bands_at_multiples_of: if set to a positive value will create a new interval list with the original
    intervals broken up at integer multiples of this value. Set to 0 to NOT break up intervals
    :param scatter_img: image to use for the job
    :param memory: job memory
    :param out_dir: output directory
    :return:
    """

    # break the calling interval list into sub-intervals
    docker_image = scatter_img if scatter_img else 'us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.3-1564508330'

    scatter_list = b.new_job(name='scatter-interval-list')

    scatter_list.image(docker_image)
    scatter_list.cpu(1)
    scatter_list.memory(f'{memory}Gi')
    scatter_list.command('mkdir /scatter_intervals')
    scatter_list.command(f'java -Xms1g -jar /usr/gitc/picard.jar \
      IntervalListTools \
      SCATTER_COUNT={scatter_count} \
      SUBDIVISION_MODE=BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
      UNIQUE=true \
      SORT=true \
      BREAK_BANDS_AT_MULTIPLES_OF={break_bands_at_multiples_of} \
      INPUT={interval_list_file} \
      OUTPUT=/scatter_intervals')
    scatter_list.command('''
    cat > my_script.py <<EOF
import sys
import os
import glob

intervals = sorted(glob.glob('/scatter_intervals/*/*.interval_list'))
for i, interval in enumerate(intervals):
      (directory, filename) = os.path.split(interval)
      newName = os.path.join(directory, str(i + 1) + filename)
      os.rename(interval, newName)
EOF
python3 my_script.py
    ''')
    scatter_list.command(f'mv /scatter_intervals {scatter_list.outfiles}')
    b.write_output(scatter_list.outfiles, f'{out_dir}/scatter-intervals')

    # We return the `scatter_list` Job object that can be used in downstream jobs.
    return scatter_list
