#!/usr/bin/env python3

__author__ = 'Lindo Nkambule'

import hail as hl


def bytes_to_gb(in_file: str):
    """
    Convert the size from bytes to GiB
    :param in_file: path to file, str
    :return: file size in GiB
    """

    file_info = hl.utils.hadoop_stat(in_file)
    size_bytes = file_info['size_bytes']
    size_gigs = size_bytes / (1024 * 1024 * 1024)

    return size_gigs
