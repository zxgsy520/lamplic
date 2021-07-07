#!/usr/bin/env python
#coding:utf-8

import os
import re
import sys
import logging
import argparse

LOG = logging.getLogger(__name__)

__version__ = "1.1.0"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []


def read_cluster(file):

    clus_id = ""
    repr_id = ""
    otus = []

    for line in open(file):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith(">"):
            if len(otus)>=1 and repr_id!='':
                yield clus_id, repr_id, otus
                clus_id = repr_id = ""
                otus = []
            clus_id = line.strip('>').replace(' ', '')
            continue
        line = line.split('>')[-1]

        if '%' not in line:
            repr_id = line.split('...')[0]
        otus.append(line.split('...')[0])
    if len(otus)>=1 and repr_id!='':
        yield clus_id, repr_id, otus


def recover_sample(otus):

    samples = []

    for i in otus:
        sample, seq_id = i.split('_', 1)
        if sample in samples:
            continue
        samples.append(sample)

    return samples


def sort_cluster(file):

    data = {}
    samples = []

    for clus_id, repr_id, otus in read_cluster(file):
        data[clus_id] = [len(otus), repr_id, otus]
        samples += recover_sample(otus)
    samples = sorted(list(set(samples)))

    return data, samples


def stat_otu(otus):

    data = {}

    for i in otus:
        sample, seq_id = i.split('_', 1)
        if sample not in data:
            data[sample] = 0
        data[sample] += 1

    return data


def sort_otus(samples, otu_dict):

    otus = []

    for i in samples:
        otu = 0
        if i in otu_dict:
            otu = otu_dict[i]
        otus.append('%.2f' % otu)

    return otus


def cluster2otu(file, min_clus=3):

    fm = open('otu_map.tsv', 'w')
    fo = open('otu_tab.tsv', 'w')
    data, samples = sort_cluster(file)

    fo.write('#OTU ID\t%s\n' % '\t'.join(samples))
    print('#Clustered Number\tOTU ID\tSeed Seq\tCluster ID')
    n = 0
    for line in sorted(data.items(), key=lambda x:x[1][0], reverse=True):
        n += 1
        clus_id = line[0]
        otu_id = "OTU_%s" % n
        line = line[1]
        otus = sorted(line[-1])

        fm.write('%s\t%s\n' % (otu_id, '\t'.join(otus)))
        print('%s\t%s\t%s\t%s' % (line[0], otu_id, line[1], clus_id))

        if min_clus > line[0]:
            continue
        otus = stat_otu(otus)
        otus = sort_otus(samples, otus)
        fo.write('%s\t%s\n' % (otu_id, '\t'.join(otus)))

    fm.close()
    fo.close()


def add_hlep_args(parser):

    parser.add_argument('input', metavar='FILE', type=str,
        help='Input the cluster result of cd-hit.')
    parser.add_argument('-mc', '--min_clus', metavar='INT', type=int, default=3,
        help='Set the minimum number of clustering sequences to filter, default=3.')

    return parser


def main():

    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
    description='''
name:
    cluster2otu.py The cluster result of cd-hit can be processed by otu.

attention:
    cluster2otu id.fa.clstr >otu_id.list
Output file：
    otu_map.tsv #sequence id corresponding to ont.
    otu_tab.tsv #OTU abundance distribution table.
version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    cluster2otu(args.input, args.min_clus)


if __name__ == "__main__":

    main()
