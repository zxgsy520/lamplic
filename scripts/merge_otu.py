#!/usr/bin/env python
#coding:utf-8

import os
import re
import sys
import logging
import argparse

from collections import OrderedDict


LOG = logging.getLogger(__name__)

__version__ = "1.0.1"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []


def read_tsv(file, sep="\t"):

    for line in open(file):
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        yield line.split(sep)


def read_otu_tax(file):

    r = {}
    tax2otu = {}
    tax_dict = {}

    for line in read_tsv(file):
        r[line[0]] = line[1]
        if "s__" not in line[1]:
            continue
        tax_dict[line[0]] = line[1::]
        if line[1] not in tax2otu:
            tax2otu[line[1]] = []
        tax2otu[line[1]].append(line[0])

    return r, tax2otu, tax_dict


def get_sample(otus):

    samples = []

    for i in otus:
        sample, seq_id = i.split('_', 1)
        if sample in samples:
            continue
        samples.append(sample)

    return samples


def read_otu_map(file):

    r = OrderedDict()
    otuids = []
    samples = []

    for line in read_tsv(file):
        r[line[0]] = line[1::]
        otuids.append(line[0])
        samples += get_sample(line[1::])
    samples = sorted(list(set(samples)))

    return r, otuids, samples


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
        otus.append(otu)

    return otus


def format_list(otus):

    r = []

    for i in otus:
        r.append('%.2f' % i)

    return r


def merge_otu(tax, otu, no_merge):

    otu2tax, tax2otu, tax_dict = read_otu_tax(tax)
    otu, otuids, samples = read_otu_map(otu)

    for i in otuids:
        if i not in otu2tax:
            continue
        taxid = otu2tax[i]
        if taxid not in tax2otu:
            continue
        for j in tax2otu[taxid]:
            if (i == j) or (i not in otu):
                continue
            if no_merge:
                continue
            otu[i] += otu[j]
            jvalue = otu.pop(j)
            LOG.info("%s merges into %s" % (j, i))

    fm = open('otu_map_final.tsv', 'w')
    fo = open('otu_tab_final.tsv', 'w')
    ft = open('otu_map_final.tax', 'w')
    fo.write('# Constructed from biom file\n')
    fo.write('#OTU ID\t%s\ttaxonomy\n' % '\t'.join(samples))

    for i in otu:
        print("%s\t%s" % (i, "\t".join(sorted(otu[i]))))
        if i not in otu2tax:
            continue
        tax = otu2tax[i].rstrip(';')
        if "s__" not in tax:
            continue
        otus = stat_otu(otu[i])
        otus = sort_otus(samples, otus)
        if sum(otus) <= len(samples)/2.5 or sum(otus)<5:
            continue
        tax = tax.replace(';', '; ')
        otus = format_list(otus)
        fo.write("%s\t%s\t%s\n" % (i, '\t'.join(otus), tax))
        fm.write("%s\t%s\n" % (i, "\t".join(sorted(otu[i]))))
        ft.write("%s\t%s\n" % (i, "\t".join(tax_dict[i])))
    fo.close()
    fm.close()
    ft.close()

    return 0


def add_hlep_args(parser):

    parser.add_argument('input', metavar='FILE', type=str,
        help='Input the out result of the cluster.')
    parser.add_argument('-tax', '--taxonomy', metavar='FILE', type=str, required=True,
        help='Input the result of species classification annotation')
    parser.add_argument("--no_merge", action="store_true",
        help="Does not merge otu by species.")

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
    merge_otu.py Combine clustering results based on species annotation results.

attention:
    merge_otu.py otu_map.tsv -tax otu_tax.txt >otu_map_new.tsv
    merge_otu.py otu_map.tsv -tax otu_tax.txt --no_merge >otu_map_new.tsv
version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    merge_otu(args.taxonomy, args.input, args.no_merge)


if __name__ == "__main__":

    main()
