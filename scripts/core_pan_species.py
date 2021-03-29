#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import logging
import argparse
import matplotlib
matplotlib.use('Agg')

from matplotlib import pyplot as plt
import numpy as np
from itertools import combinations
from scipy.optimize import curve_fit

LOG = logging.getLogger(__name__)

__version__ = "1.1.0"
__author__ = ("Xingguo Zhang",)
__email__ = "113178210@qq.com"
__all__ = []


def read_tsv(file, sep=None):

    LOG.info("reading message from %r" % file)
    if file.endswith(".gz"):
        fp = gzip.open(file)
    else:
        fp = open(file)

    for line in fp:
        if isinstance(line, bytes):
            line = line.decode('utf-8')
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        yield line.split(sep)

    fp.close()


def read_group(file):

    g2s = {}
    s2g = {}
    for line in read_tsv(file, "\t"):
        if line[1] not in g2s:
            g2s[line[1]] = set()
        g2s[line[1]].add(line[0])
        s2g[line[0]] = line[1]

    return g2s, s2g


def read_anno(file, s2g):

    r = {}

    for line in read_tsv(file, "\t"):
        try:       
            sample, seqid = line[0].split("_", 1)
        except:
            print("\t".join(line))
            continue

        group = s2g[sample]
        if group not in r:
            r[group] = {}
        if line[1] not in r[group]:
            r[group][line[1]] = set()
        r[group][line[1]].add(sample)

    return r


def all_true(list1, list2):

    jc = True
    for i in list1:
        if i not in list2:
            jc = False
            break
    return jc


def only_true(list1, list2):

    jc = False
    for i in list1:
        if i in list2:
            jc = True
            break
    return jc


def stat_include_species(data, samples):

    core = 0
    pan = 0

    for line in data.values():
        line = list(line)
        if all_true(samples, line):
            core += 1
            pan += 1
            continue
        if only_true(samples, line):
            pan += 1

    return core, pan


def simulation_sample(data, samples, number):

    cores = []
    pans = []
    n = 0

    for sample in combinations(samples, number):
        n += 1
        sample = list(sample)
        core, pan = stat_include_species(data, sample)
        if core <=50 :
            LOG.info("Abnormal number of otu cores in the sample(%s):%s " % (core, "\t".join(sample)))

        cores.append(core)
        pans.append(pan)
        if n >= 40:
            break

    return sum(cores)*1.0/len(cores), sum(pans)*1.0/len(pans)


def stat_core_pan(spec_dict, samples):

    cores = []
    pans = []

    for i in range(1, len(samples)+1, 1):
        core, pan = simulation_sample(spec_dict, samples, i)
        cores.append(core)
        pans.append(pan)

    return cores, pans


def list_float2str(nlist):

    temp = []
    for i in nlist:
        try:
            n = "{0:.2f}".format(i)
        except:
            n = i
        temp.append(n)
    return temp


def plot_dual_axis(data, prefix):

    lcolor = ["#FFA000", "#5C4314", "#FFC042", "#FFE4AB", "#FFCB52"]
    COLOR = ['#008209', '#26972D', '#210672', '#422C83', '#A68A00', '#BFA730',
    '#FF7673', '#FF4540', '#FF0700', '#FFD500', '#3B14AF', '#00C90D', '#67E46F',
    '#39E444', '#886ED7', '#6C48D7', '#FFE873', '#FFDF40', '#A60400', '#BF3330']
    fig = plt.figure(figsize=(9, 5.6))
    ax = fig.add_subplot(111)

    plt.grid(True, which='minor', linestyle='--', axis='both', lw=1.1, color='#E5C700', alpha=0.3)
    plt.grid(True, which='major', linestyle='--', axis='both', lw=1.5, color='#E2BDD5', alpha=0.3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_linewidth(2)
    ax.spines['left'].set_linewidth(2)
    #ax.spines['left'].set_visible(False)
    ax.tick_params(axis='both', which='both', color='#212121', length=5, width=1.5, direction='out')

    n = 0
    for i in data:
        x = range(1, len(data[i])+1, 1)
        ax.plot(x, data[i], "-", color=lcolor[n], lw=2, label=i)
        n += 1
        if n >= len(lcolor):
            break
 
    font1 = {'family': 'Times New Roman', 'weight': 'normal', 'color': '#212121', 'size': 16}
    ax.set_xlabel('Numbers of samples', font1)
    ax.set_ylabel('Numbers of species', font1)
    plt.legend(loc='upper right', bbox_to_anchor=(1.1,1.0), handlelength=1.0, frameon=False)
    plt.savefig('%s.species.png' % prefix, dpi=700)
    plt.savefig('%s.species.pdf' % prefix)


def core_pan_species(anno, group, prefix):

    g2s, s2g = read_group(group)
    data = read_anno(anno, s2g)
    dcore = {}
    dpan = {}

    for i in data:
        spec_dict = data[i]
        samples = list(g2s[i])
        cores, pans = stat_core_pan(spec_dict, samples)
        dcore[i] = cores
        dpan[i] = pans
        print("Group %s\tPan species\t%s" % (i, "\t".join(list_float2str(pans))))
        print("Group %s\tCore species\t%s" % (i, "\t".join(list_float2str(cores))))
    plot_dual_axis(dcore, "%s.core" % prefix)
    plot_dual_axis(dpan, "%s.pan" % prefix)

    return r


def add_hlep_args(parser):

    parser.add_argument('anno', metavar='FILE', type=str,
        help='Input species annotation file.')
    parser.add_argument('-g', '--group', metavar='FILE', type=str, required=True,
        help='Input group file.')
    parser.add_argument('-p', '--prefix', metavar='FILE', type=str, default='out',
        help='Output file prefix, default=out.')

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
    core_pan_species.py: Calculate the core and pan otu of metagenomics.
attention:
    core_pan_species.py annotation.xls -g group.list >stat_pan_gene.tsv
    core_pan_species.py annotation.xls -g group.list -p name >stat_pan_gene.tsv
version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    core_pan_species(args.anno, args.group, args.prefix)


if __name__ == "__main__":

    main()
