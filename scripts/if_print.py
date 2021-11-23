#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import logging
import argparse


LOG = logging.getLogger(__name__)

__version__ = "1.2.0"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []


def read_tsv(file, sep=None):

    LOG.info("reading message from %r" % file)

    for line in open(file):
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        yield line.split(sep)


def read_prints(prints):

    r = []

    for i in prints.strip().split(","):
        r.append(int(i)-1)

    return r


def get_line(line, prints):

    r = []

    if prints:
        for i in read_prints(prints):
            r.append(line[i])
    else:
        r = line

    return r


def if_print(file, site=1, mins=0, prints=""):

    site = site-1

    for line in read_tsv(file):
        try:
            values = float(line[site])
        except:
            continue

        if values <= mins:
            continue
        line = get_line(line, prints)
        print("\t".join(line))

    return 0


def add_help_args(parser):

    parser.add_argument("input",  metavar='FILE', type=str,
        help="")
    parser.add_argument("--site", metavar='INT', type=int, default=1,
        help="Input the judgment column, default=1")
    parser.add_argument('--mins', metavar='FLOAT', type=float, default=0,
        help='Set the minimum value of judgment, default=0')
    parser.add_argument('-p', '--prints', metavar='STR', type=str, default="",
        help='Output column, default=0')

    return parser


def main():

    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""


version: %s
contact:  %s <%s>\
    """ % (__version__, " ".join(__author__), __email__))

    parser = add_help_args(parser)
    args = parser.parse_args()
    if_print(args.input, args.site, args.mins, args.prints)


if __name__ == "__main__":
    main()
