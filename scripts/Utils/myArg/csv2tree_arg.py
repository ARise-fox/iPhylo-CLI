import argparse
import os
import sys


def add_arguments(parser):
    parser.add_argument('-f', '--file', dest='file', type=str,
                        help='csv table, each column should represent a taxon')
    parser.add_argument('-o', '--prefix', dest='prefix', type=str, default=os.path.abspath(
        os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), 'iphylo_files')),
                        help='output directory')
    parser.add_argument('-fn', '--fname', dest='out_name', type=str, default='iPHYLO_Tree', help='Output file name')
    parser.add_argument('-fg', '--fill-gap', dest='fill_gap', action='store_true', default=False,
                        help='fill the classification gap (some missing classificaltion level) in each row')
    parser.add_argument('-xh', '--header', dest='has_header', action='store_true', default=False,
                        help='remove input file\'s head row')
    parser.add_argument('-bl', '--branch-length', dest='branch_length', action='store_true', default=False,
                        help='add branch length 1')
    return parser


def get_parser():
    parser = argparse.ArgumentParser(
        description='Use this script to construct a hierarchical tree from a .csv spreadsheet.')
    parser = add_arguments(parser)
    return parser
