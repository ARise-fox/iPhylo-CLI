import argparse
import os
import sys


def add_arguments(parser):

    parser.add_argument('-f', '--file', dest="file", type=str, help='Parameter: file path, need to be a .txt file, each line '
                                                             'in txt is SMILES for a chemical compound \n'
                                                             'Example: -file smiles.txt ')
    # -subtree is replaced by -input xx|subtree

    parser.add_argument('-o', '--prefix', type=str, default=os.path.abspath(
        os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), 'iphylo_files')),
                        help='Output file prefix, default current path')
    parser.add_argument('-fn', '--fname', type=str, default='NP_Tree', help='Output file name')
    parser.add_argument('-bl', '--branch_length', dest='branch_length', action='store_true', default=False,
                        help='Choose whether you need '
                             'branch length for the tree.'
                             'Default Flase')
    parser.add_argument('-interrupt', action='store_true', default=False, help='Interrupt the tree at specified '
                                                                               'taxonomic level. You need to input '
                                                                               'the level parameter to make it work. \n'
                                                                               'Example: -file '
                                                                               '"smiles.txt" '
                                                                               '-interrupt --pathway')

    parser.add_argument('--pathway', action='store_true', default=False, help='interrupt at pathway level')
    parser.add_argument('--superclass', action='store_true', default=False, help='interrupt at superclass level')
    parser.add_argument('-x', '--threads', dest='threads', type=int, default=3,
                        help='threads number')

    return parser

def get_parser():
    """
    使用命令行操作ichem
    -i：输入多inchikeys，直接生成树
        --f 从文件输入，要求一行一个id
    -interrupt：在某一分类层次过滤
        --Super: Superclass, --Class: Class, --Sub: subclass
    """
    parser = argparse.ArgumentParser(
        description='Build a taxonomic tree by isomeric SMILES of compounds. '
                    'The classification system NPClassifier is used. '
                    'You are advised to wrap the command in double quotation marks. The generated tree is written to '
                    'the path "newick.txt" in'
                    'the current directory by default, you can change the output directory by command \'-o\' ')
    parser = add_arguments(parser)
    return parser
