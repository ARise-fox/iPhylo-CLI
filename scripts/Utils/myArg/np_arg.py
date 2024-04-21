import argparse
import os
import sys


def add_arguments(parser):
    parser.add_argument('-i', '--input', dest='items', type=str,
                        help='Generate a tree by InChI, InChIKeys or isomeric SMILES, separate the entries '
                             'with commas in English format, mixed input supported. \n'
                             'Example: -input HPDXGYHKPNONBN-NSCNQCMKSA-N, KYMWLJVXAUXWSY-USSVWMPFSA-N, '
                             'LCSIRHSEQLDODB-XSWRHYEBSA-N, LERTZPRWBVXJDZ-WSXCIZDTSA-N, LIYNAWWFJZSTGA-MYWJTRDOSA-N',
                        default="")
    parser.add_argument('-f', '--file', dest="file", type=str,
                        help='Parameter: file path, need to be a .txt file, each line '
                             'should be InChI, InChIKeys or isomeric SMILES. \n'
                             'Example: -file nature_products.txt ')
    # -subtree is replaced by -input xx|subtree

    parser.add_argument('-o', '--prefix', type=str, default=os.path.abspath(
        os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), 'iphylo_files')),
                        help='Output file prefix, default current path')
    parser.add_argument('-fn', '--fname', type=str, default='iPHYLO_NP_Tree', help='Output file name')
    parser.add_argument('-bl', '--branch_length', dest='branch_length', action='store_true', default=False,
                        help='Choose whether you need '
                             'branch length for the tree.'
                             'Default Flase')
    parser.add_argument('-interrupt', action='store_true', default=False, help='Interrupt the tree at specified '
                                                                               'taxonomic level. You need to input '
                                                                               'the level parameter to make it work. \n'
                                                                               '-interrupt --superclass')

    parser.add_argument('--pathway', action='store_true', default=False, help=' interrupt at Pathway')
    parser.add_argument('--superclass', action='store_true', default=False, help='interrupt at Superclass')

    return parser


def get_parser():
    """
    使用命令行操作ichem
    -i：输入多InChIKey，用英文逗号隔开，直接生成树
        -f 从文件输入，要求一行一个id
    -subtree：输入某一化学分类，生成子树
    -interrupt：在某一分类层次过滤
       --superclass: superclass, --pathway: pathway,
    """
    parser = argparse.ArgumentParser(
        description='Build a phylogenetic tree by InChI, InChIKeys or isomeric SMILES of compounds.'
                    'The classification system NPClassifier is used. '
                    ' You are advised to wrap '
                    'the command in double quotation marks. The generated tree is written to the path "newick.txt" in '
                    'the current directory by default, you can change the output directory by command \'-o\' ')
    parser = add_arguments(parser)

    return parser
