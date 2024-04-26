import argparse
import os
import sys


def add_arguments(parser):
    parser.add_argument('-i', '--input', dest='items', type=str,
                        help='The input item. The input format MUST be InChIKey, separated by commas.',
                        default="")
    parser.add_argument('-f', '--file', dest="file", type=str, help='Input with files. The file needs to be a .txt file, with each line containing an InChIKey. Example: -f example/inchikeys.txt.')
    # -subtree is replaced by -input xx|subtree

    parser.add_argument('-o', '--prefix', type=str, default=os.path.abspath(
        os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), 'iphylo_files')),
                        help=' Output file prefix, defaults to the current path.')
    parser.add_argument('-fn', '--fname', type=str, default='iPHYLO_Tree', help='Output file name.')
    parser.add_argument('-bl', '--branch_length', dest='branch_length', action='store_true', default=False,
                        help='Choose whether to include branch length for the tree. Default: False.')
    parser.add_argument('-interrupt', action='store_true', default=False, help=' Interrupt the tree at a specified taxonomic level. You must specify the level parameter for this to work, e.g., -interrupt --Super')

    parser.add_argument('--Super', action='store_true', default=False, help='Interrupt the tree at the Superclass level.')
    parser.add_argument('--Class', action='store_true', default=False, help='Interrupt the tree at the Class level.')
    parser.add_argument('--Sub', action='store_true', default=False, help='Interrupt the tree at the Subclass level.')

    # parser.add_argument('-c', '--color', dest='color', type=str, default=None,
    #                     help='color of tree pdf plot')

    parser.add_argument('-x', '--threads', dest='threads', type=int, default=3,
                        help='Threads number, default 3.')

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
        description='Build a taxonomic tree by inchikey of compounds. You are advised to wrap '
                    'the command in double quotation marks. The generated tree is written to the path "newick.txt" in '
                    'the current directory by default, you can change the output directory by command \'-o\' ')
    parser = add_arguments(parser)
    return parser
