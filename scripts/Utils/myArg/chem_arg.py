import argparse
import os
import sys


def add_arguments(parser):
    parser.add_argument('-input', dest='items', type=str,
                        help='Generate a tree by taxid or names, separate the entries '
                             'with commas in English format, mixed input supported. \n'
                             'Example: -input Homo sapiens,Mus musculus,9031,7227,562',
                        default="")
    parser.add_argument('-file', dest="file", type=str, help='Parameter: file path, need to be a .txt file, each line '
                                                             'in txt is a taxid or name. \n'
                                                             'Example: -file species_taxid_for_tree.txt ')
    # -subtree is replaced by -input xx|subtree

    parser.add_argument('-prefix', type=str, default=os.path.abspath(
        os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), 'iphylo_files')),
                        help='Output file prefix, defult current path')
    parser.add_argument('-fname', type=str, default='iPHYLO_Tree', help='Output file name')
    parser.add_argument('-bl', '--branch_length', dest='branch_length', action='store_true', default=False,
                        help='Choose whether you need '
                             'branch length for the tree.'
                             'Default Flase')
    parser.add_argument('-interrupt', action='store_true', default=False, help='Interrupt the tree at specified '
                                                                               'taxonomic level. You need to input '
                                                                               'the level parameter to make it work. \n'
                                                                               'Example: -file '
                                                                               '"species_for_tree_short.txt" '
                                                                               '-interrupt --g')

    parser.add_argument('--Super', action='store_true', default=False, help='Superclass')
    parser.add_argument('--Class', action='store_true', default=False, help='Class')
    parser.add_argument('--Sub', action='store_true', default=False, help='Subclass')

    # parser.add_argument('-c', '--color', dest='color', type=str, default=None,
    #                     help='color of tree pdf plot')

    return parser


def get_parser():
    """
    使用命令行操作ichem
    -i：输入多taxid或name，用英文逗号隔开，直接生成树
        --f 从文件输入，要求一行一个id
    -subtree：输入某一id或name，生成子树
    -interrupt：在某一分类层次过滤
        --p: phylum; --c: class;  --o: order; --f: family; --g: genus; --s: species
    """
    parser = argparse.ArgumentParser(
        description='Build a phylogenetic tree by inchikey of compounds. You are advised to wrap '
                    'the command in double quotation marks. The generated tree is written to the path "newick.txt" in '
                    'the current directory by default, you can change the output directory by command \'-out\' ')
    parser = add_arguments(parser)

    return parser
