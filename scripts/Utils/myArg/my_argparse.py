import argparse
import os
import sys


def add_arguments(parser):
    parser.add_argument('-i', '--input', dest='items', type=str,
                        help='Generate a tree by taxid or names, separate the entries '
                             'with commas in English format, mixed input supported. \n'
                             'Example: -input \"Homo sapiens,Mus musculus,9031,7227,562\"',
                        default="")
    parser.add_argument('-f', '--file', dest="file", type=str, help='Parameter: file path, need to be a .txt file, '
                                                                    'each line'
                                                                    'in txt is a taxid or name. \n'
                                                                    'Example: -file species_taxid_for_tree.txt ')
    # -subtree is replaced by -input xx|subtree
    parser.add_argument('--subtree', type=str, help='Draw subtree of a certain taxon. Parameter: taxid or name. '
                                                    'Example: -subtree Mammalia, also you can use \' -input xx|subtree '
                                                    '\' instead.')
    # parser.add_argument('-o', '--out', dest='out_path', type=str, help='Output the tree file to the specified file path, '
    #                                                             'support newick and phyloxml format. \n '
    #                                                             'Example: -i 9606,10090,9031,7227,562 -out '
    #                                                             '"C:\mytest.txt"',
    #                     default='newick.txt')
    parser.add_argument('-o', '--prefix', type=str, default=os.path.abspath(
        os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), 'iphylo_files')),
                        help='Output file directory, default \'iphylo_files/\' in the current project')

    parser.add_argument('-fn', '--fname', type=str, default='iPHYLO_Tree', help='Output file name')
    parser.add_argument('-bl', '--branch_length', dest='branch_length', action='store_true', default=False,
                        help='Boolean, choose whether you need '
                             'branch length for the tree.'
                             'Default False')
    parser.add_argument('-interrupt', action='store_true', default=False, help='Interrupt the tree at specified '
                                                                               'taxonomic level. You need to input '
                                                                               'the level parameter to make it work. \n'
                                                                               'Example: -file '
                                                                               '"species_for_tree_short.txt" '
                                                                               '-interrupt --g')
    parser.add_argument('--p', action='store_true', default=False, help='Phylum level')
    parser.add_argument('--c', action='store_true', default=False, help='Class level')
    parser.add_argument('--o', action='store_true', default=False, help='Order level')
    parser.add_argument('--f', action='store_true', default=False, help='Family level')
    parser.add_argument('--g', action='store_true', default=False, help='Genus level')
    parser.add_argument('--s', action='store_true', default=False, help='Species level')
    #
    # parser.add_argument('-c', '--color', dest='color', type=str, default=None,
    #                     help='color of tree pdf plot')
    return parser


def get_parser():
    """
    使用命令行操作
    -i：输入多taxid或name，用英文逗号隔开，直接生成树
        --f 从文件输入，要求一行一个id
    -subtree：输入某一id或name，生成子树
    -interrupt：在某一分类层次过滤
        --p: phylum; --c: class;  --o: order; --f: family; --g: genus; --s: species
    """
    parser = argparse.ArgumentParser(
        description='Build a phylogenetic tree by taxid or name in NCBI taxonomy database. You are advised to wrap '
                    'the command in double quotation marks. The generated tree is written to the path "newick.txt" in '
                    'the current directory by default, you can change the output directory by command \'-o\' ')
    parser = add_arguments(parser)
    return parser
