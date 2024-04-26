import argparse
import sys
from scripts import phylotree, NPtree, chemtree, csv2tree, chemonline, NPonline
from scripts.Utils.myArg import my_argparse
from scripts.Utils.myArg import chem_arg
from scripts.Utils.myArg import chem_arg_online
from scripts.Utils.myArg import csv2tree_arg
from scripts.Utils.myArg import np_arg_online
from scripts.Utils.myArg import np_arg


def main():
    parser = argparse.ArgumentParser(description='iPhylo Command Line Interface')
    subparsers = parser.add_subparsers(dest='command')

    # phylotree command
    phylotree_parser = subparsers.add_parser('phylotree', help='Run phylotree.py to generate a biological taxonomic '
                                                               'tree using the NCBI taxonomic classification.',
                                             description='Generate a tree by taxid or names, separate the entries '
                                                         'with commas in English format, mixed input supported. \n'
                                                         'Example: -input \"Homo sapiens,Mus musculus,9031,7227,562\"')
    phylotree_parser = my_argparse.add_arguments(phylotree_parser)
    # correctly get parser, tested
    phylotree_parser.set_defaults(func=phylotree)
    # print(phylotree_parser.parse_args())

    # chemtree command
    chemtree_parser = subparsers.add_parser('chemtree', help=' Run chemtree.py to generate a chemical taxonomic tree using the ChemOnt classification system based on the local database.',
                                            description='This command performs the taxonomic analysis of '
                                                        'chemical compounds by querying ClassyFire locally.')
    chemtree_parser = chem_arg.add_arguments(chemtree_parser)
    chemtree_parser.set_defaults(func=chemtree)

    # chemonline command
    chemonline_parser = subparsers.add_parser('chemonline', help='Run chemonline.py to generate a chemical taxonomic tree using the ChemOnt classification system based on the online API. This module extends the taxonomy to include compounds beyond those in the local database.',
                                              description='This command performs the taxonomic analysis of '
                                                          'chemical compounds by querying ClassyFire online. \n'
                                                          'The input format MUST be InChIKey of compounds. \n'
                                                          ' You are advised to wrap the command in double quotation marks.')
    chemonline_parser = chem_arg_online.add_arguments(chemonline_parser)
    chemonline_parser.set_defaults(func=chemonline)

    # csv2tree command
    csv2tree_parser = subparsers.add_parser('csv2tree', help='Run csv2tree.py to construct a hierarchical tree from a .csv spreadsheet.',
                                            description='This command constructs a hierarchical'
                                                        ' tree from a .csv spreadsheet.')
    csv2tree_parser = csv2tree_arg.add_arguments(csv2tree_parser)
    csv2tree_parser.set_defaults(func=csv2tree)


    # NPtree command
    NPtree_parser = subparsers.add_parser('NPtree', help='Run NPtree.py to generate a chemical taxonomic tree using the local NPClassifier database.',
                                          description='Build a chemical taxonomic tree by InChI, InChIKeys or '
                                                      'isomeric SMILES of compounds. '
                                                      'The classification system NPClassifier is used, but performed '
                                                      'locally. '
                                                      'You are advised to wrap the command in double quotation marks.')
    NPtree_parser = np_arg.add_arguments(NPtree_parser)
    NPtree_parser.set_defaults(func=NPtree)

    # NPonline command
    NPonline_parser = subparsers.add_parser('NPonline', help='Run NPonline.py to generate a chemical taxonomic tree using the NPClassifier classification system based on the online API. This module extends the taxonomy to include compounds beyond those in the local database.',
                                            description='Build a taxonomic tree by isomeric SMILES of compounds.'
                                                        'This command performs the taxonomic analysis by querying '
                                                        'NPClassifier online. '
                                                        'You are advised to wrap the command in double quotation marks.'
                                            )
    NPonline_parser = np_arg_online.add_arguments(NPonline_parser)
    NPonline_parser.set_defaults(func=NPonline)


    args = parser.parse_args()

    if args.command == 'phylotree':
        phylotree.main(args)
    elif args.command == 'chemtree':
        chemtree.main(args)
    elif args.command == 'chemonline':
        chemonline.main(args)
    elif args.command == 'csv2tree':
        csv2tree.main(args)
    elif args.command == 'NPonline':
        NPonline.main(args)
    elif args.command == 'NPtree':
        NPtree.main(args)
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == '__main__':
    main()
