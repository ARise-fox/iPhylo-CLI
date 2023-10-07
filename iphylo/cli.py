import argparse
import sys
from scripts import phylotree
from scripts import chemtree
from scripts import csv2tree
from scripts import chemonline
from scripts.Utils.myArg import my_argparse
from scripts.Utils.myArg import chem_arg
from scripts.Utils.myArg import chem_arg_online
from scripts.Utils.myArg import csv2tree_arg


def main():
    parser = argparse.ArgumentParser(description='iPhylo Command Line Interface')
    subparsers = parser.add_subparsers(dest='command')

    # phylotree command
    phylotree_parser = subparsers.add_parser('phylotree', help='Run phylotree.py')
    phylotree_parser = my_argparse.add_arguments(phylotree_parser)
    # correctly get parser, tested
    phylotree_parser.set_defaults(func=phylotree)
    # print(phylotree_parser.parse_args())

    # chemtree command
    chemtree_parser = subparsers.add_parser('chemtree', help='Run chemtree.py')
    chemtree_parser = chem_arg.add_arguments(chemtree_parser)
    chemtree_parser.set_defaults(func=chemtree)

    # chemonline command
    chemonline_parser = subparsers.add_parser('chemonline', help='Run chemonline.py')
    chemonline_parser = chem_arg_online.add_arguments(chemonline_parser)
    chemonline_parser.set_defaults(func=chemonline)

    # csv2tree command
    csv2tree_parser = subparsers.add_parser('csv2tree', help='Run csv2tree.py')
    csv2tree_parser = csv2tree_arg.add_arguments(csv2tree_parser)
    csv2tree_parser.set_defaults(func=csv2tree)

    args = parser.parse_args()

    if args.command == 'phylotree':
        phylotree.main(args)
    elif args.command == 'chemtree':
        chemtree.main(args)
    elif args.command == 'chemonline':
        chemonline.main(args)
    elif args.command == 'csv2tree':
        csv2tree.main(args)
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == '__main__':
    main()
