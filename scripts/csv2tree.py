import os
import re

import numpy as np
import pandas as pd
from Bio import Phylo

from scripts.Utils import newick
from scripts.Utils.myArg import csv2tree_arg


def check_parent_node(df):
    parent_dict = {}

    for i in range(1, len(df.columns)):
        for j in range(1, len(df)):
            for k in range(j + 1, len(df)):
                if df.iloc[j, i] == df.iloc[k, i] and df.iloc[j, i - 1] != df.iloc[k, i - 1]:
                    child_node = df.iloc[j, i]
                    parent_node_1 = df.iloc[j, i - 1]
                    parent_node_2 = df.iloc[k, i - 1]
                    if child_node not in parent_dict:
                        parent_dict[child_node] = set([parent_node_1, parent_node_2])
                    else:
                        parent_dict[child_node].add(parent_node_1)
                        parent_dict[child_node].add(parent_node_2)

    if parent_dict:
        print("\nThere are same child nodes with different parent nodes:\n")

    # 打印字典内容
    for child_node, parent_nodes in parent_dict.items():
        print(f"Child Node: {child_node}")
        print(f"Parent Nodes: {', '.join(parent_nodes)}")
        print()


def replace_branch_length(fp):
    f = open(fp, 'r')
    all_lines = f.readlines()
    f.close()
    with open(fp, 'w+') as f:
        for line in all_lines:
            a = re.sub(':0.00000', ':1.00000', line)
            f.writelines(a)


def replace_nan(row):
    if pd.isna(row).all():  # 如果一行只有NaN则删除这一行
        return None
    previous_value = None
    last_non_nan_index = row.last_valid_index()  # 找到末位的非NaN值的索引
    row_values = row.tolist()[:last_non_nan_index + 1]  # 删除结尾最后的None
    row_values.reverse()

    for i in range(len(row_values)):
        if row_values[i] is None:
            if previous_value is not None:
                row_values[i] = previous_value + '+'
                previous_value = row_values[i]
        else:
            previous_value = row_values[i]

    # print(row_values)
    row_values.reverse()
    return row_values


def csv_to_list(file_path, fill_gap=False, has_header=True):
    try:
        df = pd.read_csv(file_path, header=None if not has_header else 'infer')
        replace_dict = {
            ' ': '_',
            '\(': '_',
            '\)': '_',
            '\'': '_'
        }
        df = df.replace(replace_dict, regex=True)

        # print(df)
    except Exception as e:
        print(e)
        print("read csv error")
        return

    check_parent_node(df)
    if has_header:
        df = df.iloc[1:]  # 去掉头部行



    df = df.fillna(np.nan)
    df = df.replace({np.nan: None})

    data = df.values.tolist()


    if fill_gap:
        df = pd.DataFrame(data)

        # 替换NaN值并删除行
        df = df.apply(replace_nan, axis=1).dropna()
        data = df.values.tolist()

    data = [[item for item in row if item is not None] for row in data]

    return data


def main(args):
    # parser = csv2tree_arg.get_parser()
    #  = parser.parse_args()
    file_path = args.file
    fill_gap = args.fill_gap
    has_header = args.has_header
    output = args.out
    has_branch_length = args.branch_length

    if not (output.endswith(".txt") or output.endswith(".nwk")):
        print("output file must end with '.txt' or '.nwk'")
        exit()
    else:
        output_nex = output[:-3] + "nex"

    data = csv_to_list(file_path, fill_gap, has_header)
    # print(data)
    try:
        newick.showtree(data, output)
        print("generate tree success")
        if has_branch_length:
            # 转为nex添加分支长度0
            Phylo.convert(output, "newick", output_nex, "nexus")
            # 0改1
            replace_branch_length(output_nex)
            # nex转nwk
            Phylo.convert(output_nex, "nexus", output, "newick")
            # 删除nex
            os.remove(output_nex)

    except:
        print("convert to tree error")
