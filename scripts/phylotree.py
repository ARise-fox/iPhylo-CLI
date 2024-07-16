# from ete3 import Tree, TreeStyle
import os
import sqlite3
import re
import time
from Bio import Phylo
from pathlib import Path
import csv
from alive_progress import alive_bar
from scripts.Utils.download_db import download_db_file
# import eteUtils.render_test
# from eteUtils import tree_style
from scripts.Utils.myArg.my_argparse import get_parser

"""
DB Link: 111
"""
# 不在数据库里的id
import sys

from scripts.Utils import newick
from scripts.Utils.eteUtils import myete

null_id = []
# 不在数据库里的name，可能以后会处理，用list存放
null_name = []
# 重复名列表
multi_names = []
# 是否有层次打断flag
interrupt_flag = 0

"""
SQL语句
sqlite3占位符语法
"""
# 查询输入id的分类信息
id_lineage_tabel = "id2lineage"
sql = f"""select _superkingdom, _kingdom, _phylum, _class, _order, _family, _genus, _species from {id_lineage_tabel} where id = :id """
# sql_find_id_and_sciname_by_name = """select taxid, sciname from name_id_sci where name = :name """
# print(sql)
sql_find_sciname = """select name from alter_name_nopk_flag where taxid = :taxid and type_flag = 1 """
sql_find_id_and_flag_by_name = """select taxid, type_flag from alter_name_nopk_flag where name = :name """
# 查询名字是否有重复
sql_check_multi_name = """select count(*) from multi_taxid WHERE name = :name"""
# wrong database!!!!
# sql_find_id_and_flag_by_name = """select taxid, type_flag from name_id_flag where name = :name """
# sql_find_sciname = """select name from name_id_flag where taxid = :taxid and type_flag = 1 """
# # 用不上的两条sql
# sql_find_sciname_by_id = """select sciname from id2name where id = :id """
sql_find_id_by_name = """select taxid from alter_name_nopk_flag where name = :name """

sql_name_by_id = f"""select sciname from {id_lineage_tabel} where id = :id """
sql_rank_by_name = f"""select * from {id_lineage_tabel} where sciname = :name """
sql_rank_by_id = f"""select * from {id_lineage_tabel} where id = :name """

def write_anno(data, filename):
    with open(filename, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)

            # 写入列名
            writer.writerow(['id', 'taxid', 'domain', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'])

            # 写入数据
            for row in data:
                last_item = row[-1]  # 获取每行的最后一个元素
                row_data = [last_item]  # 将最后一个元素作为 id 写入第一列
                taxid = row[0]  # 获取每行的最后一个元素
                row_data.append(taxid) # 将最后一个元素作为 id 写入第一列

                for item in row[1:]:  # 再遍历一遍所有元素, 除了第一列taxid
                    if item.startswith('d__'):
                        row_data.append(item[3:])
                    elif item.startswith('k__'):
                        row_data.append(item[3:])
                    elif item.startswith('p__'):
                        row_data.append(item[3:])
                    elif item.startswith('c__'):
                        row_data.append(item[3:])
                    elif item.startswith('o__'):
                        row_data.append(item[3:])
                    elif item.startswith('f__'):
                        row_data.append(item[3:])
                    elif item.startswith('g__'):
                        row_data.append(item[3:])
                    elif item.startswith('s__'):
                        row_data.append(item[3:])
                writer.writerow(row_data)

def connect_db():
    """
    初始化数据库连接，打印cursor确保连接成功
    """
    # 数据库路径
    db_file = 'data/iphylo.db'
    if not os.path.exists(db_file):
        user_input = input(f"Necessary database resources (954M) not found. Do you want to download it? (y/n): ")
        if user_input.lower() == 'y':
            url = 'https://iphylo.net/resource/iphylo_db'
            download_db_file(url, db_file)
        else:
            print("Exiting iPhylo command-line.")
            exit(0)  # Exit the program gracefully

    if not os.path.exists(db_file):
        print("Database connect failed!")
        exit(1)
    try:
        # 连接数据库
        conn = sqlite3.connect(db_file)
        # 创建游标
        cursor = conn.cursor()
    except Exception as e:
        print("Database connect failed!", str(e))
        exit(1)
    else:
        print("database connect success!")
    return conn, cursor


# 解决branch length长度为0的问题,设置成1
def replace_branch_length(fp):
    f = open(fp, 'r')
    all_lines = f.readlines()
    f.close()
    with open(fp, 'w+') as f:
        for line in all_lines:
            a = re.sub(':0.00000', ':1.00000', line)
            f.writelines(a)


# 去掉branch length
def cut_branch_length(fp):
    f = open(fp, 'r')
    all_lines = f.readlines()
    f.close()
    with open(fp, 'w+') as f:
        for line in all_lines:
            a = re.sub(':0.00000', '', line)
            f.writelines(a)


def commit_close_db(conn, cursor):
    # 关闭游标
    cursor.close()
    # 提交事务
    conn.commit()
    # 关闭连接
    conn.close()


# 从name_id_flag表查2次， 生成替换名
def get_id_and_sciname_list_by_id_or_name(in_list, cursor):
    id_list = []
    name_sciname_map = {}
    bar_length = len(in_list)
    print(f"Searching for {bar_length} items… ")
    with alive_bar(bar_length) as bar:
        for item in in_list:
            # print(item)
            if item.isdigit():
                id_list.append(item)
            else:
                multi = cursor.execute(sql_check_multi_name, {"name": item}).fetchall()
                if int(multi[0][0]) > 0:  # 该名字重复
                    # print(multi)
                    multi_names.append(item)
                # _start = time.perf_counter()
                cursor.execute(sql_find_id_and_flag_by_name, {"name": item})
                out_tuple = cursor.fetchall()
                if out_tuple:
                    _id = out_tuple[0][0]
                    if out_tuple[0][1] == 0:
                        # flag 表示非标准名
                        sci_tuple = cursor.execute(sql_find_sciname, {"taxid": _id}).fetchall()
                        try:
                            _sciname = sci_tuple[0][0]
                        except:
                            print(sci_tuple)
                            print(_id)
                        name_sciname_map[item] = _sciname
                    id_list.append(_id)
                    # _end = time.perf_counter()
                    # print('Running time: %s Seconds' % (_end - _start))
                else:
                    null_name.append(item)
                    continue
            bar()

    return id_list, name_sciname_map


# 不生成替换名的表
def get_id_list_by_id_or_name(in_list, cursor):
    id_list = []
    # name_sciname_map = {}
    for item in in_list:
        print(item)
        if item.isdigit():
            id_list.append(item)
        else:
            multi = cursor.execute(sql_check_multi_name, {"name": item}).fetchall()
            if int(multi[0][0]) > 0:  # 该名字重复
                # print(multi)
                multi_names.append(item)
                # _start = time.perf_counter()
                # _start = time.perf_counter()
            cursor.execute(sql_find_id_by_name, {"name": item})
            out_tuple = cursor.fetchall()
            if out_tuple:
                _id = out_tuple[0][0]
                id_list.append(_id)
                # _end = time.perf_counter()
                # print('Running time: %s Seconds' % (_end - _start))
            else:
                null_name.append(item)
                continue
    return id_list


# 只从name_id_sci表查一次
# def get_id_and_sciname_list_by_id_or_name(in_list, cursor):
#     id_list = []
#     name_sciname_map = {}
#     for item in in_list:
#         print(item)
#         if item.isdigit():
#             id_list.append(item)
#         else:
#             # _start = time.perf_counter()
#             cursor.execute(sql_find_id_by_name, {"name": item})
#             out_tuple = cursor.fetchall()
#             if out_tuple:
#                 _id = out_tuple[0][0]
#                 sciname_tuple = cursor.execute(sql_find_sciname_by_id, {"id": _id}).fetchall()
#                 _sciname = sciname_tuple[0][0]
#                 if item != _sciname:
#                     name_sciname_map[item] = _sciname
#                 id_list.append(_id)
#                 # _end = time.perf_counter()
#                 # print('Running time: %s Seconds' % (_end - _start))
#             else:
#                 null_name.append(item)
#                 continue
#
#     return id_list, name_sciname_map


# def get_id_and_sciname_list_by_id_or_name(in_list, cursor):
#     id_list = []
#     name_sciname_map = {}
#     for item in in_list:
#         if item.isdigit():
#             id_list.append(item)
#         else:
#             print(item)
#             # _start = time.perf_counter()
#             cursor.execute(sql_find_id_and_sciname_by_name, {"name": item})
#             out_tuple = cursor.fetchall()
#             if out_tuple:
#                 _id = out_tuple[0][0]
#                 _sciname = out_tuple[0][1]
#                 if item != _sciname:
#                     name_sciname_map[item] = _sciname
#                 id_list.append(_id)
#                 # _end = time.perf_counter()
#                 # print('Running time: %s Seconds' % (_end - _start))
#             else:
#                 null_name.append(item)
#                 continue
#
#     return id_list, name_sciname_map
"""
需求：如果输入item是标准lineage则使用标准7级，如果不是，在标准lineage末尾加上item‘s name
"""


def id_list_to_newick(id_list, cursor, *interrupt_level, anno_name, path, csv_write):
    """
    根据id_list 转换成 newick格式
    :param id_list:
    :return:
    """
    # lin_list: 一个二维数组，访问第1个物种的第2个rank应该用lin_list[0][1]
    lin_list = []
    anno_2d_list = []
    # 记录scientific name 和 leaf node 的 mapping
    name_res_mapping = {}  # 用于记录name和res[-1]的映射关系

    # 如果需要在某level打断，找到该level对应的index
    if interrupt_flag == 1:
        level_list = ["p", "c", "o", "f", "g", "s"]

        for idx, level in enumerate(level_list):
            # print(interrupt_level)
            if interrupt_level[0] == level:
                level_idx = idx + 3
    id_list_rd = set(id_list)
    for item in id_list_rd:  # 去重

        int_item = int(item)
        if int_item == 1:
            continue
        cursor.execute(sql, {"id": int_item})
        out_tuple = cursor.fetchall()
        # if int_item == 730:
        #     print(out_tuple[0])
        if not out_tuple:
            global null_id
            null_id.append(item)
            continue
        full_list = list(out_tuple[0])
        if full_list[0] is None:
            print(f"Id {int_item} has no classification.")
            continue
        # name: 输入item sciname
        name = cursor.execute(sql_find_sciname, {"taxid": int_item}).fetchone()[0]

        # print(full_list)
        # 写CSV文件
        csv_list = list(full_list)
        csv_list.insert(0, name)
        csv_list.insert(0, item)
        csv_write.writerow(csv_list)
        # cursor.execute(sql_name_by_id, {"id": int(item)})
        # name = cursor.fetchone()[0]
        # print(name)
        # 获取列表最后一个非None元素， index: 最后一个非None索引/空字符串

        index = next(len(full_list) - i for i, j in enumerate(reversed(full_list), 1) if j is not None)
        # 如果sciname（从数据库获取）该索引的值，在index后插入该sciname
        if name != str(full_list[index])[3:]: # 3：过滤前3个字符如K__
            full_list.insert(index + 1, name)
        # 如果有打断
        if interrupt_flag == 1:
            full_list = full_list[:level_idx]
        # print(f'full list:{full_list}')
        # # 过滤none
        # if 'o__Pasteurellales' in full_list:
        #     print(full_list)
        #     print(int_item)
        res = list(filter(None, full_list))
        lin_list.append(res)
        # 用于注释 anno
        anno_list = res[:]
        anno_list.insert(0, item)
        anno_2d_list.append(anno_list)

        # 记录name和res[-1]的映射关系
        if res:
            name_res_mapping[name] = res[-1]

    # 二维列表去重，需要先转化为tuple
    if interrupt_flag == 1:
        lin_list = list(set([tuple(t) for t in lin_list]))
        lin_list = [list(v) for v in lin_list]
    # print(f'items:{lin_list}')
    write_anno(anno_2d_list, anno_name)
    # print(lin_list)
    # 生成树文件
    newick.showtree(data=lin_list, path=path)
    return name_res_mapping


def get_rank(item, cursor, isdigit):
    rank_list = ["_superkingdom", "_kingdom", "_phylum", "_class", "_order", "_family", "_genus", "_species"]
    rank_query_list = ["d__", "k__", "p__", "c__", "o__", "f__", "g__", "s__"]
    if isdigit:
        cursor.execute(sql_rank_by_id, {"name": item})
    else:
        cursor.execute(sql_rank_by_name, {"name": item})
    query_list = cursor.fetchone()
    rank_name = ""  # _phylum
    rank_query = ""  # p__Ascomycota

    if query_list is not None:
        sciname = query_list[1]
        for idx, item in enumerate(query_list[2:]):
            if rank_query_list[idx] + sciname == item:
                rank_name = rank_list[idx]
                rank_query = item
                break

    return rank_name, rank_query


def get_id_list_by_rank(rank_name, rank_query, cursor):
    """

    :param rank_name: _phylum
    :param rank_query: p__Ascomycota
    :param cursor:
    :return:
    """
    sql_find_id_by_rank = "SELECT id from id2lineage WHERE " + str(rank_name) + " = '" + str(rank_query) + "'"
    # print(sql_find_id_by_rank)
    cursor.execute(sql_find_id_by_rank)
    res = cursor.fetchall()
    id_list = []
    for _id in list(res):
        id_list.append(_id[0])

    # print(id_list)
    return id_list


def draw_ascii_tree(prefix_fname):
    """
    :param fp: tree file path from param -out
    """
    # 画ascii树
    # 输入xml格式
    # 判断路径是否存在
    fp = prefix_fname + '.txt'
    f_ascii = prefix_fname + '_ascii_tree.txt'
    if os.path.exists(fp):
        tree = Phylo.read(fp, "newick")
        tree.rooted = True
        with open(f_ascii, 'w') as fd:
            Phylo.draw_ascii(tree, file=fd)
            # print(f'Tree and ASCII tree is saved to {my_prefix_fname}')
            print(f'Tree and ASCII tree is saved to: {fp}')
    else:
        print('ascii_tree load failed')
        # tree = Phylo.read('newick.txt', "newick")
        # tree.rooted = True
        # with open('acsii_tree.txt', 'w') as fd:
        #     Phylo.draw_ascii(tree, file=fd)


def replace_name_record(prefix_fname, map):
    """
        :param prefix_fname: tree file path + name
        :param map: replace name map
        """
    f_map = prefix_fname + '_replaced_name.txt'  # 输出文件路径+文件名
    with open(f_map, 'w') as f:
        for name, sciname in map.items():
            f.write(name + '\t')
            f.write('->')
            f.write('\t' + sciname)
            f.write('\n')


# def convert_file(out_path):
#     """
#         处理输出文件的格式
#         :param out_path: 输出文件的路径
#     """
#     # 先在输出路径生成txt文件
#     out_dir = out_path[:-4]
#     txt_path = out_dir + ".txt"
#     # 处理输出文件结尾指定格式
#     # 输出nwk格式
#     if out_path.endswith(".nwk"):
#         f_replace = out_dir + '.nwk'
#         Phylo.convert(txt_path, "newick", f_replace, "newick")
#         if args.branch_length:
#             replace_branch_length(f_replace)
#         else:
#             cut_branch_length(f_replace)
#         # os.remove(txt_path)
#     # 输出txt格式
#     elif out_path.endswith(".txt"):
#         # 先把txt转成newick格式
#         f_nwk = out_dir + '.nwk'
#         Phylo.convert(txt_path, "newick", f_nwk, "newick")
#         if args.branch_length:
#             replace_branch_length(f_nwk)
#         else:
#             cut_branch_length(f_nwk)
#         os.remove(txt_path)  # 删除.txt
#         os.rename(f_nwk, out_path)
#     # 输出nexus格式
#     elif out_path.endswith(".nex"):
#         f_replace = out_dir + '.nex'
#         Phylo.convert(txt_path, "newick", f_replace, "nexus")
#         if args.branch_length:
#             replace_branch_length(f_replace)
#         else:
#             cut_branch_length(f_replace)
#         # os.remove(txt_path)
#     # 输出phyloxml格式
#     elif out_path.endswith(".xml"):
#         f_nex = out_dir + '.nex'
#         f_xml = out_dir + '.xml'
#         # 如果需要分支长度，先转化成nex格式添加分支长度
#         if args.branch_length:
#             Phylo.convert(txt_path, "newick", f_nex, "nexus")
#             replace_branch_length(f_nex)
#             Phylo.convert(f_nex, "nexus", f_xml, "phyloxml")
#             os.remove(f_nex)
#         else:
#             Phylo.convert(txt_path, "newick", f_xml, "phyloxml")
#         # os.remove(txt_path)  # 删除txt

def convert_file(prefix_fname, args):
    """
        处理输出文件的格式
        :param args: argument 对象
        :param prefix_fname: 输出文件的路径+文件命名
    """
    # 先在输出路径生成txt文件
    txt_path = prefix_fname + ".txt"
    # 处理输出文件结尾指定格式
    # 输出nwk格式
    f_nwk = prefix_fname + '.nwk'
    # Phylo.convert('F:/PythonProgram/my_phylo/files/test14.txt', "newick",
    # 'F:/PythonProgram/my_phylo/files/test14.nwk', "newick")

    Phylo.convert(txt_path, "newick", f_nwk, "newick")

    if args.branch_length:
        replace_branch_length(f_nwk)
    else:
        cut_branch_length(f_nwk)
        # os.remove(txt_path)
    # 输出txt格式  复制nwk重命名txt
    f_txt = prefix_fname + '.txt'
    f1 = open(f_nwk)
    f2 = open(f_txt, 'w')
    f2.write(f1.read())
    f1.close()
    f2.close()
    # 输出nexus格式
    f_nex = prefix_fname + '.nex'
    Phylo.convert(f_nwk, "newick", f_nex, "nexus")
    if args.branch_length:
        replace_branch_length(f_nex)
    else:
        cut_branch_length(f_nex)
        # os.remove(txt_path)
    # 输出phyloxml格式
    f_xml = prefix_fname + '.xml'
    # 如果需要分支长度，先转化成nex格式添加分支长度
    Phylo.convert(f_nex, "nexus", f_xml, "phyloxml")


def filter_subtree(in_list):
    # sub_list: 需要被子树处理的分类元，例如xxx|subtree，yyy|subtree，sublist存储为[xxx,yyy]
    sub_list = []
    # no_sub_list: 可以直接分析谱系的分类元，无需生成子树
    no_sub_list = []
    for i in in_list:
        match_obj = re.match(".+\|subtree$", i)
        if match_obj:
            # 包含|subtree
            sub_list.append(match_obj.group())
        else:
            no_sub_list.append(i)
    return sub_list, no_sub_list


"""
将subtree|前的id或名称转换成对应id，加入total_id列表
"""


def find_sub_id(item_list, cursor):
    sub_id = []
    # item_list:['11|subtree', '22|subtree', '33|subtree']
    bar_length = len(item_list)
    print(f"Searching for {bar_length} items for subtree… ")
    with alive_bar(bar_length) as bar:
        for raw_name in item_list:
            sciname = raw_name.split("|subtree")[0]
            isdigit = False
            if sciname.isdigit():
                isdigit = True
            rank_name, rank_query = get_rank(sciname, cursor, isdigit)
            if rank_name == "" or rank_query == "":
                # not find 存在亚种、亚门等非标准分类情况
                print(f"{sciname} not find in db!")
                continue
            id_list = get_id_list_by_rank(rank_name, rank_query, cursor)
            sub_id.extend(id_list)
            bar()
    return sub_id


# def open_scv(name):
#     # 新建scv文件
#     csv_name = name + '_items.csv'
#     f = open(csv_name, 'w', encoding='utf-8', newline='')
#     #  2.基于文件对象构建csv写入对象
#     csv_write = csv.writer(f)
#     #  3.构建列表头
#     csv_write.writerow(['input_term', 'domain', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'])
#     return f, csv_write

def write_no_match(null_id, null_name, multi_names, prefix):
    fname = os.path.join(prefix, 'no_match_result.txt')
    with open(fname, 'w') as f:
        if len(null_id) > 0:
            f.write(f"------These {len(null_id)} taxids are not in taxonomy database-----\n")
            for item in null_id:
                f.write(item + '\n')
            f.write('\n\n\n')
        if len(null_name) > 0:
            f.write(f"------These {len(null_name)} names are not in taxonomy database-----\n")
            for item in null_name:
                f.write(item + '\n')
            f.write('\n\n\n')
        if len(multi_names) != 0:
            print(
                f"------These {len(multi_names)} names these names are ambiguous because they have multiple taxids, you re advised to check and input the specific taxid! ------\n")
            for item in multi_names:
                f.write(item + '\n')


def process_line(line, in_list):
    if line.endswith("|subtree"):
        # 不能用\W直接替换全部特殊字符，不然｜subtree无法识别
        line_replaced = re.sub(r'_+', '_', re.sub(r'[:\s\(\)]+', '_', line))
    else:
        line_replaced = re.sub(r'_+', '_', re.sub(r'\W', '_', line))
    line_delete_before_after = re.sub('^_|_$', '', line_replaced)
    in_list.append(line_delete_before_after)
    return in_list




def main(args):
    global interrupt_flag
    start = time.perf_counter()
    # 建立数据库连接
    conn, cursor = connect_db()

    item_input = str(args.items)
    # 输入的字符串冒号和空白字符替换为下划线，与数据库格式匹配
    sub_name_or_id = args.subtree
    # in_str = input("Enter taxid or names, separate the entries with commas, or enter xx|subtree for subtree : ")
    in_list = []
    total_id = []


    # 写入csv


    if args.interrupt:
        interrupt_flag = 1
        if args.p:
            interrupt_level = 'p'
        if args.c:
            interrupt_level = 'c'
        if args.o:
            interrupt_level = 'o'
        if args.f:
            interrupt_level = 'f'
        if args.g:
            interrupt_level = 'g'
        if args.s:
            interrupt_level = 's'

    if sub_name_or_id:
        #  找出idlist
        # sciname = in_str.split("|")[0]
        sciname = sub_name_or_id
        isdigit = False
        if sciname.isdigit():
            isdigit = True
        rank_name, rank_query = get_rank(sciname, cursor, isdigit)
        if rank_name == "" or rank_query == "":
            # not find 存在亚种、亚门等非标准分类情况
            cursor.close()  # 先关闭游标
            conn.commit()
            conn.close()  # 再关闭数据库连接
            print(f"{sciname} is not standard taxonomic unit for input")
            sys.exit(0)
        id_list = get_id_list_by_rank(rank_name, rank_query, cursor)
        all_id = id_list
    elif args.file:
        # 输入形式是txt
        # command: -file
        try:
            with open(args.file, 'r') as f:
                in_list = []
                for line in f.readlines():
                    line = str(line.strip())
                    in_list = process_line(line, in_list)
                sub_list, no_sub_list = filter_subtree(in_list)
                print("no_sub_list", no_sub_list)
                # 生成替换name文件
                no_sub_id, name_sciname_map = get_id_and_sciname_list_by_id_or_name(no_sub_list, cursor)
                sub_id = find_sub_id(item_list=sub_list, cursor=cursor)
                all_id = list(set(no_sub_id + sub_id))
        except Exception as e:
            print("file does not exist!", e)
            exit(0)
    else:
        # command: -input
        lines = item_input.split(",")
        in_list = []
        for line in lines:
            in_list = process_line(line, in_list)
        sub_list, no_sub_list = filter_subtree(in_list)
            # 生成替换name文件
        # name_sciname_map: 生成替换name文件
        # id_list: 纯id和纯sciname转化为id的列表
        # sub_list: ['11|subtree', '33|subtree']
        # no_sub_list: ['Homo_sapiens', 'Mus_musculus', 'Gallus_gallus', 'Drosophila_melanogaster', 'Escherichia_coli']

        no_sub_id, name_sciname_map = get_id_and_sciname_list_by_id_or_name(no_sub_list, cursor)
        sub_id = find_sub_id(sub_list, cursor=cursor)

        # print(no_sub_id)
        # print(f"sub_id:{sub_id}")
        # sub_id:[11, 593907, 33, 483219, 675525, 1334629]
        # 不生成替换name
        # id_list = get_id_list_by_id_or_name(in_list, cursor)
        # all_id: id_list 和 sublist合并去重，即混合输入subtree和普通id混合后全部id
        all_id = list(set(no_sub_id + sub_id))

    """
    数据库连接后操作
    """

    """
    确定文件输出路径
    """
    my_prefix = Path(args.prefix)  # 输出文件路径前缀
    my_fname = args.fname
    if my_prefix.is_dir():  # 路径存在
        pass
    else:
        my_prefix = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), 'iphylo_files'))
        # not_exist_msg = f"The path '{args.prefix}' does not exist!"
        print(f"The path '{args.prefix}' does not exist!")

    # 指定路径下创建一级文件夹与fname同名用于存放此次运行生成的数据
    my_prefix = my_prefix / Path(my_fname)
    isExists = os.path.exists(my_prefix)
    # 判断结果
    if not isExists:
        os.makedirs(my_prefix)

    # my_prefix_name: 用户通过arg指定的输出路径和文件名称拼接，不包含后缀
    # my_prefix_name 更名为 outpath_fname
    outpath_fname = os.path.join(my_prefix, my_fname)
    txt_path = outpath_fname + ".txt"

    # print("outpath_fname",outpath_fname)
    # 新建scv文件
    """
    根据输入新建SCV文件
    """
    csv_name = outpath_fname + '_items.csv'
    f_csv = open(csv_name, 'w', encoding='utf-8', newline='')
    #  2.基于文件对象构建csv写入对象
    csv_write = csv.writer(f_csv)
    #  3.构建列表头
    csv_write.writerow(['taxid', 'input_term', 'domain', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'])
    anno_name = outpath_fname + '_items_for_anno.csv'

    """
    是否打断，画树
    """
    if interrupt_flag == 1:
        # id_list_to_newick(all_id, cursor, interrupt_level, anno_name=anno_name, path=txt_path, csv_write=csv_write)
        try:
            sciname_leaf_map = id_list_to_newick(all_id, cursor, interrupt_level, anno_name=anno_name, path=txt_path, csv_write=csv_write)
        except Exception:
            print('Parameter for interrupting needed! Seek help by -h, --help')
    else:
        sciname_leaf_map = id_list_to_newick(all_id, cursor, anno_name, anno_name=anno_name, path=txt_path, csv_write=csv_write)

    map_csv_path = outpath_fname + "_name_map.csv"
    # 写入sciname_leaf_map到CSV文件
    with open(map_csv_path, 'w', newline='') as csvfile:
        fieldnames = ['Input Scientific Name', 'Leaf Node in Tree']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for name, leaf in sciname_leaf_map.items():
            writer.writerow({'Input Scientific Name': name, 'Leaf Node in Tree': leaf})

    # 生成其他格式的树
    try:
        convert_file(prefix_fname=outpath_fname, args=args)
    except Exception as e:
        print("convert fail: ", e)

    print("Generate Tree Success!")

    try:
        # 画ascii树
        draw_ascii_tree(outpath_fname)
    except Exception:
        print("Ascii Tree Failed to Load. Maybe because the tree is too large.")

    """
    show tree structure in pdf
    """
    try:
        pdf_name = outpath_fname + '_structure.pdf'
        myete.render_pdf(txt_path, pdf_name)
    except Exception:
        print("Tree in PDF Failed to Load. You may check out the newick format.")




    # os.remove(txt_path)

    if len(null_id) != 0:
        print(f"------These {len(null_id)} taxids are not in taxonomy database------\n{null_id}")
    if len(null_name) != 0:
        # print("------These names are not in taxonomy database------\n", null_name)
        print(f"------These {len(null_name)} names are not in taxonomy database------\n{null_name}")
    if len(multi_names) != 0:
        print(
            f"------These {len(multi_names)} names these names are ambiguous because they have multiple taxids, you re advised to check and input the specific taxid! ------\n{multi_names}")


    # 将上述信息写入文件，文件名no_match_result.txt
    write_no_match(null_id, null_name, multi_names, my_prefix)

    try:
        print(f"{len(name_sciname_map)} names are replaced!")
        replace_name_record(prefix_fname=outpath_fname, map=name_sciname_map)
    except:
        print("no alter map")
    # print(name_sciname_map)
    # 关闭数据库连接，提交事务
    commit_close_db(conn, cursor)
    #  关闭csv文件
    f_csv.close()
    #
    # # 画图R脚本路径
    # r_plot_script = os.path.join(os.getcwd(), 'RScript/phylo/drawTree.R')
    # # 画图
    #
    # color = args.color
    # if args.color is not None:
    #     plot.draw_chem_plot(r_plot_script, txt_path, my_prefix, fname=my_fname, csv=csv_name, mycolor=color)
    #     # print(args)
    # else:
    #     plot.draw_chem_plot(r_plot_script, txt_path, my_prefix, fname=my_fname, csv=csv_name)

    # 打印程序运行时间
    end = time.perf_counter()
    print('Running time: %s Seconds' % (end - start))

