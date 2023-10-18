# from ete3 import Tree, TreeStyle
import os
import sqlite3
import re
import sys
import time
from Bio import Phylo
from pathlib import Path
import csv
from alive_progress import alive_bar

from scripts.Utils.download_db import download_db_file
# import eteUtils.render_test
# from eteUtils import tree_style
from scripts.Utils.myArg import chem_arg
import pandas as pd
from scripts import phylotree as pt

"""
DB Link: 111
"""
# 不在数据库里的id

from scripts.Utils import newick
from scripts.Utils.eteUtils import myete

"""
SQL语句
sqlite3占位符语法
"""
# 输入为id查询化合物分类
chem_tabel = "inchikey_info"
sql_find_classify_by_id = f"""select kingdom, superclass, class, subclass from inchikey_info where id = :id """
# sql_find_id_and_sciname_by_name = """select taxid, sciname from name_id_sci where name = :name """
sql_find_id_by_inchikey = """select id from inchikey_info where inchikey = :inchikey """
sql_find_id_by_inchi = """select id from inchikey_info where inchi = :inchi """
sql_find_id_by_smiles = """select id from inchikey_info where SMILES = :smiles """
sql_all_info = """select * from inchikey_info where id = :id """
sql_table_info = """PRAGMA table_info(inchikey_info);"""

sql_rank_by_name = """select rank from t_rank where name = :name """


def connect_db():
    """
    初始化数据库连接，打印cursor确保连接成功
    """
    # 数据库路径
    db_file = 'data/ichem.db'
    if not os.path.exists(db_file):
        user_input = input(f"Necessary database resources (587.3 MB) not found. Do you want to download it? (y/n): ")
        if user_input.lower() == 'y':
            url = 'http://222.190.246.202:35000/resource/ichem_db'
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
def get_id_by_id_or_inchikey(in_list, cursor, bar):
    id_list = []
    none_inchikey = []
    for item in in_list:
        # print(item)
        if item.isdigit():
            id_list.append(item)
        else:
            # 不考虑inchikey重复
            cursor.execute(sql_find_id_by_inchikey, {"inchikey": item})
            out_tuple = cursor.fetchall()
            if out_tuple:  # 查询到inchikey对应的id
                _id = out_tuple[0][0]
                id_list.append(_id)
                # _end = time.perf_counter()
                # print('Running time: %s Seconds' % (_end - _start))
            else:
                # 查询的inchikey不在数据库中
                none_inchikey.append(item)
                continue
        bar()
    return id_list, none_inchikey

    """
    需求：在分类末尾加入id
    """


def merge_lists(lists):
    result = {}
    for lst in lists:
        key = tuple(lst[:-1])
        value = str(lst[-1])
        if key in result:
            result[key].append(value)
        else:
            result[key] = [value]

    return [[*k, '_'.join(v)] for k, v in result.items()]


def id_list_to_newick(id_list, cursor, *interrupt_level, path, csv_write):
    """
    根据id_list 转换成 newick格式
    :param id_list: 化合物对应数据库中id
    :return:
    """
    # lin_list: 一个二维数组，访问第1个物种的第2个rank应该用lin_list[0][1]
    lin_list = []
    none_classify_id = []

    id_list_unique = set(id_list)  # 去重
    for item in id_list_unique:
        int_item = int(item)
        cursor.execute(sql_find_classify_by_id, {"id": int_item})
        out_tuple = cursor.fetchall()
        if not out_tuple:
            continue
        full_list = list(out_tuple[0])

        # print(f'full list origin:{full_list}')
        # 写CSV文件 all info
        cursor.execute(sql_all_info, {"id": int_item})
        all_info = cursor.fetchall()[0]
        csv_write.writerow(all_info)
        # cursor.execute(sql_name_by_id, {"id": int(item)})
        # name = cursor.fetchone()[0]
        # print(name)

        if interrupt_level:
            # 如果需要在某level打断，找到该level对应的index
            # level_list = ["super", "class"]
            level_list = ["superclass", "class", "subclass"]
            for idx, level in enumerate(level_list):
                if interrupt_level[0] == level:
                    # 如果打断为superclass， level_idx = 2， full_list[:2] = [kingdom, superlcass]
                    level_idx = idx + 2
            full_list = full_list[:level_idx]
        # print(f'full list:{full_list}')
        # 过滤none
        res = list(filter(None, full_list))

        if len(res) > 0:
            # 过滤掉分类全是None的情况
            res.append(str(int_item))
            lin_list.append(res)
        else:
            none_classify_id.append(int_item)

    # 二维列表去重，需要先转化为tuple
    if interrupt_level:
        lin_list = list(set([tuple(t) for t in lin_list]))
        lin_list = [list(v) for v in lin_list]
    # print(f'items:{lin_list}')

    merged_lin_list = merge_lists(lin_list)
    newick.showtree(data=lin_list, path=path)  # id分支不合并
    # newick.showtree(data=merged_lin_list, path=path)  # id分支合并
    return none_classify_id
    # 生成树文件


def get_rank(item, cursor):
    """

    :param item: chemical name without |subtree, eg: Glycerolipids
    :param cursor:  database cursor
    :return: rank name: class
    """
    cursor.execute(sql_rank_by_name, {"name": item})
    result = cursor.fetchone()
    if result is not None:
        rank_name = result[0]
        return rank_name
    else:
        return None


def get_id_list_by_rank(rank_name, chemical_name, cursor):
    """

    :param rank_name: class
    :param chemical_name: Glycerolipids
    :param cursor:
    :return:
    """
    sql_find_id_by_rank = "SELECT id from inchikey_info WHERE " + str(rank_name) + " = '" + str(chemical_name) + "'"
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


#
# def filter_subtree(in_list):
#     # sub_list: 需要被子树处理的分类元，例如xxx|subtree，yyy|subtree，sublist存储为[xxx,yyy]
#     sub_list = []
#     # no_sub_list: 可以直接分析谱系的分类元，无需生成子树
#     no_sub_list = []
#     for i in in_list:
#         match_obj = re.match(".+\|subtree$", i)
#         if match_obj:
#             # 包含|subtree
#             sub_list.append(match_obj.group())
#         else:
#             no_sub_list.append(i)
#     return sub_list, no_sub_list


def process_name(name):
    name_replaced = re.sub(r'_+', '_', re.sub(r'\W', '_', name))
    name_delete_before_after = re.sub('^_|_$', '', name_replaced)
    return name_delete_before_after


"""
将subtree|前的id或名称转换成对应id，加入total_id列表
"""


def find_sub_id(item_list, cursor, bar):
    sub_id = []
    subtree_wrong_name = []
    # item_list:['11|subtree', '22|subtree', '33|subtree']
    for raw_name in item_list:
        chemical_name = raw_name.split("|subtree")[0]
        chemical_name = process_name(chemical_name)
        rank_name = get_rank(chemical_name, cursor)
        # rank name: class, chemical_name: Glycerolipids
        if rank_name:
            id_list = get_id_list_by_rank(rank_name, chemical_name, cursor)
            sub_id.extend(id_list)
        else:
            subtree_wrong_name.append(chemical_name)
        bar()
    return sub_id, subtree_wrong_name


# def open_scv(name):
#     # 新建scv文件
#     csv_name = name + '_items.csv'
#     f = open(csv_name, 'w', encoding='utf-8', newline='')
#     #  2.基于文件对象构建csv写入对象
#     csv_write = csv.writer(f)
#     #  3.构建列表头
#     csv_write.writerow(['input_term', 'domain', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'])
#     return f, csv_write


# 判断一个字符串是否是inchikey
def is_valid_inchikey(inchikey):
    inchikey_pattern = re.compile(r'^[A-Z]{14}-[A-Z]{10}-[A-Z]$')
    return re.match(inchikey_pattern, inchikey) is not None


# 判断是否为有效inchi
def is_valid_inchi(inchi):
    inchi_pattern = re.compile(r'^InChI=1S?/[A-Za-z0-9]+')
    return re.match(inchi_pattern, inchi) is not None


def get_id_from_file(file, cursor):
    compound_no_result = []  # 数据库中查询不到的化合物
    _id = []  # 按行查询到的化合物对应id
    sub_and_no_sub_list = []
    with open(file, 'r') as f:
        lines = f.readlines()
        num_lines = len(lines)
        with alive_bar(num_lines) as bar:
            for line in lines:
                line = str(line.strip())  # 去除换行符
                sub_and_no_sub_list.append(line)
            sub_list, no_sub_list = pt.filter_subtree(sub_and_no_sub_list)
            # print(f"sublist: {sub_list}\n nosublist: {no_sub_list}")
            sub_id, subtree_wrong_name = find_sub_id(sub_list, cursor, bar)

            for no_sub in no_sub_list:
                try:
                    if is_valid_inchikey(no_sub):  # 行内容为inchikey
                        cursor.execute(sql_find_id_by_inchikey, {"inchikey": no_sub})
                        id = cursor.fetchone()[0]
                    elif is_valid_inchi(no_sub):
                        cursor.execute(sql_find_id_by_inchi, {"inchi": no_sub})
                        id = cursor.fetchone()[0]
                    else:
                        cursor.execute(sql_find_id_by_smiles, {"smiles": no_sub})
                        id = cursor.fetchone()[0]
                except:
                    compound_no_result.append(no_sub)
                else:
                    _id.append(id)
                bar()
    id_unique = list(set(_id + sub_id))
    return id_unique, compound_no_result, subtree_wrong_name


def get_table_column(cursor):
    cursor.execute(sql_table_info)
    info = cursor.fetchall()
    column_name = [tuple[1] for tuple in info]
    return column_name


def no_match_output(none_inchikey, none_classify_id, subtree_wrong_name, prefix, csv):
    fname = os.path.join(prefix, 'no_match_result.txt')
    with open(fname, 'w') as f:
        if len(none_inchikey) > 0:
            f.write(f"------These InChIKeys could not be found in the database------\n")
            for inchikey in none_inchikey:
                f.write(inchikey + '\n')
        if len(subtree_wrong_name) > 0:
            f.write(f"------These subtree could not be found in the database------\n")
            for name in subtree_wrong_name:
                f.write(name + '\n')
        if len(none_classify_id) > 0:
            try:
                # 读取csv文件
                data = pd.read_csv(csv)
                # 筛选出对应的行
                filtered_data = data[data['id'].isin(none_classify_id)]

                # 提取inchikey列并转换为列表  错误，不一定具有inchikey

                none_classified_chemicals_list = filtered_data.to_string()
                print(none_classified_chemicals_list)
                f.write('\n\n\n------No classification information found for these chemicals------\n')
                f.write(none_classified_chemicals_list)

            except FileNotFoundError:
                f.write('\n\n\n------No classification information found for these IDs------\n')
                for nc_inchikey in none_classify_id:
                    f.write(nc_inchikey + '\n')


# if __name__ == '__main__':
# 计时
def main(args):
    start = time.perf_counter()
    # 建立数据库连接
    conn, cursor = connect_db()
    # parser = chem_arg.get_parser()
    # args = parser.parse_args()

    # 输入的字符串冒号和空白字符替换为下划线，与数据库格式匹配
    # in_str = input("Enter taxid or names, separate the entries with commas, or enter xx|subtree for subtree : ")
    in_list = []
    total_id = []
    # 画图R脚本路径
    # r_plot_script = os.path.join(os.getcwd(), '../RScript/chem/chemplot.R')

    if args.file:
        # 输入形式是txt
        # command: -file
        try:
            all_id, none_inchikey, subtree_wrong_name = get_id_from_file(args.file, cursor)
        except Exception as e:
            print(e)
            print("file does not exist!")
            exit(0)
    elif args.items:

        item_input = str(args.items)
        # command: -input
        # 只输入inchikey
        in_list = item_input.split(",")
        bar_length = len(in_list)
        with alive_bar(bar_length) as bar:
            sub_list, no_sub_list = pt.filter_subtree(in_list)
            sub_id, subtree_wrong_name = find_sub_id(sub_list, cursor, bar)
            # print(f"sublist: {sub_list}\n nosublist: {no_sub_list}")
            # sublist: ['PWYBSBFLFZQRPP-UHIUAQIISA-N|subtree']
            # nosublist: ['WRXZBGOZXSMKHW-SGOGUMGMSA-N']

            """ 
            name_sciname_map: 生成替换name文件
            id_list: 纯id和纯sciname转化为id的列表
            sub_list: ['11|subtree', '33|subtree']
            no_sub_list: ['Homo_sapiens', 'Mus_musculus', 'Gallus_gallus', 'Drosophila_melanogaster', 'Escherichia_coli']
            """
            no_sub_id, none_inchikey = get_id_by_id_or_inchikey(no_sub_list, cursor, bar)
            # sub_id = find_sub_id(sub_list, cursor=cursor)
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
    column_name = get_table_column(cursor)
    # column_name: type list
    #  3.构建列表头
    csv_write.writerow(column_name)

    if args.interrupt:
        interrupt_flag = 1
        if args.Super:
            interrupt_level = 'superclass'
        elif args.Class:
            interrupt_level = 'class'
        elif args.Sub:
            interrupt_level = 'subclass'
    else:
        interrupt_flag = 0

    """
    是否打断，画树
    """
    if interrupt_flag == 1:
        try:
            none_classify_id = id_list_to_newick(all_id, cursor, interrupt_level, path=txt_path, csv_write=csv_write)
        except Exception:
            print('Parameter for interrupting needed! Seek help by -h, --help')
    else:
        none_classify_id = id_list_to_newick(all_id, cursor, path=txt_path, csv_write=csv_write)

    # # 打印树叶节点数
    # try:
    #     leaf_num = myete.get_leaf_num(txt_path)
    #     print(f"{leaf_num} leaves in tree")
    # except Exception as e:
    #     print(e)

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

    # try:
    #     leaf_num = myete.get_leaf_num(txt_path)
    #     print(f"{leaf_num} leaves in tree")
    # except Exception as e:
    #     print(e)

    """
    show tree structure in pdf
    """
    try:
        pdf_name = outpath_fname + '_structure.pdf'
        myete.render_pdf(txt_path, pdf_name)
    except Exception:
        print("Tree in PDF Failed to Load. You may check out the newick format.")

    # os.remove(txt_path)
    try:
        if len(none_inchikey) != 0:
            print(f"------These {len(none_inchikey)} inchikeys are not in database------\n{none_inchikey}")
    except:
        pass

    try:
        if len(subtree_wrong_name) != 0:
            print(f"------These {len(subtree_wrong_name)} subtrees are not in database------\n{subtree_wrong_name}")
    except:
        pass

    if len(none_classify_id) != 0:
        # print("------These names are not in taxonomy database------\n", null_name)
        print(f"------These {len(none_classify_id)} inchikeys are not classified currently------\n{none_classify_id}")
    # if len(compound_not_find) != 0:
    #     print(f"------These {len(compound_not_find)} items are not recorded as InChI, InChIKey, or SMILES------\n{compound_not_find}")
    # # 关闭数据库连接，提交事务
    # commit_close_db(conn, cursor)
    # #  关闭csv文件
    f_csv.close()
    no_match_output(none_inchikey, none_classify_id, subtree_wrong_name, my_prefix, csv_name)
    # color = args.color
    #
    #
    # #  画图，先注释掉吧，最后再说
    # if args.color is not None:
    #     plot.draw_chem_plot(r_plot_script, txt_path, my_prefix, fname=my_fname, csv=csv_name, mycolor=color)
    #     # print(args)
    # else:
    #     plot.draw_chem_plot(r_plot_script, txt_path, my_prefix, fname=my_fname, csv=csv_name)
    #

    # 打印程序运行时间
    end = time.perf_counter()
    print('Running time: %s Seconds' % (end - start))
