# from ete3 import Tree, TreeStyle
import os
import queue
import sqlite3
import re
import sys
import threading
import time
import random

import requests
from Bio import Phylo
from pathlib import Path
import csv

# import eteUtils.render_test
# from eteUtils import tree_style
from scripts.Utils.myArg import chem_arg_online
import pandas as pd
from alive_progress import alive_bar

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


class myThread(threading.Thread):
    def __init__(self, name, q, error_q, result_q, inchikey_not_find, bar):
        threading.Thread.__init__(self)
        self.name = name
        self.q = q
        self.result_q = result_q
        self.error_q = error_q
        self.inchikey_not_find = inchikey_not_find
        self.bar = bar

    def run(self):
        # print("Starting " + self.name)
        while True:
            try:
                crawl(self.name, self.q, self.error_q, self.result_q, self.inchikey_not_find)
                self.bar()

                # self.q.join()
                # while self.result_q.qsize():  # 返回队列的大小
                #     res = self.result_q.get()
                #     data = [item.values() for item in res]
                #     print(data)

            except:
                break


def connect_db():
    """
    初始化数据库连接，打印cursor确保连接成功
    """
    try:
        # 连接数据库
        conn = sqlite3.connect('data/iphylo.db')
        # 创建游标
        cursor = conn.cursor()
    except Exception:
        print("------database connect failed!------")
    else:
        print("------database connect success!------")
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
def get_id_by_id_or_inchikey(in_list, cursor):
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

    return id_list, none_inchikey

    """
    需求：在分类末尾加入id
    """


def merge_lists(lists):
    """
    接受一个列表的列表 lists，将具有相同前缀的列表合并成一个新的列表，新列表中的最后一个元素是原来每个列表的最后一个元素的字符串形式以 "_" 连接的结果。
    :param lists:二维列表
    :return: lists
    """
    result = {}
    for lst in lists:
        key = tuple(lst[:-1])
        value = str(lst[-1])
        if key in result:
            result[key].append(value)
        else:
            result[key] = [value]

    return [[*k, '_'.join(v)] for k, v in result.items()]


def id_list_to_newick(classify_lis, *interrupt_level, path, csv_write):
    lin_list = []
    none_classify_id = []
    for idx, item in enumerate(classify_lis):

        # 每条分类，list
        classy_info = item[0:-1]  # 除掉最后一个，即分类信息
        inchikey = item[-1]
        _id = idx + 1
        row = [_id] + [inchikey] + classy_info
        csv_write.writerow(row)

        # "注意要写入csv"

        # cursor.execute(sql_name_by_id, {"id": int(item)})
        # name = cursor.fetchone()[0]
        # print(name)
        if interrupt_level:
            # 如果需要在某level打断，找到该level对应的index
            level_list = ["superclass", "class", "subclass"]
            for idx, level in enumerate(level_list):
                if interrupt_level[0] == level:
                    # 如果打断为superclass， level_idx = 2， full_list[:2] = [kingdom, superlcass]
                    level_idx = idx + 2
            classy_info = classy_info[:level_idx]
        #      print(classy_info)
        # 'Organic_compounds', 'Lipids_and_lipid_like_molecules']
        # ['Organic_compounds', 'Organic_oxygen_compounds']
        # ['Organic_compounds', 'Organic_oxygen_compounds']
        # ['Organic_compounds', 'Organic_oxygen_compounds']
        #
        # print(f'full list:{full_list}')
        # 过滤none
        full_list = list(filter(None, classy_info))

        if len(classy_info) > 0:
            # 过滤掉分类全是None的情况
            full_list.append(str(_id))  # inchikey
            lin_list.append(full_list)
        else:
            # these inchikey found no classification
            none_classify_id.append(item[-1])

    # 二维列表去重，需要先转化为tuple
    if interrupt_level:
        lin_list = list(set([tuple(t) for t in lin_list]))
        lin_list = [list(v) for v in lin_list]
    # print(f'items:{lin_list}')

    # 合并分类相同，inchikey不同的分支
    # merged_lin_list = merge_lists(lin_list)
    newick.showtree(data=lin_list, path=path)  # id分支不合并
    # newick.showtree(data=merged_lin_list, path=path)  # id分支合并
    return none_classify_id
    # 生成树文件


# def get_rank(item, cursor, isdigit):
#     rank_list = ["_superkingdom", "_kingdom", "_phylum", "_class", "_order", "_family", "_genus", "_species"]
#     rank_query_list = ["d__", "k__", "p__", "c__", "o__", "f__", "g__", "s__"]
#     if isdigit:
#         cursor.execute(sql_rank_by_id, {"name": item})
#     else:
#         cursor.execute(sql_rank_by_name, {"name": item})
#     query_list = cursor.fetchone()
#     rank_name = ""  # _phylum
#     rank_query = ""  # p__Ascomycota
#
#     if query_list is not None:
#         sciname = query_list[1]
#         for idx, item in enumerate(query_list[2:]):
#             if rank_query_list[idx] + sciname == item:
#                 rank_name = rank_list[idx]
#                 rank_query = item
#                 break
#
#     return rank_name, rank_query


def get_id_list_by_rank(rank_name, rank_query, cursor):
    """

    :param rank_name: _phylum
    :param rank_query: p__Ascomycota
    :param cursor:
    :return:
    """
    sql_find_id_by_rank = "SELECT id from id2lineage_copy1 WHERE " + str(rank_name) + " = '" + str(rank_query) + "'"
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
            print(f'Tree and ASCII tree is saved with path and name prefix: {prefix_fname}')
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


"""
将subtree|前的id或名称转换成对应id，加入total_id列表
"""


# def find_sub_id(item_list, cursor):
#     sub_id = []
#     # item_list:['11|subtree', '22|subtree', '33|subtree']
#     for raw_name in item_list:
#         sciname = raw_name.split("|")[0]
#         isdigit = False
#         if sciname.isdigit():
#             isdigit = True
#         rank_name, rank_query = get_rank(sciname, cursor, isdigit)
#         if rank_name == "" or rank_query == "":
#             # not find 存在亚种、亚门等非标准分类情况
#             print(f"{sciname} not find in db!")
#             continue
#         id_list = get_id_list_by_rank(rank_name, rank_query, cursor)
#         sub_id.extend(id_list)
#     return sub_id


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
    with open(file, 'r') as f:
        for line in f.readlines():
            line = str(line.strip())  # 去除换行符
            try:
                if is_valid_inchikey(line):  # 行内容为inchikey
                    cursor.execute(sql_find_id_by_inchikey, {"inchikey": line})
                    id = cursor.fetchone()[0]
                elif is_valid_inchi(line):
                    cursor.execute(sql_find_id_by_inchi, {"inchi": line})
                    id = cursor.fetchone()[0]
                else:
                    cursor.execute(sql_find_id_by_smiles, {"smiles": line})
                    id = cursor.fetchone()[0]
            except:
                compound_no_result.append(line)
            else:
                _id.append(id)
    id_unique = list(set(_id))
    return id_unique, compound_no_result


def get_table_column(cursor):
    cursor.execute(sql_table_info)
    info = cursor.fetchall()
    column_name = [tuple[1] for tuple in info]
    return column_name


def no_match_output(none_inchikey, none_classify_id, prefix, csv):
    fname = os.path.join(prefix, 'no_match_result.txt')
    with open(fname, 'w') as f:
        if len(none_inchikey) > 0:
            f.write(f"------These InChIKeys could not be found in the database------\n")
            for inchikey in none_inchikey:
                f.write(inchikey + '\n')
        if len(none_classify_id) > 0:
            try:
                # 读取csv文件
                data = pd.read_csv(csv)
                # 筛选出对应的行
                filtered_data = data[data['id'].isin(none_classify_id)]
                # 提取inchikey列并转换为列表
                none_classified_inchikey_list = filtered_data['inchikey'].tolist()
                f.write('\n\n\n------No classification information found for these InChIKeys------\n')
                for nc_inchikey in none_classified_inchikey_list:
                    f.write(nc_inchikey + '\n')
            except FileNotFoundError:
                f.write('\n\n\n------No classification information found for these IDs------\n')
                for nc_inchikey in none_classify_id:
                    f.write(nc_inchikey + '\n')


def q_work(inchikey_lis, thread_num):
    res_list = []
    inchikey_not_find = []  # ClassyFire中不存在分类的
    equal_times = 2
    cnt = 0
    workQueue = queue.Queue(len(inchikey_lis))
    for inchikey in inchikey_lis:
        workQueue.put(inchikey)
    equal_list = []  # 判断error_q equal_times次相等
    while True:
        with alive_bar(workQueue.qsize()) as bar:
            cnt += 1
            resq = queue.Queue()
            error_q = queue.Queue()  # 网络错没查到的
            threads = []
            print(f'Search round {cnt},  {workQueue.qsize()} items to process')
            for i in range(1, thread_num + 1):
                # 创建n个新线程
                thread = myThread("Thread-" + str(i), q=workQueue, error_q=error_q, result_q=resq, inchikey_not_find=inchikey_not_find, bar=bar)
                # 判断两次循环error_q长度不变则结束循环，输出查不到的inchikey到txt

                # 开启新线程
                thread.start()
                # 添加新线程到线程列表
                threads.append(thread)

            # 等待所有线程完成
            for thread in threads:
                thread.join()

            # print("Exiting Main Thread")

            # 输出列表
            while resq.qsize():  # 返回队列的大小
                res = (resq.get())  # 将结果存入数据库中
                res_list.append(res)
            workQueue = error_q

            if len(equal_list) == equal_times:
                if all(x == equal_list[0] for x in equal_list):
                    break
                equal_list[cnt % equal_times] = error_q.qsize()  # 取余循环覆写
            else:
                equal_list.append(error_q.qsize())

    return res_list, inchikey_not_find


def replace_non_alphanumeric(s):
    if s is not None:
        if isinstance(s, str):
            return re.sub(r'\W+', '_', s)
        else:
            return re.sub(r'\W+', '_', s['name'])
    else:
        return None


def crawl(name, q, error_q, result_q, inchikey_no_classify):
    inchikey = q.get(timeout=2)

    # url_feihn = "https://cfb.fiehnlab.ucdavis.edu/entities/{}.json".format(inchikey)
    url_whishart = 'http://classyfire.wishartlab.com/entities/' + inchikey + '.json'

    # 创建UA池
    # 浏览器代理
    agents = [
        "Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:17.0; Baiduspider-ads) Gecko/17.0 Firefox/17.0",
        "Mozilla/5.0 (Linux; U; Android 2.3.6; en-us; Nexus S Build/GRK39F) AppleWebKit/533.1 (KHTML, like Gecko) Version/4.0 Mobile Safari/533.1",
        "Avant Browser/1.2.789rel1 (http://www.avantbrowser.com)",
        "Mozilla/5.0 (Windows; U; Windows NT 6.1; en-US) AppleWebKit/532.5 (KHTML, like Gecko) Chrome/4.0.249.0 Safari/532.5",
        "Mozilla/5.0 (Windows; U; Windows NT 5.2; en-US) AppleWebKit/532.9 (KHTML, like Gecko) Chrome/5.0.310.0 Safari/532.9",
        "Mozilla/5.0 (Windows; U; Windows NT 5.1; en-US) AppleWebKit/534.7 (KHTML, like Gecko) Chrome/7.0.514.0 Safari/534.7",
        "Mozilla/5.0 (Windows; U; Windows NT 6.0; en-US) AppleWebKit/534.14 (KHTML, like Gecko) Chrome/9.0.601.0 Safari/534.14",
        "Mozilla/5.0 (Windows; U; Windows NT 6.1; en-US) AppleWebKit/534.14 (KHTML, like Gecko) Chrome/10.0.601.0 Safari/534.14",
        "Mozilla/5.0 (Windows; U; Windows NT 6.1; en-US) AppleWebKit/534.20 (KHTML, like Gecko) Chrome/11.0.672.2 Safari/534.20",
        "Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/534.27 (KHTML, like Gecko) Chrome/12.0.712.0 Safari/534.27",
        "Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/535.1 (KHTML, like Gecko) Chrome/13.0.782.24 Safari/535.1",
        "Mozilla/5.0 (Windows NT 6.0) AppleWebKit/535.2 (KHTML, like Gecko) Chrome/15.0.874.120 Safari/535.2",
        "Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/535.7 (KHTML, like Gecko) Chrome/16.0.912.36 Safari/535.7",
        "Mozilla/5.0 (Windows; U; Windows NT 6.0 x64; en-US; rv:1.9pre) Gecko/2008072421 Minefield/3.0.2pre",
        "Mozilla/5.0 (Windows; U; Windows NT 5.1; zh-CN; rv:1.9b4) Gecko/2008030317 Firefox/3.0b4",
        "Mozilla/5.0 (Windows; U; Windows NT 5.1; en-US; rv:1.9.0.10) Gecko/2009042316 Firefox/3.0.10",
        "Mozilla/5.0 (Windows; U; Windows NT 6.0; en-GB; rv:1.9.0.11) Gecko/2009060215 Firefox/3.0.11 (.NET CLR 3.5.30729)",
        "Mozilla/5.0 (Windows; U; Windows NT 6.0; en-US; rv:1.9.1.6) Gecko/20091201 Firefox/3.5.6 GTB5",
        "Mozilla/5.0 (Windows; U; Windows NT 5.1; tr; rv:1.9.2.8) Gecko/20100722 Firefox/3.6.8 ( .NET CLR 3.5.30729; .NET4.0E)",
        "Mozilla/5.0 (Windows; U; MSIE 6.0; Windows NT 5.1; SV1; .NET CLR 2.0.50727; BIDUBrowser 7.6)",
        "Mozilla/5.0 (Windows NT 6.3; WOW64; Trident/7.0; rv:11.0) like Gecko",
        "Mozilla/5.0 (Windows NT 6.3; WOW64; rv:46.0) Gecko/20100101 Firefox/46.0",
        "Mozilla/5.0 (Windows NT 6.3; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/45.0.2454.99 Safari/537.36",
        "Mozilla/5.0 (Windows NT 6.3; Win64; x64; Trident/7.0; Touch; LCJB; rv:11.0) like Gecko",
        "Mozilla/5.0 (Windows NT 6.1; rv:2.0.1) Gecko/20100101 Firefox/4.0.1",
        "Mozilla/5.0 (Windows NT 6.1; Win64; x64; rv:2.0.1) Gecko/20100101 Firefox/4.0.1",
        "Mozilla/5.0 (Windows NT 5.1; rv:5.0) Gecko/20100101 Firefox/5.0",
        "Mozilla/5.0 (Windows NT 6.1; WOW64; rv:6.0a2) Gecko/20110622 Firefox/6.0a2",
        "Mozilla/5.0 (Windows NT 6.1; WOW64; rv:7.0.1) Gecko/20100101 Firefox/7.0.1",
        "Mozilla/5.0 (Windows NT 6.1; WOW64; rv:2.0b4pre) Gecko/20100815 Minefield/4.0b4pre",
        "Mozilla/4.0 (compatible; MSIE 5.5; Windows NT 5.0 )",
        "Mozilla/4.0 (compatible; MSIE 5.5; Windows 98; Win 9x 4.90)",
        "Mozilla/5.0 (Windows; U; Windows XP) Gecko MultiZilla/1.6.1.0a",
        "Mozilla/2.02E (Win95; U)",
        "Mozilla/3.01Gold (Win95; I)",
        "Mozilla/4.8 [en] (Windows NT 5.1; U)",
        "Mozilla/5.0 (Windows; U; Win98; en-US; rv:1.4) Gecko Netscape/7.1 (ax)",
        "HTC_Dream Mozilla/5.0 (Linux; U; Android 1.5; en-ca; Build/CUPCAKE) AppleWebKit/528.5  (KHTML, like Gecko) Version/3.1.2 Mobile Safari/525.20.1",
        "Mozilla/5.0 (hp-tablet; Linux; hpwOS/3.0.2; U; de-DE) AppleWebKit/534.6 (KHTML, like Gecko) wOSBrowser/234.40.1 Safari/534.6 TouchPad/1.0",
        "Mozilla/5.0 (Linux; U; Android 1.5; en-us; sdk Build/CUPCAKE) AppleWebkit/528.5  (KHTML, like Gecko) Version/3.1.2 Mobile Safari/525.20.1",
        "Mozilla/5.0 (Linux; U; Android 2.1; en-us; Nexus One Build/ERD62) AppleWebKit/530.17 (KHTML, like Gecko) Version/4.0 Mobile Safari/530.17",
        "Mozilla/5.0 (Linux; U; Android 2.2; en-us; Nexus One Build/FRF91) AppleWebKit/533.1 (KHTML, like Gecko) Version/4.0 Mobile Safari/533.1",
        "Mozilla/5.0 (Linux; U; Android 1.5; en-us; htc_bahamas Build/CRB17) AppleWebKit/528.5  (KHTML, like Gecko) Version/3.1.2 Mobile Safari/525.20.1",
        "Mozilla/5.0 (Linux; U; Android 2.1-update1; de-de; HTC Desire 1.19.161.5 Build/ERE27) AppleWebKit/530.17 (KHTML, like Gecko) Version/4.0 Mobile Safari/530.17",
        "Mozilla/5.0 (Linux; U; Android 2.2; en-us; Sprint APA9292KT Build/FRF91) AppleWebKit/533.1 (KHTML, like Gecko) Version/4.0 Mobile Safari/533.1",
        "Mozilla/5.0 (Linux; U; Android 1.5; de-ch; HTC Hero Build/CUPCAKE) AppleWebKit/528.5  (KHTML, like Gecko) Version/3.1.2 Mobile Safari/525.20.1",
        "Mozilla/5.0 (Linux; U; Android 2.2; en-us; ADR6300 Build/FRF91) AppleWebKit/533.1 (KHTML, like Gecko) Version/4.0 Mobile Safari/533.1",
        "Mozilla/5.0 (Linux; U; Android 2.1; en-us; HTC Legend Build/cupcake) AppleWebKit/530.17 (KHTML, like Gecko) Version/4.0 Mobile Safari/530.17",
        "Mozilla/5.0 (Linux; U; Android 1.5; de-de; HTC Magic Build/PLAT-RC33) AppleWebKit/528.5  (KHTML, like Gecko) Version/3.1.2 Mobile Safari/525.20.1 FirePHP/0.3",
        "Mozilla/5.0 (Linux; U; Android 1.6; en-us; HTC_TATTOO_A3288 Build/DRC79) AppleWebKit/528.5  (KHTML, like Gecko) Version/3.1.2 Mobile Safari/525.20.1",
        "Mozilla/5.0 (Linux; U; Android 1.0; en-us; dream) AppleWebKit/525.10  (KHTML, like Gecko) Version/3.0.4 Mobile Safari/523.12.2",
        "Mozilla/5.0 (Linux; U; Android 1.5; en-us; T-Mobile G1 Build/CRB43) AppleWebKit/528.5  (KHTML, like Gecko) Version/3.1.2 Mobile Safari 525.20.1",
        "Mozilla/5.0 (Linux; U; Android 1.5; en-gb; T-Mobile_G2_Touch Build/CUPCAKE) AppleWebKit/528.5  (KHTML, like Gecko) Version/3.1.2 Mobile Safari/525.20.1",
        "Mozilla/5.0 (Linux; U; Android 2.0; en-us; Droid Build/ESD20) AppleWebKit/530.17 (KHTML, like Gecko) Version/4.0 Mobile Safari/530.17",
        "Mozilla/5.0 (Linux; U; Android 2.2; en-us; Droid Build/FRG22D) AppleWebKit/533.1 (KHTML, like Gecko) Version/4.0 Mobile Safari/533.1",
        "Mozilla/5.0 (Linux; U; Android 2.0; en-us; Milestone Build/ SHOLS_U2_01.03.1) AppleWebKit/530.17 (KHTML, like Gecko) Version/4.0 Mobile Safari/530.17",
        "Mozilla/5.0 (Linux; U; Android 2.0.1; de-de; Milestone Build/SHOLS_U2_01.14.0) AppleWebKit/530.17 (KHTML, like Gecko) Version/4.0 Mobile Safari/530.17",
        "Mozilla/5.0 (Linux; U; Android 3.0; en-us; Xoom Build/HRI39) AppleWebKit/525.10  (KHTML, like Gecko) Version/3.0.4 Mobile Safari/523.12.2",
        "Mozilla/5.0 (Linux; U; Android 0.5; en-us) AppleWebKit/522  (KHTML, like Gecko) Safari/419.3",
        "Mozilla/5.0 (Linux; U; Android 1.1; en-gb; dream) AppleWebKit/525.10  (KHTML, like Gecko) Version/3.0.4 Mobile Safari/523.12.2",
        "Mozilla/5.0 (Linux; U; Android 2.0; en-us; Droid Build/ESD20) AppleWebKit/530.17 (KHTML, like Gecko) Version/4.0 Mobile Safari/530.17",
        "Mozilla/5.0 (Linux; U; Android 2.1; en-us; Nexus One Build/ERD62) AppleWebKit/530.17 (KHTML, like Gecko) Version/4.0 Mobile Safari/530.17",
        "Mozilla/5.0 (Linux; U; Android 2.2; en-us; Sprint APA9292KT Build/FRF91) AppleWebKit/533.1 (KHTML, like Gecko) Version/4.0 Mobile Safari/533.1",
        "Mozilla/5.0 (Linux; U; Android 2.2; en-us; ADR6300 Build/FRF91) AppleWebKit/533.1 (KHTML, like Gecko) Version/4.0 Mobile Safari/533.1",
        "Mozilla/5.0 (Linux; U; Android 2.2; en-ca; GT-P1000M Build/FROYO) AppleWebKit/533.1 (KHTML, like Gecko) Version/4.0 Mobile Safari/533.1",
        "Mozilla/5.0 (Linux; U; Android 3.0.1; fr-fr; A500 Build/HRI66) AppleWebKit/534.13 (KHTML, like Gecko) Version/4.0 Safari/534.13",
        "Mozilla/5.0 (Linux; U; Android 3.0; en-us; Xoom Build/HRI39) AppleWebKit/525.10  (KHTML, like Gecko) Version/3.0.4 Mobile Safari/523.12.2",
        "Mozilla/5.0 (Linux; U; Android 1.6; es-es; SonyEricssonX10i Build/R1FA016) AppleWebKit/528.5  (KHTML, like Gecko) Version/3.1.2 Mobile Safari/525.20.1",
        "Mozilla/5.0 (Linux; U; Android 1.6; en-us; SonyEricssonX10i Build/R1AA056) AppleWebKit/528.5  (KHTML, like Gecko) Version/3.1.2 Mobile Safari/525.20.1"
    ]
    try:
        while True:
            # # 发请求时，随机从池中选择UA
            header = {
                # ip代理池获取时用的url的headers
                'User-Agent': random.choice(agents),
                "Content-Type": "application/json"
            }

            requests.adapters.DEFAULT_RETRIES = 3  # 增加重连次数
            s = requests.session()  # session
            s.keep_alive = False
            s.headers = header
            r = s.get(url_whishart, timeout=5, verify=False)
            time.sleep(random.uniform(0, 0.5))

            code = r.status_code  # 状态码

            # 用列表存储，最终: [kingdom, superclass, class, subclass]
            item_result = []
            # item_result['table'] = table
            if code == 429:
                # Wait for some time before making another request
                sleep_time = int(r.headers.get('Retry-After', 30))
                # print(f'Received status code 429. Sleeping for {sleep_time} seconds...')
                # print(url_whishart)
                # print(r.text)
                # print(f'cid: {cid}, URL: {url}')
                time.sleep(sleep_time)

                # break loop
            elif code == 200:
                # print(f'{name}: {inchikey}')
                url_res = r.json()
                _kingdom = url_res['kingdom']
                _superclass = url_res['superclass']
                _class = url_res['class']
                _subclass = url_res['subclass']
                parent_level_1 = None
                parent_level_2 = None
                if url_res['direct_parent'] not in [_kingdom, _superclass, _class, _subclass]:  # 存在比subclass更深的分类
                    parent_level_1 = url_res['direct_parent']
                    if len(url_res['intermediate_nodes']) == 1:
                        parent_level_1 = url_res['intermediate_nodes'][0]  # 只有一个值也是列表
                        parent_level_2 = url_res['direct_parent']
                    elif len(url_res['intermediate_nodes']) >= 2:
                        parent_level_1 = url_res['intermediate_nodes'][0]
                        parent_level_2 = url_res['intermediate_nodes'][1]
                # if len(url_res['intermediate_nodes']) > 0:
                #     # 存在内部level
                #     if len(url_res['intermediate_nodes']) == 1:
                #         parent_level_1 = url_res['intermediate_nodes']
                #         parent_level_2 = url_res['direct_parent']
                #     else:  # 有两个以上内部节点
                #         parent_level_1 = url_res['intermediate_nodes'][0]
                #         parent_level_2 = url_res['intermediate_nodes'][1]
                # else:
                #     if url_res['direct_parent'] == _kingdom or url_res['direct_parent'] == _superclass or url_res['direct_parent'] == _class or url_res['direct_parent'] == _subclass:
                #         parent_level_1 = None
                #         parent_level_2 = None
                #     else:
                #         parent_level_1 = url_res['direct_parent']
                #         parent_level_2 = None

                # r.raise_for_status()
                # page_text = r.text
                #
                # html = etree.HTML(page_text)
                # kingdom = html.xpath(
                #     "/html/body[@class='entities-c  show-a']/main/div[@class='main_card'][2]/div[@class='card_minor'][1]/p[@class='text_card']//text()")  # 所有li的a节点
                # superclass = html.xpath(
                #     "/html/body[@class='entities-c  show-a']/main/div[@class='main_card'][2]/div[@class='card_minor'][2]/p[@class='text_card'][2]//text()")
                # _class = html.xpath(
                #     "/html/body[@class='entities-c  show-a']/main/div[@class='main_card'][2]/div[@class='card_minor'][3]/p[@class='text_card'][2]//text()")
                # subclass = html.xpath(
                #     "/html/body[@class='entities-c  show-a']/main/div[@class='main_card'][2]/div[@class='card_minor'][4]/p[@class='text_card'][2]//text()")
                # k = kingdom[0] if len(kingdom) > 0 else None
                # sup = superclass[0] if len(superclass) > 0 else None
                # c = _class[0] if len(_class) > 0 else None
                # sub = subclass[0] if len(subclass) > 0 else None
                item_result.append(replace_non_alphanumeric(_kingdom))
                item_result.append(replace_non_alphanumeric(_superclass))
                item_result.append(replace_non_alphanumeric(_class))
                item_result.append(replace_non_alphanumeric(_subclass))
                item_result.append(replace_non_alphanumeric(parent_level_1))
                item_result.append(replace_non_alphanumeric(parent_level_2))

                item_result.append(inchikey)

                break
            else:
                if r.text == '{"status":404,"error":"Not Found"}':
                    inchikey_no_classify.append(inchikey)
                else:
                    error_q.put(inchikey)
                break

        if len(item_result) > 0:
            result_q.put(item_result)  # put是向结果集 队列里添加元素

    except (requests.exceptions.ProxyError, requests.exceptions.ConnectionError, requests.exceptions.ReadTimeout):
        error_q.put(inchikey)

    except Exception as e:
        print(name, "Error: ", e, url_whishart)


def classify(inchikey_lis, thread_num):
    # URLs = ['https://cfb.fiehnlab.ucdavis.edu/entities/' + i + '.json' for i in inchikey_lis]
    res, not_find_list = q_work(inchikey_lis, thread_num)

    # 获取分类列表
    return res, not_find_list


# if __name__ == '__main__':
#     with open('./test/test_inchikeys.txt', 'r') as f:
#         inchikey_lis = f.read().splitlines()
#
#     print(inchikey_lis)
#
#     res = classify(inchikey_lis, 20)
#     print(res)

# if __name__ == '__main__':
    # 计时ter()

def main(args):
    # 建立数据库连接
    # parser = chem_arg_online.get_parser()
    # args = parser.parse_args()

    # 输入的字符串冒号和空白字符替换为下划线，与数据库格式匹配
    # in_str = input("Enter taxid or names, separate the entries with commas, or enter xx|subtree for subtree : ")
    in_list = []
    total_id = []
    # 画图R脚本路径
    # r_plot_script = os.path.join(os.getcwd(), '../RScript/chem/chemplot.R')

    if args.file:
        # 读取文件中inchikey
        try:
            f = open(args.file, 'r')
            inchikeys = [i.strip('\n') for i in f]
            f.close()
            inchikey_list = inchikeys
        except Exception as e:
            print(e)
            print("file does not exist!")
            exit(0)
    elif args.items:
        item_input = str(args.items)
        # command: -input
        # 只输入inchikey
        inchikey_list = item_input.split(",")

    threads_num = int(args.threads)
    # 二维列表
    res_lis, error_list = classify(inchikey_lis=inchikey_list, thread_num=threads_num)

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
    column_name = ['id', 'inchikey', 'kingdom', 'superclass', 'class', 'subclass', 'parent_level_1', 'parent_level_2']
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
            none_classify_id = id_list_to_newick(res_lis, interrupt_level, path=txt_path, csv_write=csv_write)
        except Exception:
            print('Parameter for interrupting needed! Seek help by -h, --help')
    else:
        none_classify_id = id_list_to_newick(res_lis, path=txt_path, csv_write=csv_write)

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
        # print("------These names are not in taxonomy database------\n", null_name)
        print(f"------These {len(error_list)} inchikeys are not found in classyfire------\n{error_list}")
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
    no_match_output(error_list, none_classify_id, my_prefix, csv_name)
    # color = args.color
    #
    # #  画图，先注释掉吧，最后再说
    # if args.color is not None:
    #     plot.draw_chem_plot(r_plot_script, txt_path, my_prefix, fname=my_fname, csv=csv_name, mycolor=color)
    #     # print(args)
    # else:
    #     plot.draw_chem_plot(r_plot_script, txt_path, my_prefix, fname=my_fname, csv=csv_name)

    # 打印程序运行时间
    # end = time.perf_counter()
    # print('Running time: %s Seconds' % (end - start))
