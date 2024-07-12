# from ete3 import Tree, TreeStyle
import os
import queue
import sqlite3
import re
import sys
import threading
import time
import random
import traceback

import requests
from Bio import Phylo
from pathlib import Path
import csv

# import eteUtils.render_test
# from eteUtils import tree_style
from scripts.Utils.myArg import chem_arg_online
import pandas as pd
from alive_progress import alive_bar

from scripts.Utils import newick
from scripts.Utils.eteUtils import myete


class myThread(threading.Thread):
    def __init__(self, name, q, error_q, result_q, smiles_not_found, bar, agents, csv_writer):
        threading.Thread.__init__(self)
        self.name = name
        self.q = q
        self.result_q = result_q
        self.error_q = error_q
        self.smiles_not_found = smiles_not_found
        self.bar = bar
        self.agents = agents
        self.csv_writer = csv_writer
        self.session = requests.Session()

    def run(self):
        # print("Starting " + self.name)
        while True:
            if self.q.empty():
                break
            else:
                try:
                    crawl(self.name, self.q, self.error_q, self.result_q, self.smiles_not_found, self.agents,
                          self.session,
                          self.csv_writer)
                    self.bar()
                except Exception as e:
                    print(f"Error in {self.name}: {str(e)}")
                    # traceback.print_exc()  # 打印异常的详细信息和堆栈跟踪
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


def id_list_to_newick(classify_lis, *interrupt_level, path):
    lin_list = []
    none_classify_id = []
    for idx, item in enumerate(classify_lis):
        # classy_info:分类信息
        classy_info = item[0:-1]
        _id = item[-1]

        if interrupt_level:
            # 如果需要在某level打断，找到该level对应的index
            level_list = ["pathway", "superclass"]
            level_idx = -1
            for _idx, level in enumerate(level_list):
                if interrupt_level[0] == level:
                    # 如果打断为pathway， level_idx = 1， full_list[:1] = [Fatty_acids]
                    level_idx = _idx + 1
            classy_info = classy_info[:level_idx]

        full_list = list(filter(None, classy_info))

        if len(classy_info) > 0:
            # 过滤掉分类全是None的情况
            full_list.append(str(_id))
            lin_list.append(full_list)
        else:
            # these inchikey found no classification
            none_classify_id.append(item[-1])

    # 二维列表去重，需要先转化为tuple
    if interrupt_level:
        lin_list = list(set([tuple(t) for t in lin_list]))
        lin_list = [list(v) for v in lin_list]

    newick.showtree(data=lin_list, path=path)  # id分支不合并
    return none_classify_id
    # 生成树文件


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


# 判断一个字符串是否是inchikey
def is_valid_inchikey(inchikey):
    inchikey_pattern = re.compile(r'^[A-Z]{14}-[A-Z]{10}-[A-Z]$')
    return re.match(inchikey_pattern, inchikey) is not None


# 判断是否为有效inchi
def is_valid_inchi(inchi):
    inchi_pattern = re.compile(r'^InChI=1S?/[A-Za-z0-9]+')
    return re.match(inchi_pattern, inchi) is not None


def no_match_output(none_inchikey, none_classify_id, prefix, csv):
    fname = os.path.join(prefix, 'no_match_result.txt')
    with open(fname, 'w') as f:
        if len(none_inchikey) > 0:
            unique_inchikeys = list(set(none_inchikey))
            f.write(f"------These SMILES could not be found in the database------\n")
            for inchikey in unique_inchikeys:
                f.write(inchikey + '\n')
        if len(none_classify_id) > 0:
            try:
                # 读取csv文件
                data = pd.read_csv(csv)
                # 筛选出对应的行
                filtered_data = data[data['id'].isin(none_classify_id)]
                # 提取inchikey列并转换为列表
                none_classified_inchikey_list = filtered_data['inchikey'].tolist()
                f.write('\n\n\n------No classification information found for these SMILES------\n')
                for nc_inchikey in none_classified_inchikey_list:
                    f.write(nc_inchikey + '\n')
            except FileNotFoundError:
                f.write('\n\n\n------No classification information found for these IDs------\n')
                for nc_inchikey in none_classify_id:
                    f.write(nc_inchikey + '\n')


def q_work(workQueue, thread_num, agents, csv_writer):
    res_list = []
    smiles_not_found = []  # NP_Classifier中不存在分类的
    equal_times = 2
    cnt = 0

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
                # print(str(i), list(workQueue.queue))
                thread = myThread("Thread-" + str(i), q=workQueue, error_q=error_q, result_q=resq,
                                  smiles_not_found=smiles_not_found, bar=bar, agents=agents, csv_writer=csv_writer)

                # 开启新线程
                thread.start()
                # 添加新线程到线程列表
                threads.append(thread)

            # 等待所有线程完成
            for thread in threads:
                thread.join()

            print("Exiting Main Thread")

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

    return res_list, smiles_not_found


# 替换空格为下划线，空白为none
def replace_non_alphanumeric(s, smiles, idx, errq):
    if s is None:
        return None
    if isinstance(s, list):
        if len(s) > 0:
            s = s[0]
            if s:  # 确保s不为空
                processed = re.sub(r'\W+', '_', s)
                return processed
        else:
            return None
    else:
        errq.put((idx, smiles))


def crawl(name, q, error_q, result_q, smiles_not_found, agents, session, csv_writer):
    # print(list(q.queue))
    idx, smiles = q.get()
    url = f"https://npclassifier.gnps2.org/classify?smiles={smiles}"
    try:
        while True:
            headers = {'User-Agent': random.choice(agents)}
            response = session.get(url, headers=headers, timeout=3)
            time.sleep(random.uniform(0, 0.5))
            code = response.status_code  # 状态码

            # 用列表存储，最终: [pathway, superclass, class]
            item_result = []

            # break loop
            if code == 200:
                data = response.json()

                _pathway = replace_non_alphanumeric(data['pathway_results'], idx, smiles, error_q)
                _superclass = replace_non_alphanumeric(data['superclass_results'], idx, smiles, error_q)
                _class = replace_non_alphanumeric(data['class_results'], smiles, idx, error_q)

                # 树的分支：_pathway, _superclass, _class, idx
                item_result.append(_pathway)
                item_result.append(_superclass)
                item_result.append(_class)
                item_result.append(idx + 1)

                csv_writer.writerow([idx + 1, smiles, _pathway, _superclass, _class])
                #  ['Fatty_acids', 'Sphingolipids', 'Neutral_glycosphingolipids', 17],

                break
            else:
                smiles_not_found.append(smiles)
                error_q.put((idx, smiles))
                # print(list(error_q.queue))
                break

        if len(item_result) > 0:
            result_q.put(item_result)  # put是向结果集 队列里添加元素

    except Exception as e:
        print(name, "Error: ", e, url)


def np_classify(workQueue, thread_num, agents, csv_writer):
    res, not_find_list = q_work(workQueue, thread_num, agents, csv_writer)
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

# 将化合物名称及其编号加入队列
def enqueue_compounds(compounds, queue):
    for idx, compound in enumerate(compounds):
        # 创建元组，包括编号和化合物名称
        queue.put((idx, compound))


def main(args):
    user_agents = [
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

    # 线程数
    threads_num = int(args.threads)
    # 输出文件路径前缀
    my_prefix = Path(args.prefix)
    # 输出文件名
    my_fname = args.fname
    # 输入文件
    finput = args.file
    # 输入SMILES和编号队列
    compound_queue = queue.Queue()

    # 读取文件中的smiles
    try:
        f = open(finput, 'r')
        smiles_list = [i.strip('\n') for i in f]
        enqueue_compounds(smiles_list, compound_queue)
        f.close()
    except Exception as e:
        print(e)
        print("file does not exist!")
        exit(0)

    # 解析输出路径

    if my_prefix.is_dir():
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

    # outpath_fname: 用户通过arg指定的输出路径和文件名称拼接，不包含后缀
    outpath_fname = os.path.join(my_prefix, my_fname)
    txt_path = outpath_fname + ".txt"

    """
    根据输入新建SCV文件
    """
    csv_name = outpath_fname + '_items.csv'

    scv_file = open(csv_name, 'w', encoding='utf-8', newline='')
    writer = csv.writer(scv_file)
    column_name = ['id', 'SMILES', 'pathway', 'superclass', 'class']
    #  3.构建列表头
    writer.writerow(column_name)
    """
    开始分类
    """
    # res_list: 结果列表 二维列表
    # error_list: 错误列表 二维列表
    res_list, error_list = np_classify(workQueue=compound_queue, thread_num=threads_num, agents=user_agents,
                                       csv_writer=writer)

    if args.interrupt:
        interrupt_flag = 1
        if args.pathway:
            interrupt_level = 'pathway'
        elif args.superclass:
            interrupt_level = 'superclass'

    else:
        interrupt_flag = 0

    """
    是否打断，画树
    """
    if interrupt_flag == 1:
        try:
            none_classify_id = id_list_to_newick(res_list, interrupt_level, path=txt_path)
        except Exception as e:
            print('Parameter for interrupting needed! Seek help by -h, --help \n', str(e))
    else:
        none_classify_id = id_list_to_newick(res_list, path=txt_path)

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
    try:
        if len(error_list) > 0:
            uniq_error_list = list(set(error_list))
            print(f"------These {len(uniq_error_list)} SMILES are not found.------\n{uniq_error_list}")
    except:
        pass
    if len(none_classify_id) != 0:
        # print("------These names are not in taxonomy database------\n", null_name)
        print(f"------These {len(none_classify_id)} SMILES are not classified currently------\n{none_classify_id}")

    no_match_output(error_list, none_classify_id, my_prefix, csv_name)
    scv_file.close()
