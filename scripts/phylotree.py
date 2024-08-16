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
import sys
from scripts.Utils import newick
from scripts.Utils.eteUtils import myete

# Global variables for storing IDs and names that are not in the database
null_id = []
null_name = []
multi_names = []
interrupt_flag = 0

# SQL queries for fetching data from the database
id_lineage_tabel = "id2lineage"
sql = f"""select _superkingdom, _kingdom, _phylum, _class, _order, _family, _genus, _species from {id_lineage_tabel} where id = :id """
sql_find_sciname = """select name from alter_name_nopk_flag where taxid = :taxid and type_flag = 1 """
sql_find_id_and_flag_by_name = """select taxid, type_flag from alter_name_nopk_flag where name = :name """
sql_check_multi_name = """select count(*) from multi_taxid WHERE name = :name"""
sql_find_id_by_name = """select taxid from alter_name_nopk_flag where name = :name """
sql_name_by_id = f"""select sciname from {id_lineage_tabel} where id = :id """
sql_rank_by_name = f"""select * from {id_lineage_tabel} where sciname = :name """
sql_rank_by_id = f"""select * from {id_lineage_tabel} where id = :name """


def write_anno(data, filename):
    """
       Writes the annotation data to a CSV file.

       Parameters:
       data (list): The data to be written to the file.
       filename (str): The name of the file to write to.
       """
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)

        # write row name
        writer.writerow(['id', 'taxid', 'domain', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'])

        # write in data
        for row in data:
            last_item = row[-1]
            row_data = [last_item]  # the last item is the id, write it to the first column
            taxid = row[0]  # get the taxid, which is the first element in the row
            row_data.append(taxid)  # write the taxid to the second column

            for item in row[1:]:  # iterate over the rest of the items in the row, except the first
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
    Connects to the database and returns the connection and cursor.

    Returns:
    conn (sqlite3.Connection): The connection to the database.
    cursor (sqlite3.Cursor): The cursor for the database.
    """
    # database directory
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


# make a newick file, replace the branch length to 1
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
                if int(multi[0][0]) > 0:  # this name is ambiguous
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
            if int(multi[0][0]) > 0:  # this name is ambiguous
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


"""
Requirement: If the input item is a standard lineage, use the standard 7 levels. 
If not, append the item's name at the end of the standard lineage.
"""


def id_list_to_newick(id_list, cursor, *interrupt_level, anno_name, path, csv_write):
    """
    make a newick file base on the id_list
    :param id_list: a list of taxid
    :return: name_res_mapping, a dict, key is the input name, value is the leaf node in the tree
    """
    # lin_list: a two-dimensional array. To access the second rank of the first species, you should use lin_list[0][1].
    lin_list = []
    anno_2d_list = []
    # 记录scientific name 和 leaf node 的 mapping
    name_res_mapping = {}  # record the mapping between name and res[-1]

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
        # name: input item: sciname
        name = cursor.execute(sql_find_sciname, {"taxid": int_item}).fetchone()[0]

        # write in csv
        csv_list = list(full_list)
        csv_list.insert(0, name)
        csv_list.insert(0, item)
        csv_write.writerow(csv_list)

        # get the last not None element in the list， index: the last not None element's index

        index = next(len(full_list) - i for i, j in enumerate(reversed(full_list), 1) if j is not None)
        if name != str(full_list[index])[3:]:  # 3：filter the first 3 character, eg: K__
            full_list.insert(index + 1, name)

        if interrupt_flag == 1:
            full_list = full_list[:level_idx]

        res = list(filter(None, full_list))
        lin_list.append(res)
        # 用于注释 anno
        anno_list = res[:]
        anno_list.insert(0, item)
        anno_2d_list.append(anno_list)

        if res:
            name_res_mapping[name] = res[-1]

    # remove duplicate
    if interrupt_flag == 1:
        lin_list = list(set([tuple(t) for t in lin_list]))
        lin_list = [list(v) for v in lin_list]

    write_anno(anno_2d_list, anno_name)

    # make newick file
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
    f_map = prefix_fname + '_replaced_name.txt'  # output filepath+filename
    with open(f_map, 'w') as f:
        for name, sciname in map.items():
            f.write(name + '\t')
            f.write('->')
            f.write('\t' + sciname)
            f.write('\n')


def convert_file(prefix_fname, args):
    """
    :param prefix_fname: tree file path + name
    :param args: args from command line
    """
    # first generate txt file
    txt_path = prefix_fname + ".txt"
    # handle the defined suffix

    f_nwk = prefix_fname + '.nwk'

    Phylo.convert(txt_path, "newick", f_nwk, "newick")

    if args.branch_length:
        replace_branch_length(f_nwk)
    else:
        cut_branch_length(f_nwk)
        # os.remove(txt_path)
    # Output txt format copy nwk rename txt
    f_txt = prefix_fname + '.txt'
    f1 = open(f_nwk)
    f2 = open(f_txt, 'w')
    f2.write(f1.read())
    f1.close()
    f2.close()
    # output nexus format
    f_nex = prefix_fname + '.nex'
    Phylo.convert(f_nwk, "newick", f_nex, "nexus")
    if args.branch_length:
        replace_branch_length(f_nex)
    else:
        cut_branch_length(f_nex)
        # os.remove(txt_path)
    # ouput phyloxml format
    f_xml = prefix_fname + '.xml'
    # If need a branch length, convert it to nex format and add the branch length first
    Phylo.convert(f_nex, "nexus", f_xml, "phyloxml")


def filter_subtree(in_list):
    # sub_list: subtree input need to be processed，eg: xxx|subtree，yyy|subtree，sublist is stored as[xxx,yyy]
    sub_list = []
    # no_sub_list: taxa which don't need to be processed as subtree
    no_sub_list = []
    for i in in_list:
        match_obj = re.match(".+\|subtree$", i)
        if match_obj:
            sub_list.append(match_obj.group())
        else:
            no_sub_list.append(i)
    return sub_list, no_sub_list


def find_sub_id(item_list, cursor):
    """
    Convert the id or name of subtree| to the corresponding id and add it to the total_id list
    """

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
                # not find: non-standard classification such as subspecies, subphyla etc.
                print(f"{sciname} not find in db!")
                continue
            id_list = get_id_list_by_rank(rank_name, rank_query, cursor)
            sub_id.extend(id_list)
            bar()
    return sub_id


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
    # connect to database
    conn, cursor = connect_db()
    item_input = str(args.items)
    sub_name_or_id = args.subtree  # get subtree from args

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
        sciname = sub_name_or_id
        isdigit = False
        if sciname.isdigit():
            isdigit = True
        rank_name, rank_query = get_rank(sciname, cursor, isdigit)
        if rank_name == "" or rank_query == "":
            # not find:  other non-standard classification
            cursor.close()
            conn.commit()
            conn.close()
            print(f"{sciname} is not standard taxonomic unit for input")
            sys.exit(0)
        id_list = get_id_list_by_rank(rank_name, rank_query, cursor)
        all_id = id_list
    elif args.file:
        try:
            with open(args.file, 'r') as f:
                in_list = []
                for line in f.readlines():
                    line = str(line.strip())
                    in_list = process_line(line, in_list)
                sub_list, no_sub_list = filter_subtree(in_list)
                # generate no_sub_id and name_sciname_map
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

        no_sub_id, name_sciname_map = get_id_and_sciname_list_by_id_or_name(no_sub_list, cursor)
        sub_id = find_sub_id(sub_list, cursor=cursor)

        all_id = list(set(no_sub_id + sub_id))

    """
    get output path
    """
    my_prefix = Path(args.prefix)
    my_fname = args.fname
    if my_prefix.is_dir():
        pass
    else:
        my_prefix = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), 'iphylo_files'))
        # not_exist_msg = f"The path '{args.prefix}' does not exist!"
        print(f"The path '{args.prefix}' does not exist!")

    # Create a level-1 folder in the specified path with the same name as fname to store the generated data
    my_prefix = my_prefix / Path(my_fname)
    isExists = os.path.exists(my_prefix)

    if not isExists:
        os.makedirs(my_prefix)

    # my_prefix_name: The output path specified by the user concatenated with the file name, without suffixes
    # my_prefix_name renamed to outpath_fname
    outpath_fname = os.path.join(my_prefix, my_fname)
    txt_path = outpath_fname + ".txt"

    csv_name = outpath_fname + '_items.csv'
    f_csv = open(csv_name, 'w', encoding='utf-8', newline='')
    csv_write = csv.writer(f_csv)
    csv_write.writerow(
        ['taxid', 'input_term', 'domain', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'])
    anno_name = outpath_fname + '_items_for_anno.csv'

    if interrupt_flag == 1:
        # id_list_to_newick(all_id, cursor, interrupt_level, anno_name=anno_name, path=txt_path, csv_write=csv_write)
        try:
            sciname_leaf_map = id_list_to_newick(all_id, cursor, interrupt_level, anno_name=anno_name, path=txt_path,
                                                 csv_write=csv_write)
        except Exception:
            print('Parameter for interrupting needed! Seek help by -h, --help')
    else:
        sciname_leaf_map = id_list_to_newick(all_id, cursor, anno_name, anno_name=anno_name, path=txt_path,
                                             csv_write=csv_write)

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

    commit_close_db(conn, cursor)
    f_csv.close()
    #
    # # plot r script path
    # r_plot_script = os.path.join(os.getcwd(), 'RScript/phylo/drawTree.R')
    # # plot
    #
    # color = args.color
    # if args.color is not None:
    #     plot.draw_chem_plot(r_plot_script, txt_path, my_prefix, fname=my_fname, csv=csv_name, mycolor=color)
    #     # print(args)
    # else:
    #     plot.draw_chem_plot(r_plot_script, txt_path, my_prefix, fname=my_fname, csv=csv_name)

    # print run time
    end = time.perf_counter()
    print('Running time: %s Seconds' % (end - start))
