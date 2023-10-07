import os
from queue import Queue
from Bio import Phylo


class Node:
    def __init__(self, parent=None, value=None):
        """
        self.children = [] 如果使用函数默认值[]，会导致所有node的children 这个list指向内存同一个区域
        :param parent:
        :param value:
        """
        self.parent = parent
        self.children = []
        self.value = value
        self.children_dict = {}

    def add_child(self, node):
        self.children.append(node)
        self.children_dict[node.value] = node

    def __str__(self):
        return self.value if self.value is not None else "None"


class Tree:
    def __init__(self, root):
        self.root = root
        self.num = 1

    def insert(self, chain):
        """
        插入一条数据链
        :param chain:
        :return:
        """
        node_num = len(chain)
        self.build(self.root, chain, 0, node_num)

    def build(self, current_node, chain, idx, node_num):
        if idx == node_num:
            return
        node_value = chain[idx]
        # if node_value is None:
        #     node_value = ""
        next_node = None
        if node_value in current_node.children_dict:
            """
            有这个节点
            """
            next_node = current_node.children_dict[node_value]
        else:
            """
            没有这个结点
            """
            new_node = Node(parent=current_node, value=node_value)
            current_node.add_child(new_node)
            # current_node.children.append(new_node)
            # current_node.children_dict[node_value] = new_node
            next_node = new_node
            self.num += 1
        self.build(next_node, chain, idx + 1, node_num)

    def __str__(self):
        return f"total {self.num} nodes in tree"

    def level_traversal(self, root):
        print("level_traversal..........")
        q = Queue()
        prev_level = 0  # 记录上一层
        q.put((root, 0))
        while not q.empty():
            node, cur_level = q.get()
            if cur_level != prev_level:
                print('')
            print(str(node) + f'({node.parent})', end='\t\t')
            prev_level = cur_level
            for child in node.children:
                q.put((child, cur_level + 1))

    def post_traversal(self, node):
        for child in node.children:
            self.post_traversal(child)
        print(str(node) + f'({node.parent})')

    def print_newick_post(self, root, path):
        stack = []
        newick_str = ""
        prev_level = -1
        stack.append((root, 0, True))
        while len(stack) != 0:
            node, level, tag = stack.pop()
            if tag is True:
                tag = False
                stack.append((node, level, tag))
                # print(f'{node} children: {node.children}')
                for child in node.children[::-1]:
                    """
                    从最右边的孩子开始入栈
                    """
                    stack.append((child, level + 1, True))
            else:
                if level > prev_level:
                    """
                    进入了更深的层，要打左括号
                    """
                    num = level - prev_level
                    newick_str += '(' * num
                    # print('(' * num, end='')
                elif level < prev_level:
                    """
                    level < prev_level, 说明prev_level层所有节点处理完了，回到了上一层，要打右括号
                    此时 level - prev_level 一定为1，因为回到上一层，只能回一层
                    """
                    newick_str += ')'
                    # print(')', end='')
                newick_str += str(node)
                # print(node, end='')
                newick_str += ','
                # print(',', end='')
                prev_level = level

        # 补根节点 最后面的右括号v
        newick_str += ')'
        # print(')', end='')
        newick_str = newick_str.replace(',)', ')')
        # print(newick_str)
        cut_newick = newick_str[1:-1] + ';\n'
        # 打印tree
        # print(cut_newick)
        try:
            with open(path, 'w') as fp:
                fp.write(cut_newick)
        except Exception:
            print("File path does not exist! \nTree is written to /newick.txt in current directory.")
            cur_path = os.getcwd()
            defult_path = os.path.join(cur_path, 'newick.txt')
            with open(defult_path, 'w') as dfp:
                dfp.write(cut_newick)
            # with open('acsii_tree.txt', 'w') as fd:
            #     Phylo.draw_ascii(tree, file=fd)

    def print_newick_level(self, root):
        q = Queue()
        newick_str = ""
        prev_level = 0
        num = 0
        q.put((root, 0))
        while not q.empty():
            node, cur_level = q.get()
            if cur_level != prev_level or cur_level == 0:
                # print('(', end='')
                newick_str += '('
                num += 1
            else:
                # print(',', end='')
                newick_str += ','
            # print(node, end='')
            newick_str += str(node)
            prev_level = cur_level
            for child in node.children:
                q.put((child, cur_level + 1))
        # print(')' * num, end='')
        # print()
        newick_str += ')' * num
        # print(newick_str)


def test():
    a = []
    a.reverse()
    for i in a:
        print(1)


def showtree(data, path):
    data.reverse()

    root = Node(value="")
    # print(root)
    tree = Tree(root)
    for chain in data:
        tree.insert(chain)
    # print(tree)
    # tree.level_traversal(tree.root)
    # tree.post_traversal(tree.root)
    # tree.print_newick_level(root)
    tree.print_newick_post(root, path)
#
# #
# if __name__ == '__main__':
# #     # data = [[None, 'Proteobacteria', 'Betaproteobacteria', 'Neisseriales', 'Neisseriaceae', 'Neisseria', 'Neisseria meningitidis', 'Neisseria meningitidis serogroup B'],
# #     #         [None, 'Proteobacteria', 'Gammaproteobacteria', 'Enterobacterales', 'Enterobacteriaceae', 'Salmonella', 'Salmonella enterica', 'Salmonella enterica subsp. enterica serovar Abortusequi'],
# #     #         [None, 'Firmicutes', None, None, None, None, None],
# #     #         [None, 'Firmicutes', 'Bacilli', 'Lactobacillales', 'Lactobacillaceae', 'Leuconostoc', 'Leuconostoc fallax']]
# #
#     data = [
#         ['F', 'A'],
#         ['F', 'A', 'G', None],
#         ['F', 'A', 'H', None, 'M'],
#         ['F', 'B'],
#         ['F', 'E', 'C'],
#         ['F', 'E', 'D'],
#     ]
# #
# #     """
# #     是否翻转，取决于输入数据和要求
# #     """
# #     data.reverse()
# #
#     root = Node(value="")
#     print(root)
#     tree = Tree(root)
#     for chain in data:
#         tree.insert(chain)
#     print(tree)
# #     # tree.level_traversal(tree.root)
# #     # tree.post_traversal(tree.root)
# #     # tree.print_newick_level(root)
#     tree.print_newick_post(root, "./test3.txt")
