from ete3 import Tree, TreeStyle


def set_pdf():
    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.show_branch_length = False
    ts.show_branch_support = False
    return ts


def render_pdf(path, outpath):
    # 生成可视化pdf
    t = Tree(newick=path, format=1)
    pdf_ts = set_pdf()
    t.render(outpath, tree_style=pdf_ts)
    print((f"{len(t)} leaves in tree"))




def get_leaf_num(path):
    t = Tree(newick=path, format=1)
    return len(t)
