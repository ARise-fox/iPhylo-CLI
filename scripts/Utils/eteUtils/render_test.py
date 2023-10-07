from ete3 import Tree
import tree_style

def render_pdf(path, outpath):
    # 生成可视化pdf
    t = Tree(newick=path, format=1)
    pdf_ts = tree_style.set_pdf()
    t.render(outpath, tree_style=pdf_ts)
# t = Tree(newick='D:/桌面2/test/1k_Tree.txt', format=1)


render_pdf(path='../files/sub_22.txt', outpath='D:/桌面2/test/sub_22.pdf')
