from ete3 import TreeStyle


def set_ts():
    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.show_branch_length = True
    ts.show_branch_support = False
    ts.mode = "c"
    ts.arc_start = -180  # 0 degrees = 3 o'clock
    ts.arc_span = 360
    # ts.optimal_scale_level = "full"
    # ts.scale = 80
    return ts

def set_pdf():
    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.show_branch_length = True
    ts.show_branch_support = False
    return ts
