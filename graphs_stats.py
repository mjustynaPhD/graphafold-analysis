import os

import numpy as np
import pandas as pd

from _metrics import filter_to_noncanon, read_ct, read_gt, to_bin

GT_PATH = "/home/mjustyna/data/graphafold_data/casp/"


def main():
    # Get the ground truth data
    gt_files = os.listdir(GT_PATH)
    all_gt_files = sorted([f for f in gt_files if f.endswith('.amt')])
    can_gt_files = sorted([f for f in gt_files if f.endswith('.cmt')])

    gts = {}
    print("PDB ID, Unpaired, paired, nodes")
    for all_g, can_g in zip(all_gt_files, can_gt_files):
        # print(all_g, can_g)
        all_g_path = os.path.join(GT_PATH, all_g)
        can_g_path = os.path.join(GT_PATH, can_g)

        all_g_vec, leng_all = read_gt(all_g_path)
        can_g_vec, leng_can = read_gt(can_g_path)
        assert len(all_g_vec) >= len(can_g_vec)
        assert leng_all == leng_can, f'{all_g} {can_g} {leng_all} {leng_can}'
        # get difference between all_g_vec and can_g_vec
        # print(f"{all_g.replace('.amt', '')} {len(can_g_vec)//2} {(len(all_g_vec)-len(can_g_vec))//2}")
        
        all_paired = np.zeros(leng_all)
        all_paired[all_g_vec] = 1
        all_paired[can_g_vec] = 1
        print(f"{all_g.split('.')[0]} {int(leng_all - np.sum(all_paired))} {int(np.sum(all_paired))} {leng_all}")
if __name__ == "__main__":
    main()