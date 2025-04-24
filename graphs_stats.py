import os

import numpy as np
import pandas as pd

from _metrics import filter_to_noncanon, read_ct, read_gt, to_bin

# GT_PATH = "/home/mjustyna/data/graphafold_data/casp/"
GT_PATH = "large"

def get_uniq_chains(idx_path):
    with open(idx_path, 'r') as f:
        lines = f.readlines()
    lines = [line.strip().split(',')[1].split('.')[0] for line in lines]
    return list(set(lines))

def main():
    # Get the ground truth data
    gt_files = os.listdir(GT_PATH)
    all_gt_files = sorted([f for f in gt_files if f.endswith('.amt')])
    can_gt_files = sorted([f for f in gt_files if f.endswith('.cmt')])

    gts = {}
    print("PDB ID, Length, Unpaired, paired, nodes, canonical, non-canonical, backbone")
    for all_g, can_g in zip(all_gt_files, can_gt_files):
        # print(all_g, can_g)
        all_g_path = os.path.join(GT_PATH, all_g)
        can_g_path = os.path.join(GT_PATH, can_g)

        all_g_vec, leng_all = read_gt(all_g_path)
        can_g_vec, leng_can = read_gt(can_g_path)
        uniq_chains = get_uniq_chains(all_g_path.replace('.amt', '.idx'))
        assert len(all_g_vec) >= len(can_g_vec)
        assert leng_all == leng_can, f'{all_g} {can_g} {leng_all} {leng_can}'
        # get difference between all_g_vec and can_g_vec
        # print(f"{all_g.replace('.amt', '')} {len(can_g_vec)//2} {(len(all_g_vec)-len(can_g_vec))//2}")
        name = all_g.split('.')[0]
        name = name.replace('PDB_0000', '')[:4]
        all_paired = np.zeros(leng_all)
        all_paired[all_g_vec] = 1
        all_paired[can_g_vec] = 1
        print(f"{name} {leng_all} {int(leng_all - np.sum(all_paired))} {int(np.sum(all_paired))} {leng_all} {len(can_g_vec)//2} {(len(all_g_vec)-len(can_g_vec))//2} {leng_all - len(uniq_chains)}")

if __name__ == "__main__":
    main()