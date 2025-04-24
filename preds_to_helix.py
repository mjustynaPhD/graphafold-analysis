import os

import numpy as np
import pandas as pd

from metrics_as_matrix import filter_to_noncanon, read_ct, get_gt

preds_methods = ['sincfold', 'ufold', 'spotrna']
GT_PATH = "/home/mjustyna/data/graphafold_data/casp/"

type_categories = {
    1: 'PredictedGoodNonCanonical',
    2: 'PredictedBadNonCanonical',
    3: 'NotPredictedNonCanonical'
}

def read_gt(file):
    with open(file, 'r') as f:
        lines = f.readlines()

    lines = [[int(i) for i in l1.strip().split(',')] for l1 in lines]
    arr = np.array(lines)
    x,y = np.where(arr > 1)
    nc_pairs = np.array([[int(i), int(j)] for i, j in zip(x, y)])
    return nc_pairs

def get_gt_nc_pairs():
    gt_files = os.listdir(GT_PATH)
    all_gt_files = sorted([f for f in gt_files if f.endswith('.amt')])
    gts = {}
    for all_g in all_gt_files:
        all_g_path = os.path.join(GT_PATH, all_g)

        nc_pairs = read_gt(all_g_path)
        gts[all_g.split('.')[0]] = nc_pairs
    
    return gts

def matrix_to_pairs(matrix):
    pairs = []
    for i, row in enumerate(matrix):
        for j, col in enumerate(row):
            if col > 0:
                pairs.append([i, j])
    return np.array(pairs)


def to_helix(gt_pairs, pred_pairs):
    helix = []
    list_of_pred_x = pred_pairs[:, 0]
    pairs_added = []
    for x, y in gt_pairs:
        p = f'{y},{x}'
        if p in pairs_added:
            continue
        pairs_added.append(f'{x},{y}')

        if x in list_of_pred_x or y in list_of_pred_x:
            # get the index of pred_x in gt_pairs
            index = np.where(list_of_pred_x == x)
            if len(index[0]) == 0:
                index = np.where(list_of_pred_x == y)
            index = index[0][0]  # get the first index
            # get the corresponding y value from gt_pairs
            
            pred_pair = pred_pairs[index]
            if all(np.array([x, y]) == pred_pair):
                helix.append([x, y, 1, 1.0])
            else:
                helix.append([x, y, 2, 1.0])
        else:
            helix.append([x, y, 3, 1.0])
    helix = np.array(helix)
    out_df = pd.DataFrame(helix, columns=['x', 'y', 'type', 'leng'])
    out_df['type'] = out_df['type'].map(type_categories)
    # convert start and end to int
    out_df['x'] = out_df['x'].astype(int)
    out_df['y'] = out_df['y'].astype(int)
    return out_df

def main():
    # Get the ground truth data
    gts = get_gt(flatten=False)

    # print('GTs:', gts.keys())
    # # create a dataframe for each method
    for method in preds_methods:
        print(method)
        results = sorted(os.listdir(f'{method}_results'))
        results = sorted([f for f in results if f.endswith('.ct')])
        
        for res in results:
            # print(res)
            res_path = os.path.join(f'{method}_results', res)
            mat, seq = read_ct(res_path)
            
            mat = filter_to_noncanon(mat, seq)
            gt_pairs = matrix_to_pairs(gts[res.split('.')[0]])
            pred_pairs = matrix_to_pairs(mat)


            # assert len(bps) %2 == 0, f'len(bps) {len(bps)} is not even'
            # print(bps)
            if len(pred_pairs) > 0:
                helix_df = to_helix(gt_pairs, pred_pairs)
                # increas 'x' and 'y' by 1
                helix_df['x'] = helix_df['x'] + 1
                helix_df['y'] = helix_df['y'] + 1
                # save to csv. No index, no header
                helix_df.to_csv(f'{method}_results/{res.split(".")[0]}.helix', index=False, header=False)


    pass

if __name__ == "__main__":
    main()