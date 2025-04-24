import os
import numpy as np
import pandas as pd
from sklearn.metrics import precision_score, recall_score, f1_score

GT_PATH = "/home/mjustyna/data/graphafold_data/casp/"
preds_methods = ['sincfold', 'ufold', 'spotrna']

# Funkcja do obliczania metryk
def calculate_metrics(predicted, labels):
    
    if len(np.unique(labels)) == 1:
        precision = recall = f1 = accuracy = 1.0 if np.all(labels == predicted) else 0.0
    else:
        precision = precision_score(labels, predicted, zero_division=0)
        recall = recall_score(labels, predicted, zero_division=0)
        f1 = f1_score(labels, predicted, zero_division=0)
        accuracy = np.mean(predicted == labels)
        tp = np.sum((predicted == 1) & (labels == 1))
        fn = np.sum((predicted == 0) & (labels == 1))
        fp = np.sum((predicted == 1) & (labels == 0))
        ppv = tp / (tp + fp) if (tp + fp) > 0 else 0
        tpr = tp / (tp + fn) if (tp + fn) > 0 else 0
        inf_metric = np.sqrt(ppv * tpr) if (ppv * tpr) > 0 else 0
    return round(accuracy, 2), round(precision, 2), round(recall,2), round(f1,2), round(inf_metric, 2), tp, fp, fn

def read_ct(file):
    with open(file, 'r') as f:
        lines = f.readlines()
    lines = [l.split() for l in lines[1:] if len(l.split()) == 6]
    # pairs = [[int(i) for i in l1] for l1 in pairs]
    lines = np.array(lines)
    seq = lines[:, 1]
    pairs = lines[:, 4]
    # print('Pairs:', pairs)
    pairs = np.array([int(i) for i in pairs])
    return pairs, seq

def filter_to_noncanon(pairs, seq, return_bps=False):
    # Filter pairs to only include non-canonical pairs (A-U, G-C and G-U)
    filtered_pairs = []
    bps = []
    assert len(pairs) == len(seq)
    assert max(pairs) <= len(seq), f'Max pair:{max(pairs)}\nsequence length {len(seq)}'
    
    for i, p in enumerate(pairs):
        if p == 0:
            filtered_pairs.append(0)
            continue
        
        # print(f'Pair: {i+1} {p} {pairs[p-1]}')
        nt1 = seq[pairs[p-1]-1]
        nt2 = seq[p-1]
        bp=f'{nt1}-{nt2}'
        # print(f'Base pair: {bp}')
        if bp not in ['A-U', 'U-A', 'G-C', 'C-G']: # We treat wobble as non-canonical: 'G-U', 'U-G'
            # print(f'Pair: {i+1} {p} {pairs[p-1]}\t{nt1}-{nt2}')
            filtered_pairs.append(p)
            bps.append([p, pairs[p-1]])
        else:
            filtered_pairs.append(0)
    assert len(filtered_pairs) == len(pairs)
    
    if return_bps:
        return np.array(filtered_pairs), np.array(bps)
    else:
        return np.array(filtered_pairs)

def to_bin(arr):
    # Convert to binary
    arr[arr > 0] = 1
    return arr

def read_gt(file):
    with open(file, 'r') as f:
        lines = f.readlines()

    lines = [[int(i) for i in l1.strip().split(',')] for l1 in lines]
    arr = np.array(lines)

    return np.where(arr > 0)[0], len(arr[0])

def get_vector_diff(v1, v2):
    # v1 = [0, 1, 2, 3, 4, 4, 5]
    # v2 = [0, 1, 2, 4, 5]
    # output: [3, 4]
    output = []
    vj_index = 0
    for vi in v1:
        vj_index = min(len(v2)-1, vj_index)
        if vi == v2[vj_index]:
            vj_index+=1
        else:
            output.append(vi)
    return np.array(output)


def get_gt():
    gt_files = os.listdir(GT_PATH)
    all_gt_files = sorted([f for f in gt_files if f.endswith('.amt')])
    can_gt_files = sorted([f for f in gt_files if f.endswith('.cmt')])

    gts = {}
    for all_g, can_g in zip(all_gt_files, can_gt_files):
        print(all_g, can_g)
        all_g_path = os.path.join(GT_PATH, all_g)
        can_g_path = os.path.join(GT_PATH, can_g)

        all_g_vec, leng_all = read_gt(all_g_path)
        can_g_vec, leng_can = read_gt(can_g_path)
        print("all g", all_g_vec)
        print(can_g_vec)
        print("Lengs: ", len(all_g_vec), len(can_g_vec), len(all_g_vec)-len(can_g_vec))
        assert len(all_g_vec) >= len(can_g_vec)
        assert leng_all == leng_can, f'{all_g} {can_g} {leng_all} {leng_can}'
        # get difference between all_g_vec and can_g_vec
        # diff = np.setdiff1d(all_g_vec, can_g_vec)
        diff = get_vector_diff(all_g_vec, can_g_vec)
        print(diff)
        print(len(np.unique(diff)))
        print("Leng diff:", len(diff), len(all_g_vec)-len(can_g_vec))
        
        
        nc = np.zeros(leng_all)
        nc[diff] = 1
        gts[all_g.split('.')[0]] = nc
    
    return gts

def main():
    
    gts = get_gt()

    print('GTs:', gts.keys())
    # create a dataframe for each method
    for method in preds_methods:
        print(method)
        results = sorted(os.listdir(f'{method}_results'))
        results = sorted([f for f in results if f.endswith('.ct')])
        method_df = pd.DataFrame(columns=['PDB', 'Precision', 'Recall', 'F1', 'INF', 'TP', 'FP', 'FN'])
        for res in results:
            # print(res)
            res_path = os.path.join(f'{method}_results', res)
            pairs, seq = read_ct(res_path)
            
            pairs = filter_to_noncanon(pairs, seq)
            bin_pairs = to_bin(pairs)
            accuracy, precision, recall, f1, inf, tp, fp, fn = calculate_metrics(bin_pairs, gts[res.split('.')[0]])
            # print(f'Accuracy: {accuracy}, Precision: {precision}, Recall: {recall}, F1: {f1}')
            row = pd.DataFrame({'PDB': res.split('.')[0], 'Precision': precision, 'Recall': recall, 'F1': f1, 'INF':inf, 'TP':tp, 'FP':fp, 'FN': fn}, index=[0])
            method_df = pd.concat([method_df, row])
        method_df.to_csv(f'{method}_results.csv', index=False)
    

if __name__ == "__main__":
    main()
