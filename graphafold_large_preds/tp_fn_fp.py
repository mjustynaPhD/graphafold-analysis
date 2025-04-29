from pathlib import Path
from collections import Counter


# results_dir = Path("graphafold_large_preds/helix_outputs_large")
results_dir = Path("helix_outputs_large")

def get_tp_fn_fp(helix_file):
    with open(helix_file, "r") as f:
        lines = f.readlines()
    res = [l.split(',')[2] for l in lines]
    indeces = [[int(l.split(',')[0]), int(l.split(',')[1])] for l in lines]
    res_uniq = []
    
    
    name = helix_file.stem
    name = name.replace('.amt', "").replace('PDB_0000', "")
    # print(name)
    added = {}
    for ind, r in zip(indeces, res):
        i, j = ind
        uid = f'{i}-{j}'
        if uid not in added:
            added[f'{j}-{i}'] = 1
            added[uid] = 1
            res_uniq.append(r)
    # use counter to count the number of each type
    counts = Counter(res_uniq)
    print(name[:4], counts['PredictedGoodNonCanonical'], counts['NotPredictedNonCanonical'], counts['PredictedBadNonCanonical'])

def main():
    # iterate over helix files
    for helix_file in Path.iterdir(results_dir):
        get_tp_fn_fp(helix_file)
    pass

if __name__ == "__main__":
    main()