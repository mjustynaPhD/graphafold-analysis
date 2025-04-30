from pathlib import Path
import numpy as np

elements_dir = Path("structural-elements/casp")
idxs = Path("/home/mjustyna/data/graphafold_data/casp")
preds_methods = ['sincfold', 'ufold', 'spotrna', 'graphafold_casp']
save_distribution_path = Path("distributions")

def load_elements(path: Path) -> dict:
    """Load structural elements from the file."""
    elements = {
        "stems": [],
        "loops": [],
        "ss": [],
    }

    with Path.open(path / "elements.txt", "r") as f:
        lines = f.readlines()

    elements["stems"] = [l for l in lines if l.startswith("Stem")]
    elements["loops"] = [l for l in lines if l.startswith("Loop")]
    elements["ss"] = [l for l in lines if l.startswith("Single strand")]

    return elements

def parse_elements(elements: dict) -> tuple[dict, dict, dict]:
    """Parse structural elements into a set of base pairs ranges."""
    stems = parse(elements["stems"], "Stem")
    loops = parse(elements["loops"], "Loop")
    ss = parse(elements["ss"], "Single strand")
    return stems, loops, ss

def parse(elems: list, element:str = "Stem") -> dict:
    """Parse specific elements into a set of base pairs ranges."""
    elems = [s.replace(element, "").strip() for s in elems]
    elems = [s.split() for s in elems]
    elems_dict = {}
    for e in elems:
        if element == "Single strand" and len(e) == 7:
            e = e[1:]
        e_id = e[0]
        arr_e = np.array(e[1:])
        arr_e = arr_e.reshape((-1, 5))
        elems_dict[e_id] = arr_e[:, :2]
    return elems_dict

def parse_idx(idx_file: Path) -> dict:
    """Parse idx file into a dictionary."""
    idx_dict = {}
    idx_seq_dict = {}
    with Path.open(idx_file, "r") as f:
        lines = f.readlines()

    idxs = [l.split(',') for l in lines]

    for i, l in idxs:
        idx_dict[int(i)] = l[:2] + l[3:].strip()
        idx_seq_dict[int(i)] = l[2]

    return idx_dict, idx_seq_dict

def map_idx_reverse(idx_dict: dict) -> dict:
    """Map idx dictionary to reverse."""
    idx_dict_reverse = {}
    for k, v in idx_dict.items():
        if v in idx_dict_reverse:
            raise ValueError(f"Duplicate value found in idx_dict: {v}")
        idx_dict_reverse[v] = k
    return idx_dict_reverse

def preds_from_helix(helix_path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Load and parse predictions from the helix file."""
    if not helix_path.exists():
        print(f"Warning: Helix file does not exist: {helix_path}")
        return {}, {}
    preds = {}
    with Path.open(helix_path, "r") as f:
        lines = f.readlines()
    lines = [l.split(",") for l in lines]
    lines = [l[:-1] for l in lines]
    preds = [l for l in lines if l[2].startswith("PredictedGoodNonCanonical") or l[2].startswith("PredictedBadNonCanonical")]  # noqa: E501
    gts = [l for l in lines if l[2].startswith("PredictedGoodNonCanonical") or l[2].startswith("NotPredictedNonCanonical")] # noqa: E501
    preds = [l[:2] for l in preds]
    gts = [l[:2] for l in gts]
    preds = [[int(i) for i in l] for l in preds]
    gts = [[int(i) for i in l] for l in gts]
    preds = np.array(preds)
    gts = np.array(gts)
    # if helix_path startswith "graphafold" then add 1 to all values
    if str(helix_path).startswith("graphafold"):
        preds += 1
        gts += 1
    return preds, gts

def elements_to_indeces(elements: dict, idx_rev: dict) -> dict:
    """Convert elements like A.C10 to integer indeces using idx reverse mapping."""
    out = {}
    for k, v in elements.items():
        out[k] = []
        for i, j in v:
            if i in idx_rev and j in idx_rev:
                out[k].append([idx_rev[i], idx_rev[j]])
            else:
                print(f"Index {i} or {j} not found in idx reverse mapping.")

    out = {k: np.array(v) for k, v in out.items()}
    return out

def detect_regions(stems: dict,
                   loops: dict,
                   ss: dict,
                   pairs: np.ndarray) -> dict:
    """Detect regions in the pairs based on stems, loops, and ss."""
    #hierarchy: stems, loops, ss (if any remaining)

    pass

def pairs_distribution(preds:np.ndarray, gts:np.ndarray, idx:dict) -> tuple[dict, dict]:
    """Calculate the distribution of pairs in ground-truth and predictions."""
    gt_distr = distribution(gts, idx)
    preds_distr = distribution(preds, idx)
    return preds_distr, gt_distr

def distribution(arr:np.ndarray, idx:dict) -> dict:
    """Calculate the distribution of pairs in the array."""
    distr = {}
    for i, j in arr:
        pair = sorted([idx[i], idx[j]]) # get nucleotides
        pair = "-".join(pair)
        if pair in distr:
            distr[pair] += 1
        else:
            distr[pair] = 1
    return distr

def save_distribution(distr:dict, path:Path, overwrite:bool=True) -> None:
    """Save the distribution to a CSV file."""
    if not overwrite and path.exists():
        return
    if not path.parent.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
    with Path.open(path, "w") as f:
        for k, v in distr.items():
            f.write(f"{k},{v}\n")
    print(f"Saved distribution to {path}")


def main() -> None:
    """Load and parse structural elements from the directory."""
    files = Path.iterdir(elements_dir)
    for f in files:
        print(f.name)
        elements = load_elements(f)
        stems, loops, ss = parse_elements(elements)
        # add .idx extension to the file name
        idx_file = idxs / f"{f.name}.idx"
        idx_parsed, idx_seq_parsed = parse_idx(idx_file)
        idx_parsed_reverse = map_idx_reverse(idx_parsed)
        stems = elements_to_indeces(stems, idx_parsed_reverse)
        loops = elements_to_indeces(loops, idx_parsed_reverse)
        ss = elements_to_indeces(ss, idx_parsed_reverse)
        for m in preds_methods:
            helix_path = Path(f"{m}_results/{f.name}.helix")
            preds, gts = preds_from_helix(helix_path)
            preds_dit, gt_dist = pairs_distribution(preds, gts, idx_seq_parsed)
            save_distribution(preds_dit, save_distribution_path / f"{m}_{f.name}.csv")
            save_distribution(gt_dist, save_distribution_path / f"gt_{f.name}.csv", overwrite=False)
            # do something with the preds
    

if __name__ == "__main__":
    main()