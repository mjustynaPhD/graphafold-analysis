import os

from Bio import PDB

desc_dir = "descs"

def main():
    files = os.listdir(desc_dir)
    files = sorted([f for f in files if f.endswith('.pdb')])
    for f in files:
        # load with biopython and check number of residues in file
        parser = PDB.PDBParser()
        structure = parser.get_structure(f, os.path.join(desc_dir, f))
        model = structure[0]
        num_residues = len(list(model.get_residues()))
        # get the number of residues in the first chain
        chain = model.get_list()[0]
        num_residues_chain = len(list(chain.get_residues()))
        print(f, num_residues, num_residues_chain)

if __name__ == "__main__":
    main()

