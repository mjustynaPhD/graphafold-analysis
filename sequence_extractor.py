import os

dir_path = "large/"
out_dir = "large_sequences/"
files = os.listdir(dir_path)
files = [f for f in files if f.endswith('.idx')]
print(files)

os.makedirs(out_dir, exist_ok=True)
for f in files:
    with open(dir_path + f, 'r') as file:
        lines = file.readlines()
        lines = [l.strip() for l in lines]
        lines = [l.split('.')[1][0] for l in lines]
        seq = "".join(lines)
        seq = f">{f.split('.')[0]}\n{seq}\n"
        print(seq)
    with open(out_dir + f.split('.')[0] + ".fasta", 'w') as out:
        out.write(seq)