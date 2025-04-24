#!/bin/bash



DATA="/home/mjustyna/data/graphafold_data/sequences"
source .venv/bin/activate

for sequence in $(ls $DATA/*.fasta); do
    echo "Processing $sequence"
    # Extract the base name of the file (without path and extension)
    base_name=$(basename "$sequence" .fasta)
    
    # echo $DATA/$sequence
    # Run the Python script with the sequence file as input
    sincFold pred $sequence -o my_results/$base_name
done