#!/bin/bash

# Define base directories
base_input_dir="../../raw_data/"
base_output_dir="../../cellbender"

# Define samples and fpr values
samples=("A9_2" "A10_2" "A11_2" "A12_2" "B1_2" "B2_2") 
fpr_values=(0.01 0.03 0.05 0.07 0.1)

# Run cellbender with specified model type
run_cellbender_full() {
    local input_dir=$1
    local sample=$2
    local fpr=$3

    local input_file="$input_dir/$sample/raw_feature_bc_matrix"
    local output_dir="$base_output_dir/${fpr}_full/$sample"

    mkdir -p $output_dir
    cd $output_dir

    echo "Running full model for $input_file with fpr $fpr..."

    cellbender remove-background \
        --input $input_file \
        --output $output_dir/${sample}.h5 \
        --cuda \
        --fpr $fpr \
        --model full
}

# Loop through each sample and fpr value
for sample in "${samples[@]}"; do
    for fpr in "${fpr_values[@]}"; do
        run_cellbender_full "$base_input_dir" "$sample" "$fpr"
    done
done
