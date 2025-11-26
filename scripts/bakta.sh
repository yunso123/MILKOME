#!/bin/bash

# Load modules 
module load tools
module load bakta/1.9.3

# path
BTdb="/home/projects/dtu_00032/db/bakta/db"
BTpath="/home/projects/dtu_00032/analysis/milk-cohort/meta/bakta/non-stool"
cd /home/projects/dtu_00032/analysis/milk-cohort || exit 1
mkdir -p $BTpath

# Map ids
samples="E100079703_L01_UDB-100 E100079703_L01_UDB-101"
IFS=' ' read -r -a samples_array <<< "$samples"

declare -A sample_map

for i in "${!samples_array[@]}"; do
    old_name="${samples_array[$i]}"
    new_name="META$(printf "%04d" $((i+1)))"
    sample_map["$old_name"]="$new_name"
    echo "Old name: $old_name -> New name: $new_name"
done

# Create a mapping file
mapping_file="/home/projects/dtu_00032/analysis/milk-cohort/meta/sample_mapping.txt"
: > "$mapping_file"
for old_name in "${!sample_map[@]}"; do
    echo "$old_name -> ${sample_map[$old_name]}" >> "$mapping_file"
done

# Run Bakta
for old_name in "${!sample_map[@]}"; do
    N="${sample_map[$old_name]}"
    echo 
    echo "Running bakta for $old_name"
    echo
    bakta --threads 40 \
          --db $BTdb \
          --compliant \
          --meta \
          --verbose \
          --output $BTpath/$N \
          --locus-tag $N \
          --prefix $N \
          "data/megahit/${old_name}/final.contigs.fa"
done
