#!/bin/bash

datasets="kosarak mushroom retail"
epsilons="0.9 0.5 0.2"
lengths="3 5 7"

kosarak_th=9900
mushroom_th=81
retail_th=881

program="./dph"

data_dir="data"
dataset_dir="datasets"
freq_dir="frequencies"
dataset_ext="dat"

output_dir="outputs"
output_subdir="outputs.3"

# "short" for O(2^n) generation
# "long" for true O(3^n) generation
rule_gen="long"

# comment if there is no limit on item generation (RULE_LEN_LIMIT 0)
item_gen="limited"

# "xy" - displacement
# "delta" - delta
# "pow" - new
# "d" - distance
quality="xy"

# min conf parametere
conf=0.9

# epsilon share for the itemset counts
share=0.1

for dataset in ${datasets}; do
    for epsilon in ${epsilons}; do
        for length in ${lengths}; do
            for item in ${data_dir}/${freq_dir}/${dataset}/*; do
                limit=$((${dataset}_th))
                $program ${data_dir}/${dataset_dir}/${dataset}.${dataset_ext} ${item} ${conf} ${epsilon} ${share} ${length} ${limit} > ${output_dir}/${output_subdir}/${item##*/}_${conf}_${epsilon}_${share}_${length}_${limit}_${rule_gen}_${item_gen}_${quality}
            done
        done &
    done
done
