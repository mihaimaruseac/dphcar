#!/bin/bash

seeds="42 132 243 314 4562"

run_on_db () {
    db=$1
    rmax=$2
    nis="$3"
    bfs="$4"

    for ni in ${nis}; do
        for bf in ${bfs}; do
            for seed in ${seeds}; do
                echo ${db} ${rmax} ${ni} ${bf} ${seed}
            done
        done
    done | xargs -n5 -P7 ./do_run_tree.sh
}

#run_on_db mushroom 5 "60 70" "64"
#run_on_db mushroom 5 "110 120" "1 2 4 8 16"
#run_on_db pumsb_star 3 "90 100 150 200 250 300" "64"
#run_on_db pumsb_star 5 "30 40" "16 32 64"
run_on_db pumsb_star 5 "60 70 80 90 100 150 200 250 300" "8 16 32"
