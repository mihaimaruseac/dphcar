#!/bin/bash

seeds="42 132 243 314 4562"

run_on_db () {
    db=$1
    maxni=$2
    rmax=$3
    bfs="$4"

    for ni in `seq 70 10 ${maxni}`; do
        for bf in ${bfs}; do
            for exp in dph_*; do
                for seed in ${seeds}; do
                    echo ${exp} ${db} ${rmax} ${ni} ${bf} ${seed}
                done
            done
        done
    done | xargs -n6 -P7 ./do_run_tree.sh
}

#run_on_db mushroom 120 3 "1 2 4 8 16 32 64"
#run_on_db mushroom 80 5 "1 2 4 8 16 32 64"
#run_on_db pumsb_star 3000 3 "1 2 4 8 16 32 64"
run_on_db pumsb_star 90 5 "1 2 4 8 16 32 64"
#run_on_db retail 100 3 "1 2 4 8 16 32 64"
#run_on_db retail 50 5 "1 2 4 8"
