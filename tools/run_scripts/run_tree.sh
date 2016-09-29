#!/bin/bash

seeds="42 132 243 314 4562"
dbs="mushroom pumsb_star"
epss="0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0"
rmaxs="3 5"
nis="50 60"
bfs="2 4 8"

for db in ${dbs}; do
    for eps in ${epss}; do
        for rmax in ${rmaxs}; do
            for ni in ${nis}; do
                for bf in ${bfs}; do
                    for seed in ${seeds}; do
                        echo ${db} ${eps} ${rmax} ${ni} ${bf} ${seed}
                    done
                done
            done
        done
    done
done | xargs -n6 -P7 ./do_run_tree.sh
