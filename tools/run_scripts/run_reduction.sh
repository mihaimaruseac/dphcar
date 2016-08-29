#!/bin/bash

rmaxes="3 5"
rfs="10 20 30 40 50 60 70 80 90 100 125 150 175 200"
seeds="42 132 243 314 4562"

run_on_db () {
    df=$1
    maxni=$2

    for rmax in ${rmaxes}; do
        for ni in `seq 10 10 ${maxni}`; do
            for rf in ${rfs}; do
                for seed in ${seeds}; do
                    echo ${df} ${rmax} ${ni} ${rf} ${seed}
                done
            done
        done
    done | xargs -n5 -P7 ./run_one_reduction.sh
}

run_on_db mushroom 120
run_on_db retail 1000
run_on_db pumsb_star 600
