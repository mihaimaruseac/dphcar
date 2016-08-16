#!/bin/bash

datafiles="mushroom retail pumsb_star"
rmaxes="3 5"
nis="10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 170 180 190 200"
rfs="10 100 1000 10000"
seeds="42 132 243 314 4562"

for df in ${datafiles}; do
    for rmax in ${rmaxes}; do
        for ni in ${nis}; do
            for rf in ${rfs}; do
                for seed in ${seeds}; do
                    echo ${df} ${rmax} ${ni} ${rf} ${seed}
                done
            done
        done
    done
done | xargs -n5 -P7 ./run_one_reduction.sh
