#!/bin/bash

#datafiles="mushroom retail pumsb_star"
#datafiles="pumsb_star" # TODO: need to run retail 5
datafiles="retail"
#rmaxes="3 5"
rmaxes="5"
nis=`seq 10 10 2500`
rfs="10 20 50 80 100 200 500 1000 5000 10000"
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
