#!/bin/bash

datafile=$1
rmax=$2
ni=$3

seeds="42 132 243 314 4562"
outdir=icde_2016_outputs/private/phase_tran/10rules

for seed in ${seeds}; do
    (./dph p inputs/data/datasets/${datafile}.dat\
        0.5 0.01 0.5 ${rmax} 10 ${ni} ${seed} | tee\
        ${outdir}/${datafile}_ni${ni}_rmax${rmax}_seed${seed}.txt)&
done
