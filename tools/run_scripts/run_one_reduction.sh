#!/bin/bash

datafile=$1
rmax=$2
ni=$3
rf=$4
seed=$5
outdir=icde_2016_outputs/private/reduction_factor_2/10rules

./dph p inputs/data/datasets/${datafile}.dat 0.5 0.01 0.5 ${rmax} 10\
    ${ni} ${rf} ${seed} >\
    ${outdir}/${datafile}_ni${ni}_rf${rf}_rmax${rmax}_seed${seed}.txt
