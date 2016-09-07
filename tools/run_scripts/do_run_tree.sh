#!/bin/bash

db=$1
rmax=$2
ni=$3
bf=$4
seed=$5
outdir=icde_2016_outputs/private/tree3

./dph inputs/data/datasets/${db}.dat 0.5 0.01 0.5\
    ${rmax} ${ni} ${bf} ${seed} \
    > ${outdir}/${db}_rmax${rmax}_ni${ni}_bf${bf}_seed${seed}.txt
