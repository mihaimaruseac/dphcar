#!/bin/bash

db=$1
eps=$2
rmax=$3
ni=$4
bf=$5
seed=$6
outdir=icde_2016_outputs/private/tree3

./dph inputs/data/datasets/${db}.dat\
    inputs/data/datasets/dph-recall/${db}.dat_${rmax}_${ni}\
    ${eps} 0.01 0.5 ${rmax} ${ni} ${bf} ${seed} \
    > ${outdir}/${db}_eps${eps}_rmax${rmax}_ni${ni}_bf${bf}_seed${seed}.txt
