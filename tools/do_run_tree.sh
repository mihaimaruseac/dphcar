#!/bin/bash

exp=$1
db=$2
rmax=$3
ni=$4
bf=$5
seed=$6
outdir=icde_2016_outputs/private/tree

./$exp inputs/data/datasets/${db}.dat 0.5 0.01 0.5 ${rmax} ${ni} ${bf}\
    ${seed} > ${outdir}/${exp}_${db}_rmax${rmax}_ni${ni}_bf${bf}_seed${seed}.txt
