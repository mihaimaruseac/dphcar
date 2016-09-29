#!/bin/sh

datafile="$1"
datakfile="$2"
databases="mushroom pumsb"
lengths="3 5"

declare -A bf_keys=(
    [mushroom_3]="70 80 90 100 110"
    [mushroom_5]="70 80 90 100 110"
    [pumsb_3]="70 80 90 100 150 200"
    [pumsb_5]="70 80 90 100 150 200"
    )
declare -A Zs=(
    [mushroom_3]=300000
    [mushroom_5]=4.5e7
    [pumsb_3]=1.4e6
    [pumsb_5]=1.2e8
    )

#terminal='post eps enhanced font "Helvetica,28"'
#plotext='eps'
terminal='png'
plotext='png'

# plot_star ${database} ${rmax} ${keytitle} ${keycolumn} ${keys} $X ${xlabel} ${xcolumn} $Y ${ylabel} ${ycolumn}
plot_star () {
    database=$1
    rmax=$2
    keytitle=$3
    keycolumn=$4
    keys=$5

    X=$6
    xlabel=$7
    xcolumn=$8
    Y=$9
    ylabel=${10}
    ycolumn=${11}

    eps=0.5 #fixed
    c0=0.5 #fixed

    output=${Y}_vs_${X}_${database}_${rmax}.${plotext}

    cat << END | gnuplot
set xlabel "${xlabel}"
#set xrange [0:32] #  3
#set yrange [0:60000]
#set xrange [0:16] # mushroom 5
#set yrange [0:2e6]
set ylabel "${ylabel}"
set key top left
set terminal ${terminal}
set output "${output}"
key="${keys}"
plot \
    for [i=1:words(key)] \
        "<awk -F, '{if(index(\$1, \"${database}\") && \$2 == ${eps} && \$3 == ${c0}\
                    && \$4 == ${rmax} && \$${keycolumn} == ".word(key, i).")\
                        {print \$${xcolumn}, \$${ycolumn}}}' ${datafile}"\
            w lp ps 2 lt 1 pt (2+2*i) title "${keytitle}=".word(key, i)
END
}

plot_R_vs_bf () {
    database=$1
    rmax=$2
    keys="$3"
    plot_star ${database} ${rmax} n_i 5 "${keys}" bf "Branching factor" 6 R "Rules extracted" 7
}

plot_R50_vs_bf () {
    database=$1
    rmax=$2
    keys="$3"
    plot_star ${database} ${rmax} n_i 5 "${keys}" bf "Branching factor" 6 R50 "Good rules extracted" 8
}

plot_R_vs_ni () {
    database=$1
    rmax=$2
    keys="$3"
    plot_star ${database} ${rmax} b 6 "${keys}" ni "Number of items" 5 R "Rules extracted" 7
}

plot_R50_vs_ni () {
    database=$1
    rmax=$2
    keys="$3"
    plot_star ${database} ${rmax} b 6 "${keys}" ni "Number of items" 5 R50 "Good rules extracted" 8
}

# plot_R50_v_R ${database} ${rmax} $z
plot_R50_v_R () {
    database=$1
    rmax=$2
    z=$3

    eps=0.5 #fixed
    c0=0.5 #fixed

    output=R50_vs_R_${database}_${rmax}.${plotext}

    cat << END | gnuplot
set xlabel "Rules extracted"
set ylabel "Good rules extracted"
set key top left
set terminal ${terminal}
set output "${output}"
plot \
    "<awk -F, '{if(index(\$1, \"${database}\") && \$2 == ${eps} && \$3 == ${c0} && \$4 == ${rmax})\
                        {print \$7, \$8}}' ${datafile} | sort -nk1"\
            w lp ps 2 lt 1 pt 4,\
    [0:$z] x w l lt 2
END
}

# plot_codaspy ${database} ${rmax}
plot_codaspy () {
    database=$1
    rmax=$2
    z=$3

    eps=0.5 #fixed
    c0=0.5 #fixed

    output=cmp_CODASPY_${database}_${rmax}.${plotext}

    cat << END | gnuplot
set xlabel "Rules extracted"
set xrange [0:256]
set ylabel "Precision"
set yrange [-.1:1.1]
set key bottom right
set terminal ${terminal}
set output "${output}"
plot \
    "<awk -F, '{if(index(\$1, \"${database}\") && \$2 == ${eps} && \$3 == ${c0}\
                    && \$4 == ${rmax} && \$5 == 1)\
                        {print \$6, \$7}}' ${datakfile}"\
            w lp ps 2 lt 1 pt 5 title "HCR, {/Symbol a}=1%"
END
}

for db in ${databases}; do
    for rmax in ${lengths}; do
        key=${db}_${rmax}
        # plot R vs *
        plot_R_vs_bf ${db} ${rmax} "${bf_keys[${key}]}"
        plot_R50_vs_bf ${db} ${rmax} "${bf_keys[${key}]}"
        plot_R_vs_ni ${db} ${rmax} "1 2 4 8 16 32 64"
        plot_R50_vs_ni ${db} ${rmax} "1 2 4 8 16 32 64"
        # plot R50 vs R
        plot_R50_v_R ${db} ${rmax} ${Zs[${key}]}
        # compare with CODASPY
        plot_codaspy ${db} ${rmax}
    done
done
