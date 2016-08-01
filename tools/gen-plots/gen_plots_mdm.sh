#!/bin/bash

csvfile="data.csv"
datasets="bk10 bk11 bk14 mushroom pumsb"
gdatasets="bk10\|bk11\|bk14\|mushroom"
lengths="3 4 5 6 7"
ks="10 25 50 100"
ngcsvfile="ngrams_data.csv"
ngdatasets="bk11 mushroom"
nglengths="3"

plot_one () {
    infile=${1}
    data=${2}
    length=${3}
    k=${4}
    file=`mktemp`

    grep "${data},${length},${k}," ${infile} > ${file}
    c=`wc -l ${file} | cut -d\  -f1`

    if [ $c -eq 10 ]; then
        cat << END | gnuplot
set datafile separator ","
set xrange [0:1.1]
set xlabel "{/Symbol e}"
set yrange [-0.1:1.1]
set ylabel "Precision"
set key top left
set terminal png #post eps enhanced font "Helvetica,28"
set output "single-${infile}-${data}-${length}-${k}rules.png" #eps"
plot "${file}" using 4:6 w lp ps 2 lt 1 pt 2 title "conf 90%",\
     "${file}" using 4:8 w lp ps 2 lt 2 pt 4 title "conf 70%",\
     "${file}" using 4:10 w lp ps 2 lt 3 pt 6 title "conf 50%"
END
    fi

    rm ${file}
}

# plot_one csvfile dataset length k
for dataset in ${datasets}; do
    for len in ${lengths}; do
        for k in ${ks}; do
            plot_one ${csvfile} ${dataset} ${len} $k
        done
    done
done
for dataset in ${ngdatasets}; do
    for len in ${nglengths}; do
        for k in ${ks}; do
            plot_one ${ngcsvfile} ${dataset} ${len} $k
        done
    done
done

plot_baseline () {
    data=${1}
    length=${2}
    k=${3}
    file=`mktemp`

    # output is our method and then, after empty line, ngram method
    grep "${data},${length},${k}," ${csvfile} >> ${file}
    echo "" >> ${file}
    grep "${data},${length},${k}," ${ngcsvfile} >> ${file}
    c=`wc -l ${file} | cut -d\  -f1`

    if [ $c -eq 21 ]; then
        cat << END | gnuplot
set datafile separator ","
set xrange [0:1.1]
set xlabel "{/Symbol e}"
set yrange [-0.1:1.1]
set ylabel "Precision"
set key top left
set terminal png #post eps enhanced font "Helvetica,28"
set output "baseline-${data}-${length}-${k}rules.png" #eps"
plot\
    "${file}" using  4:8 every :::0::1 w lp ps 2 lt 1 pt 2 title "HCSR - 70%",\
    "${file}" using  4:9 every :::0::1 w lp ps 2 lt 1 pt 4 title "HCSR - 60%",\
    "${file}" using 4:10 every :::0::1 w lp ps 2 lt 1 pt 6 title "HCSR - 50%",\
    "${file}" using  4:8 every :::1::2 w lp ps 2 lt 2 pt 2 title "ngram - 70%",\
    "${file}" using  4:9 every :::1::2 w lp ps 2 lt 2 pt 4 title "ngram - 60%",\
    "${file}" using 4:10 every :::1::2 w lp ps 2 lt 2 pt 6 title "ngram - 50%"
END
    fi

    rm ${file}
}

# plot comparison with the baseline (ngrams)
for dataset in ${ngdatasets}; do
    for len in ${nglengths}; do
        for k in ${ks}; do
            plot_baseline ${dataset} ${len} $k
        done
    done
done

plot_vs_k () {
    length=${1}
    conf=${2}
    confCol=${3}
    epsilon="0\.5"
    file=`mktemp`

    grep "[a-z0-9A-Z]*,${length},[0-9]*,${epsilon}," ${csvfile} | grep ${gdatasets} | sort -t',' -nk3 >> ${file}
    c=`wc -l ${file} | cut -d\  -f1`

    if [ $c -eq 16 ]; then
        cat << END | gnuplot
set datafile separator ","
set logscale x 2
set xrange [8:128]
set xlabel "Rule count (k)"
set yrange [-0.1:1.1]
set ylabel "Precision"
set key top left
set terminal png #post eps enhanced font "Helvetica,28"
set output "k-${length}-${conf}.png" #eps"
plot\
    "${file}" using 3:${confCol} every 4::0 w lp ps 2 lt 1 pt 2 title "bk10 - ${conf}%",\
    "${file}" using 3:${confCol} every 4::1 w lp ps 2 lt 1 pt 4 title "bk11 - ${conf}%",\
    "${file}" using 3:${confCol} every 4::2 w lp ps 2 lt 1 pt 6 title "bk14 - ${conf}%",\
    "${file}" using 3:${confCol} every 4::3 w lp ps 2 lt 1 pt 8 title "mushroom - ${conf}%"
END
    fi

    rm ${file}
}

## plot precision vs k
for len in ${lengths}; do
    plot_vs_k ${len} 50 10
    plot_vs_k ${len} 60 9
    plot_vs_k ${len} 70 8
done
