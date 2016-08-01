#!/bin/bash

datafiles="mushroom pumsb retail"
datafilesnames="Mushroom Pumsb Retail"
alphas="f-1000 f-2000 p-1 p-2 p-5 p-10"
nmaxs="3 5 7"
ks="10 20 50 100 200"
epsilons="03 05 07 09" # remove the dot, LaTeX doesn't like it

# plot precision based on epsilon, varying all other params
plot-prec-eps() {
    column=$1
    conf=$2
    cat << END | gnuplot
set terminal post eps enhanced font "Helvetica,28"
set yrange [-0.1:1.1]
set xrange [0:1.1]
set xlabel "{/Symbol e}"
set ylabel "Precision"
set nokey
#set label "  Retail k = 10" at 0.05, 1.17 point ps 3 pt 8
#set label "  Retail k = 50" at 0.6, 1.17 point ps 3 pt 9
start = 0
do for [alpha = 1:words('${alphas}')] {
    do for [nmax = 1:words('${nmaxs}')] {
        # output filename
        set output 'precision-${conf}-'.word('${alphas}', alpha).\
            '-'.word('${nmaxs}', nmax).'.eps'
        # plot
        plot for [i=1:words('${datafiles}')]\
            word('${datafiles}',i).".csv" using 1:${column}\
                every ::start+1::start+10 w lp ps 2 lt 1 pt 2+2*i\
                title word('${datafiles}',i),\
             for [i=1:words('${datafiles}')]\
            word('${datafiles}',i).".csv" using 1:${column}\
                every ::60+start+1::60+start+10 w lp ps 2 lt 1 pt 3+2*i\
                title word('${datafiles}',i).'2'
        start = start + 10
    }
}
END
}

## 90% confidence to 50%
#plot-prec-eps 5 90
#plot-prec-eps 6 80
#plot-prec-eps 7 70
#plot-prec-eps 8 60
#plot-prec-eps 9 50

# plot precision based on alpha, varying all other params
plot-prec-alpha() {
    ifile=$1
    datafile=${ifile%%.csv}
    cat << END | gnuplot
set terminal post eps enhanced font "Helvetica,28"
set yrange [-0.1:1.1]
set xrange [0:1.1]
set xlabel "{/Symbol e}"
set ylabel "Precision"
set nokey
#set label "  {/Symbol a} = 10%, c = 70%" at 0.05, 1.17 point ps 3 pt 8
#set label "  {/Symbol a} = 10%, c = 50%" at 0.6, 1.17 point ps 3 pt 9
start = 301
do for [k = 1:words('${ks}')] {
    do for [nmax = 1:words('${nmaxs}')] {
        # output filename
        set output 'alpha-${datafile}-'.word('${nmaxs}', nmax).\
            '-'.word('${ks}', k).'.eps'
        # plot
        plot \
            '${ifile}' using 1:7 every ::start::start+9\
                w lp ps 2 lt 1 pt 4\
                title '{/Symbol a} = 1%, c = 70%',\
            '${ifile}' using 1:9 every ::start::start+9\
                w lp ps 2 lt 1 pt 5\
                title '{/Symbol a} = 1%, c = 50%',\
            '${ifile}' using 1:7 every ::start+300::start+300+9\
                w lp ps 2 lt 1 pt 6\
                title '{/Symbol a} = 5%, c = 70%',\
            '${ifile}' using 1:9 every ::start+300::start+300+9\
                w lp ps 2 lt 1 pt 7\
                title '{/Symbol a} = 5%, c = 50%',\
            '${ifile}' using 1:7 every ::start+450::start+450+9\
                w lp ps 4 lt 1 pt 8\
                title '{/Symbol a} = 10%, c = 70%',\
            '${ifile}' using 1:9 every ::start+450::start+450+9\
                w lp ps 2 lt 1 pt 9\
                title '{/Symbol a} = 10%, c = 50%'
        start = start + 10
    }
}
END
}

for i in ${datafiles}; do
    plot-prec-alpha ${i}.csv
done

# plot precision based on nmax, varying all other params
plot-prec-nmax() {
    column=$1
    conf=$2
    cat << END | gnuplot
set terminal post eps enhanced font "Helvetica,28"
set yrange [-0.1:1.1]
set xrange [2:8]
set xlabel "{n_{max}}"
set ylabel "Precision"
set key bottom right
start = 0
do for [alpha = 1:words('${alphas}')] {
    do for [k = 1:words('${ks}')] {
        do for [epsilon = 1:words('${epsilons}')] {
            # output filename
            set output 'nmax-${conf}-'.word('${alphas}', alpha).\
                '-'.word('${ks}', k).\
                '-'.word('${epsilons}', epsilon).'.eps'
            # plot
            plot for [i=1:words('${datafiles}')]\
                word('${datafiles}',i).".csv" using 2:${column}\
                    every 10::start+2*epsilon+1::start+2*epsilon+30 w lp ps 2 lt 1 pt 2+2*i\
                    title word('${datafiles}',i)
        }
        start = start + 30
    }
}
END
}

# 90% confidence
#plot-prec-nmax 5 90
#plot-prec-nmax 6 80
#plot-prec-nmax 7 70
#plot-prec-nmax 8 60
#plot-prec-nmax 9 50

# plot precision based on number of rules, varying all other params
plot-prec-k() {
    column=$1
    conf=$2
    cat << END | gnuplot
set terminal post eps enhanced font "Helvetica,28"
set yrange [-0.1:1.1]
set logscale x 2
set xrange [8:256]
set xlabel "Rule count (k)"
set ylabel "Precision"
set key bottom left
start = 0
do for [alpha = 1:words('${alphas}')] {
do for [nmax = 1:words('${nmaxs}')] {
do for [epsilon = 1:words('${epsilons}')] {
            # output filename
            set output 'k-${conf}-'.word('${alphas}', alpha).\
                '-'.word('${nmaxs}', nmax).\
                '-'.word('${epsilons}', epsilon).'.eps'
            # plot
            plot for [i=1:words('${datafiles}')]\
                word('${datafiles}',i).".csv" using 3:${column}\
                    every 30::start+(nmax-1)*10 + 2*epsilon+1::start+(nmax-1)*10 + 2*epsilon+150 w lp ps 2 lt 1 pt 2+2*i\
                    title word('${datafilesnames}',i)
        }
    }
    start = start + 150
}
END
}

# 90% confidence
#plot-prec-k 5 90
#plot-prec-k 6 80
plot-prec-k 7 70
#plot-prec-k 8 60
plot-prec-k 9 50

# plot baseline (PrivBasis) comparison
titles="HCR-0.9 PB-0.9 HCR-0.8 PB-0.8 HCR-0.7 PB-0.7 HCR-0.6 PB-0.6 HCR-0.5 PB-0.5"
pts="- - - - 4 5 6 7 8 9"
plot-baseline () {
    ifile=$1
    datafile=${ifile%%.csv}
    cat << END | gnuplot
set terminal post eps enhanced font "Helvetica,28"
set output "baseline-${datafile}.eps"
set yrange [-0.1:1.1]
set xrange [0:1.1]
set xlabel "{/Symbol e}"
set ylabel "Precision"
set key bottom right
plot for [i=5:words('${titles}')] '${ifile}' using 1:i+1\
    w lp ps 2 lt 1 pt word('${pts}', i)\
    title word('${titles}', i)
END
}

#for i in pb-*.csv; do
#    plot-baseline $i
#done

# plot expanded comparison (1: precision, 2: rules)
titles="3 5 7"
pts="4 6 8"
plot-expanded-1 () {
    ifile=$1
    datafile=${ifile%%.csv}
    column=$2
    conf=$3
    cat << END | gnuplot
set terminal post eps enhanced font "Helvetica,28"
set output "${datafile}-${conf}-precision.eps"
set yrange [-0.1:1.1]
set xrange [0:1.1]
set xlabel "{/Symbol e}"
set ylabel "Precision"
set key bottom right
plot for [i=1:words('${titles}')] '${ifile}' using 1:${column}\
    every :::i::i\
    w lp ps 2 lt 1 pt word('${pts}', i)\
    title 'length '.word('${titles}', i)
END
}

plot-expanded-2 () {
    ifile=$1
    datafile=${ifile%%.csv}
    cat << END | gnuplot
set terminal post eps enhanced font "Helvetica,28"
set output "${datafile}-rules.eps"
set yrange [0:12]
set ytics 1
set xrange [0:1.1]
set xlabel "{/Symbol e}"
set ylabel "Ratio of extra rules"
set key top right
plot for [i=1:words('${titles}')] '${ifile}' using 1:3 \
    every :::i::i\
    w lp ps 2 lt 1 pt word('${pts}', i)\
    title 'length '.word('${titles}', i)
END
}

do_plot-1 () {
    column=$1
    conf=$2
    for i in expanded-*.csv; do
        plot-expanded-1 $i $column $conf
    done
}

#do_plot-1 4 90
#do_plot-1 5 80
#do_plot-1 6 70
#do_plot-1 7 60
#do_plot-1 8 50
#
#for i in expanded*.csv; do
#    plot-expanded-2 $i
#done
