#!/bin/bash

plot () {
    k=${1}

    cat << END | gnuplot
set xrange [0:1.1]
set xlabel "{/Symbol e}"
set yrange [0:1.1]
set ylabel "Precision"
set key top left
set terminal post eps enhanced font "Helvetica,28"
set output "k${k}.eps"
plot \
    "<awk -F, '{if(\$11 == ${k} && \$9 == 3 && \$10 == 0){print \$6,\$19}}' new_results.csv"\
        w lp ps 2 lt 1 pt 4 title "HCR 3",\
    "<awk -F, '{if(\$11 == ${k} && \$9 == 5 && \$10 == 0){print \$6,\$19}}' new_results.csv"\
        w lp ps 2 lt 1 pt 6 title "HCR 5",\
    "<awk -F, '{if(\$11 == ${k} && \$9 == 3 && \$10 == 3){print \$6,\$19}}' new_results.csv"\
        w lp ps 2 lt 1 pt 5 title "GRF 3",\
    "<awk -F, '{if(\$11 == ${k} && \$9 == 5 && \$10 == 4){print \$6,\$19}}' new_results.csv"\
        w lp ps 2 lt 1 pt 7 title "GRF 5"
END
}

plot 10
plot 50
plot 100
