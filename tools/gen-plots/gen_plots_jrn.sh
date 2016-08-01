#!/bin/bash

databases="mushroom pumsb retail"
lengths="3 5 7"
declare -a confs=(0.9 0.8 0.7 0.6 0.5)

terminal='post eps enhanced font "Helvetica,28"'
plotext="eps"

#--

# plot_one_heatmap ${database} ${ix} - plot heatmap for database at conf threshold ${confs[ix]}
plot_one_heatmap() {
    datafile=$1
    i=$2
    conf=${confs[i]}

    output=heatmap_${datafile}_${conf}.${plotext}
    datafile=${datafile}.csv

    cat << END | gnuplot
set xrange [-0.5:9.5]
set yrange [-0.5:9.5]
set xlabel "#bins"
set ylabel "#layers"
set terminal ${terminal}
set output "${output}"
set view map
plot "${datafile}" index $i matrix w image
END
}

# plot_heatmap ${database} - plot all heatmaps for database
plot_heatmap() {
    for i in `seq 0 4`; do
        plot_one_heatmap $1 $i
    done
}

# plot_heatmaps - plot all heatmaps
plot_heatmaps() {
    for db in ${databases}; do
        for l in ${lengths}; do
            plot_heatmap ${db}_${l}
        done
    done
}

#--

# plot_ext_bins_data ${database} - plot precision/bins w/ c_0={..}, l_max={..}
plot_ext_bins_data() {
    database=$1

    l1=3
    l2=5
    inf1=${database}_${l1}.csv
    inf2=${database}_${l2}.csv

    c1i=4
    c2i=2
    c1=${confs[c1i]}
    c2=${confs[c2i]}
    oli1=$((12 * c1i + 2))
    oli2=$((12 * c2i + 2))

    output=ext_bins_${database}.${plotext}

    case ${database} in
        pumsb*)
            legend="bottom right"
            ;;
        *)
            legend="top right"
            ;;
    esac

    cat << END | gnuplot
set xrange [0:11]
set yrange [-0.1:1.1]
set xlabel "Bin count"
set ylabel "Precision"
set terminal ${terminal}
set output "${output}"
set key ${legend}
plot \
    "<awk 'NR == ${oli1} { for (i = 1; i <= NF; i++) print i,(\$i/100)}' ${inf1}"\
        w lp ps 2 lt 1 pt 4 title "c_0=${c1}, l_{max}=${l1}",\
    "<awk 'NR == ${oli2} { for (i = 1; i <= NF; i++) print i,(\$i/100)}' ${inf1}"\
        w lp ps 2 lt 1 pt 6 title "c_0=${c2}, l_{max}=${l1}",\
    "<awk 'NR == ${oli1} { for (i = 1; i <= NF; i++) print i,(\$i/100)}' ${inf2}"\
        w lp ps 2 lt 1 pt 5 title "c_0=${c1}, l_{max}=${l2}",\
    "<awk 'NR == ${oli2} { for (i = 1; i <= NF; i++) print i,(\$i/100)}' ${inf2}"\
        w lp ps 2 lt 1 pt 7 title "c_0=${c2}, l_{max}=${l2}"
END
}

# plot_ext_bins - plot precision/bins w/ c_0={...}, l_max={...}
plot_ext_bins() {
    # note that we don't plot all combinations
    for d in mushroom pumsb; do
        plot_ext_bins_data $d
    done
}

# plot_db_vs_bins_ix ${database} ${ix} - plot p vs bins for database at conf threshold ${confs[ix]}
plot_db_vs_bins_ix() {
    datafile=$1
    i=$2
    conf=${confs[i]}

    output=bins_${datafile}_${conf}.${plotext}
    title=${datafile}
    datafile=${datafile}.csv
    outputline=$((12 * i + 2))

    cat << END | gnuplot
set xrange [0:11]
set yrange [-10:110]
set xlabel "#bins"
set ylabel "Rules extracted"
set terminal ${terminal}
set output "${output}"
plot \
    "<awk 'NR == ${outputline} { for (i = 1; i <= NF; i++) print i,\$i}' ${datafile}"\
        w lp ps 2 lt 1 pt 4 title "${title}"
END
}

# plot_db_vs_bins ${database} - plot p vs bins for database
plot_db_vs_bins() {
    for i in `seq 0 4`; do
        plot_db_vs_bins_ix $1 $i
    done
}

# plot_db_bins ${database} - plot p vs bins for all databases
plot_db_bins() {
    for db in ${databases}; do
        for l in ${lengths}; do
            plot_db_vs_bins ${db}_${l}
        done
    done
}

#--

# plot_ext_layers_data ${database} - plot precision/layers w/ c_0={...}, bins={...}
plot_ext_layers_data() {
    database=$1
    datafile=${database}.csv

    c1i=4
    c2i=2
    c1=${confs[c1i]}
    c2=${confs[c2i]}
    oli1=$((12 * c1i + 2))
    oli2=$((12 * c2i + 2))
    ole1=$((oli1 + 9))
    ole2=$((oli2 + 9))

    output=ext_layers_${database}.${plotext}

    columns="5 10"

    case ${database} in
        pumsb*)
            legend="bottom right"
            ;;
        *)
            legend="top right"
            ;;
    esac

    cat << END | gnuplot
set xrange [0:11]
set yrange [-0.1:1.1]
set xlabel "Layer count"
set ylabel "Precision"
set terminal ${terminal}
set output "${output}"
set key ${legend} Left
columns="${columns}"
plot\
    for [i=1:words(columns)] \
        "<awk 'NR == ${oli1}, NR == ${ole1} {print (1 + NR - ${oli1}), ($".word(columns, i)."/100)}' ${datafile}"\
            w lp ps 2 lt 1 pt (2+2*i) title "c_0=${c1}, ".sprintf("%3s", word(columns,i))." bins",\
    for [i=1:words(columns)] \
        "<awk 'NR == ${oli2}, NR == ${ole2} {print (1 + NR - ${oli2}), ($".word(columns, i)."/100)}' ${datafile}"\
            w lp ps 2 lt 1 pt (3+2*i) title "c_0=${c2}, ".sprintf("%3s", word(columns,i))." bins"
END
}

# plot_ext_layers - plot precision/layers w/ c_0={...}, bins={...}
plot_ext_layers() {
    # note that we don't plot all combinations
    for db in ${databases}; do
        for l in ${lengths}; do
            plot_ext_layers_data ${db}_${l}
        done
    done
}
# plot_db_vs_layers_ix ${database} ${ix} - plot p vs layers for database at conf threshold ${confs[ix]}
plot_db_vs_layers_ix() {
    datafile=$1
    i=$2
    conf=${confs[i]}

    output=layers_${datafile}_${conf}.${plotext}
    title=${datafile}
    datafile=${datafile}.csv

    stline=$((12 * i + 2))
    enline=$((stline + 10))

    columns="5 10"

    cat << END | gnuplot
set xrange [0:11]
set yrange [-10:110]
set xlabel "#layers"
set ylabel "Rules extracted"
set terminal ${terminal}
set output "${output}"
columns="${columns}"
plot for [i=1:words(columns)] \
    "<awk 'NR == ${stline}, NR == ${enline} {print (1 + NR - ${stline}), $".word(columns, i)."}' ${datafile}"\
        w lp ps 2 lt 1 pt (2+2*i) title word(columns,i)." bins"
END
    }

# plot_db_vs_layers ${database} - plot p vs layers for database
plot_db_vs_layers() {
    for i in `seq 0 4`; do
        plot_db_vs_layers_ix $1 $i
    done
}

# plot_db_layers ${database} - plot p vs layers for all databases
plot_db_layers() {
    for db in ${databases}; do
        for l in ${lengths}; do
            plot_db_vs_layers ${db}_${l}
        done
    done
}

#--

# main: comment to not redo some plots
#plot_heatmaps
#plot_db_bins
#plot_db_layers
plot_ext_bins
plot_ext_layers
