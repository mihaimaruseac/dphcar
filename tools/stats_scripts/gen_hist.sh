#!/bin/bash

file=$1

process() {
    echo
    echo "Histogram:"
    for d in `seq 9 -1 0`; do
        echo '>0.'${d}00000: `grep -e '->' ${file} | sort | uniq | grep "0\.$d[^ ]*\$" | wc -l`;
    done
}

process ${file} >> ${file}
