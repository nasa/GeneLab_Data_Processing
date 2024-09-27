#!/usr/bin/env bash
set -e

ls benchmarks/ > benchmark-filenames.tmp

head -n 1 benchmarks/$( head -n 1 benchmark-filenames.tmp ) > benchmark-header.tmp

paste <( printf "process" ) benchmark-header.tmp > building-tab.tmp

for file in $(cat benchmark-filenames.tmp)
do

    cat <( paste <( echo ${file} | sed 's/-benchmarks.tsv//' ) <( tail -n +2 benchmarks/${file} ) ) >> building-tab.tmp

done

mv building-tab.tmp all-benchmarks.tsv
rm -rf benchmark-filenames.tmp benchmark-header.tmp
