#!/usr/bin/env bash
set -e

# only built for use on N288 cluster

# example usage: bash clean-paths.sh <input-file> <workflow-dir>

# making sure by chance we are not overwriting a wanted file called 't'

if [ -s t ]; then
    printf "\n    This simple program temporarily writes to a file called 't'\n"
    printf "    Since that exists already here, we are not going to continue.\n\n"
    exit
fi


ROOT_DIR=$(echo $2 | awk '{N=split($0,a,"/"); for(i=0; i < N-1; i++) printf "%s/", a[i]}' | sed 's|//|/|')

 
sed -E 's|.*/GLDS_Datasets/(.+)|\1|g'  ${1} \
    | sed -E 's|.+/miniconda.+/envs/[^/]*/||g' \
    | sed -E 's|/[^ ]*/GLDS-|GLDS-|g' \
    | sed -E 's|/[a-z]{6}/[^ ]*|<path-removed-for-security-purposes>|g' \
    | sed -E "s|${ROOT_DIR}||g" > t && mv t ${1}