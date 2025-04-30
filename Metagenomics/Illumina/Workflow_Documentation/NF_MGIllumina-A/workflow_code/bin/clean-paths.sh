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

FILE=$1
ROOT_DIR=$(echo $2 | awk '{N=split($0,a,"/"); for(i=0; i < N-1; i++) printf "%s/", a[i]}' | sed 's|//|/|')


# Remove path in paired end runsheet
if [ `awk 'NR==1{print}' ${FILE} | grep -c reverse` -gt 0 ]; then

      awk 'BEGIN{FS=OFS=","} NR==1{print} NR>1{split($2, f, "/");split($3, r, "/"); print $1,f[length(f)],r[length(r)],$4}' ${FILE} > temp && mv temp ${FILE}

# Remove path in single end runsheet
elif [ `awk 'NR==1{print}' ${FILE} | grep -c forward` -gt 0 ]; then


     awk 'BEGIN{FS=OFS=","} NR==1{print} NR>1{split($2, f, "/"); print $1,f[length(f)],$3}' ${FILE} > temp && mv temp ${FILE}

fi
 
sed -E 's|.*/GLDS_Datasets/(.+)|\1|g' ${FILE} \
    | sed -E 's|.+/miniconda.+/envs/[^/]*/||g' \
    | sed -E 's|/[^ ]*/GLDS-|GLDS-|g' \
    | sed -E 's|/[a-z]{6}/[^ ]*|<path-removed-for-security-purposes>|g' \
    | sed -E "s|${ROOT_DIR}||g" > t && mv t  ${FILE}
