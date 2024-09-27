#!/usr/bin/env bash

# Addresses issue: https://github.com/nextflow-io/nextflow/issues/1210

CONFILE=${1:-nextflow.config}
OUTDIR=${2:-./singularity}

if [ ! -e $CONFILE ]; then
        echo "$CONFILE does not exist"
        exit
fi

TMPFILE=`mktemp`

CURDIR=$(pwd)

mkdir -p $OUTDIR

cat ${CONFILE}|grep 'container'|perl -lane 'if ( $_=~/container\s*\=\s*\"(\S+)\"/ ) { $_=~/container\s*\=\s*\"(\S+)\"/; print $1 unless ( $1=~/^\s*$/ or $1=~/\.sif/ or $1=~/\.img/ ) ; }' > $TMPFILE

cd ${OUTDIR}

while IFS= read -r line; do
        name=$line
        name=${name/:/-}
        name=${name//\//-}
        echo $name
        singularity pull ${name}.img docker://$line
done < $TMPFILE

cd $CURDIR