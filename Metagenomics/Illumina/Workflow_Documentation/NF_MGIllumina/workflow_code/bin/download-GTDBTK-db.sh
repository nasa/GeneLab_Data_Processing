#!/usr/bin/env bash

set -e

echo "Downloading the GTDB-Tk database to ${GTDBTK_DATA_PATH}..."

# GTDBTK_DB_PATH is defined in build.sh, store the db there


DB_URL=$1
TAR_FILE=$(basename ${DB_URL})


downloadFile=true

while ${downloadFile};do

    wget --timeout=3600 --tries=0 --continue ${DB_URL} -P ${GTDBTK_DATA_PATH} && downloadFile=false

done

tar xvzf ${GTDBTK_DATA_PATH}/${TAR_FILE} -C ${GTDBTK_DATA_PATH} --strip 1
rm ${GTDBTK_DATA_PATH}/${TAR_FILE}

echo "GTDB-Tk database has been successfully downloaded."

exit 0
