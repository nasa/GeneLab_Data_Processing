#!/usr/bin/env bash
set -e
# only built for use on N288 cluster
# example usage: bash clean-paths.sh <input-file>
# making sure by chance we are not overwriting a wanted file called 't'

if [ -s t ]; then
printf "\n This simple program temporarily writes to a file called 't'\n"
printf " Since that exists already here, we are not going to continue.\n\n"
exit
fi


sed 's|/global/data/Data_Processing/Metagenomics_Datasets/GLDS_Datasets/||g' ${1} \
| sed 's|/global/data/Data_Processing/Amplicon_Datasets/GLDS_Datasets/||g' \
| sed 's|/global/data/Data_Processing/Metagenomics_Datasets/||g' \
| sed 's|/global/data/Data_Processing/Amplicon_Datasets/||g' \
| sed 's|/global/smf/miniconda38_admin/envs/[^/]*/||g' \
| sed 's|/[^ ]*/GLDS-|GLDS-|g' \
| sed 's|/global/[^ ]*|<path-removed-for-security-purposes>|g' > t && mv t ${1}