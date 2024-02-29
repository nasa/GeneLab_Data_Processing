#!/bin/bash
set -u

software_versions_file="software_versions_GLmicroarray.md"

# Read the markdown table
while read -r line; do
    # Extract program, version, and link
    program=$(echo "$line" | awk -F'|' '{gsub(/^[[:blank:]]+|[[:blank:]]+$/,"",$1); print $1}')
    version=$(echo "$line" | awk -F'|' '{gsub(/^[[:blank:]]+|[[:blank:]]+$/,"",$2); print $2}')

    # Skip the header row and rows without version information
    if [[ $program != "Program" && $version != "Version" && ! -z $version ]]; then
        # Replace invalid characters in program name with underscores
        sanitized_program=$(echo "$program" | tr -cd '[:alnum:]_')

        # Create environment variable name
        env_var_name="${sanitized_program}_VERSION"

        # Set the environment variable
        export "$env_var_name=$version"
    fi
done < <(sed -n '/|/p' "$software_versions_file" | sed 's/^ *|//;s/|$//')

# Print the extracted versions
env | grep "_VERSION"

# Determine mapped sections
source meta.sh

# List of organisms
organism_list=("Mus musculus" "Homo sapiens" "Rattus norvegicus" "Drosophila melanogaster" "Danio rerio")

# Check the value of 'organism' variable and set 'GENE_MAPPING_STEP' accordingly
if [[ $organism == "Arabidopsis thaliana" ]]; then
    GENE_MAPPING_STEP="Ensembl gene ID mappings were retrieved for each probeset using the Plants Ensembl database ftp server (plants.ensembl.org, release 54). "
elif [[ " ${organism_list[*]} " == *"${organism//\"/}"* ]]; then
    GENE_MAPPING_STEP="Ensembl gene ID mappings were retrieved for each probeset using biomaRt (version ${biomaRt_VERSION}), Ensembl database (ensembl.org, release 107). "
else
    GENE_MAPPING_STEP="TBD"
fi

# Check the value of 'organism' variable and set 'GENE_MAPPING_STEP' accordingly
if [[ $organism == "Arabidopsis thaliana" ]]; then
    GENE_ANNOTATION_DB="org.At.tair.db"
elif [[ $organism == "Homo sapiens" ]]; then
    GENE_ANNOTATION_DB="org.Hs.eg.db"
elif [[ $organism == "Mus musculus" ]]; then
    GENE_ANNOTATION_DB="org.Mm.eg.db"
elif [[ $organism == "Rattus norvegicus" ]]; then
    GENE_ANNOTATION_DB="org.Rn.eg.db"
elif [[ $organism == "Drosophila melanogaster" ]]; then
    GENE_ANNOTATION_DB="org.Dm.eg.db"
elif [[ $organism == "Caenorhabditis elegans" ]]; then
    GENE_ANNOTATION_DB="org.Ce.eg.db"
else
    GENE_ANNOTATION_DB="TBD"
fi

# Read the template file
template="Data were processed as described in GL-DPPD-7114 (https://github.com/nasa/GeneLab_Data_Processing/blob/master/Microarray/Affymetrix/Pipeline_GL-DPPD-7114_Versions/GL-DPPD-7114.md) using NF_MAAffymetrix version 1.0.3 (https://github.com/nasa/GeneLab_Data_Processing/tree/NF_MAAffymetrix_1.0.3/Microarray/Affymetrix/Workflow_Documentation/NF_MAAffymetrix).  In short, a RunSheet containing raw data file location and processing metadata from the study's *ISA.zip file was generated using dp_tools (version ${dp_tools_VERSION}). The raw array data files were loaded into R (version ${R_VERSION}) using oligo (version ${oligo_VERSION}). Raw data quality assurance density plot, pseudo images, MA plots, and boxplots were generated using oligo (version ${oligo_VERSION}). The raw probe level intensity data was background corrected and normalized across arrays via the oligo (version ${oligo_VERSION}) quantile method. Normalized probe level data quality assurance density plot, pseudo images, MA plots, and boxplots were generated using oligo (version ${oligo_VERSION}).  Normalized probe level data was summarized to the probeset level using the oligo (version ${oligo_VERSION}) RMA method. ${GENE_MAPPING_STEP} Differential expression analysis was performed in R (version ${R_VERSION}) using limma (version ${limma_VERSION}); all groups were compared pairwise for each probeset to generate a moderated t-statistic and associated p- and adjusted p-value. Gene annotations were assigned for every probeset that mapped to exactly one Ensembl gene ID using the custom annotation tables generated in-house as detailed in GL-DPPD-7110 (https://github.com/nasa/GeneLab_Data_Processing/blob/GL_RefAnnotTable_1.0.0/GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110/GL-DPPD-7110.md), with STRINGdb (version 2.8.4), PANTHER.db (version 1.0.11), and ${GENE_ANNOTATION_DB} (version 3.15.0)."

# Output the filled template
echo "$template" > PROTOCOL_GLmicroarray.txt