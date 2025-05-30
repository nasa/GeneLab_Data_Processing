#!/usr/bin/env python
"""
This script generates a protocol text file for GeneLab RNA-seq data processing.
It reads software versions from a YAML file and incorporates other parameters.
"""

import argparse
import yaml
import os
import sys
from datetime import datetime
import pandas as pd
import re

def parse_args():
    parser = argparse.ArgumentParser(description='Generate protocol file for GeneLab RNA-seq pipeline')
    parser.add_argument('--mode', choices=['standard', 'microbes'], default='standard',
                        help='Processing mode (standard or microbes)')
    parser.add_argument('--outdir', required=True,
                        help='Output directory for the protocol file')
    parser.add_argument('--software_table', required=True,
                        help='Path to YAML file containing software versions')
    parser.add_argument('--assay_suffix', required=True,
                        help='Suffix for the assay type')
    parser.add_argument('--paired_end', required=True,
                        help='Boolean indicating paired-end sequencing')
    parser.add_argument('--strandedness', required=True,
                        help='Strandedness information')
    parser.add_argument('--has_ercc', required=True,
                        help='Boolean indicating presence of ERCC spike-ins')
    parser.add_argument('--workflow_version', default='unknown',
                        help='Version of the NF_RCP workflow manifest')
    parser.add_argument('--organism', required=False,
                        help='Organism name for the reference genome')
    parser.add_argument('--reference_source', required=False,
                        help='Source of the reference genome')
    parser.add_argument('--reference_version', required=False,
                        help='Version of the reference genome')
    parser.add_argument('--reference_fasta', required=True,
                        help='Path to the reference genome FASTA file')
    parser.add_argument('--reference_gtf', required=True,
                        help='Path to the reference genome GTF file')
    parser.add_argument('--runsheet', required=False,
                        help='Path to the runsheet CSV file')
    return parser.parse_args()

def read_software_versions(yaml_file):
    try:
        with open(yaml_file, 'r') as f:
            return yaml.safe_load(f)
    except Exception as e:
        sys.stderr.write(f"Error reading software versions file: {e}\n")
        sys.exit(1)

def generate_protocol_content(args, software_versions):
    # Get current date
    current_date = datetime.now().strftime("%Y-%m-%d")
    
    # Create header
    header = f"# GeneLab RNA-Seq Pipeline Protocol{args.assay_suffix}\n"
    header += f"# Date: {current_date}\n\n"
    
    # Add appropriate protocol reference based on mode
    if args.mode == 'microbes':
        header += "Data were processed as described in GL-DPPD-7115 "
        header += "(https://github.com/nasa/GeneLab_Data_Processing/blob/master/RNAseq/Pipeline_GL-DPPD-7115_Versions/GL-DPPD-7115.md), "
    else:
        header += "Data were processed as described in GL-DPPD-7101-G "
        header += "(https://github.com/nasa/GeneLab_Data_Processing/blob/master/RNAseq/Pipeline_GL-DPPD-7101_Versions/GL-DPPD-7101-G.md), "
    
    header += f"using NF_RCP version {args.workflow_version} "
    header += f"(https://github.com/nasa/GeneLab_Data_Processing/tree/NF_RCP-{args.workflow_version}/RNAseq/Workflow_Documentation/NF_RCP). "
    
    # Start building the description as a single paragraph
    # Add processing description with software versions
    trim_galore_version = software_versions.get('TrimGalore!', 'unknown')
    cutadapt_version = software_versions.get('Cutadapt', 'unknown')
    fastqc_version = software_versions.get('FastQC', 'unknown')
    multiqc_version = software_versions.get('MultiQC', 'unknown')
    
    description = f"In short, raw fastq files were filtered using Trim Galore! (version {trim_galore_version}) powered by Cutadapt (version {cutadapt_version}). "
    description += f"Trimmed fastq file quality was evaluated with FastQC (version {fastqc_version}), and MultiQC (version {multiqc_version}) was used to generate MultiQC reports. "
    
    # Add reference description based on mode, reference source, and ERCC status
    reference_description = ""
    
    # Format organism name if provided - replace underscores with spaces and title case
    organism_name = ""
    organism_name_italics = ""  # For use in text (with underscores for Jira)
    if hasattr(args, 'organism') and args.organism:
        organism_name = args.organism.replace('_', ' ').title()
        organism_name_italics = f"_{organism_name}_"  # Surrounded by underscores for Jira italics
    
    # Get reference FASTA and GTF basenames
    ref_fasta_name = ""
    ref_gtf_name = ""
    genome_assembly = ""
    if args.reference_fasta and args.reference_fasta != "null":
        ref_fasta_name = os.path.basename(args.reference_fasta)
        
        # Handle different reference naming conventions
        if args.reference_source and "ensembl" in args.reference_source.lower():
            # For Ensembl: extract after first dot (e.g., Bacillus_subtilis.ASM904v1.dna.toplevel.fa)
            if '.' in ref_fasta_name:
                parts = ref_fasta_name.split('.')
                if len(parts) > 1:
                    genome_assembly = parts[1]
        elif args.reference_source and "ncbi" in args.reference_source.lower():
            # For NCBI: extract only the assembly name (e.g., ASM746v2 from GCF_000007465.2_ASM746v2_genomic.fna.gz)
            if '_genomic' in ref_fasta_name:
                # Get the part right before _genomic
                parts = ref_fasta_name.split('_genomic')[0].split('_')
                # The assembly name is typically the last part before _genomic (after the GCF_accession_)
                if len(parts) > 2:
                    # Skip the GCF part and accession, take the assembly name
                    genome_assembly = parts[-1]
            # Fallback to old method if _genomic is not in the name
            elif '_' in ref_fasta_name:
                parts = ref_fasta_name.split('_')
                if len(parts) > 2:  # Should have at least 3 parts
                    # The assembly name is usually the second part after the GCF_accession
                    genome_assembly = parts[1]
                    # If there are more parts before "genomic", include them in assembly name
                    for i in range(2, len(parts)):
                        if "genomic" in parts[i]:
                            break
                        genome_assembly += "_" + parts[i]
        else:
            # For other sources, try to extract anything that looks like an assembly version
            if '_' in ref_fasta_name and '.' in ref_fasta_name:
                # First try the Ensembl pattern (after first dot)
                parts = ref_fasta_name.split('.')
                if len(parts) > 1:
                    genome_assembly = parts[1]
                
                # If that didn't work, try the NCBI pattern (after first underscore)
                if not genome_assembly:
                    parts = ref_fasta_name.split('_')
                    if len(parts) > 1:
                        genome_assembly = parts[1]
    
    if args.reference_gtf and args.reference_gtf != "null":
        ref_gtf_name = os.path.basename(args.reference_gtf)
    
    # Get software versions
    star_version = software_versions.get('STAR', 'unknown')
    rsem_version = software_versions.get('RSEM', 'unknown')
    bowtie2_version = software_versions.get('Bowtie2', 'unknown')
    
    # Check if reference source contains "ensembl" or "ncbi"
    is_ensembl = False
    is_ncbi = False
    ref_source_formatted = ""
    if args.reference_source:
        if "ensembl" in args.reference_source.lower():
            is_ensembl = True
            # Format ensembl source nicely (e.g., ensembl_bacteria -> Ensembl Bacteria)
            ref_source_formatted = args.reference_source.replace('_', ' ').title()
        elif "ncbi" in args.reference_source.lower():
            is_ncbi = True
            # Format NCBI source nicely (e.g., ncbi_refseq -> NCBI RefSeq)
            ref_source_formatted = args.reference_source.replace('_', ' ').upper()
    
    # Generate reference description
    if args.mode == 'microbes':
        # For microbes, only mention Bowtie 2
        reference_description = f"{organism_name_italics} Bowtie 2 reference was built using Bowtie 2 (version {bowtie2_version})"
        
        # Add source and version information
        if is_ensembl and args.reference_version:
            reference_description += f", {ref_source_formatted} release {args.reference_version}"
        elif is_ncbi and args.reference_version:
            reference_description += f", {ref_source_formatted} version {args.reference_version}"
        
        # Add assembly information if available
        if genome_assembly:
            reference_description += f", genome assembly {genome_assembly}"
            
        reference_description += f" ({ref_fasta_name}) and the following gtf annotation file: {ref_gtf_name}. "
    else:
        # For standard mode, mention STAR and RSEM
        reference_description = f"{organism_name_italics} STAR and RSEM references were built using STAR (version {star_version}) and RSEM (version {rsem_version}), respectively"
        
        # Add source and version information
        if is_ensembl and args.reference_version:
            reference_description += f", {ref_source_formatted} release {args.reference_version}"
        elif is_ncbi and args.reference_version:
            reference_description += f", {ref_source_formatted} version {args.reference_version}"
            
        # Add assembly information if available
        if genome_assembly:
            reference_description += f", genome assembly {genome_assembly}"
        
        # Handle ERCC differently based on has_ercc
        if args.has_ercc.lower() == "true":
            reference_description += f" ({ref_fasta_name}) concatenated with ERCC92.fa from ThermoFisher (https://assets.thermofisher.com/TFS-Assets/LSG/manuals/ERCC92.zip), and the following gtf annotation file: {ref_gtf_name} concatenated with the ERCC92.gtf file from ThermoFisher (https://assets.thermofisher.com/TFS-Assets/LSG/manuals/ERCC92.zip). "
        else:
            reference_description += f" ({ref_fasta_name}) and the following gtf annotation file: {ref_gtf_name}. "
    
    # Add reference description to the protocol
    description += reference_description
    
    # Add alignment info based on mode and ERCC presence
    if args.mode == 'microbes':
        # For microbes workflow
        if args.has_ercc.lower() == "true":
            description += f"Trimmed reads were aligned to the {organism_name_italics} + ERCC Bowtie 2 reference with Bowtie 2 (version {bowtie2_version}) and alignment log files were compiled with MultiQC (version {multiqc_version}). "
        else:
            description += f"Trimmed reads were aligned to the {organism_name_italics} Bowtie 2 reference with Bowtie 2 (version {bowtie2_version}) and alignment log files were compiled with MultiQC (version {multiqc_version}). "
    else:
        # For standard workflow
        if args.has_ercc.lower() == "true":
            description += f"Trimmed reads were aligned to the {organism_name_italics} + ERCC STAR reference with STAR (version {star_version}) and alignment log files were compiled with MultiQC (version {multiqc_version}). "
        else:
            description += f"Trimmed reads were aligned to the {organism_name_italics} STAR reference with STAR (version {star_version}) and alignment log files were compiled with MultiQC (version {multiqc_version}). "
    
    # Add RSeQC assessment sentence
    infer_exp_version = software_versions.get('infer_experiment.py', 'unknown')
    gene_body_version = software_versions.get('geneBody_coverage.py', 'unknown')
    inner_dist_version = software_versions.get('inner_distance.py', 'unknown')
    read_dist_version = software_versions.get('read_distribution.py', 'unknown')
    
    description += f"Aligned reads were assessed for strandedness using RSeQC infer_experiment.py (version {infer_exp_version}), gene body coverage using RSeQC geneBody_coverage.py (version {gene_body_version}), "
    
    # Only include inner distance for paired-end data
    if args.paired_end.lower() == "true":
        description += f"inner distance was assessed using RSeQC inner_distance.py (version {inner_dist_version}), "
    
    description += f"and distribution of RNA features using RSeQC read_distribution.py (version {read_dist_version}). "
    
    # Add read quantification with appropriate strandedness setting
    rsem_version = software_versions.get('RSEM', 'unknown')
    featurecounts_version = software_versions.get('Subread', 'unknown')
    
    # Define strandedness mappings
    rsem_strandedness_map = {
        "sense": "forward",
        "antisense": "reverse",
        "unstranded": "none"
    }
    
    featurecounts_strandedness_map = {
        "unstranded": "0",
        "sense": "1",
        "antisense": "2"
    }
    
    # Get the appropriate strandedness setting based on the tool and strandedness value
    rsem_strandedness = rsem_strandedness_map.get(args.strandedness, "none")
    featurecounts_strandedness = featurecounts_strandedness_map.get(args.strandedness, "0")
    
    if args.mode == 'microbes':
        # For microbes workflow (featureCounts)
        description += f"Aligned reads from all samples were quantified using featureCounts from the Subread package (version {featurecounts_version}), "
        description += f"with isStrandSpecific set to {featurecounts_strandedness} "
        description += f"and featureCounts log files were compiled with MultiQC (version {multiqc_version}). "
        
        # Add rRNA removal sentence (for microbes workflow)
        id_type = "Ensembl IDs" if args.reference_source and "ensembl" in args.reference_source.lower() else "gene IDs"
        description += f"To create a set of rRNA-removed quantification data, {id_type} mapping to ribosomal RNA (rRNA) were removed from the featureCounts gene count data. Both the original set of featureCounts quantification data and the rRNA-removed data were subject to differential gene expression analysis as follows. "
    else:
        # For standard workflow (RSEM)
        description += f"Aligned reads from all samples were quantified using RSEM (version {rsem_version}), "
        description += f"with strandedness set to {rsem_strandedness} "
        description += f"and RSEM count log files were compiled with MultiQC (version {multiqc_version}). "
        
        # Add rRNA removal sentence for standard workflow
        id_type = "Ensembl IDs" if args.reference_source and "ensembl" in args.reference_source.lower() else "gene IDs"
        description += f"To create a set of rRNA-removed quantification data, {id_type} mapping to ribosomal RNA (rRNA) were removed from the RSEM gene count data. Both the original set of RSEM quantification data and the rRNA-removed data were subject to differential gene expression analysis as follows. "
    
    # Add runsheet generation and normalization sentence
    dp_tools_version = software_versions.get('dp_tools', 'unknown')
    r_version = software_versions.get('R', 'unknown')
    tximport_version = software_versions.get('tximport', 'unknown')
    deseq2_version = software_versions.get('DESeq2', 'unknown')
    
    if args.mode == 'microbes':
        # For microbes workflow (omit tximport)
        description += f"A runsheet was generated with dp_tools (version {dp_tools_version}) and the runsheet and quantification data were imported to R (version {r_version}) and normalized with DESeq2 (version {deseq2_version}) median of ratios method. "
    else:
        # For standard workflow (include tximport)
        if args.has_ercc.lower() == "true":
            description += f"A runsheet was generated with dp_tools (version {dp_tools_version}) and the runsheet and quantification data were imported to R (version {r_version}) with tximport (version {tximport_version}), ERCC genes were removed, and the non-ERCC genes were normalized with DESeq2 (version {deseq2_version}) median of ratios method. "
        else:
            description += f"A runsheet was generated with dp_tools (version {dp_tools_version}) and the runsheet and quantification data were imported to R (version {r_version}) with tximport (version {tximport_version}) and normalized with DESeq2 (version {deseq2_version}) median of ratios method. "
    
    # Parse runsheet for technical replicate handling
    tech_rep_sentence = ""
    if hasattr(args, 'runsheet') and args.runsheet and os.path.exists(args.runsheet):
        try:
            runsheet_df = pd.read_csv(args.runsheet)
            # Check if Source Name column exists
            if 'Source Name' in runsheet_df.columns:
                # Count occurrences of each Source Name
                source_counts = runsheet_df['Source Name'].value_counts()
                n_reps = list(source_counts.values())
                unique_n = set(n_reps)
                
                if all(x == 1 for x in n_reps):
                    # No technical replicates at all
                    tech_rep_sentence = ""
                elif len(unique_n) == 1 and list(unique_n)[0] > 1:
                    # All samples have the same number of tech reps
                    tech_rep_sentence = ("Counts from all technical replicates for each sample were summed using DESeq2's collapseReplicates function. "
                                         "These collapsed counts were then used for count normalization and differential expression analysis. ")
                elif len(unique_n) > 1 and min(unique_n) > 1:
                    # All samples have tech reps, but unequal number
                    min_reps = min(unique_n)
                    tech_rep_sentence = (f"For each sample, counts from the first {min_reps} technical replicates were summed using DESeq2's collapseReplicates function. "
                                         "These collapsed counts were then used for count normalization and differential expression analysis. ")
                else:
                    # Some samples have tech reps, some don't (mixed scenario)
                    tech_rep_sentence = ("For samples with technical replicates, only the first replicate was used for count normalization and differential expression analysis. ")
            # Fallback to old method if Source Name column doesn't exist
            elif 'Sample Name' in runsheet_df.columns:
                # Remove whitespace and NA
                sample_names = runsheet_df['Sample Name'].dropna().astype(str).tolist()
                # Remove trailing/leading whitespace
                sample_names = [s.strip() for s in sample_names]
                # Remove empty
                sample_names = [s for s in sample_names if s]
                # Find base names (remove _techrepN if present)
                base_names = [re.sub(r'_techrep\d+$', '', s) for s in sample_names]
                from collections import Counter
                base_counts = Counter(base_names)
                n_reps = list(base_counts.values())
                unique_n = set(n_reps)
                if all(x == 1 for x in n_reps):
                    # No technical replicates at all
                    tech_rep_sentence = ""
                elif len(unique_n) == 1 and list(unique_n)[0] > 1:
                    # All samples have the same number of tech reps
                    tech_rep_sentence = ("Counts from all technical replicates for each sample were summed using DESeq2's collapseReplicates function. "
                                         "These collapsed counts were then used for count normalization and differential expression analysis. ")
                elif len(unique_n) > 1 and min(unique_n) > 1:
                    # All samples have tech reps, but unequal number
                    min_reps = min(unique_n)
                    tech_rep_sentence = (f"For each sample, counts from the first {min_reps} technical replicates were summed using DESeq2's collapseReplicates function. "
                                         "These collapsed counts were then used for count normalization and differential expression analysis. ")
                else:
                    # Some samples have tech reps, some don't
                    tech_rep_sentence = ("For samples with technical replicates, only the first replicate was used for count normalization and differential expression analysis. ")
        except Exception as e:
            tech_rep_sentence = ""
    # If no runsheet, leave tech_rep_sentence as empty
    
    # Add normalization and differential expression analysis sentence
    description += "Normalized gene counts were subject to differential expression analysis. "
    
    # Add tech rep sentence
    if tech_rep_sentence:
        description += tech_rep_sentence
        
    # Add differential expression analysis sentence
    description += f"Differential expression analysis was performed in R (version {r_version}) using DESeq2 (version {deseq2_version}); all groups were compared pairwise using the Wald test and the likelihood ratio test was used to generate the F statistic p-value. "
    
    # Add gene annotations section
    # Define versions for annotation packages
    stringdb_version = "2.16.4"
    pantherdb_version = "1.0.12"
    
    # Define organism to annotation package mapping using scientific names
    organism_annotation_packages = {
        "arabidopsis_thaliana": ("org.At.tair.db", "3.19.1"),
        "caenorhabditis_elegans": ("org.Ce.eg.db", "3.19.1"),
        "drosophila_melanogaster": ("org.Dm.eg.db", "3.19.1"),
        "danio_rerio": ("org.Dr.eg.db", "3.19.1"),
        "homo_sapiens": ("org.Hs.eg.db", "3.19.1"),
        "mus_musculus": ("org.Mm.eg.db", "3.19.1"),
        "rattus_norvegicus": ("org.Rn.eg.db", "3.19.1"),
        "saccharomyces_cerevisiae": ("org.Sc.sgd.db", "3.19.1")
    }
    
    # List of organisms that use custom annotation packages via AnnotationForge
    organisms_with_custom_annotations = [
        "Bacillus subtilis subsp. subtilis str. 168",
        "Brachypodium distachyon",
        "Escherichia coli str. K-12 substr. MG1655",
        "Oryzias latipes",
        "Salmonella enterica subsp. enterica serovar Typhimurium str. LT2"
    ]
    
    # Organisms that don't use STRING database
    organisms_without_string = [
        "Pseudomonas aeruginosa UCBPP-PA14",
        "Staphylococcus aureus MRSA252"
    ]
    
    # Organisms that don't use PANTHER database
    organisms_without_panther = [
        "Caenorhabditis elegans",
        "Lactobacillus acidophilus NCFM",
        "Mycobacterium marinum M",
        "Oryza sativa Japonica",
        "Pseudomonas aeruginosa UCBPP-PA14",
        "Serratia liquefaciens ATCC 27592",
        "Staphylococcus aureus MRSA252",
        "Streptococcus mutans UA159",
        "Vibrio fischeri ES114"
    ]
    
    # Format the organism name for lookup
    organism_formatted = ""
    if hasattr(args, 'organism') and args.organism:
        organism_formatted = args.organism.replace(' ', '_').replace('-', '_').lower()
    
    # Build gene annotations sentence
    annotation_sources = []
    # Only include STRINGdb if not in organisms_without_string
    if organism_formatted not in organisms_without_string:
        annotation_sources.append(f"STRINGdb (version {stringdb_version})")
    # Only include PANTHER.db if not in organisms_without_panther
    if organism_formatted not in organisms_without_panther:
        annotation_sources.append(f"PANTHER.db (version {pantherdb_version})")
    # Custom annotation package
    custom_pkg = None
    if hasattr(args, 'organism') and args.organism and args.organism in organisms_with_custom_annotations:
        custom_pkg = "a custom annotation package generated in-house using AnnotationForge"
        annotation_sources.append(custom_pkg)
    # org.*.eg.db package
    if organism_formatted and organism_formatted in organism_annotation_packages:
        package_name, package_version = organism_annotation_packages[organism_formatted]
        annotation_sources.append(f"{package_name} (version {package_version})")
    # Build the sentence
    base_text = ("Gene annotations were assigned using the custom annotation tables generated in-house as detailed in GL-DPPD-7110-A "
                 "(https://github.com/nasa/GeneLab_Data_Processing/blob/GL_RefAnnotTable-A_1.1.0/GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A.md)")
    if annotation_sources:
        if len(annotation_sources) == 1:
            gene_annotations_text = f"{base_text}, with {annotation_sources[0]}."
        else:
            gene_annotations_text = f"{base_text}, with {', '.join(annotation_sources[:-1])}, and {annotation_sources[-1]}."
    else:
        gene_annotations_text = f"{base_text}."
    description += gene_annotations_text
    
    # Add ERCC assessment sentence if ERCC spike-ins were used
    if args.has_ercc.lower() == "true":
        description += f" Assessment of ERCC spiked genes was performed as detailed in the ERCC_analysis{args.assay_suffix}.html report."
    
    # Add paragraph break after the main description for separation from config section
    description += "\n\n"
    
    # Create configuration section
    config = "## Configuration\n\n"
    config += f"- Mode: {args.mode}\n"
    config += f"- Paired-end: {args.paired_end}\n"
    config += f"- Strandedness: {args.strandedness}\n"
    config += f"- ERCC spike-ins: {args.has_ercc}\n"
    
    # Add reference information if provided
    if hasattr(args, 'organism') and args.organism:
        config += f"- Organism: {organism_name}\n"
    if hasattr(args, 'reference_source') and args.reference_source:
        config += f"- Reference source: {args.reference_source}\n"
    if hasattr(args, 'reference_version') and args.reference_version:
        config += f"- Reference version: {args.reference_version}\n"
    
    # For reference_fasta and reference_gtf, check if they're "null" (from Nextflow)
    if args.reference_fasta and args.reference_fasta != "null":
        config += f"- Reference FASTA: {os.path.basename(args.reference_fasta)}\n"
    if args.reference_gtf and args.reference_gtf != "null":
        config += f"- Reference GTF: {os.path.basename(args.reference_gtf)}\n"
    
    config += "\n"
    
    # Create software versions section
    sw_section = "## Software Versions\n\n"
    for software, version in software_versions.items():
        sw_section += f"- {software}: {version}\n"
    
    # Combine all sections
    content = header + description + config + sw_section
    
    return content

def main():
    args = parse_args()
    
    # Read software versions from YAML file
    software_versions = read_software_versions(args.software_table)
    
    # Generate protocol content
    protocol_content = generate_protocol_content(args, software_versions)
    
    # Write to output file
    output_file = os.path.join(args.outdir, f"protocol{args.assay_suffix}.txt")
    try:
        with open(output_file, 'w') as f:
            f.write(protocol_content)
        print(f"Protocol file generated successfully: {output_file}")
    except Exception as e:
        sys.stderr.write(f"Error writing protocol file: {e}\n")
        sys.exit(1)

if __name__ == "__main__":
    main()
