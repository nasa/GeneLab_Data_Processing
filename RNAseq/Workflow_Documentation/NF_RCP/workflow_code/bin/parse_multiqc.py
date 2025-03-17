#!/usr/bin/env python
from statistics import mean, median
from io import TextIOWrapper
from zipfile import ZipFile
from glob import glob
import pandas as pd
import numpy as np
import requests
import argparse
import json
import csv
import sys
import re
import os


def main(osd_num, paired_end, assay_suffix, mode):

    osd_num = osd_num.split('-')[1]

    # Create the multiqc_data list with conditionally selected parsers based on mode
    multiqc_data = [
        parse_isa(),
        parse_fastqc('raw', assay_suffix),
        parse_fastqc('trimmed', assay_suffix)
    ]
    
    # Add the appropriate parsers based on mode
    if mode == 'microbes':
        multiqc_data.append(parse_bowtie2(assay_suffix))
        #multiqc_data.append(parse_featurecounts(assay_suffix))
        multiqc_data.append(parse_featurecounts(assay_suffix))
    else:
        multiqc_data.append(parse_star(assay_suffix))
        multiqc_data.append(parse_rsem(assay_suffix))
    
    # Add the remaining common parsers
    multiqc_data.extend([
        parse_genebody_cov(assay_suffix),
        parse_infer_exp(assay_suffix),
        parse_read_dist(assay_suffix),
        get_genecount(assay_suffix, mode)
    ])

    if paired_end:
        multiqc_data.append(parse_inner_dist(assay_suffix))

    samples = set([s for ss in multiqc_data for s in ss])

    metadata = get_metadata(osd_num)

    fieldnames = [
        'osd_num', 'sample', 'organism', 'tissue', 'sequencing_instrument', 'library_selection', 'library_layout', 'strandedness', 'read_depth', 'read_length', 'rrna_contamination', 'rin', 'organism_part', 'cell_line', 'cell_type', 'secondary_organism', 'strain', 'animal_source', 'seed_source', 'source_accession', 'mix',

        # Gene count
        'gene_detected_gt10', 'gene_total', 'gene_detected_gt10_pct',

        # FastQC (raw)
        'raw_total_sequences_f', 'raw_avg_sequence_length_f', 'raw_median_sequence_length_f', 'raw_quality_score_mean_f', 'raw_quality_score_median_f', 'raw_percent_duplicates_f',
        'raw_percent_gc_f', 'raw_gc_min_1pct_f', 'raw_gc_max_1pct_f', 'raw_gc_auc_25pct_f', 'raw_gc_auc_50pct_f', 'raw_gc_auc_75pct_f', 'raw_n_content_sum_f',
        'raw_total_sequences_r', 'raw_avg_sequence_length_r', 'raw_median_sequence_length_r', 'raw_quality_score_mean_r', 'raw_quality_score_median_r', 'raw_percent_duplicates_r',
        'raw_percent_gc_r', 'raw_gc_min_1pct_r', 'raw_gc_max_1pct_r', 'raw_gc_auc_25pct_r', 'raw_gc_auc_50pct_r', 'raw_gc_auc_75pct_r', 'raw_n_content_sum_r',

        # FastQC (trimmed)
        'trimmed_total_sequences_f', 'trimmed_avg_sequence_length_f', 'trimmed_median_sequence_length_f', 'trimmed_quality_score_mean_f', 'trimmed_quality_score_median_f', 'trimmed_percent_duplicates_f',
        'trimmed_percent_gc_f', 'trimmed_gc_min_1pct_f', 'trimmed_gc_max_1pct_f', 'trimmed_gc_auc_25pct_f', 'trimmed_gc_auc_50pct_f', 'trimmed_gc_auc_75pct_f', 'trimmed_n_content_sum_f',
        'trimmed_total_sequences_r', 'trimmed_avg_sequence_length_r', 'trimmed_median_sequence_length_r', 'trimmed_quality_score_mean_r', 'trimmed_quality_score_median_r', 'trimmed_percent_duplicates_r',
        'trimmed_percent_gc_r', 'trimmed_gc_min_1pct_r', 'trimmed_gc_max_1pct_r', 'trimmed_gc_auc_25pct_r', 'trimmed_gc_auc_50pct_r', 'trimmed_gc_auc_75pct_r', 'trimmed_n_content_sum_r',

        # STAR
        'uniquely_mapped_percent', 'multimapped_percent', 'multimapped_toomany_percent', 'unmapped_tooshort_percent', 'unmapped_other_percent',

        # Bowtie2
        'total_reads', 'overall_alignment_rate', 'aligned_none', 'aligned_one', 'aligned_multi',

        # FeatureCounts
        'total_count', 'num_assigned', 'pct_assigned', 'num_unassigned_nofeatures', 'num_unassigned_ambiguity', 'pct_unassigned_nofeatures', 'pct_unassigned_ambiguity',
        
        # RSeQC
        'mean_genebody_cov_5_20', 'mean_genebody_cov_40_60', 'mean_genebody_cov_80_95', 'ratio_genebody_cov_3_to_5',
        'pct_sense', 'pct_antisense', 'pct_undetermined',
        'peak_inner_dist', 'peak_inner_dist_pct_reads',
        'cds_exons_pct', '5_utr_exons_pct', '3_utr_exons_pct', 'introns_pct', 'tss_up_1kb_pct', 'tss_up_1kb_5kb_pct', 'tss_up_5kb_10kb_pct', 'tes_down_1kb_pct', 'tss_down_1kb_5kb_pct', 'tss_down_5kb_10kb_pct', 'other_intergenic_pct',

        # RSEM
        'num_uniquely_aligned', 'pct_uniquely_aligned', 'pct_multi_aligned', 'pct_filtered', 'pct_unalignable'
    ]

    output_filename = f'qc_metrics{assay_suffix}.csv'
    with open(output_filename, mode='w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for sample in samples:
            fields = {k: v for d in [i[sample] for i in multiqc_data if sample in i] for k, v in d.items()}

            writer.writerow({'osd_num': 'OSD-' + osd_num, 'sample': sample, **metadata, **fields})


def get_metadata(osd_num):
    r = requests.get('https://osdr.nasa.gov/osdr/data/osd/meta/' + osd_num)
    data = {}

    # Source accession
    comments = r.json()['study']['OSD-' + osd_num]['studies'][0]['comments']
    source_accession = [c for c in comments if c['name'] == 'Data Source Accession'][0]

    if source_accession and source_accession['value']:
        data['source_accession'] = source_accession['value']

    return data


def parse_isa():
    try:
        with ZipFile(glob('*ISA.zip')[0]) as z:
            with z.open([i for i in z.namelist() if i.startswith('s_')][0]) as f:
                s = list(csv.DictReader(TextIOWrapper(f, 'utf-8'), delimiter='\t'))

            with z.open([i for i in z.namelist() if i.startswith('a_') and 'rna-seq' in i.lower().replace('_', '-') and 'mirna' not in i.lower()][0]) as f:
                a = list(csv.DictReader(TextIOWrapper(f, 'utf-8'), delimiter='\t'))

        data = {}

        sample_fields = {
            'characteristics[organism]': 'organism',
            'characteristics[material type]': 'tissue',
            'material type': 'tissue',
            'characteristics[organism part]': 'organism_part',
            'factor value[organism part]': 'organism_part',
            'characteristics[cell line]': 'cell_line',
            'characteristics[cell line,http://purl.obolibrary.org/obo/clo_0000031,obi]': 'cell_line',
            'characteristics[cell type]': 'cell_type',
            'characteristics [cell type]': 'cell_type',
            'characteristics[secondary organism]': 'secondary_organism',
            'characteristics[strain]': 'strain',
            'characteristics[animal source]': 'animal_source',
            'characteristics[seed source]': 'seed_source'
        }
        assay_fields = {
            'parameter value[sequencing instrument]': 'sequencing_instrument',
            'parameter value[library selection]': 'library_selection',
            'parameter value[library layout]': 'library_layout',
            'parameter value[strandedness]': 'strandedness',
            'parameter value[stranded]': 'strandedness',
            'parameter value[read depth]': 'read_depth',
            'parameter value[read length]': 'read_length',
            'parameter value[rrna contamination]': 'rrna_contamination',
            'parameter value[rrna estimation]': 'rrna_contamination',
            'parameter value[qa score]': 'rin',
            'parameter value[spike-in mix number]': 'mix'
        }

        for row in a:
            data[row['Sample Name'].strip()] = {assay_fields[k.lower()]:v.strip() for k, v in row.items() if k.lower() in assay_fields}

        for row in s:
            if row['Sample Name'].strip() not in data:  # samples could be in other assays
                continue

            sample_data = {sample_fields[k.lower()]:v.strip() for k, v in row.items() if k.lower() in sample_fields}
            data[row['Sample Name'].strip()] = {**data[row['Sample Name'].strip()], **sample_data}

        return data
    except (FileNotFoundError, KeyError, IndexError, ValueError, Exception) as e:
        print(f"WARNING: Error processing ISA data: {str(e)}")
        return {}


def parse_fastqc(prefix, assay_suffix):
    try:
        with open(f'{prefix}_multiqc{assay_suffix}_data/multiqc_data.json') as f:
            j = json.loads(f.read())

        # Group the samples by base name for paired end data
        sample_groups = {}
        for sample in j['report_general_stats_data'][-1].keys():
            # Handle various naming patterns
            if ' Read 1' in sample:
                base_name = sample.replace(' Read 1', '')
                if base_name not in sample_groups:
                    sample_groups[base_name] = {'f': None, 'r': None}
                sample_groups[base_name]['f'] = sample
            elif ' Read 2' in sample:
                base_name = sample.replace(' Read 2', '')
                if base_name not in sample_groups:
                    sample_groups[base_name] = {'f': None, 'r': None}
                sample_groups[base_name]['r'] = sample
            elif '_R1' in sample:
                base_name = sample.replace('_R1', '')
                if base_name not in sample_groups:
                    sample_groups[base_name] = {'f': None, 'r': None}
                sample_groups[base_name]['f'] = sample
            elif '_R2' in sample:
                base_name = sample.replace('_R2', '')
                if base_name not in sample_groups:
                    sample_groups[base_name] = {'f': None, 'r': None}
                sample_groups[base_name]['r'] = sample
            else:
                # For single-end or non-paired samples
                base_name = sample
                if base_name not in sample_groups:
                    sample_groups[base_name] = {'f': None, 'r': None}
                sample_groups[base_name]['f'] = sample

        data = {}
        # Process each sample group
        for base_name, reads in sample_groups.items():
            data[base_name] = {}
            
            # Process forward read
            if reads['f']:
                for k, v in j['report_general_stats_data'][-1][reads['f']].items():
                    if k != 'percent_fails':
                        data[base_name][prefix + '_' + k + '_f'] = v
                        
            # Process reverse read
            if reads['r']:
                for k, v in j['report_general_stats_data'][-1][reads['r']].items():
                    if k != 'percent_fails':
                        data[base_name][prefix + '_' + k + '_r'] = v

        # Process other stats sections (quality, GC, etc)
        for section, suffix in [
            ('fastqc_per_base_sequence_quality_plot', 'quality_score'),
            ('fastqc_per_sequence_gc_content_plot', 'gc'),
            ('fastqc_per_base_n_content_plot', 'n_content')
        ]:
            if section in j['report_plot_data']:
                for data_item in j['report_plot_data'][section]['datasets'][0]['lines']:
                    sample = data_item['name']
                    
                    # Determine if it's forward or reverse read
                    read_suffix = '_f'  # Default to forward
                    base_name = sample
                    
                    if ' Read 2' in sample:
                        read_suffix = '_r'
                        base_name = sample.replace(' Read 2', '')
                    elif ' Read 1' in sample:
                        base_name = sample.replace(' Read 1', '')
                    elif '_R2' in sample:
                        read_suffix = '_r'
                        base_name = sample.replace('_R2', '')
                    elif '_R1' in sample:
                        base_name = sample.replace('_R1', '')
                        
                    # Skip if we don't have this sample
                    if base_name not in data:
                        continue
                        
                    # Process based on the section
                    if suffix == 'quality_score':
                        data[base_name][prefix + '_quality_score_mean' + read_suffix] = mean([i[1] for i in data_item['pairs']])
                        data[base_name][prefix + '_quality_score_median' + read_suffix] = median([i[1] for i in data_item['pairs']])
                    elif suffix == 'gc':
                        gc_data_1pct = [i[0] for i in data_item['pairs'] if i[1] >= 1]
                        if gc_data_1pct:
                            data[base_name][prefix + '_gc_min_1pct' + read_suffix] = gc_data_1pct[0]
                            data[base_name][prefix + '_gc_max_1pct' + read_suffix] = gc_data_1pct[-1]
                            
                            gc_data_cum = list(np.cumsum([i[1] for i in data_item['pairs']]))
                            data[base_name][prefix + '_gc_auc_25pct' + read_suffix] = list(i >= 25 for i in gc_data_cum).index(True)
                            data[base_name][prefix + '_gc_auc_50pct' + read_suffix] = list(i >= 50 for i in gc_data_cum).index(True)
                            data[base_name][prefix + '_gc_auc_75pct' + read_suffix] = list(i >= 75 for i in gc_data_cum).index(True)
                    elif suffix == 'n_content':
                        data[base_name][prefix + '_n_content_sum' + read_suffix] = sum([i[1] for i in data_item['pairs']])

        return data
    except (FileNotFoundError, KeyError, IndexError, json.JSONDecodeError, ValueError):
        print(f"WARNING: Could not process {prefix} FastQC data")
        return {}


def parse_star(assay_suffix):
    try:
        with open(f'align_multiqc{assay_suffix}_data/multiqc_data.json') as f:
            j = json.loads(f.read())

        data = {}

        align_fields = ['uniquely_mapped_percent', 'multimapped_percent', 'multimapped_toomany_percent', 'unmapped_tooshort_percent', 'unmapped_other_percent']

        for sample in j['report_general_stats_data'][0].keys():
            data[sample] = {k:v for k, v in j['report_general_stats_data'][0][sample].items() if k in align_fields}

        return data
    except (FileNotFoundError, KeyError, IndexError):
        return {}


def parse_bowtie2(assay_suffix):
    """
    Parse Bowtie2 alignment statistics from MultiQC report
    """
    try:
        # Open the data file
        with open(f'align_multiqc{assay_suffix}_data/multiqc_data.json') as f:
            data = json.load(f)
        
        # Create a dictionary to hold sample data
        samples_data = {}
        
        # Extract the stats for each sample
        for section in data['report_general_stats_data']:
            for sample, stats in section.items():
                # Clean up sample name to remove read identifiers
                base_sample = re.sub(r'_R[12]$', '', sample)
                
                if base_sample not in samples_data:
                    samples_data[base_sample] = {}
                
                # Report these 5 key metrics regardless of whether data is paired-end or single-end
                if 'total_reads' in stats:
                    samples_data[base_sample]['total_reads'] = stats['total_reads']
                
                if 'overall_alignment_rate' in stats:
                    samples_data[base_sample]['overall_alignment_rate'] = stats['overall_alignment_rate']
                
                # Handle aligned none, one, and multi for both paired and unpaired data
                # For paired data
                if 'paired_aligned_none' in stats:
                    samples_data[base_sample]['aligned_none'] = stats['paired_aligned_none']
                    samples_data[base_sample]['aligned_one'] = stats['paired_aligned_one']
                    samples_data[base_sample]['aligned_multi'] = stats['paired_aligned_multi']
                
                # For unpaired/single-end data
                elif 'unpaired_aligned_none' in stats:
                    samples_data[base_sample]['aligned_none'] = stats['unpaired_aligned_none']
                    samples_data[base_sample]['aligned_one'] = stats['unpaired_aligned_one']
                    samples_data[base_sample]['aligned_multi'] = stats['unpaired_aligned_multi']
        
        return samples_data
    
    except FileNotFoundError:
        print(f"WARNING: Could not find Bowtie2 MultiQC data file")
        return {}
    except Exception as e:
        print(f"ERROR parsing Bowtie2 data: {str(e)}")
        return {}


def parse_genebody_cov(assay_suffix):
    try:
        with open(f'geneBody_cov_multiqc{assay_suffix}_data/multiqc_data.json') as f:
            j = json.loads(f.read())

        data = {}

        for cov_data in j['report_plot_data']['rseqc_gene_body_coverage_plot']['datasets'][0]['lines']:
            sample = cov_data['name']

            data[sample] = {
                'mean_genebody_cov_5_20': mean([i[1] for i in cov_data['pairs'] if 5 <= i[0] <= 20]),
                'mean_genebody_cov_40_60': mean([i[1] for i in cov_data['pairs'] if 40 <= i[0] <= 60]),
                'mean_genebody_cov_80_95': mean([i[1] for i in cov_data['pairs'] if 80 <= i[0] <= 95])
            }

            # Handle division by zero
            if data[sample]['mean_genebody_cov_5_20'] == 0:
                data[sample]['ratio_genebody_cov_3_to_5'] = None
            else:
                data[sample]['ratio_genebody_cov_3_to_5'] = data[sample]['mean_genebody_cov_80_95'] / data[sample]['mean_genebody_cov_5_20']

        return data
    except (FileNotFoundError, KeyError, IndexError):
        return {}


def parse_infer_exp(assay_suffix):
    try:
        with open(f'infer_exp_multiqc{assay_suffix}_data/multiqc_data.json') as f:
            j = json.loads(f.read())

        key_dict = {
            'se_sense': 'pct_sense',
            'se_antisense': 'pct_antisense',
            'pe_sense': 'pct_sense',
            'pe_antisense': 'pct_antisense',
            'failed': 'pct_undetermined'
        }

        data = {s.replace('_infer_expt', ''):{key_dict[k]:v * 100 for k, v in d.items() if k in key_dict} for s, d in j['report_saved_raw_data']['multiqc_rseqc_infer_experiment'].items()}

        return data
    except (FileNotFoundError, KeyError, IndexError, json.JSONDecodeError):
        print(f"WARNING: Could not process infer experiment data")
        return {}


def parse_inner_dist(assay_suffix):
    try:
        with open(f'inner_dist_multiqc{assay_suffix}_data/multiqc_data.json') as f:
            j = json.loads(f.read())

        data = {}

        for dist_data in j['report_plot_data']['rseqc_inner_distance_plot']['datasets'][1]['lines']:
            sample = dist_data['name']

            max_dist = sorted(dist_data['pairs'], key=lambda i: i[1], reverse=True)[0]
            data[sample] = {
                'peak_inner_dist': max_dist[0],
                'peak_inner_dist_pct_reads': max_dist[1]
            }

        return data
    except (FileNotFoundError, KeyError, IndexError, json.JSONDecodeError):
        print(f"WARNING: Could not process inner distance data")
        return {}


def parse_read_dist(assay_suffix):
    try:
        with open(f'read_dist_multiqc{assay_suffix}_data/multiqc_data.json') as f:
            j = json.loads(f.read())

        data = {}

        read_fields = ['cds_exons_tag_pct', '5_utr_exons_tag_pct', '3_utr_exons_tag_pct', 'introns_tag_pct', 'tss_up_1kb_tag_pct', 'tes_down_1kb_tag_pct', 'other_intergenic_tag_pct']

        for sample, read_data in j['report_saved_raw_data']['multiqc_rseqc_read_distribution'].items():
            data[sample] = {k:v for k, v in read_data.items() if k in read_fields}

            data[sample]['tss_up_1kb_5kb_pct'] = read_data['tss_up_5kb_tag_pct'] - read_data['tss_up_1kb_tag_pct']
            data[sample]['tss_up_5kb_10kb_pct'] = read_data['tss_up_10kb_tag_pct'] - read_data['tss_up_5kb_tag_pct']
            data[sample]['tss_down_1kb_5kb_pct'] = read_data['tes_down_5kb_tag_pct'] - read_data['tes_down_1kb_tag_pct']
            data[sample]['tss_down_5kb_10kb_pct'] = read_data['tes_down_10kb_tag_pct'] - read_data['tes_down_5kb_tag_pct']

        data = {s.replace('_read_dist', ''):{k.replace('_tag', ''):v for k, v in d.items()} for s, d in data.items()}

        return data
    except (FileNotFoundError, KeyError, IndexError, json.JSONDecodeError):
        print(f"WARNING: Could not process read distribution data")
        return {}


def parse_rsem(assay_suffix):
    try:
        with open(f'RSEM_count_multiqc{assay_suffix}_data/multiqc_data.json') as f:
            j = json.loads(f.read())

        data = {}

        for sample, count_data in j['report_saved_raw_data']['multiqc_rsem'].items():
            total_reads = count_data['Unique'] + count_data['Multi'] + count_data['Filtered'] + count_data['Unalignable']

            data[sample] = {
                'num_uniquely_aligned': count_data['Unique'],
                'pct_uniquely_aligned': count_data['Unique'] / total_reads * 100,
                'pct_multi_aligned': count_data['Multi'] / total_reads * 100,
                'pct_filtered': count_data['Filtered'] / total_reads * 100,
                'pct_unalignable': count_data['Unalignable'] / total_reads * 100
            }

        return data
    except (FileNotFoundError, KeyError, IndexError, json.JSONDecodeError, ZeroDivisionError):
        print(f"WARNING: Could not process RSEM data")
        return {}


def parse_featurecounts(assay_suffix):
    try:
        with open(f'featureCounts_multiqc{assay_suffix}_data/multiqc_data.json') as f:
            j = json.loads(f.read())

        data = {}

        for sample, count_data in j['report_saved_raw_data']['multiqc_featurecounts'].items():
            data[sample] = {
                'total_count': count_data['Total'],
                'num_assigned': count_data['Assigned'],
                'pct_assigned': count_data['percent_assigned'],
                'num_unassigned_nofeatures': count_data['Unassigned_NoFeatures'],
                'num_unassigned_ambiguity': count_data['Unassigned_Ambiguity'],
                'pct_unassigned_nofeatures': count_data['Unassigned_NoFeatures'] / count_data['Total'] * 100 if count_data['Total'] > 0 else 0,
                'pct_unassigned_ambiguity': count_data['Unassigned_Ambiguity'] / count_data['Total'] * 100 if count_data['Total'] > 0 else 0
            }

        return data
    except (FileNotFoundError, KeyError, IndexError, json.JSONDecodeError):
        print(f"WARNING: Could not process FeatureCounts data")
        return {}


def get_genecount(assay_suffix, mode):
    try:
        # For microbes mode, prioritize FeatureCounts, otherwise prioritize RSEM
        if mode == 'microbes':
            if os.path.exists(f'FeatureCounts_Unnormalized_Counts{assay_suffix}.csv'):
                count_file = f'FeatureCounts_Unnormalized_Counts{assay_suffix}.csv'
            elif os.path.exists(f'RSEM_Unnormalized_Counts{assay_suffix}.csv'):
                count_file = f'RSEM_Unnormalized_Counts{assay_suffix}.csv'
            else:
                print(f"WARNING: Could not find any count files for gene count statistics")
                return {}  # No counts file found
        else:
            # Default behavior - RSEM first, then FeatureCounts
            if os.path.exists(f'RSEM_Unnormalized_Counts{assay_suffix}.csv'):
                count_file = f'RSEM_Unnormalized_Counts{assay_suffix}.csv'
            elif os.path.exists(f'FeatureCounts_Unnormalized_Counts{assay_suffix}.csv'):
                count_file = f'FeatureCounts_Unnormalized_Counts{assay_suffix}.csv'
            else:
                print(f"WARNING: Could not find any count files for gene count statistics")
                return {}  # No counts file found
            
        df = pd.read_csv(count_file, index_col=0)
        df = df[~df.index.str.contains('^ERCC-')]

        data = (df > 10).sum().to_frame(name='gene_detected_gt10')
        data['gene_total'] = df.shape[0]
        data['gene_detected_gt10_pct'] = data.gene_detected_gt10 / data.gene_total * 100

        return data.to_dict(orient='index')
    except (FileNotFoundError, pd.errors.EmptyDataError, ValueError, Exception) as e:
        print(f"WARNING: Error processing gene count data: {str(e)}")
        return {}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--osd-num')
    parser.add_argument('--paired', action='store_true')
    parser.add_argument('--assay_suffix', default='_GLbulkRNAseq')
    parser.add_argument('--mode', default='default')
    args = parser.parse_args()

    main(args.osd_num, args.paired, args.assay_suffix, args.mode)
