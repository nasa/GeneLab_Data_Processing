#!/usr/bin/env python

import json
import os
import zipfile
import re
import sys
import shutil

def clean_path(path):
    """Extract just the filename from the full path"""
    return os.path.basename(path)

def clean_log_line(line):
    """Clean paths in a log line while preserving the log format"""
    path_pattern = r'(?<=: )(/[^\s\]]+)'
    return re.sub(path_pattern, lambda m: clean_path(m.group(1)), line)

def detect_module_from_data(data):
    """Detect which MultiQC module was used based on the data structure"""
    if 'report_data_sources' in data:
        modules = list(data['report_data_sources'].keys())
        if modules:
            return modules[0]
    return None

def clean_multiqc_data(input_dir, output_dir):
    """Clean paths in MultiQC data files"""
    # Clean multiqc_data.json
    json_path = os.path.join(input_dir, 'multiqc_data.json')
    if os.path.exists(json_path):
        with open(json_path, 'r') as f:
            data = json.load(f)
        
        module = detect_module_from_data(data)
        
        if 'report_data_sources' in data:
            for module_name in data['report_data_sources']:
                module_data = data['report_data_sources'][module_name]
                
                # Handle both formats: with all_sections and direct module structure
                if 'all_sections' in module_data:
                    for sample in module_data['all_sections']:
                        old_path = module_data['all_sections'][sample]
                        module_data['all_sections'][sample] = clean_path(old_path)
                else:
                    # Handle direct module structure (like RSeQC)
                    for section in module_data:
                        if isinstance(module_data[section], dict):
                            for sample in module_data[section]:
                                old_path = module_data[section][sample]
                                module_data[section][sample] = clean_path(old_path)
        
        if 'config_output_dir' in data:
            data['config_output_dir'] = clean_path(data['config_output_dir'])
        if 'config_analysis_dir_abs' in data:
            data['config_analysis_dir_abs'] = [clean_path(p) for p in data['config_analysis_dir_abs']]
        if 'config_script_path' in data:
            data['config_script_path'] = clean_path(data['config_script_path'])
        
        if module:
            data['config_module'] = module
        
        with open(os.path.join(output_dir, 'multiqc_data.json'), 'w') as f:
            json.dump(data, f, indent=4)
    
    # Clean multiqc_sources.txt
    sources_path = os.path.join(input_dir, 'multiqc_sources.txt')
    if os.path.exists(sources_path):
        with open(sources_path, 'r') as f:
            lines = f.readlines()
        
        cleaned_lines = []
        for line in lines:
            if line.strip():
                parts = line.split('\t')
                if len(parts) >= 4:
                    parts[3] = clean_path(parts[3])
                cleaned_lines.append('\t'.join(parts))
        
        with open(os.path.join(output_dir, 'multiqc_sources.txt'), 'w') as f:
            f.writelines(cleaned_lines)

    # Clean multiqc.log
    log_path = os.path.join(input_dir, 'multiqc.log')
    if os.path.exists(log_path):
        with open(log_path, 'r') as f:
            lines = f.readlines()
        
        cleaned_lines = [clean_log_line(line) for line in lines]
        
        with open(os.path.join(output_dir, 'multiqc.log'), 'w') as f:
            f.writelines(cleaned_lines)

def process_input(input_path, output_dir):
    """Process the MultiQC data directory"""
    temp_dir = "temp_multiqc"
    os.makedirs(temp_dir, exist_ok=True)
    
    # Copy input directory contents directly to temp
    for item in os.listdir(input_path):
        s = os.path.join(input_path, item)
        d = os.path.join(temp_dir, item)
        if os.path.isfile(s):
            shutil.copy2(s, d)
        elif os.path.isdir(s):
            shutil.copytree(s, d)
    
    # Clean the data
    clean_multiqc_data(temp_dir, temp_dir)
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Create new zip file
    output_name = os.path.basename(input_path)
    if not output_name.endswith('.zip'):
        output_name += '.zip'
    output_zip = os.path.join(output_dir, output_name)
    
    # Get the directory name (remove .zip if present)
    dir_name = output_name[:-4] if output_name.endswith('.zip') else output_name
    
    with zipfile.ZipFile(output_zip, 'w', zipfile.ZIP_DEFLATED) as zip_ref:
        for root, _, files in os.walk(temp_dir):
            for file in files:
                file_path = os.path.join(root, file)
                arcname = os.path.join(dir_name, os.path.relpath(file_path, temp_dir))
                zip_ref.write(file_path, arcname)
    
    # Clean up
    shutil.rmtree(temp_dir)
    return output_zip

def main():
    if len(sys.argv) != 3:
        print("Usage: clean_multiqc_paths.py <multiqc_data_dir> <output_dir>")
        sys.exit(1)
    
    input_path = sys.argv[1]
    output_dir = sys.argv[2]
    
    process_input(input_path, output_dir)

if __name__ == "__main__":
    main()
