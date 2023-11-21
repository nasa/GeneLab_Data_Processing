import os
import shutil
import yaml

def copy_file(src, dest):
    try:
        shutil.copy(src, dest)
    except Exception as e:
        print(f"Error copying {src} to {dest}: {e}")

def copy_directory(src, dest):
    try:
        shutil.copytree(src, dest, dirs_exist_ok=True)
    except Exception as e:
        print(f"Error copying {src} directory to {dest}: {e}")

def main(config, sample_IDs_file):
    info_out_dir = config["info_out_dir"]
    output_prefix = config.get("output_prefix", "")  # Get the output_prefix, default to empty string if not found
    os.makedirs(info_out_dir, exist_ok=True)
    os.makedirs(os.path.join(info_out_dir, "benchmarks"), exist_ok=True)

    # Files to copy
    files_to_copy = [
        ("config.yaml", os.path.join(info_out_dir, "config.yaml")),
        (sample_IDs_file, os.path.join(info_out_dir, os.path.basename(sample_IDs_file))),
        (config["runsheet"], os.path.join(info_out_dir, os.path.basename(config["runsheet"]))),
        ("R-processing.log", os.path.join(info_out_dir, "R-processing.log")),
        ("R-visualizations.log", os.path.join(info_out_dir, "R-visualizations.log")),
        ("all-benchmarks.tsv", os.path.join(info_out_dir,"all-benchmarks.tsv")),
        ("Snakefile", os.path.join(info_out_dir, "Snakefile"))
    ]

    # Optional ISA archive
    if config.get("isa_archive") and os.path.isfile(config["isa_archive"]):
        files_to_copy.append((config["isa_archive"], os.path.join(info_out_dir, os.path.basename(config["isa_archive"]))))

    # Directories to copy
    directories_to_copy = [
        ("benchmarks", os.path.join(info_out_dir, "benchmarks")),
        ("envs", os.path.join(info_out_dir, "envs")),
        ("scripts", os.path.join(info_out_dir, "scripts")),
        ("config", os.path.join(info_out_dir, "config"))
    ]

    # Copy directories
    for src, dest in directories_to_copy:
        copy_directory(src, dest)

    # Copy files
    for src, dest in files_to_copy:
        copy_file(src, dest)


if __name__ == "__main__":
    with open('config.yaml') as f:
        config = yaml.safe_load(f)
    sample_IDs_file = config['sample_info_file']

    main(config, sample_IDs_file)