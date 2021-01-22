__author__ = 'Jonathan Rubin'

import os, sys, subprocess, config, metadata_process,shlex
from difflib import SequenceMatcher

#Finds similarity between two strings a, b. Not currently used but here in case it's needed (JDR 7/5/18)
def similar(a, b):
    return SequenceMatcher(None, a, b).ratio()

#Function that simply unzips/decompresses an archive based on its extension. If a 'tar' or 'zip' archive,
#also remove it to save space.
def unzip(out_file_path,rawdata_out):
    extension=out_file_path.split('.')[-1]
    if 'zip' in extension:
        unzip_command = ["unzip", "-o", "-qq",out_file_path,"-d",rawdata_out]
        subprocess.call(unzip_command)
        remove_command = ["rm", out_file_path]
        subprocess.call(remove_command)
    if 'gz' in extension:
        gunzip_command = ["gunzip", "-f",out_file_path]
        subprocess.call(gunzip_command)
    if 'tar' in extension:
        untar_command = ["tar", "-xf", out_file_path, "-C", rawdata_out]
        subprocess.call(untar_command)
        remove_command = ["rm", out_file_path]
        subprocess.call(remove_command)


#Copies most recent metadata and all raw data and then unzips, and removes compressed files
def copy(rawdata_directory):
    #Find name of GLDS number
    GLDS = os.path.basename(os.path.dirname(rawdata_directory))
    rawdata_out = os.path.join(config.outdir,GLDS,config.microarray_out)

    #Make appropriate output directory
    if not os.path.exists(rawdata_out):
        os.makedirs(rawdata_out)


    #Now search the microarray folder for raw data files (this part could be done in a smarter way...)
    for file1 in os.listdir(rawdata_directory):
        out_file_path = os.path.join(rawdata_out,file1)

        #Command for copying the raw files to desired output
        cp_command = ["cp", os.path.join(rawdata_directory,file1),out_file_path]

        #md5sum command to check original files
        config.get_md5sum(cp_command[1],'original',action='copy')

        #Then execute the copy command to copy raw files to output directory
        subprocess.call(cp_command)

        #md5sum command to check copied files
        config.get_md5sum(cp_command[2],'new')

        #Once copied, unzip/untar/gunzip compressed directories (if there are any)
        unzip(out_file_path,rawdata_out)


    #Sometimes compressed files spit out more compressed files so loop through the files once again to catch those and uncompress them
    for file2 in os.listdir(rawdata_out):
        out_file_path = os.path.join(rawdata_out,file2)
        unzip(out_file_path,rawdata_out)

    #Sometimes compressed files spit out more compressed files so loop through the files once again to catch those and uncompress them
    for file2 in os.listdir(rawdata_out):
        out_file_path = os.path.join(rawdata_out,file2)
        unzip(out_file_path,rawdata_out)

#This function renames all raw data files in a specified GLDS_path and copies them into a 'raw_files' directory
def rename(GLDS_path):
    #First get all the proper paths according to specifications
    metadata_out = os.path.join(GLDS_path,'metadata')
    rawdata_out = os.path.join(GLDS_path,config.microarray_out)
    GLDS = os.path.basename(GLDS_path)
    assay_dict = metadata_process.read_assay(metadata_out)
    final_rawdata_out = os.path.join(rawdata_out,'raw_files')
    extension = ''

    #Make the 'raw_files' directory if it doesn't already exist
    if not os.path.isdir(final_rawdata_out):
        os.makedirs(final_rawdata_out)

    #Loop through the raw data files 
    for filename in os.listdir(rawdata_out):
        if not os.path.isdir(os.path.join(rawdata_out,filename)):

            #Boolean to detect whether the first column corresponds well to the filenames
            sample_in_first_column = False
            sample_in_other_column = False
            extension = filename.split('.')[-1]

            #If the filename is an annotation type, don't include 'raw' in filename
            if '.adf.' in filename:
                new_filename = filename.replace('_','-').replace('(','-').replace(')','-').replace(' ','-').replace(GLDS,'').replace('microarray','').replace('--','-').replace('.adf','-adf').strip('-').split('.')[0]
                move_command = ["mv", "'"+os.path.join(rawdata_out,filename)+"'", os.path.join(final_rawdata_out,GLDS+'_'+new_filename+'_microarray_annotation.adf.'+extension)]
                new_md5sum_file = os.path.join(final_rawdata_out,GLDS+'_'+new_filename+'_microarray_annotation.adf.'+extension)
            elif 'GPL' in filename:
                new_filename = filename.replace('_','-').replace('(','-').replace(')','-').replace(' ','-').replace(GLDS,'').replace('microarray','').replace('--','-').replace('.adf','-adf').strip('-').split('.')[0]
                move_command = ["mv", "'"+os.path.join(rawdata_out,filename)+"'", os.path.join(final_rawdata_out,GLDS+'_'+new_filename+'_microarray_annotation.'+extension)]
                new_md5sum_file = os.path.join(final_rawdata_out,GLDS+'_'+new_filename+'_microarray_annotation.'+extension)
                config.GPL=True
            else:
                for key in assay_dict:
                    #If the first column does correspond well, then assume the first column is the sample name and rename accordingly
                    if key in filename:
                        sample_in_first_column = True
                        sample_name = key.replace(' ','-').replace('_','-').replace('(','-').replace(')','-').replace(GLDS,'').replace('microarray','').replace('--','-').strip('-')
                        move_command = ["mv", "'"+os.path.join(rawdata_out,filename)+"'", os.path.join(final_rawdata_out,GLDS+'_'+sample_name+'_microarray_raw.'+extension)]
                        new_md5sum_file = os.path.join(final_rawdata_out,GLDS+'_'+sample_name+'_microarray_raw.'+extension)

                #If the first column does not correspond to any filenames, look in other columns. Still rename as first column without special characters
                if not sample_in_first_column:
                    for key in assay_dict:
                        for item in assay_dict[key]:
                            if item == filename and item != '':
                                sample_in_other_column = True
                                new_filename = key.replace('_','-').replace('(','-').replace(')','-').replace(' ','-').replace(GLDS,'').replace('microarray','').replace('--','-').strip('-')
                                move_command = ["mv","'"+os.path.join(rawdata_out,filename)+"'", os.path.join(final_rawdata_out,GLDS+'_'+new_filename+'_microarray_raw.'+extension)]
                                new_md5sum_file = os.path.join(final_rawdata_out,GLDS+'_'+new_filename+'_microarray_raw.'+extension)

                #If the filename isn't in the metadata, just remove special characters and append appropriate information. Also, don't consider it a 'raw' file
                if not sample_in_first_column and not sample_in_other_column:
                    new_filename = filename.replace('_','-').replace('(','-').replace(')','-').replace(' ','-').replace(GLDS,'').replace('microarray','').replace('--','-').strip('-').split('.')[0]
                    move_command = ["mv","'"+os.path.join(rawdata_out,filename)+"'", os.path.join(final_rawdata_out,GLDS+'_'+new_filename+'_microarray_other.'+extension)]
                    new_md5sum_file = os.path.join(final_rawdata_out,GLDS+'_'+new_filename+'_microarray_other.'+extension)

            #Execute the command if the file was in metadata - catch whether the file already exists and don't output an error
            try:
                config.get_md5sum(move_command[1].strip("'"),'original',action='rename')
                with open(os.devnull,'w') as FNULL:
                    subprocess.check_call(' '.join(move_command),shell=True,stdout=FNULL, stderr=subprocess.STDOUT)
                config.get_md5sum(new_md5sum_file,'new')
            except subprocess.CalledProcessError:
                config.md5sum['new'].append(('Move Error',' '.join(move_command)))

    #Add appropriate columns and filenames to the assay file in ISA metadata
    if len(extension) != 0:
        metadata_process.modify_assay(metadata_out,GLDS,extension)

#A function to simply detect the array type. Outputs into the 'QC_reporting' directory. Assumes the input is already within 'raw_files'. This means that the rename
#function has already been called.
def detect_array(GLDS_path):
    GLDS = os.path.basename(GLDS_path)
    rawdata_out = os.path.join(GLDS_path,config.microarray_out)
    R_script = os.path.join(config.R_dir,'affyNormQC.R')
    R_command = ["Rscript", R_script, 
                    "--arrayInfoOnly=TRUE",
                    "--outDir="+rawdata_out, 
                    "--QCDir="+os.path.join(rawdata_out,'QC_reporting'), 
                    "-i", os.path.join(rawdata_out,'raw_files'),
                    "--GLDS="+GLDS]
    subprocess.call(R_command)

    array_info = os.path.join(rawdata_out,'QC_reporting','summary_report',GLDS+'_arrayInfo.txt')
    if os.path.exists(array_info):
        with open(array_info) as F:
            array = F.readline().strip('\n')
            arraytype = F.readline().strip('\n')
            if 'Pae_G1a' in arraytype:
                array = 'Pae_G1a'
            elif 'PrimeView' in arraytype:
                array = 'PrimeView'
    else:
        twocolor = False
        for infile in os.listdir(os.path.join(rawdata_out,'raw_files')):
            twocolor = config.detect_2channel(os.path.join(rawdata_out,'raw_files',infile))
        if twocolor:
            array='TwoColor'
        else:
            array='Agilent'

    return array

#Function for normalizing and performing QC on agilent arrays
def sChAgilNormQC(rawdata_out,GLDS):
    R_script = os.path.join(config.R_dir,'sChAgilNormQC.R')
    if config.GPL:
        R_command = ["Rscript", R_script,
                        "-i", os.path.join(rawdata_out,'raw_files'),
                        "-o", os.path.join(rawdata_out,'processed_data',GLDS+"_microarray_normalized"),
                        "-t", 'txt',
                        "--QCDir="+os.path.join(rawdata_out,'QC_reporting'),
                        "--QCpackage=R",
                        "--GLDS="+GLDS]
    else:
        R_command = ["Rscript", R_script,
                        "-i", os.path.join(rawdata_out,'raw_files'),
                        "-o", os.path.join(rawdata_out,'processed_data',GLDS+"_microarray_normalized"),
                        "-t", 'txt',
                        "--pullIDs=true",
                        "--QCDir="+os.path.join(rawdata_out,'QC_reporting'),
                        "--QCpackage=R",
                        "--GLDS="+GLDS]
    subprocess.call(R_command)

#Function for normalizing and performing QC on agilent arrays
def dChAgilNormQC(rawdata_out,GLDS):
    R_script = os.path.join(config.R_dir,'dChAgilNormQC.R')
    if config.GPL:
        R_command = ["Rscript", R_script,
                        "-i", os.path.join(rawdata_out,'raw_files'),
                        "-o", os.path.join(rawdata_out,'processed_data',GLDS+"_microarray_normalized"),
                        "--QCDir="+os.path.join(rawdata_out,'QC_reporting'),
                        "--QCpackage=R",
                        "--GLDS="+GLDS]
    else:
        R_command = ["Rscript", R_script,
                        "-i", os.path.join(rawdata_out,'raw_files'),
                        "-o", os.path.join(rawdata_out,'processed_data',GLDS+"_microarray_normalized"),
                        "--pullIDs=true",
                        "--QCDir="+os.path.join(rawdata_out,'QC_reporting'),
                        "--QCpackage=R",
                        "--GLDS="+GLDS]
    subprocess.call(R_command)

#This function simply runs the R script affyNormQC.R specifying the correct inputs
def qc_and_normalize(rawdata_out,GLDS):
    R_script = os.path.join(config.R_dir,'affyNormQC.R')
    R_command = ["Rscript", R_script, 
                    "-n", "rma", 
                    "-o", GLDS+"_microarray_normalized",
                    "--outDir="+os.path.join(rawdata_out, 'processed_data'), 
                    "--QCDir="+os.path.join(rawdata_out,'QC_reporting'), 
                    "-i", os.path.join(rawdata_out,'raw_files'),
                    "--outType=txt", 
                    "--outputData=TRUE",
                    "--QCoutput=TRUE", 
                    "--QCpackage=R",
                    "--GLDS="+GLDS]
    subprocess.call(R_command)
    if not GLDS+'_microarray_normalized.txt' in os.listdir(os.path.join(rawdata_out,'processed_data')):
        print "Warning: Normalized expression file missing, some processing steps may have failed"

def annotateTwoColor(rawdata_out,GLDS):
    R_script = os.path.join(config.R_dir,'annotateAgilent.R')
    normalized_expression = os.path.join(rawdata_out,'processed_data',GLDS+"_microarray_normalized.rda")
    R_command = ["Rscript", R_script,
                "-i", normalized_expression,
                "--gplDir="+os.path.join(rawdata_out,'raw_files'),
                "-o", os.path.join(rawdata_out,'processed_data',GLDS+"_microarray_normalized-annotated"),
                "-g", 'BSU',
                "-p", 'FeatureNum',
                "-t", 'txt',
                "--QCDir=" + os.path.join(rawdata_out,'QC_reporting'),
                "--GLDS="+GLDS]
    subprocess.call(R_command)

def annotateAgilent(rawdata_out,GLDS):
    R_script = os.path.join(config.R_dir,'annotateAgilent.R')
    normalized_expression = os.path.join(rawdata_out,'processed_data',GLDS+"_microarray_normalized.txt")
    R_command = ["Rscript", R_script,
                "-i", normalized_expression,
                "--gplDir="+os.path.join(rawdata_out,'raw_files'),
                "-o", os.path.join(rawdata_out,'processed_data',GLDS+"_microarray_normalized-annotated"),
                "-g", 'detect',
                "-p", 'detect',
                "-t", 'txt',
                "--QCDir=" + os.path.join(rawdata_out,'QC_reporting'),
                "--GLDS="+GLDS]
    subprocess.call(R_command)


def annotatePae_G1a(rawdata_out,GLDS):
    R_script = os.path.join(config.R_dir,'annotateAgilent.R')
    normalized_expression = os.path.join(rawdata_out,'processed_data',GLDS+"_microarray_normalized.txt")
    R_command = ["Rscript", "--no-save", "--no-restore", R_script,
                    "-i", normalized_expression,
                    "--annotation="+os.path.join(config.R_dir,'Pae_G1a_GPL.txt'),
                    "-o", os.path.join(rawdata_out,'processed_data',GLDS+"_microarray_normalized-annotated"),
                    "-g", 'Representative.Public.ID',
                    "-p", 'Probe.Set.ID',
                    "-t", 'txt',
                    "--QCDir=" + os.path.join(rawdata_out,'QC_reporting'),
                    "--GLDS="+GLDS]
    subprocess.call(R_command)

def annotatePrimeView(rawdata_out,GLDS):
    R_script = os.path.join(config.R_dir,'annotateAgilent.R')
    normalized_expression = os.path.join(rawdata_out,'processed_data',GLDS+"_microarray_normalized.txt")
    R_command = ["Rscript", "--no-save", "--no-restore", R_script,
                    "-i", normalized_expression,
                    "--annotation="+os.path.join(config.R_dir,'GLDS-125_microarray_GPL15207-14509.txt'),
                    "-o", os.path.join(rawdata_out,'processed_data',GLDS+"_microarray_normalized-annotated"),
                    "-g", 'RefSeq.Transcript.ID',
                    "-t", 'txt',
                    "--QCDir=" + os.path.join(rawdata_out,'QC_reporting'),
                    "--GLDS="+GLDS]
    subprocess.call(R_command)


def annotate(rawdata_out,GLDS):
    R_script = os.path.join(config.R_dir,'annotateProbes.R')
    normalized_expression = os.path.join(rawdata_out,'processed_data',GLDS+"_microarray_normalized.txt")
    array_info = os.path.join(rawdata_out,'QC_reporting','summary_report',GLDS+"_arrayInfo.txt")
    R_command = ["Rscript", R_script, 
                    "-i", normalized_expression,
                    "-a", array_info, 
                    "-o", os.path.join(rawdata_out,'processed_data',GLDS+"_microarray_normalized-annotated"),
                    "--QCDir=" + os.path.join(rawdata_out,'QC_reporting')]
    subprocess.call(R_command)

    
#This function runs limma differential script limmaDiffExp.R
def limma_differential(rawdata_out,metadata_out,GLDS):
    condition1,condition2,pval_cut = config.visualize.split(',')
    limma_script = os.path.join(config.R_dir,'limmaDiffExp.R')
    if GLDS + "_microarray_normalized_annotated.txt" in os.listdir(rawdata_out):
        d_path = os.path.join(rawdata_out,GLDS + "_microarray_normalized-annotated.txt")
    elif GLDS + "_microarray_normalized.txt" in os.listdir(rawdata_out):
        print "Warning: No annotated expression file detected, using unannotated file instead - Labelling will be done with probe IDs and will need to be converted to gene names manually."
        d_path = os.path.join(rawdata_out,GLDS + "_microarray_normalized.txt")
    else:
        print "Error: No expression count file detected, exiting..."
        sys.exit(1)
    limma_differential_command = ["Rscript", "--vanilla", limma_script, 
                                    "-d", d_path, 
                                    "-i", metadata_out, 
                                    "--group1=" + condition1, 
                                    "--group2=" + condition2, 
                                    "-o", os.path.join(rawdata_out,GLDS + "_microarray_DGE.txt")]
    subprocess.call(limma_differential_command)



