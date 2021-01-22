__author__ = 'Jonathan Rubin'

from stat import ST_MTIME
import os, sys, time, subprocess, config

def clean(metadata_directory):
    #Path to the directory (absolute)
    dirpath = metadata_directory

    #Get all entries in the directory w/ stats
    entries = (os.path.join(dirpath, fn) for fn in os.listdir(dirpath))
    entries = ((os.stat(path), path) for path in entries)

    #Insert creation date so we can get the last modified metadata
    entries = ((stat[ST_MTIME], path) for stat, path in entries)

    #Find name of GLDS number
    GLDS = os.path.basename(os.path.dirname(metadata_directory))
    #metadata_out is the path to the output metadata
    metadata_out = os.path.join(config.outdir,GLDS,'metadata')

    #Make appropriate output directory
    if not os.path.exists(metadata_out):
        os.makedirs(metadata_out)

    #Get last modified metadata zip file, copy to the output directory, unzip it, remove the zipped directory, and finally bring all files within folders to the top metadata directory
    i = 0
    for cdate, path in sorted(entries,reverse=True):
        if 'zip' in path and i == 0:
            metadata_zip = os.path.join(metadata_directory,os.path.basename(path))
            zip_filename = os.path.basename(metadata_zip)

            #Check md5sum of original zip file
            config.get_md5sum(metadata_zip,'original',action='copy')

            #Copy the last modified metadata
            cp_command = ["cp","-r",metadata_zip,metadata_out]
            #Unzip it into the metadata_out directory
            unzip_command = ["unzip", "-o", "-qq", os.path.join(metadata_out,zip_filename), "-d", metadata_out]
            #Remove the .zip compressed file to avoid confusion and save space
            remove_zip_command = ["rm",os.path.join(metadata_out,zip_filename)]

            #Execute copy command
            subprocess.call(cp_command)
            subprocess.call(unzip_command)

            #Verify md5sum for 'new' file
            config.get_md5sum(os.path.join(metadata_out,zip_filename),'new')

            #Execute unzipping and zip removal commands
            subprocess.call(remove_zip_command)

            i += 1

    #Loop through the metadata_out directory in case the unzipping produces a folder. If so, mv contents of folder up one directory and remove folder
    for filename in os.listdir(metadata_out):
        if os.path.isdir(os.path.join(metadata_out,filename)):
            move_command = ["mv", os.path.join(metadata_out,filename,"*"), metadata_out]

            #This is needed because the subprocess command cannot inherently deal with wildcards...
            shell_move_command =  ' '.join(move_command)

            remove_folder_command = ["rm", "-r",os.path.join(metadata_out,filename)]
            subprocess.call(shell_move_command, shell=True)
            subprocess.call(remove_folder_command)

    #Rename all metadata files to a standard naming convention
    for filename in os.listdir(metadata_out):
        isa = filename.split('_')[0]
        newfilename = isa + '_' + GLDS + '_microarray_metadata.txt'
        if not os.path.exists(os.path.join(metadata_out,newfilename)):
            config.get_md5sum(os.path.join(metadata_out,filename),'original',action='rename')
            move_command = ["mv","'"+os.path.join(metadata_out,filename)+"'",os.path.join(metadata_out,newfilename)]
            try:
                with open(os.devnull,'w') as FNULL:
                    subprocess.check_call(' '.join(move_command), shell=True,stdout=FNULL, stderr=subprocess.STDOUT)
                config.get_md5sum(os.path.join(metadata_out,newfilename),'new')
            except subprocess.CalledProcessError:
                config.md5sum['new'].append(('Move Error',' '.join(move_command)))

    #Modify the investigation file to account for sample and assay renaming
    modify_i(GLDS,os.path.join(metadata_out,'i_' + GLDS + '_microarray_metadata.txt'))



#Modifies investigation file to correctly detect sample and assay files that have been renamed
def modify_i(GLDS,i_file):
    #First read in all lines in investigation file
    with open(i_file,'r') as infile:
        lines = infile.readlines()

    #Then rewrite the same lines simply substituting in the correct filenames
    with open(i_file,'w') as outfile:
        for line in lines:
            if 'Study File Name' in line:
                line = 'Study File Name\t"s_'+GLDS+'_microarray_metadata.txt"\n'
            if 'Study Assay File Name' in line:
                line = 'Study Assay File Name\t"a_'+GLDS+'_microarray_metadata.txt"\n'
            outfile.write(line)


#Creates an assay dictionary which is basically just the metadata (specifically the assay file in the isa metadata)
#where the key is the first column (assumed to be sample name) and the value is the rest of the columns. This is 
#mainly used to rename files.
def read_assay(metadata_out):
    #Loop through metadata files, find the assay file (starts with 'a_')
    for filename in os.listdir(metadata_out):
        if 'a_' in filename[:2] and config.microarray_out in filename:
            assay_file = os.path.join(metadata_out,filename)

    #Create an assay dictionary where the key is the name of the sample file
    assay_dict = dict()
    try:
        with open(assay_file) as F:
            F.readline()
            for line in F:
                linelist = [item.strip('"') for item in line.strip('\n').split('\t')]
                assay_dict[linelist[0]] = linelist[1:]

        return assay_dict
    except:
        print "Warining: No assay file found in ISA metadata. Files will be renamed without considering metadata."
        return assay_dict

def modify_assay(metadata_out,GLDS,extension):
    for filename in os.listdir(metadata_out):
        if 'a_' in filename[:2]:
            assay_file = os.path.join(metadata_out,filename)

    #Create an assay dictionary where the key is the name of the sample file
    new_assay_file = list()
    try:
        with open(assay_file) as F:
            new_assay_file.append(F.readline().replace('\r','').replace('\n','').replace('^M','')+'\t'+'\t'.join(['"Protocol REF"','"Parameter Value[Raw File]"',
                '"Term Source REF"','"Term Accession Number"','"Parameter Value[Processed Data]"',
                '"Term Source REF"','"Term Accession Number"']))
            for line in F:
                linelist = line.strip('\n').split('\t')
                basefilename = linelist[0].split('.')[0].replace('_','-').replace('(','-').replace(')','-').replace(' ','-').replace(GLDS,'').replace('microarray','').replace('--','-').strip('-').strip('"')
                raw_filename = GLDS + '_' + basefilename + '_microarray_raw.'+extension
                new_assay_file.append(line.replace('\r','').replace('\n','').replace('^M','')+'\t'+'\t'.join(['"GeneLab data processing protocol"',
                    '"'+raw_filename+'"','""', '""',
                    '"'+GLDS+'/microarray/processed_data/"','""', '""']))
        with open(assay_file,'w') as outfile:
            outfile.write('\r\n'.join(new_assay_file))
    except:
        print "Warning: No assay file found in ISA metadata. Proceeding without metadata modifications."

#Creates a .txt file with all md5sum output
def create_md5sum_out(rawdata_out,GLDS):
    with open(os.path.join(rawdata_out,'raw_files',GLDS+'_md5sum.txt'),'w') as outfile:
        outfile.write('#Action\tCheck\tOriginalFile,md5sum -> NewFile,md5sum\n')
        for original,new in zip(config.md5sum['original'],config.md5sum['new']):
            if original[0] != 'remove':
                check = original[2] == new[1]
                if not check:
                    print "Warning: " + original[0] + " action for file " + original[1] + " did not pass md5sum check."
                outfile.write(original[0] + '\t' + str(check) + '\t' + ','.join(original[1:])+' -> '+','.join(new)+'\n')
            else:
                outfile.write(original[0] + '\tN/A\t' + ','.join(original[1:])+' -> '+','.join(new)+'\n')

if __name__ == "__main__":
    metadata_out = '/Users/jonathanrubin/Google_Drive/NASA/home/processed_GLDS/GLDS-4/metadata/'
    GLDS = 'GLDS-4'
    extension = 'CEL'
    outdir = '/Users/jonathanrubin/Google_Drive/NASA/home/assay_test/GLDS-4/metadata/'
    modify_assay(metadata_out,GLDS,extension,outdir)


