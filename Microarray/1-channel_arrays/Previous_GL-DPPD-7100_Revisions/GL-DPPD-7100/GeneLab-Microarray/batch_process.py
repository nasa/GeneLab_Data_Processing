__author__ = 'Jonathan Rubin'

import os, config, metadata_process, rawdata_process

def run(batch_file):
    compatible_arrays = ['Affymetrix','Pae_G1a','Agilent','PrimeView','TwoColor']
    batch_list = list()
    with open(batch_file) as F:
        parent_dir = F.readline().strip('\n').split('=')[1]
        header = F.readline()
        for line in F:
            linelist = line.strip('\n').split()
            if len(linelist) != 5:
                print "Error, batch file line not formatted properly: " + line + " skipping..."
            else:
                batch_list.append(linelist)


    for i in range(len(batch_list)):
        if 'False' in batch_list[i]:
            GLDS, copy, array, norm_qc, annotate = batch_list[i]
            GLDS_path = os.path.join(config.outdir,GLDS)
            config.microarray_out='microarray'
            rawdata_out = os.path.join(config.outdir,GLDS,config.microarray_out)
            metadata_out = os.path.join(config.outdir,GLDS,'metadata')

            #Copy module, copies and unzips both metadata and raw data. If precise directories are not found,
            #that GLDS is skipped.
            if copy == 'False':
                print "Copying files for " + GLDS + "..."
                config.md5sum = {"original": [], "new": []}

                #Process metadata
                metadata_in = os.path.join(parent_dir,GLDS,'metadata')
                if os.path.isdir(metadata_in) and 'GSE' not in GLDS:
                    metadata_process.clean(metadata_in)
                else:
                    print "metadata directory within " + GLDS + " not found, skipping..."
                    copy, array, norm_qc, annotate = ['Skipped' for j in range(4)]
                    batch_list[i] = [GLDS, copy, array, norm_qc, annotate]

                #Copy rawdata and metadata into output
                rawdata_in = os.path.join(parent_dir,GLDS,'microarray')
                config.GPL = False
                if os.path.isdir(rawdata_in):
                    rawdata_process.copy(rawdata_in)
                    rawdata_process.rename(os.path.join(config.outdir,GLDS))
                    metadata_process.create_md5sum_out(rawdata_out,GLDS)
                    batch_list[i][1] = 'True'
                elif os.path.isdir(os.path.join(parent_dir,GLDS,'micoarray')):
                    rawdata_in = os.path.join(parent_dir,GLDS,'micoarray')
                    rawdata_process.copy(rawdata_in)
                    rawdata_process.rename(os.path.join(config.outdir,GLDS))
                    metadata_process.create_md5sum_out(rawdata_out,GLDS)
                    batch_list[i][1] = 'True'
                elif os.path.isdir(os.path.join(parent_dir,GLDS,'transcriptomics')):
                    rawdata_in = os.path.join(parent_dir,GLDS,'transcriptomics')
                    rawdata_process.copy(rawdata_in)
                    rawdata_process.rename(os.path.join(config.outdir,GLDS))
                    metadata_process.create_md5sum_out(rawdata_out,GLDS)
                    batch_list[i][1] = 'True'
                elif True in ['microarray' in x for x in os.listdir(os.path.join(parent_dir,GLDS))]:
                    for folder in os.listdir(os.path.join(parent_dir,GLDS)):
                        if 'microarray' in folder:
                            config.microarray_out = folder
                            rawdata_in = os.path.join(parent_dir,GLDS,folder)
                            rawdata_out = os.path.join(config.outdir,GLDS,folder)
                            config.GPL = False
                            rawdata_process.copy(rawdata_in)
                            rawdata_process.rename(os.path.join(config.outdir,GLDS))
                            metadata_process.create_md5sum_out(rawdata_out,GLDS)
                            array = rawdata_process.detect_array(GLDS_path)
                            if array == 'Pae_G1a':
                                rawdata_process.qc_and_normalize(rawdata_out,GLDS)
                                rawdata_process.annotatePae_G1a(rawdata_out,GLDS)
                            elif array == 'PrimeView':
                                rawdata_process.qc_and_normalize(rawdata_out,GLDS)
                                rawdata_process.annotatePrimeView(rawdata_out,GLDS)
                            elif array == 'Affymetrix':
                                rawdata_process.qc_and_normalize(rawdata_out,GLDS)
                                rawdata_process.annotate(rawdata_out,GLDS)
                            elif array == 'TwoColor':
                                rawdata_process.TwoColorNormQC(rawdata_out,GLDS)
                                rawdata_process.annotateTwoColor(rawdata_out,GLDS)
                            else:
                                rawdata_process.sChAgilNormQC(rawdata_out,GLDS)
                                rawdata_process.annotateAgilent(rawdata_out,GLDS)
                            if os.path.exists(os.path.join(rawdata_out,'processed_data')):
                                for file1 in os.listdir(os.path.join(rawdata_out,'processed_data')):
                                    if 'annotated' in file1:
                                        annotate = 'True'
                    copy, norm_qc = ['Multi' for j in range(2)]
                    batch_list[i] = [GLDS, copy, array, norm_qc, annotate]
                else:
                    print "microarray directory within " + GLDS + " not found, skipping..."
                    copy, array, norm_qc, annotate = ['Skipped' for j in range(4)]
                    batch_list[i] = [GLDS, copy, array, norm_qc, annotate]

                update_batch(parent_dir,header,batch_file,batch_list)
                print "done"
            elif copy != 'True':
                print "Warning: Files were not copied for " + GLDS + ". If this was not desired, check batch file and make sure this GLDS was set to 'False'."

            #Array module, this part simply generates an arrayInfo.txt file and reads it in. If the array is part of a list of arrays that we can process then
            #continue, otherwise skip the GLDS
            if array == 'False':
                print "Detecting array type for " + GLDS + "..."
                array = rawdata_process.detect_array(GLDS_path)
                if array != 'Skipped':
                    batch_list[i][2] = array
                    if array in compatible_arrays:
                        update_batch(parent_dir,header,batch_file,batch_list)
                        print "done"
                    else:
                        print "Warning: " + GLDS + " " + array + " arrays not currently supported, skipping..."
                        norm_qc, annotate = ['Skipped' for j in range(2)]
                        batch_list[i] = batch_list[i][:2]+[array,norm_qc,annotate]
                        update_batch(parent_dir,header,batch_file,batch_list)
                else:
                    norm_qc, annotate = ['Skipped' for j in range(2)]
                    batch_list[i] = batch_list[i][:2]+[array,norm_qc,annotate]
                    update_batch(parent_dir,header,batch_file,batch_list)
            elif array != 'True' and array != 'Skipped' and array != 'Multi':
                print "Warning: Array was not detected for " + GLDS + ". If this was not desired, check batch file and make sure this GLDS was set to 'False'."

            #Performs normalization and qc pre- and post-normalization
            if norm_qc == 'False':
                print "Performing QC, normalization, and post-normalization QC on data for " + GLDS + "..."
                if array == 'Pae_G1a' or array == 'PrimeView' or array == 'Affymetrix':
                    rawdata_process.qc_and_normalize(rawdata_out,GLDS)
                elif array == 'TwoColor':
                    rawdata_process.dChAgilNormQC(rawdata_out,GLDS)
                else:
                    rawdata_process.sChAgilNormQC(rawdata_out,GLDS)
                if os.path.exists(os.path.join(rawdata_out,'processed_data')):
                    for file1 in os.listdir(os.path.join(rawdata_out,'processed_data')):
                        if 'normalized' in file1:
                            norm_qc = 'True'
                batch_list[i][3] = norm_qc
                update_batch(parent_dir,header,batch_file,batch_list)
                print "done"
            elif norm_qc != 'True' and norm_qc != 'Skipped' and norm_qc != 'Multi':
                print "Warning: QC and normalization not performed for " + GLDS + ". If this was not desired, check batch file and make sure this GLDS was set to 'False'."

            #Annotates probeIDs with gene names. Autodetection of array annotation package is attempted but if it fails then return 'Skipped'.
            if annotate == 'False':
                print "Annotating probe IDs with gene names for " + GLDS + "..."
                if array == 'Pae_G1a':
                    rawdata_process.annotatePae_G1a(rawdata_out,GLDS)
                elif array == 'PrimeView':
                    rawdata_process.annotatePrimeView(rawdata_out,GLDS)
                elif array == 'Affymetrix':
                    rawdata_process.annotate(rawdata_out,GLDS)
                elif array == 'TwoColor':
                    rawdata_process.annotateTwoColor(rawdata_out,GLDS)
                else:
                    rawdata_process.annotateAgilent(rawdata_out,GLDS)
                if os.path.exists(os.path.join(rawdata_out,'processed_data')):
                    for file1 in os.listdir(os.path.join(rawdata_out,'processed_data')):
                        if 'annotated' in file1:
                            annotate = 'True'
                batch_list[i][4] = annotate
                update_batch(parent_dir,header,batch_file,batch_list)
                print "done"

    print "done."

#This function updates the batch file
def update_batch(parent_dir,header,batch_file,batch_list):
    with open(batch_file,'w') as outfile:
        outfile.write('#Directory='+parent_dir+'\n')
        outfile.write(header)
        for linelist in batch_list:
            outfile.write('\t'.join(linelist)+'\n')
        