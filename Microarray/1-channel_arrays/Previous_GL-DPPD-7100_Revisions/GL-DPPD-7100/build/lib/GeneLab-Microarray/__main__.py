__author__ = 'Jonathan Rubin'

import os,sys,argparse

def run():
    parser = argparse.ArgumentParser(prog='GeneLab-Microarray',usage='%(prog)s [options] Output',description='Standardized processing pipeline for microarray data on GeneLab.')
    parser.add_argument('Output', help='The full path to the desired output directory.')
    parser.add_argument('-p','--process',help='Specify for process mode. If specified, give a directory to a GLDS directory to be processed.',metavar='',default=False)
    parser.add_argument('-b','--batch',help='Specify for batch processing submode (must also specify process). If specified, input the full directory to a batch.txt file (see README for format guidelines) to the process flag.'
     ,default=False,action='store_const',const=True,metavar='')
    parser.add_argument('-v','--visualize',help='Specify for visualization mode. If selected, must input a comma-separated list of factor values, an adjusted p-value cutoff, and a list of outliers (ex. --visualize flight,ground,0.1,GSM1234_GSM5678) to compare. Multiple factor values can be specified with an underscore ("_") delimiter',
        default=False,metavar='')
    parser.add_argument('-d','--diff',help='[options: limma, voom; default: limma]. Specify type of differential analysis to perform. Limma is used for microarray, Limma-voom is used for RNA-Seq',
        default='limma',metavar='')
    parser.add_argument('-c','--counts',help='Specify a counts table to be used with visualize mode. Must be tab-delimited with gene or probe IDs in the first column and sample values in other columns. Header contains sample names.',
        default=False,metavar='')
    parser.add_argument('-m','--metadata',help='Metadata in ISA tab format. Specifically the s file.',
        default=False,metavar='')
    parser.add_argument('-g','--galaxy',help='For use with Galaxy only. Same outputs as visualization mode simply formatted in a way thats compatible with Galaxy tools.',
        default=False,metavar='')


    #If user does not provide any arguments, simply display help message
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)


    #Parse user-provided arguments 
    args = parser.parse_args()
    batch = args.batch
    indir = args.process
    outdir = os.path.normpath(args.Output)
    visualize = args.visualize
    diff = args.diff
    counts = args.counts
    metadata = args.metadata
    galaxy = args.galaxy



    #Get full paths to locations within this package
    srcdir = os.path.dirname(os.path.realpath(__file__))
    wrkdir = os.getcwd()
    tempdir = os.path.join(os.path.dirname(srcdir),'temp')
    R_dir = os.path.join(srcdir,'R_scripts')


    #Write full paths to locations to a config.py file to be used by other scripts in this package
    with open(os.path.join(srcdir,'config.py'),'w') as outfile:
        outfile.write('indir = "' + str(indir) + '"\n')
        outfile.write('outdir = "' + outdir + '"\n')
        outfile.write('srcdir = "' + srcdir + '"\n')
        outfile.write('wrkdir = "' + wrkdir + '"\n')
        outfile.write('tempdir = "' + tempdir + '"\n')
        outfile.write('R_dir = "' + R_dir + '"\n')
        outfile.write('md5sum = {"original": [], "new": []}\n')
        outfile.write('batch = "' + str(batch) + '"\n')
        outfile.write('visualize = "' + str(visualize) + '"\n')
        outfile.write('microarray_out="microarray"\n')
        outfile.write('GPL=False\n')
        outfile.write("""def get_md5sum(filepath,key,action=False):
    import os, subprocess
    if not action:
        md5sum_command = ["md5sum",filepath]
        md5sum_character = subprocess.check_output(md5sum_command).split(' ')[0].encode("utf-8")
        md5sum[key].append((os.path.basename(os.path.normpath(filepath)),md5sum_character))
    else:
        md5sum_command = ["md5sum",filepath]
        md5sum_character = subprocess.check_output(md5sum_command).split(' ')[0].encode("utf-8")
        md5sum[key].append((action,os.path.basename(os.path.normpath(filepath)),md5sum_character))
def detect_2channel(infile):
    is2channel = False
    with open(infile,'r') as F:
        for line in F:
            if 'rMedianSignal' in line or 'rBGMedianSignal' in line or 'F635' in line or 'Cy5' in line or 'CH2' in line:
                is2channel = True
                return is2channel
    return is2channel\n""")


    #Either run batch module or just run the processing steps on a single dataset
    if indir != False:
        indir = os.path.normpath(indir)
        if batch:
            print "Working Directory: ", wrkdir
            print "Batch option specified.\nUsing batch file: " + indir + "\nWriting output to: " + outdir
            import batch_process
            batch_process.run(indir)
        else:
            import batch_process
            print "Working Directory: ", wrkdir
            print "Single processing mode specified\nInput directory: " + indir + "\nWriting output to: " + outdir
            GLDS = os.path.basename(indir)
            with open(os.path.join(tempdir,'temp_batch.txt'),'w') as outfile:
                outfile.write("#Directory="+os.path.dirname(indir)+"\nGLDS#\tCopied\tArrayType\tNormalize/QC\tAnnotated\n"+GLDS+"\tFalse\tFalse\tFalse\tFalse")
            batch_process.run(os.path.join(tempdir,'temp_batch.txt'))
    elif visualize != False:
        print "Visualization mode specified.\nInput GLDS directory: "+outdir
        visualize_list = visualize.split(',')
        if len(visualize_list) == 3:
            condition1,condition2,padj_cutoff = visualize_list
            outliers = 'None'
        else:
            condition1,condition2,padj_cutoff,outliers = visualize_list
        print "Comparing: " + condition1 + " vs. " + condition2 + "\nAdjusted p-value cutoff set at: " + padj_cutoff + "\nSpecified outliers: " + outliers
        GLDS = os.path.basename(outdir)
        if not counts:
            counts = os.path.join(outdir,'microarray','processed_data',GLDS+'_microarray_normalized-annotated.txt')
            if not os.path.exists(counts):
                counts = os.path.join(outdir,'microarray','processed_data',GLDS+'_microarray_normalized.txt')
                if not os.path.exists(counts):
                    print "No counts files found within specified GLDS directory. Exiting..."
                    sys.exit(1)
        if not metadata:
            metadata = os.path.join(outdir,'metadata','s_'+GLDS+'_microarray_metadata.txt')
            if not os.path.exists(metadata):
                print "No metadata files found within specified GLDS directory. Exiting..."
                sys.exit(1)
        import galaxy_mode
        html_main = os.path.join(outdir,'GeneLab-Visualize_'+condition1+'-'+condition2+'_pval-'+padj_cutoff,'results.html')
        html_folder = os.path.join(outdir,'GeneLab-Visualize_'+condition1+'-'+condition2+'_pval-'+padj_cutoff)
        galaxy_mode.run(counts,metadata,diff,condition1,condition2,padj_cutoff,outliers,html_main,html_folder)
        print "done. Output in: " + outdir
    elif galaxy != False:
        import galaxy_mode
        counts_table,metadata,diff_analysis,condition1,condition2,padj_cutoff,outliers,html_main,html_folder = galaxy.split(',_,')
        galaxy_mode.run(counts_table,metadata,diff_analysis,condition1,condition2,padj_cutoff,outliers,html_main,html_folder)
        print "done."
    else:
        print "Error: No mode selected. See help for information on how to run GeneLab-Microarray. Exiting..."
        sys.exit(1)


