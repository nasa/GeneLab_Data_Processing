indir = "False"
outdir = "."
srcdir = "/Users/jonathanrubin/Google_Drive/NASA/home/GeneLab-Microarray/GeneLab-Microarray"
wrkdir = "/private/var/folders/rr/fbksjg6n2d7_3gh_8t_jgj_r0000gn/T/tmpg3DLPi/job_working_directory/000/3/working"
tempdir = "/Users/jonathanrubin/Google_Drive/NASA/home/GeneLab-Microarray/temp"
R_dir = "/Users/jonathanrubin/Google_Drive/NASA/home/GeneLab-Microarray/GeneLab-Microarray/R_scripts"
md5sum = {"original": [], "new": []}
batch = "False"
visualize = "False"
microarray_out="microarray"
GPL=False
def get_md5sum(filepath,key,action=False):
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
            if 'rMedianSignal' in line or 'rBGMedianSignal' in line or 'F635' in line or 'Cy5' in line:
                is2channel = True
                return is2channel
    return is2channel
