__author__ = 'Jonathan Rubin'

GL_datasets = './GL_Datasets.txt'
outfile = open('./batch.txt','w')

outfile.write('#Directory=/opt/genelab-genomespace-dev_mount_point/\n')
outfile.write('GLDS#\tCopied\tArrayType\tNormalize/QC\tAnnotated\n')

with open(GL_datasets) as F:
    F.readline()
    for line in F:
        line = line.strip('\n').split('\t')
        outfile.write(line[0]+'\tFalse\tFalse\tFalse\tFalse\n')
