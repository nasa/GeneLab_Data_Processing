__author__ = 'Jonathan Rubin'


#This file is for use with the galaxy mode of GeneLab-Microarray.

#This function looks into metadata files for condition options from a specific field.
def get_conditions(dataset, field_name):
    options = list()

    with open(dataset) as F:
        header = F.readline().strip('\n').split('\t')
        indexes = list()
        i = 0
        for item in header:
            if field_name in item:
                indexes.append(i)
            i += 1
        for line in F:
            line = line.strip('\n').split('\t')
            for ind in indexes:
                options.append((line[ind].strip('"'),line[ind].strip('"'),False))
    options = list(set(options))

    return options

#This function looks into counts table txt file and checks to see whether the column names correspond to the data. If they are off by any amount,
#it assumes that the first column name was ommitted. Otherwise it assumes the first column has gene names.
def get_sample_names(dataset):
    options = list()
    with open(dataset) as F:
        header = F.readline().strip('\n').split('\t')
        testline = F.readline().strip('\n').split('\t')
        if len(testline) == len(header):
            for item in header[1:]:
                options.append((str(item),str(item),False))
        else:
            for item in header:
                options.append((str(item),str(item),False))

    return options

def test_function(value):
    options = [(str(value),str(value),False)]

    return options