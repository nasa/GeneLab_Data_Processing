# GeneLab bioinformatics processing pipeline for one-channel microarray data

> **This page holds an overview and example code of how GeneLab processes one-channel Microarray datasets. Exact processing code and GL-DPPD-7100 revision used for specific datasets are available in the [GLDS_Processing_Scripts](https://developer.nasa.gov/asaravia/GeneLab_Data_Processing/tree/master/Microarray/1-channel_arrays/GLDS_Processing_Scripts) sub-directory and are also provided with their processed data in the [GeneLab Data Systems (GLDS) repository](https://genelab-data.ndc.nasa.gov/genelab/projects).**  

---

**Date:** November 19, 2018  
**Revision:** -   
**Document Number:** GL-DPPD-7100  

**Submitted by:**  
Matt Geniza (GeneLab Data Processing Team)  
Jonathan Rubin and Daniel Mattox (GeneLab Data Processing Team, summer 2018 Interns)

**Approved by:**  
Sylvain Costes (GeneLab Project Manager)  
Marla Smithwick (GeneLab Deputy Project Manager)    
Jonathan Galazka (GeneLab Project Scientist)    
John M. Costa Sr. (GeneLab Configuration Manager)

---

## Table of Contents
1. <A href=#Installation>Installation</A>
   - <A href=#Dependencies>Dependencies</A>
        + <A href=#R>R</A>
        + <A href=#Python>Python</A>
2. <A href=#Running>Running GeneLab-Microarray</A>
   - <A href=#ProcessMode>Process Mode</A>
        + <A href=#BatchSubmode>Batch Submode</A>
   - <A href=#VisualizationMode>Visualization Mode</A>
   - <A href=#GalaxyMode>Galaxy Mode</A>
3. <A href=#Directory>Directory Structure</A>
4. <A href=#BatchFile>Batch File Format</A>
5. <A href=#RScripts>R Scripts</A>
   - <A href=#NormalizationQC>Normalization and QC</A>
        + <A href=#AffyNormQC>Affymetrix</A>
        + <A href=#sCh>Single Channel Agilent</A>
        + <A href=#dCh>Dual Channel Agilent</A>
   - <A href=#Annotation>AnnotationC</A>
        + <A href=#AffyNormQC>Affymetrix</A>
        + <A href=#Agil>Agilent/Non-Affy</A>
   - <A href=#DiffExp>Differential Expression</A>
        + <A href=#limmaDiffExp>Microarrays</A>
        + <A href=#limmaVoom>RNA-Seq</A>


<H2 id="Installation">Installation</H2>
GeneLab-Microarray writes to itself so to ensure you can easily run it, it is recommended to pip install with the --user flag. Example:

```
pip install --user GeneLab-Microarray
```

<H3 id="Dependencies">Dependencies</H3>
GeneLab-Microarray uses md5sum to check files that are copied/renamed. If you do not have this installed, and you're on MacOSX, the md5 -r option should be identical output to md5sum so within your ~/.bashrc or ~/.profile add the following:

```
alias md5sum='md5 -r'
```

<H4 id="R">R</H4>
And to finish installation, open an R session and run the following commands (note: you should only need to run the 'Primary package' installation commands, we provide the dependencies in case installing it this way doesn't work):

```
install.packages("optparse")
source("http://bioconductor.org/biocLite.R")

# Primary package:
biocLite("affy")
 ## affy dependencies:
 biocLite("zlibbioc")
 biocLite("Biobase")

# Primary package:
biocLite("affyPLM")
 ## affyPLM dependencies:
 biocLite("S4Vectors")
 biocLite("IRanges")
 biocLite("XVector")
 biocLite("Biostrings")

# Primary package:
biocLite("oligo")
 ## oligo dependencies:
 biocLite("bit")
 biocLite("ff")
 biocLite("bitops")
 biocLite("RCurl")
 biocLite("GenomicRanges")
 biocLite("matrixStats")
 biocLite("Rcpp")
 biocLite("bit64")
 biocLite("digest")
 biocLite("RSQLite")

# Primary package:
biocLite("genefilter")
 ## genefilter dependencies:
 biocLite("XML")

# Primary package:
biocLite("limma")

# Primary package:
biocLite("arrayQualityMetrics")
 ## arrayQualityMetrics dependencies:
 biocLite("hexbin")
 biocLite("jsonlite")
 biocLite("openssl")
 biocLite("stringi")
 biocLite("reshape2")
 biocLite("Cairo")
 # dependencies list incomplete



```

<H4 id="Python">Python</H4>
Python packages are installed automatically when installing GeneLab-Microarray using pip. But for documentation purposes, these are the packages (with version numbers) that are known to work:

matplotlib - v2.2.2
mpld3 - v
sklearn - 

Once the above steps are completed without error, you should be able to call GeneLab-Microarray from any directory. Try:

```
GeneLab-Microarray --help
```

<H2 id="Running">Running GeneLab-Microarray</H2>
Gene-Lab microarray has two modes: process and visualization. There are 3 possible flags to give GeneLab-Microarray (--process,--batch, and --visualize). In all cases, an output directory must be specified as the last argument (no flags) - this is a positional argument.

<H3 id="ProcessMode">Process Mode</H3>
Processing mode requires as input a path to a GLDS directory. GeneLab-Microarray will copy raw files from your specified directory into your output directory, rename them according to standard specifications, and perform QC and normalization.


Example:
```
GeneLab-Microarray --process /opt/genelab-genomespace-dev_mount_point/GLDS-4/ /opt/jdrubin/batch_out/
```

<H4 id="BatchSubmode">Batch Submode</H4>
Within the processing mode, there is a submode called batch (specified with -b/--batch). This submode will batch process all GLDS directories specified within a batch.txt file (see <A href=#BatchFile>batch file format</A>). If specified, input into process a batch.txt file instead of a GLDS directory. There are example batch files within the batch/ directory.



Example:
```
GeneLab-Microarray --batch --process /opt/jdrubin/GeneLab-Microarray/batch/batch.txt /opt/jdrubin/batch_out/
```

<H3 id="VisualizationMode">Visualization Mode</H3>
The visualization mode for GeneLab-Microarray is specified with the -v/--visualize flag and takes as input a comma separated list of factor values (multiple factor values can be specified with an underscore delimiter) followed by an adjusted p-value cutoff (for plotting purposes), and optionally a list of outliers (delimited by underscores). The visualization mode will output an html file with interactive graphs and a png file with identical graphs that can be used for publication. This mode can be run in two ways: First, if the data is structured as specified by this package, the output directory can be used as the input directory and this package will automatically look for input files in specific places. Alternatively, a user can specify a counts table (using the -c option) and a metadata 's' file (using the -m option). When using it this way, the final GLDS directory will simply be the output directory (i.e. can be any directory)


Example (data processed using GeneLab-Microarray):
```
GeneLab-Microarray --visualize 'flight,ground,0.1' /opt/jdrubin/batch_out/GLDS-4/
```

Example (counts and metadata file specified by user):
```
GeneLab-Microarray -c /opt/jdrubin/batch_out/GLDS-4/processed_data/GLDS-4_microarray_normalized-annotated.txt -m /opt/jdrubin/batch_out/GLDS-4/processed_data/s_GLDS-4_microarray_metadata.txt --visualize 'flight,ground,0.1' /opt/jdrubin/batch_out/GLDS-4/
```

The visualization mode can also do comparisons with multiple factor values which are specified with underscore delimiters.



Example:
```
GeneLab-Microarray --visualize 'flight_geneKO,flight_noKO,0.1' /opt/jdrubin/batch_out/GLDS-4/
```


Finally, visualization mode can take as inputs a list of outliers.


Example:
```
GeneLab-Microarray --visualize 'flight_geneKO,flight_noKO,0.1,GSM1234_GSM5678' /opt/jdrubin/batch_out/GLDS-4/
```

<H3 id="GalaxyMode">Galaxy Mode</H3>
Galaxy mode is a special module specifically designed for use within galaxy. This mode can be run in the command line but its output will be the exact same as the visualization mode. For completeness, details on how to run this mode on the command line are provided here. Galaxy mode takes a single list of inputs separated by a custom delimiter (`,_,`). The command within galaxy is run as follows:

```
GeneLab-Microarray --galaxy '$counts_filename',_,'$meta_filename',_,'$diff_analysis',_,'$condition1',_,'$condition2',_,$padj,_,'$outliers',_,'$html_file',_,'${html_file.extra_files_path}' .
```

These inputs are provided by the user and formatted specifically to be used within Galaxy but those same inputs could be provided within the command line (the visualization mode accomplishes this same task).

<H2 id="Directory">Directory Structure</H2>
GeneLab-Microarray expects directories to be in a specific structure. A parent directory with GLDS-# followed by two subdirectories (where one is named metadata and the other microarray) each of which contain zipped archives with either raw microarray data or ISA formatted metadata. For example:

```
GLDS-#/
|-metadata/
  |--metadata_ISA.zip
|-microarray/
  |--microarray_raw.tar
```


<H2 id="BatchFile">Batch File Format</H2>
If `-b,--batch` option is desired. In addition to calling the flag, users must submit the full path to a batch.txt file (examples and a simple script to create this batch.txt file is located in the batch subdirectory). 

Briefly, the batch.txt file expects the first line to begin with '#' followed by 'Directory=' then a full path to a directory. The rest of the file is a tab delimited txt file with 3 columns (header is required). The first column is the name of a folder within the specified Directory. All subsequent columns are booleans (True or False) and are used to keep track of the progress of processing the desired data in batch. For example:

```
#Directory=/opt/genelab-genomespace-dev_mount_point/
GLDS#     Copied    ArrayType    Normalize/QC    Annotated
GLDS-4    False     False        False           False
```

GeneLab-Microarray will overwrite the specified batch.txt file changing booleans to True or Skipped when the specific step is finished. An example of a batch.txt file can be found within the `batch/` folder

<H2 id="RScripts">R Scripts</H2>
The R scripts used to read in and process the data are also available to be run manually outside of the GeneLab-Microarray package. In general, the set default options for the individual R scripts are the options used within the wrapped package. To call the R scripts from the command line, the recommended usage is with the --no-save --no-restore flags as seen below (as called from the parent directory of the package):
```
Rscript --no-save --no-restore GeneLab-Microarray/R_scripts/RScriptHere.R 
```

<H3 id="NormalizationQC">Normalization and QC</H3>
The R scripts `affyNormQC.R`, `sChAgilNorm.R`, and `dChAgilNorm.R` can be used to read in raw microarray data, provide raw QC reports, process and normalize the data, and provide post-normalization QC reports for Affymetrix, single-channel Agilent, and two-channel Agilent microarray experiments respectively. Generally, these scripts will read in all of the raw files in the indicated input directory and return a normalized, unannotated table of log2 intensities with probes IDs in the first column and subsequent columns for each sample. Two-channel microarray experiments will not fit this format but that is discussed further below. Additionally, these scripts will generate a QC reporting directory (`./QC_reporting/` by default) with `raw_report` and `normalized_report` sub-directories. If the QC reports are generated with the base R plotting options, each sub-directory will contain a collection of QC figures saved as .png files. If the QC reports are generated with the arrayQualityMetrics tool (aqm), the respective sub-directories will contain the individual figures saved as a combination of .pdf, .svg, and .png files, with the interactive summary report accessible through the `index.html` file. At present, it is recommended to open the `index.html` file with either Chrome or Firefox as some bugs have been observed with some versions of Safari.
Each of these scripts apply as similar of normalization techniques as possible, using the "normexp" (normal plus exponential convolution) background correction method, quantile normalization, and log2 transformation.

<H4 id="AffyNormQC">Affymetrix</H4>
The most recent list of current options and set defaults can be viewed by running the following line from the parent directory:
```
Rscript --no-save --no-restore GeneLab-Microarray/R_scripts/affyNormQC.R --help
```
The only required option is the `-i/--input` to indicate a path to the directory containing the raw .CEL files. Any .CEL file in that directory will be read in and processed. The script is set such that it will use the "affy" package for older Affymetrix arrays and the "oligo" package for the newer Affymetrix arrays (primarily the ST arrays).
After the .CEL files are read in, an array info text file is created along with a `summary_report` directory within the QC report directory. The text file contains the array manufacturer on the first line and the array version on the second line.
The normalized output intensity table is saved to the specified file and file type, either a text file, binary RData file, or both.

<H4 id="sCh">Single Channel Agilent</H4>
The most recent list of current options and set defaults can be viewed by running the following line from the parent directory:
```
Rscript --no-save --no-restore GeneLab-Microarray/R_scripts/sChAgil.R --help
```
Again, the only required option is the `-i/input` option to point to a directory containing the raw microarray data files as text files. The script will try to read in any file in that directory with the "raw" semantic tag ("_raw.txt"). This script uses the "limma" package to read in and normalize the data. It assumes the data is in the typical Agilent data formatting, containing the columns "gMedianSignal", "gBGMedianSignal", and "FeatureNum". This could be adapted by altering the source type and column names to be able to read in other array types with other naming schema.
















<H2 id="affyNormQC">Affy QC and Normalization</H2>
The affyNormQC.R script can be run from any directory, but requires to be pointed to the appropriate directory containing Affymetrix .CEL file microarray data with the `-i/--input` option. It can determine the version of the array and load the appropriate packages (ie "affy" for earlier microarrays and "oligo" for the newer arrays). No inputs are required to run it, but to view the available options, simply run the line below:

```
Rscript --no-save --no-restore affyNormQC.R --help
```

Before running this script, it may be necessary to run the commented out lines immediately below the shebang in an R session to be sure all of the necessary packages are installed

```
install.packages("optparse")
source("http://bioconductor.org/biocLite.R")
biocLite("affy")
biocLite("affyPLM")
biocLite("oligo")
biocLite("arrayQualityMetrics")
```
Please note: bugs in the quality control report have been logged when opening the html file with Safari. At this point, either Chrome or Firefox are recommended for the viewing of the html QC report.
An example run with all of the options explicitly set to the default or example options:

```
Rscript --no-save --no-restore affyNormQC.R -i path/to/input/files/ -n rma -o expValues --outType=both --outputData=TRUE --arrayInfoOnly=FALSE --QCoutput=TRUE --QCDir=./QC_reporting/ --GLDS=21
```

This script can also be used to detect the Affymetrix array information only, outputting a text file containing the manufacturer and the array version and quitting before normalizing the data or performing QC. This option can be accessed by setting `--arrayInfoOnly=TRUE`. However, the array information txt file will be output in the standard mode as well.

<H2 id="annotateProbe">Probe annotation</H2>
The `annotateProbes.R` script can be used to map probe IDs to RefSeq gene IDs using annotation packages. If an array type has not been seen before, the annotation package will need to be manually loaded into the array:annotation pseudo-dictonary. The available options for the script are viewable by the following command:

```
Rscript --no-save --no-restore annotateProbes.R --help
```

The required packages to be install prior to running shown here, as well as in the script immediately below the shebang.

```
install.packages("optparse")
source("http://bioconductor.org/biocLite.R")
biocLite("genefilter")
biocLite("mogene10sttranscriptcluster.db")
```

An example call with all of the default/recommended options explicitly defined:

```
Rscript --no-save --no-restore annotateProbes.R -i path/to/normalized/data.txt -a GLDS-4_arrayInfo.txt -o annotExpValues --outType=both --dupProbes=max
```

<H2 id="limmaDiffExp">Limma Differential Expression</H2>
The `limmaDiffExp.R` script can be used to calculate changes in expression between two groups of samples for a given dataset. The available options can be examined by calling:

```
Rscript --no-save --no-restore limmaDiffExp.R --help
```

This script requires the packages: `optparse` and `limma`. If these packages are not installed, they can be by running the following line in an R session:

```
install.packages("optparse")
source("http://bioconductor.org/biocLite.R")
biocLite("limma")
```

An example call with options set for all parameters is shown below:

```
Rscript --no-save --no-restore limmaDiffExp.R -d /path/to/normalized/data.txt -i ../path/to/sample/s_metadata.txt --group1=flight_KO --group2=ground_KO -o GLDS-4_microarray_DGE.txt --rmOutlies=GSM1234_GSM1235
```

