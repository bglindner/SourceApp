      ________                                   _____                 
     /   ____/ ____  __ _________  ____  ____   /  _  \ ______ ______  
     \____  \ /  _ \|  |  \_  __ \/ ____/ __ \ /  /_\  \  __ \\  __ \ 
     /       (  |_| |  |  /|  | \\  \__\  ___//   /|   |  |_ ||  |_ |
    /________/\____/|____/ |__|   \_____\_____\___||__ |   __/|   __/ 
                                                      \|__|   |__|   

Python implementation of the Unix-based environmental monitoring tool.

SourceApp is still in active development. Check back here for updates and the eventual release.

# Description 
SourceApp is a bioinformatic workflow designed to apportion fecal signal amongst
multiple competing sources in short read metagenomes collected from impaired 
waterways. SourceApp was developed for use on Unix or Unix-like operating systems
and its implementation in Python has the same requirement.

SourceApp is designed to automate all tasks necessary to detect and quantify fecal
signal in metagenomes collected from the water environment if users will simply 
provide gzipped copies of the raw reads in FASTQ format. SourceApp requires a
database of reference genomes with known source-associations. Consequently,
a database was developed in parallel with this software through systematic 
literature review and should be downloaded for use at: <URL>

SourceApp reports metrics akin to prokaryotic cell fraction and has been bench
marked in laboratory trials. SourceApp was found to perform best when fecal signal 
is relatively fresh and when the user is able to provide supplementary reference 
genomes to bolster SourceApp's default genomic reference database. De novo 
databases can be constructed by the user using sourceapp_build.py and then 
supplied here to SourceApp in place of the default database. 

For functional population proportioning in biological wastewater treatment datasets, see "WasteApp" -- the built envrionment extension of SourceApp's working principles (https://github.com/blindner6/WasteApp). 

# Installation

The development version of SourceApp is available here for those interested in testing the software before its publication. It is anticipated that at the time of publication, we will create and maintain a Bioconda recipe for direct installation with `conda` or `mamba`. In the interim, we suggest the following for obtaining a working version of the tool:
```
# clone this repo
git clone https://github.com/blindner6/SourceApp.git

# create two environments using the provided .yml files, one supporting the 
# dependencies for "sourceapp.py" and the other for "sourceapp_build.py"
mamba create -n sourceapp
mamba env update -n sourceapp --file SourceApp/sourceapp.yml
mamba create -n sourceapp_build
mamba env update -n sourceapp_build --file SourceApp/sourceapp_build.yml

# test your installation:
micromamba activate sourceapp
python sourceapp.py -h

```
# Usage

`sourceapp.py` expects as input paired short read metagenomic data. The user should supply these reads as gzipped FASTQ files. No prior adapter trimming or QC is necessary as SourceApp will automate read trimming with `fastp` but the user can disable this step if they would like to perform it themselves with `--skip-trimming`. In addition to short reads, SourceApp needs a genomic database specifically formatted for use by the tool. This can be created by SourceApp from genomes provided by the user with `sourceapp_build.py` or a default database will be provided in a forthcoming publication (Graham et al, _in prep_). Users should specify the location of an output directory for results to be written to. 

SourceApp is primarily designed for use with a Unix-based HPC and no support is offered for deployment with alternative operating systems.

```
usage: sourceapp.py [-h] -i  -o  -d  [-l] [-r] [-q] [-t] [--use-geq] [--no-limits] [--skip-trimming]

SourceApp: Python implementation of the Unix-based environmental monitoring tool.

options:
  -h, --help            Show this help message and exit
  -i , --input-files    Comma-delimited path to forward and reverse metagenomic reads. Must 
                        be in FASTQ format and gzipped (reads.1.fastq.gz,reads.2.fastq.gz)
  -o , --output-dir     Path to the desired output directory
  -d , --sourceapp-database 
                        Path to directory containing a SourceApp formatted database. Default 
                        database available for download or produced de novo as the output directory 
                        from sourceapp_build.py
  -l , --limit-threshold 
                        Sequence breadth needed to consider a genome detected. Increasing this value 
                        will increase false negative rate. Decreasing this value will increase false 
                        positive rate (float; default 0.1)
  -r , --percent-identity 
                        Minimum BLAST-like percent identity of alignment between read and reference 
                        genome (float; default 0.95)
  -q , --query-coverage 
                        Minimum fraction of read covered by an alignment between read and reference 
                        genome (float; default 0.7)
  -t , --threads        Threads available to SourceApp and its subroutines
  --use-geq             Report results normalized to genome equivalents
  --no-limits           Disable the analytical limit of detection used in estimating sequence depth.   
                        Synonymous with -l 0
  --skip-trimming       Disable read trimming and QC

```

# Citations

SourceApp wraps several essential pieces of software developed by other teams. 
Those finding SourceApp useful are encouraged to also cite these excellent works:
    
### MicrobeCensus

Nayfach, S.; Pollard, K. S. Average Genome Size Estimation Improves Comparative Metagenomics and Sheds Light on the Functional Ecology of the Human Microbiome. Genome Biology 2015, 16 (1), 51. https://doi.org/10.1186/s13059-015-0611-7.

### BWA-MEM2

Vasimuddin, Md.; Misra, S.; Li, H.; Aluru, S. Efficient Architecture-Aware Acceleration of BWA-MEM for Multicore Systems. In 2019 IEEE International Parallel and Distributed Processing Symposium (IPDPS); 2019; pp 314–324. https://doi.org/10.1109/IPDPS.2019.00041.

### Fastp

Chen, S. Ultrafast One‐pass FASTQ Data Preprocessing, Quality Control, and Deduplication Using Fastp. iMeta 2023, 2 (2), e107. https://doi.org/10.1002/imt2.107.

### CoverM

Woodcroft, B. https://github.com/wwood/CoverM

### CheckM2

Chklovski, A.; Parks, D. H.; Woodcroft, B. J.; Tyson, G. W. CheckM2: A Rapid, Scalable and Accurate Tool for Assessing Microbial Genome Quality Using Machine Learning. Nat Methods 2023, 20 (8), 1203–1212. https://doi.org/10.1038/s41592-023-01940-w.

### dRep

Olm, M. R.; Brown, C. T.; Brooks, B.; Banfield, J. F. dRep: A Tool for Fast and Accurate Genomic Comparisons That Enables Improved Genome Recovery from Metagenomes through de-Replication. ISME J 2017, 11 (12), 2864–2868. https://doi.org/10.1038/ismej.2017.126.
