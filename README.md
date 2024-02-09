      ________                                   _____                 
     /   ____/ ____  __ _________  ____  ____   /  _  \ ______ ______  
     \____  \ /  _ \|  |  \_  __ \/ ____/ __ \ /  /_\  \  __ \\  __ \ 
     /       (  |_| |  |  /|  | \\  \__\  ___//   /|   |  |_ ||  |_ |
    /________/\____/|____/ |__|   \_____\_____\___||__ |   __/|   __/ 
                                                      \|__|   |__|   

Python implementation of the Unix-based environmental monitoring tool.

SourceApp is in active development. Check back here for updates and the eventual release.

# Description 
SourceApp is a bioinformatic workflow designed to apportion fecal signal amongst multiple competing sources in short read metagenomes collected from impaired waterways. SourceApp was developed for use on Unix or Unix-like operating systemsand its implementation in Python has the same requirement.

SourceApp is designed to automate all tasks necessary to detect and quantify fecal signal in metagenomes collected from the water environment if users will simply provide gzipped copies of the raw reads in FASTQ format. SourceApp requires adatabase of reference genomes with known source-associations. Consequently, a database was developed in parallel with this software through systematic literature review and will be made available upon completion of that work (Graham et al, _in prep_)

SourceApp reports metrics akin to prokaryotic cell fraction and has been benchmarked in laboratory trials. SourceApp was found to perform best when fecal signal is relatively fresh and when the user is able to provide supplementary reference genomes to bolster SourceApp's default genomic reference database. De novo databases can be constructed by the user using sourceapp_build.py and then supplied here to SourceApp in place of the default database. 

For functional population proportioning in biological wastewater treatment datasets, see "WasteApp" -- the built envrionment extension of SourceApp's working principles (https://github.com/blindner6/WasteApp). 

# Installation

We are currently recommending installation with `conda` / `mamba` and aim to create a Bioconda recipe for easier user installation after publication. In the interim, the following guide should allow users to set up two environments for running SourceApp. One environment is for the main pipeline of SourceApp (`sourceapp.py`) and the other is for building your own databases (`sourceapp_build.py`). 
```
# clone this repo
git clone https://github.com/blindner6/SourceApp.git

# create two environments using the provided .yml files, one supporting the 
# dependencies for "sourceapp.py" and the other for "sourceapp_build.py"
mamba create -n sourceapp --file path/to/SourceApp/dependencies/sourceapp.yml
mamba create -n sourceapp_build --file path/to/SourceApp/dependencies/sourceapp_build.yml

# pip install one additional dependency
mamba activate sourceapp
pip install MicrobeCensus-SourceApp

# test your environments:
micromamba activate sourceapp
python path/to/SourceApp/pipelines/sourceapp.py -h

micromamba activate sourceapp_build
python path/to/SourceApp/pipelines/sourceapp_build.py -h

```
# Usage

`sourceapp.py` expects as input paired short read metagenomic data. The user should supply these reads as gzipped FASTQ files. No prior adapter trimming or QC is necessary as SourceApp will automate read trimming with `fastp` but the user can disable this step with `--skip-trimming`. In addition to short reads, SourceApp needs a genomic database specifically formatted for use by the tool. This can be created by SourceApp from genomes provided by the user with `sourceapp_build.py` or a default database will be available here soon. Users should specify the location of an output directory for results to be written to. 

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

# Database construction

SourceApp performs best when users are able to supply genomes recovered from the contaminating sources the user expects to be present. `sourceapp_build.py` allows users to provide a set of genomes which they have collected/curated and create a SourceApp database for use with the main pipeline. To do this, users should gather genomes as FASTA files in an input directory and record in a tab-separated list, the name of each genome (col1) and its fecal source (col2). From these two inputs, SourceApp will output a directory containing a database which can be passed into `sourceapp.py -d`.

```
usage: sourceapp_build.py [-h] -i  -o  -s  [-a] [-t] [-q] [--remove-crx] [--no-dereplication]

SourceApp: Python implementation of the Unix-based environmental monitoring tool.

sourceapp.py requires a database directory that is properly configured:

SourceApp_db:
  database      # bwamem2 index of a concatenated FASTA file of input genomes
  gdef.tsv      # TSV, the name of the genome each contig came from in "database" ['genome', 'contig_name']
  sources.txt   # TSV, the name of each genome and its associated source ['genome', 'source']

options:
  -h, --help            show this help message and exit
  -i , --input-dir      Path to directory containing input genomes
  -o , --output-name    Name of the database to be created. SourceApp will create an output directory in the current working
                        directory containing the finished database with the provided string + "_SourceAppdb/"
  -s , --source-associations 
                        Text file describing source associations of input genomes
  -a , --ani            ANI threshold for calling genome clusters
  -t , --threads        Threads available to SourceApp
  -q , --genome-quality 
                        Aggregate quality score threshold for accepting input genomes
  --remove-crx          Remove genomes found in the same cluster but belonging to different sources
  --no-dereplication    Disable genome dereplication. This will create the database using all of the provided genomes which pass
                        quality requirements.
```
# Pre-built database
To assist users with achieving good performance from SourceApp, we have built a default database which users can download and immediately run. This database includes entries representing the following fecal sources:
```
source > n > note
bird > x > non-domesticated
cat > x > domesticated
chicken > x > domesticated
cow > x > domesticated; bos *
dog > x > domesticated
human > x >
pig > x > domesticated 
ruminant > x > mostly domesticated (goats and sheep)
septage > x > 
wastewater > x > genomes from collection system, raw influent, and secondary
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
