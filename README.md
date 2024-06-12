      ________                                   _____                 
     /   ____/ ____  __ _________  ____  ____   /  _  \ ______ ______  
     \____  \ /  _ \|  |  \_  __ \/ ____/ __ \ /  /_\  \  __ \\  __ \ 
     /       (  |_| |  |  /|  | \\  \__\  ___//   /|   |  |_ ||  |_ |
    /________/\____/|____/ |__|   \_____\_____\___||__ |   __/|   __/ 
                                                      \|__|   |__|   

Python implementation of the Unix-based environmental monitoring tool.

SourceApp is in active development but we encourage interested users to give it a try ahead of publication.

# Description 
SourceApp is a bioinformatic workflow designed to apportion fecal signal amongst multiple competing sources in short read metagenomes collected from impaired waterways. SourceApp was developed for use on Unix or Unix-like operating systems and its implementation in Python has the same requirement.

SourceApp is designed to automate all tasks necessary to detect and quantify fecal signal in metagenomes collected from the water environment if users will simply provide gzipped copies of the raw reads in FASTQ format. SourceApp requires a database of reference genomes with known source-associations. Consequently, a database was developed in parallel with this software through systematic literature review and details on obtaining it are available below (see also Lindner et al, _in prep_ which describes: https://doi.org/10.5281/zenodo.10728776).

SourceApp reports metrics akin to prokaryotic cell fraction (when `--use-geq` is called) and has been benchmarked in laboratory trials. SourceApp was only tested on fresh fecal signal and works best when the user is able to provide supplementary reference genomes to bolster SourceApp's default genomic reference database. De novo databases can be constructed by the user with `sourceapp_build.py` and then supplied here to SourceApp in place of the default database. 

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
pip install MicrobeCensus-SourceApp==1.1.2

# test your environments:
python path/to/SourceApp/pipelines/sourceapp.py -h
mamba deactivate

mamba activate sourceapp_build
python path/to/SourceApp/pipelines/sourceapp_build.py -h
```
If you have trouble solving the environment, first check to see if your channel priority is flexible and change if needed:
```
# see channel_priority: (--list in conda)
mamba config list

# change as needed: (--set in conda)
mamba config set channel_priority flexible
```
Both pipelines have a fair amount of dependencies and depending on your system architecture and its disk space limits, you may want to clean up some of the installation files before moving on:
```
mamba clean --all
```
# Usage

`sourceapp.py` expects as input paired short read metagenomic data. The user should supply these reads as gzipped FASTQ files. No prior adapter trimming or QC is necessary as SourceApp will automate read trimming with `fastp` but the user can disable this step with `--skip-trimming`. In addition to short reads, SourceApp needs a genomic database specifically formatted for use by the tool. This can be created by SourceApp from genomes provided by the user with `sourceapp_build.py` or using the default database described below. Users should specify the location of an output directory for results to be written to. 

SourceApp is primarily designed for use with a Unix-based HPC and no support is offered for deployment with alternative operating systems.

```
usage: sourceapp.py [-h] -i  -o  -d  [-l] [-f] [-r] [-q] [-t] [--use-geq] [--aggregate-human] [--no-limits] [--skip-trimming] [--drop-env]

SourceApp: Python implementation of the Unix-based environmental monitoring tool.

options:
  -h, --help            show this help message and exit
  -i , --input-files    Comma-delimited path to forward and reverse metagenomic reads. Must be in FASTQ format and compressed with gzip
  -o , --output-dir     Path to the desired output directory
  -d , --sourceapp-database 
                        Path to directory containing a SourceApp formatted database. Default database available for download or produced de novo as the output directory
                        from sourceapp_build.py
  -l , --limit-threshold 
                        Sequence breadth needed to consider a genome detected. Increasing this value will increase false negative rate. Decreasing this value will
                        increase false positive rate (float; default 0.1)
  -f , --min-frac       The minimum read or cell fraction a source must have to be considered detected and therefore apportionable (float; default 0.0001)
  -r , --percent-identity 
                        Minimum BLAST-like percent identity of alignment between read and reference genome (float; default 0.95)
  -q , --query-coverage 
                        Minimum fraction of read covered by an alignment between read and reference genome (float; default 0.7)
  -t , --threads        Threads available to SourceApp
  --use-geq             Report results normalized to genome equivalents
  --aggregate-human     Treat human signal as wastewater signal. This will result in specific human signal being treated as cross-reactive wastewater signal. It alone
                        cannot be used for attribution but it will contribute to apportioning if non-crossreactive wastewater signal was detected.
  --no-limits           Disable the analytical limit of detection used in estimating sequence depth. Synonymous with -l 0
  --skip-trimming       Disable read trimming and QC
  --drop-env            Discard environmental signal from final results. This can significantly impact apportioning results -- you probably want it on if apportioning
                        fecal contamination.

```

For example:
```
python sourceapp.py --use-geq --drop-env --aggregate-human -i raw_reads/reads.1.fastq.gz,raw_reads/reads.2.fastq.gz -o sourceapp_results -d path/to/SourceApp_db
```

# Outputs

SourceApp aims to provide information on three metrics: attribution, apportionment, and fractioning. If SourceApp was ran with the following:

```
python sourceapp.py -i path/to/reads.1.fastq.gz,path/to/reads.2.fastq.gz -o path/to/output_dir -d path/to/database_dir --use-geq --aggregate-human --drop-env
```

`sourceapp.py` will output a table of for each result (as a `.csv`) with similar formats. Users should note the usage guidelines of the tool, as the results described in these tables can vary depending on parameter choice and flags passed. 

## Attribution

The attribution table shows each source described in the database as the first column (`Source`) and that source's presence/absence in the second column (`Detection`).

```
Source,Detection
bird,0
cat,1
chicken,0
cow,0
dog,0
pig,1
ruminant,0
septage,0
wastewater,1
```

## Apportionment

The apportionment table shows each source described in the database as the first column (`Source`) and, as the second column (`Portion`) the portion of the total fecal signal SourceApp predicted that source contributed. All attributed sources will have non-zero portions. 

```
Source,Portion
bird,0.0
cat,0.600
chicken,0.0
cow,0.0
dog,0.0
pig,0.035
ruminant,0.0
septage,0.0
wastewater,0.365
```

## Fractioning

The fractioning table shows each source described in the database as the first column (`Source`) and then that source's predicted fraction within the sample (as c/c if `--use-geq` was passed) in the second column (`Fraction`). 

```
Source,Fraction
bird,0.0
cat,0.032
chicken,0.0
cow,0.0
dog,0.0
pig,0.002
ruminant,0.0
septage,0.0
wastewater,0.020
```

# Database construction

SourceApp performs best when users are able to supply genomes recovered from the contaminating sources the user expects to be present. `sourceapp_build.py` allows users to provide a set of genomes which they have collected/curated and create a SourceApp database for use with the main pipeline. To do this, users should gather genomes as FASTA files in an input directory and record in a tab-separated list, the name of each genome (col1) and its fecal source (col2). From these two inputs, SourceApp will output a directory containing a database which can be passed into `sourceapp.py -d`.

```
usage: sourceapp_build.py [-h] -i  -o  -s  [-a] [-t] [-q] [--no-dereplication] [-d] [-c]

SourceApp: Python implementation of the Unix-based environmental monitoring tool.

sourceapp.py requires a properly formatted database. This script automates its creation.

options:
  -h, --help            show this help message and exit
  -i , --input-dir      Path to directory containing input genomes (path/to/dir/*.fna)
  -o , --output-name    Path to new directory which SourceApp should build the database in. SourceApp will create the database output directory at the specified path
                        appended with + '_SourceAppdb/'
  -s , --source-associations 
                        Text file describing source associations of input genomes
  -a , --ani            ANI threshold for calling genome clusters
  -t , --threads        Threads available to SourceApp
  -q , --genome-quality 
                        Aggregate quality score threshold for accepting input genomes (float, 50 default)
  --no-dereplication    Disable genome dereplication. This will create the database using all of the provided genomes which pass quality requirements.
  -d , --checkm2_db     Path to a local installation of the CheckM2 database (.dmnd). If not passed, SourceApp assumes you have let CheckM2 install the database in the
                        default location (~/databases). See 'checkm2 databases -h' for more information.
  -c , --checkm2_info   If you've already run CheckM2 yourself, path to the quality_report.tsv output.
```
`sourceapp_build.py` expects the tab-separated input source association table to be structured like so (no header):
```
genome1.fna      sourceA
genome2.fna      sourceA
genome3.fna      sourceB
genome4.fna      sourceC
...
genomeN.fna      sourceN
```
For example:
```
python sourceapp_build.py -i path/to/genomes/dir -o path/to/output_dir -s path/to/source_associations.tsv -d path/to/CheckM2_db/*.dmnd
```

The output database created by `sourceapp_build.py` is the database input into `sourceapp.py`.


# Database structure

# Default database
To assist users with achieving good performance from SourceApp, we have built a default database which users can download and immediately use in their runs (or to augment their own genomic datasets i.e., with `sourceapp_build.py`). This database includes entries representing the following fecal sources:

| source | specific | crx | note |
|--------|---|------|
| bird | 40 | 16 | non-domesticated |
|cat | 6 | 91 | domesticated |
|chicken | 883 | 431 | domesticated |
|cow | 19 | 24 | domesticated; bos * |
|dog | 58 | 89 | domesticated |
|human | 3333  | 1151 | |
|pig | 2038 | 728 | domesticated |
|ruminant | 725 | 33 |  mostly domesticated (goats and sheep) |
|septage | 32 | 7 |  |
|wastewater | 3076 | 47 | primary and secondary |
|environmental | 6559 | 43 | freshwater |

# Citations

SourceApp wraps several essential pieces of software developed by other teams. 
Those finding SourceApp useful are encouraged to also cite these excellent works:
    
### MicrobeCensus

Nayfach, S.; Pollard, K. S. Average Genome Size Estimation Improves Comparative Metagenomics and Sheds Light on the Functional Ecology of the Human Microbiome. Genome Biology 2015, 16 (1), 51. https://doi.org/10.1186/s13059-015-0611-7.

### bwa mem

Li, H.; Durbin, R. Fast and Accurate Short Read Alignment with Burrows-Wheeler Transform. Bioinformatics 2009, 25 (14), 1754–1760. https://doi.org/10.1093/bioinformatics/btp324.

### Fastp

Chen, S. Ultrafast One‐pass FASTQ Data Preprocessing, Quality Control, and Deduplication Using Fastp. iMeta 2023, 2 (2), e107. https://doi.org/10.1002/imt2.107.

### CoverM

Woodcroft, B. https://github.com/wwood/CoverM

### CheckM2

Chklovski, A.; Parks, D. H.; Woodcroft, B. J.; Tyson, G. W. CheckM2: A Rapid, Scalable and Accurate Tool for Assessing Microbial Genome Quality Using Machine Learning. Nat Methods 2023, 20 (8), 1203–1212. https://doi.org/10.1038/s41592-023-01940-w.

### dRep

Olm, M. R.; Brown, C. T.; Brooks, B.; Banfield, J. F. dRep: A Tool for Fast and Accurate Genomic Comparisons That Enables Improved Genome Recovery from Metagenomes through de-Replication. ISME J 2017, 11 (12), 2864–2868. https://doi.org/10.1038/ismej.2017.126.
