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

For functional population proportioning in biological wastewater treatment datasets, see "WasteApp" the built envrionment extension of SourceApp (https://github.com/blindner6/WasteApp). 

# Installation

The development version of SourceApp is available here for those interested in testing the software before its publication. It is anticipated that at the time of publication, we will create and maintain a Bioconda recipe for direct installation with `conda` or `mamba`. In the interim, we suggest the following for obtaining a working version of the tool:
```
# clone this repo
git clone https://github.com/blindner6/SourceApp.git

# create two environments using the provided .yml files, one supporting the dependencies for "sourceapp.py" and the other for "sourceapp_build.py"
mamba create -n sourceapp
mamba env update -n sourceapp --file SourceApp/sourceapp.yml
mamba create -n sourceapp_build
mamba env update -n sourceapp_build --file SourceApp/sourceapp_build.yml

# when calling sourceappy.py or sourceapp_build.py, be sure to activate the environment corresponding to your task:

```
# Usage

In prep

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

