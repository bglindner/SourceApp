      _________                                   _____                 
     /   _____/ ____  __ _________  ____  ____   /  _  \ ______ ______  
     \_____  \ /  _ \|  |  \_  __ \/ ____/ __ \ /  /_\  \____ \\____ \ 
     /        (  |_| |  |  /|  | \\  \__\  ___//   /|   |  |_ ||  |_ |
    /_______  /\____/|____/ |__|   \___  \___  \___||__ |   __/|   __/ 
            \/                         \/    \/        \|__|   |__|    

Python implementation of the Unix-based environmental monitoring tool.

# Description 
SourceApp is a bioinformatic workflow designed to partition fecal signal amongst
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

SourceApp wraps several essential pieces of software developed by other teams. 
Those finding SourceApp useful are encouraged to also cite these excellent works:
    
::MicrobeCensus::
    Nayfach, S.; Pollard, K. S. Average Genome Size Estimation Improves 
Comparative Metagenomics and Sheds Light on the Functional Ecology of the Human 
Microbiome. Genome Biology 2015, 16 (1), 51. 
https://doi.org/10.1186/s13059-015-0611-7.

::BWA-MEM2::
    Vasimuddin, Md.; Misra, S.; Li, H.; Aluru, S. Efficient Architecture-Aware 
Acceleration of BWA-MEM for Multicore Systems. In 2019 IEEE International Parallel 
and Distributed Processing Symposium (IPDPS); 2019; pp 314–324. 
https://doi.org/10.1109/IPDPS.2019.00041.

::Fastp::
    Chen, S. Ultrafast One‐pass FASTQ Data Preprocessing, Quality Control, and 
Deduplication Using Fastp. iMeta 2023, 2 (2), e107. 
https://doi.org/10.1002/imt2.107.

::CoverM::
    Woodcroft, B. https://github.com/wwood/CoverM
