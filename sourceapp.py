#!/usr/bin/env python3

'''
  _________                                   _____                 
 /   _____/ ____  __ _________  ____  ____   /  _  \ ______ ______  
 \_____  \ /  _ \|  |  \_  __ \/ ____/ __ \ /  /_\  \____ \\____ \ 
 /        (  |_| |  |  /|  | \\  \__\  ___//   /|   |  |_ ||  |_ |
/_______  /\____/|____/ |__|   \___  \___  \___||__ |   __/|   __/ 
        \/                         \/    \/        \|__|   |__|    

Python implementation of the Unix-based environmental monitoring tool.



sourceapp.py -i -o -d (-p -t -m --use-geq --no-limits)

Req'd: -i <string> -o <string> -d <string>
Opt'l: -p <float> -t <int> -m <int> --use-geq --no-limits

'''
#%%
### Libraries:
import argparse
import pandas as pd
import textdistance as td


# =============================================================================
### Functions:
#%%
# 1. validate_inputs
def validate_inputs(var1, var2):
    
#%%  
# 2. accession
def accession():
    
#%%    
# 3. read_trim
def read_trim():
    
#%% 
# 4. geq_estimation
def geq_estimation():
    
#%%
# 5. read_map
def read_map():
    
#%%    
# 6. read_filter
def read_filter():
    
#%%    
# 7. summarize
def summarize():
    return output

# =============================================================================

#%%
### Workflow arguments:  
# =============================================================================
def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--input-file',
        help='Path to metagenomic reads. These are expected to be paired reads\
             in FASTQ format and compressed with gzip)',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--output-file',
        help='Desired output file (text, tab-delim)',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-d', '--sourceapp-database',
        help='Path to directory containing a SourceApp formatted database. \
            Default database available for download or produced de novo as the \
            output directory from sourceapp_build.py',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-l', '--limit-threshold',
        help='Sequence breadth needed to consider a genome detected. Increasing\
             this value will increase false negative rate. Decreasing this value will\
             increase false positive rate',
        metavar='',
        type=float,
        default=0.1
        required=False
        )
    parser.add_argument(
        '-r', '--percent-identity',
        help='Minimum BLAST-like percent identity of alignment between read and reference genome',
        metavar='',
        type=float,
        default=0.95
        required=False
        )
    parser.add_argument(
        '-q', '--query-coverage',
        help='Minimum fraction of read covered by an alignment between read and reference genome ',
        metavar='',
        type=float,
        default=0.7
        required=False
        )
    parser.add_argument(
        '-t','--threads',
        help='Threads available to SourceApp',
        metavar='',
        type=int,
        default=1,
        required=False
        )
    parser.add_argument(
        '-m','--memory',
        help='Memory available to SourceApp (in gigabytes)',
        metavar='',
        type=int,
        default=32,
        required=False
        )
    parser.add_argument(
        '--use-geq',
        help='Report results normalized to genome equivalents',
        metavar='',
        action='store_true',
        required=False
        )
    parser.add_argument(
        '--no-limits',
        help='Disable the analytical limit of detection used in estimating \
            sequence depth. Synonymous with -l 0',
        metavar='',
        action='store_true',
        required=False
        )
    args=vars(parser.parse_args())
# =============================================================================
#%%
# Workflow:
    #args['threads']
    print('Running SourceApp...')

if __name__ == "__main__":
    main()
