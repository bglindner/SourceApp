#!/usr/bin/env python3

'''
  _________                                   _____                 
 /   _____/ ____  __ _________  ____  ____   /  _  \ ______ ______  
 \_____  \ /  _ \|  |  \_  __ \/ ____/ __ \ /  /_\  \____ \\____ \ 
 /        (  |_| |  |  /|  | \\  \__\  ___//   /|   |  |_ ||  |_ |
/_______  /\____/|____/ |__|   \___  \___  \___||__ |   __/|   __/ 
        \/                         \/    \/        \|__|   |__|    

Python implementation of the Unix-based environmental monitoring tool.



sourceapp_build.py -i -o -d (-p -t -m --use-geq --no-limits)

Req'd: -i <string> -o <string> -d <string>
Opt'l:

'''

### Libraries:
import argparse
import pandas as pd
import textdistance as td

### Functions:
def translate_gene_ids(file, clusters):
    return genome

def score_synteny(g1, g2, n):

### Workflow: 
#%%
def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--input-dir',
        help='Path to directory containing genomes',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--output-dir',
        help='Path to desired output directory',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-s', '--source-associations',
        help='Text file describing source associations of input genomes',
        metavar='',
        type=str,
        required=True,
        )
    parser.add_argument(
        '-a','--ani',
        help='ANI threshold for calling genome clusters',
        metavar='',
        type=float,
        default=0.95,
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
        '--remove-crx',
        help='Remove genomes found in the same cluster but belonging to different\
             sources',
        metavar='',
        action='store_true',
        required=False
        )
    parser.add_argument(
        '--no-dereplication',
        help='Disable genome dereplication. This will create the database using all\
            of the provided genomes which pass quality requirements.',
        metavar='',
        action='store_true',
        required=False
        )
    args=vars(parser.parse_args())
#%%
   
   # collect the cluster information provided by the user
    clust = pd.read_csv(args['cluster_info'],sep='\t',names=['rep','child'])
           
    print('Running SourceApp...')

if __name__ == "__main__":
    main()
