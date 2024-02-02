#!/usr/bin/env python3

'''
SourceApp: Python implementation of the Unix-based environmental monitoring tool.

sourceapp.py requires a database directory that is properly configured:

SourceApp_db:
  database      # bwamem2 index of a concatenated FASTA file of input genomes
  gdef.tsv      # TSV, the name of the genome each contig came from in "database" ['genome', 'contig_name']
  sources.txt   # TSV, the name of each genome and its associated source ['genome', 'source']
'''
# =============================================================================
### Libraries:
import os
import sys
import argparse
import subprocess
import pandas as pd

### Core Functions:
def genome_qc(args):

    threads = args['threads']
    input_dir=args['input-dir'] # Note: what file extensions does checkm2 want? .fna? .fasta? .f*?
    output_dir=args['output-name'] + '_SourceAppdb'
    subprocess.run(["checkm2 predict --threads "+ str(threads)+" --input "+ input_dir + " --output-directory "+output_dir], shell=True, check=True)

def genome_derep(args):
    allowderep=args['no-dereplication']
    if allowderep:
        threads=args['threads']
        input_dir=args['input-dir']
        output_dir=args['output-name']+'_SourceAppdb'
        ani=args['ani']
        ginfo=1 #from checkm2
        subprocess.call(['dRep dereplicate '+ output_dir + " -g " + input_dir + " --S_ani " + ani + " --genomeInfo " + ginfo +
                        " --S_algorithm fastANI  -p " + str(threads)])
    else:
        copy_dir=1 # copy all input genomes into ...
        output_dir=1 #                    output directory that dRep would've used.
        ginfo=1
        copy_ginfo=1
        subprocess.call(['cp',copy_dir,output_dir])
        subprocess.call(['cp',ginfo,copy_ginfo])

def genome_selection(args):
    # toss out bad genomes
    quality = 1
    input_dir=1
    output_dir=1
    ginfo=1 # genome quality
    cinfo=1 # genomes sharing the same cluster but from different sources
    sinfo=1 # the source info associated with each genome

    # read ginfo into df
    df = pd.read_csv(ginfo, sep='\t',index_col=['name','comp','redun'])

    # get rid of genomes below quality threshold
    df['qual']= df['comp']-5*df['redun']
    good_df = df[df['qual']>quality]

    # add source info
    # loop that reads sinfo in the order of genomes in good_df['name'], and creates a list that it concatenates as a new column called ['source']

    # if allowing crx, just move the good genomes based on good_df into a new output_dir
    allowcrx=1
    if allowcrx:
        # if we don't care about cross-reactions between
        for row in range(len(good_df)):
            genome = good_df['name'].iloc[row]
            subprocess.call(['cp',input_dir,genome,output_dir])
    else:
        print('Removing cross reactive genomes')# if we do care about crx, then remove all genomes in the cinfo object and proceed
        ### finish me!

def build_database(args):
    cinfo=1
    sinfo=1
    input_dir=1

    ## step 1: rename all genome contigs

    ## step 2: create a list of the contig name and the genome it belongs to

    ## step 3: concatenate all renamed genomes

    input_fasta=1
    prefix="database" # users will provide output_dir so we fix prefix name so that sourceapp.py can find it.
    ## step 4: run bwa-mem2 indexing
    subprocess.call(["bwa-mem2 index", "-p", prefix, input_fasta]) # gotta pay attention to paths here

### Pipeline: 
def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--input-dir',
        help='Path to directory containing input genomes',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--output-name',
        help='Name of the database to be created. SourceApp will create an output directory in the current working directory containing the finished database with the provided string + "_SourceAppdb/"',
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
        '-q','--genome-quality',
        help='Aggregate quality score threshold for accepting input genomes',
        metavar='',
        type=float,
        default=0.95,
        required=False
        )
    parser.add_argument(
        '--remove-crx',
        help='Remove genomes found in the same cluster but belonging to different\
             sources',
        action='store_true',
        required=False
        )
    parser.add_argument(
        '--no-dereplication',
        help='Disable genome dereplication. This will create the database using all\
            of the provided genomes which pass quality requirements.',
        action='store_true',
        required=False
        )
    args=vars(parser.parse_args())

    # accession

    # read in the genomes from their input directory
    #genome_qc(args)
    # checkm2 outputs "quality_report.tsv" in the output_dir which is our ginfo:
    #file=args['output-name'] + '_SourceAppdb/quality_report.tsv'
    #df = read.csv(file, sep='\t') # we need to create a new file, comma-delimited with just header "genome,completeness,contamination" for dRep to read it properly.
    # so that means we need to strip off everything after the first 1-3 columns in the current output ("file") which we're reading in to pandas to manipulate.
    # we also need to give it a new header and replace tsv for csv. We'll do this when we write.
    #####
    ##### Start here

    # pass in all the genomes in the input directory to dRep but also pass ginfo with checkm2 quality scores
    # and dRep will take care of the rest
    genome_derep(args)

if __name__ == "__main__":
    main()
