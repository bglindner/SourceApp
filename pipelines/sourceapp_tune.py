#!/usr/bin/env python3
"""
SourceApp: Python implementation of the Unix-based environmental monitoring tool.
"""

import os
import sys
import argparse
import subprocess
import pandas as pd

### Core functions:

def read_filter(args):
    if args['output_dir'][-1] == '/': # in the event user provides trailing '/'
        bam = args['output_dir'][:-1] + '/mappings.bam'
    else:
        bam = args['output_dir'] + '/mappings.bam'
    if args['sourceapp_database'][-1] == '/': # in the event user provides trailing '/'
        gdef = args['sourceapp_database'][:-1] + '/gdef.txt' # this is a file describing which contigs belong to each genome as prepared by sourceapp_build.py
    else:
        gdef = args['sourceapp_database'] + '/gdef.txt'
    #gdef is a two col TSV with genome_name'\t'contig_name
    pid=args['percent_identity']
    qcov=args['query_coverage']
    threads=args['threads']
    trunc=args['limit_threshold']
    if args['output_dir'][-1] == '/': # in the event user provides trailing '/'
        output = args['output_dir'][:-1] + '/mappings_filtered.txt'
    else:
        output = args['output_dir'] + '/mappings_filtered.txt'
    nolimits=args['no_limits'] # if passed by the user, is true.
    usegeq=args['use_geq']
    print("Filtering read mapping results...")
    if trunc == 0 or nolimits: # if -l 0 or --no-limits was passed
        if usegeq:
            try:
                subprocess.run(["coverm genome -b "+bam+" --genome-definition "+gdef+" --min-read-percent-identity "+str(pid*100)+
                            " --min-read-aligned-percent "+str(qcov*100)+" --output-format dense -t "+str(threads)+" -m mean covered_bases variance "+
                            "-o "+output], shell=True, check=True)
            except Exception as e:
                print('Error in step 4: filtering of read mappings. Exiting . . .')
                sys.exit()
        else:
            try:
                subprocess.run(["coverm genome -b "+bam+" --genome-definition "+gdef+" --min-read-percent-identity "+str(pid*100)+
                            " --min-read-aligned-percent "+str(qcov*100)+" --output-format dense -t "+str(threads)+" -m relative_abundance "+
                            " -o "+output], shell=True, check=True)
            except Exception as e:
                print('Error in step 4: filtering of read mappings. Exiting . . .')
                sys.exit()
    else: # if -l is nonzero and --no-limits wasn't passed
        if usegeq:
            try:
                subprocess.run(["coverm genome -b "+bam+" --genome-definition "+gdef+" --min-read-percent-identity "+str(pid*100)+
                            " --min-read-aligned-percent "+str(qcov*100)+" --output-format dense -t "+str(threads)+" -m trimmed_mean covered_bases variance "+
                            " -o "+output+" --trim-min "+str(trunc)+" --trim-max "+str(100-(trunc*100))], shell=True, check=True)
            except Exception as e:
                print('Error in step 4: filtering of read mappings. Exiting . . .')
                sys.exit()
        else:
            try:
                subprocess.run(["coverm genome -b "+bam+" --genome-definition "+gdef+" --min-read-percent-identity "+str(pid*100)+
                            " --min-read-aligned-percent "+str(qcov*100)+" --output-format dense -t "+str(threads)+" -m relative_abundance "+
                            " -o "+output+" --trim-min "+str(trunc)+" --trim-max "+str(100-(trunc*100))], shell=True, check=True)
            except Exception as e:
                print('Error in step 4: filtering of read mappings. Exiting . . .')
                sys.exit()

def summarize(args):
    #produce the final dataframe and make some visuals.
    usegeq=args['use_geq']
    sourcedict = read_source(args['sourceapp_database'].replace('/','') + '/sources.txt')
    sources = list(set(sourcedict.values()))
    if args['output_dir'][-1] == '/': # in the event user provides trailing '/'
        df = pd.read_csv(args['output_dir'][:-1] + '/mappings_filtered.txt', header=0, sep='\t')
    else:
        df = pd.read_csv(args['output_dir'] + '/mappings_filtered.txt', header=0, sep='\t')
    portions=[]
    if usegeq:
        geq = get_geq(args)
        for source in sources:
            gsum=0
            glist = [key for key, val in sourcedict.items() if val == source]
            for genome in glist:
                gsum = gsum + df[df['Genome']==genome].iloc[:,1].sum()
            geq_frac = gsum/geq
            portions.append([source,geq_frac])
        sourceApp = pd.DataFrame(portions, columns=['Source','Portion'])
        if sourceApp['Portion'].sum() > 1:
            sourceApp['Portion'] = sourceApp['Portion']/sourceApp['Portion'].sum()
            print('Warning: the sum of GEQ-based relative abundances exceeds 1. Source portions have been rescaled.')
            print('It is recommended to re-run SourceApp without the --use-geq flag to examine what percentage of reads are recovered. If the value is <~90%, then GEQ-based normalization may not be robust for this dataset.')
    else:
        for source in sources: # if we don't normalize to GEQ, then we should just report DNA relabd.
            gsum=0
            glist = [key for key, val in sourcedict.items() if val == source]
            for genome in glist:
                gsum = gsum + df[df['Genome']==genome]['Relative Abundance (%)'].sum()
            portions.append([source,gsum])
        portions = pd.DataFrame(portions, header=['Source', 'Portion'])
    return portions

### Helper functions:
def read_source(file): # helper function to read in the source dictionary stored by sourceapp_build.py
    dic = {}
    with open(file) as fh:
        for line in fh:
            key, val = line.rstrip().split()
            dic[key] = val
    return dic

def get_geq(args):
    if args['output_dir'][-1] == '/': # in the event user provides trailing '/'
        file = args['output_dir'][:-1] + '/geq.txt'
    else:
        file = args['output_dir'] + '/geq.txt'
    with open(file) as fh:
        for index, line in enumerate(fh):
            if index == 12:
                censusline = line.split()
    output = float(censusline[1])
    return output

### Pipeline:
def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--input-dir',
        help='Path to an output directory from a previous "sourceapp.py" run',
        metavar='',
        type=str, 
        required=True
        )
    parser.add_argument(
        '-o', '--output-file',
        help='Path to the desired output directory',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-d', '--sourceapp-database',
        help='Path to directory containing a SourceApp formatted database. Default database available for download or produced de novo as the output directory from sourceapp_build.py',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-l', '--limit-threshold',
        help='Sequence breadth needed to consider a genome detected. Increasing this value will increase false negative rate. Decreasing this value will increase false positive rate (float; default 0.1)',
        metavar='',
        type=float,
        default=0.1,
        required=False
        )
    parser.add_argument(
        '-r', '--percent-identity',
        help='Minimum BLAST-like percent identity of alignment between read and reference genome (float; default 0.95)',
        metavar='',
        type=float,
        default=0.95,
        required=False
        )
    parser.add_argument(
        '-q', '--query-coverage',
        help='Minimum fraction of read covered by an alignment between read and reference genome (float; default 0.7)',
        metavar='',
        type=float,
        default=0.7,
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
        '--use-geq',
        help='Report results normalized to genome equivalents',
        action='store_true',
        required=False
        )
    parser.add_argument(
        '--no-limits',
        help='Disable the analytical limit of detection used in estimating sequence depth. Synonymous with -l 0',
        action='store_true',
        required=False
        )
    args=vars(parser.parse_args())

    print('Beginning step 1: filtering of read mappings')
    read_filter(args)

    print('Beginning step 2: results summarization')
    results = summarize(args)
    print(results)
    print(results.__class__)
    output_table = pd.DataFrame(results)

    if args['output_dir'][-1] == '/': # in the event user provides trailing '/'
        output_table.to_csv(args['output_dir'][-1]+'/results.csv', index=False, header=["Source","Portion"])
    else:
        output_table.to_csv(args['output_dir']+'/results.csv', index=False, header=["Source","Portion"])

    print('The following results printed to results.csv in output directory:')
    print(output_table)

    print('Thank you for using SourceApp.')

if __name__ == "__main__":
    main()
