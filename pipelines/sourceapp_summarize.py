#!/usr/bin/env python3
"""
Helper script that will re-summarize a completed SourceApp run if provided the output directory.
This will overwrite existing SourceApp results in the given directory.
"""

import os
import sys
import argparse
import subprocess
import pandas as pd

### Core functions:

def summarize(args):
    #produce the final dataframe and make some visuals.
    usegeq=args['use_geq']
    sourcedict = pd.read_csv(args['sourceapp_database'].replace('/','') + '/sources.txt',sep="\t",header=0).iloc[:,0:2].set_index('genome')['source'].to_dict()
    sources = sorted(list(set(sourcedict.values())))
    df = pd.read_csv(args['output_dir'] + '/mappings_filtered.txt', header=0, sep='\t')
    portions=[]
    if usegeq:
        geq = get_geq(args)
        for source in sources:
            gsum=0
            gcount=0
            nount=0
            glist = [key for key, val in sourcedict.items() if val == source]
            for genome in glist:
                if df[df['Genome']==genome].iloc[:,1].sum()>0:
                    gcount = gcount + 1
                ncount=ncount+1
                gsum = gsum + df[df['Genome']==genome].iloc[:,1].sum()
            geq_frac = gsum/geq
            fractions.append([source,geq_frac,gcount,ncount])
        fractions = pd.DataFrame(fractions, columns=['Source','Fraction','Detected Genomes','Total Genomes'])
        if fractions['Fraction'].sum() > 1:
            fractions['Fraction'] = fractions['Fraction']/fractions['Fraction'].sum()
            print('Warning: the sum of GEQ-based relative abundances exceeds 1. Source portions have been rescaled.', flush=True)
            print('It is recommended to re-run SourceApp without the --use-geq flag to examine what percentage of reads are recovered. If the value is <~90%, then GEQ-based normalization may not be robust for this dataset.', flush=True)
    else:
        for source in sources: # if we don't normalize to GEQ, then we should just report DNA relabd. checkM reports relative abundances as ff.fff
            gsum=0
            glist = [key for key, val in sourcedict.items() if val == source]
            for genome in glist:
                if df[df['Genome']==genome].iloc[:,1].sum()>0:
                    gcount = gcount + 1
                ncount=ncount+1
                gsum = gsum + (df[df['Genome']==genome].iloc[:,1].sum())
            fractions.append([source,gsum,gcount,ncount])
        fractions = pd.DataFrame(fractions, columns=['Source', 'Fraction','Detected Genomes','Total Genomes'])
    return fractions

def clean_output(table, thresh, args):
    thresh =  args['min_frac']
    addhum = args["aggregate_human"] # if true, add human signal to wastwater
    sources = table.index[~table.index.str.contains("_crx")]
    
    df_att = table.loc[sources]
    df_att.drop("environmental",inplace=True)
    df_att.where(df_att >= thresh, 0, inplace=True)
    df_app = df_att.copy()
    df_att.where(df_att <= 0, 1, inplace=True)
    
    for source in sources.drop("environmental"):
        if not df_app.loc[source].iloc[0] == 0:
            df_app.loc[source].iloc[0] = df_app.loc[source].iloc[0] + table.loc[source+"_crx"].iloc[0]
            
    if addhum:
        if df_app.loc["wastewater"] > 0:
            df_app.loc["wastewater"] = df_app.loc["wastewater"] + df_app.loc["human"]
        df_app.drop("human",inplace=True)
        
    df_frac = df_app.copy()
    df_app = df_app/df_app.sum()
    
    return df_app, df_att, df_frac

### Pipeline:
def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-o', '--output-dir',
        help='Path to an existing SourceApp output directory',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-f', '--min-frac',
        help='The minimum read or cell fraction a source must have to be considered detected and therefore apportionable (float; default 0.0001)',
        metavar='',
        type=float,
        default=0.0001,
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
        '--aggregate-human',
        help='Treat human signal as wastewater signal. This will result in specific human signal being treated as cross-reactive wastewater signal. It alone cannot be used for attribution but it will contribute to apportioning if non-crossreactive wastewater signal was detected.',
        action='store_true',
        required=False
        )
    args=vars(parser.parse_args())

    if args['output_dir'][-1] == '/': # in the event user provides trailing '/'
        args['output_dir'] = args['output_dir'][:-1]
        
    if os.path.isdir(args['output_dir']):
        print('SourceApp output directory found. Searching for read mapping results. . .', flush=True)
        if os.path.isfile(args['output_dir']+'/mappings_filtered.txt'):
            print('Filtered read mapping results found. Proceeding.', flush=True)
        else:
            print('Cannot locate filtered read mapping results. Did SourceApp complete successfully?', flush=True)
            print('Exiting.', flush=True)
    else:
        print('SourceApp output directory not found. Exiting.', flush=True)
        sys.exit()

    print('Summarizing run SourceApp results.', flush=True)
    output_table = summarize(args)
    output_table.set_index("Source", inplace=True)

    app, att, frac = clean_output(output_table, args)

    output_table.to_csv(args['output_dir']+'/raw_results.csv', index=True, header=["Fraction"])
    att.to_csv(args['output_dir']+'/attributions.csv', index=True, header=["Detection"])
    app.to_csv(args['output_dir']+'/apportions.csv', index=True, header=["Portion"])
    frac.to_csv(args['output_dir']+'/fractions.csv', index=True, header=["Fraction"])

    print('The following results printed to raw_results.csv in output directory:', flush=True)
    print(output_table)

    print('Thank you for using SourceApp.', flush=True)

if __name__ == "__main__":
    main()
