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
    sourcedict = pd.read_csv(args['sourceapp_database'] + '/sources.txt',sep="\t",header=0).iloc[:,0:2].set_index('genome')['source'].to_dict()
    sources = sorted(list(set(sourcedict.values())))
    df = pd.read_csv(args['output_dir'] + '/mappings_filtered.txt', header=0, sep='\t')
    fractions=[]
    if usegeq:
        geq = get_geq(args)
        for source in sources:
            gsum=0
            gcount=0
            ncount=0
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
            gcount=0
            ncount=0
            glist = [key for key, val in sourcedict.items() if val == source]
            for genome in glist:
                if df[df['Genome']==genome].iloc[:,1].sum()>0:
                    gcount = gcount + 1
                ncount=ncount+1
                gsum = gsum + (df[df['Genome']==genome].iloc[:,1].sum())
            fractions.append([source,gsum,gcount,ncount])
        fractions = pd.DataFrame(fractions, columns=['Source', 'Fraction','Detected Genomes','Total Genomes'])
    return fractions

def clean_output(table, args):
    sources = table.index[~table.index.str.contains("_crx")]
    
    df = table.loc[sources]

    dlst=[]
    tlst=[]
    flst=[]
    for source in sources:
        flst.append(table.loc[source+"_crx"]["Fraction"])
        dlst.append(table.loc[source+"_crx"]["Detected Genomes"])
        tlst.append(table.loc[source+"_crx"]["Total Genomes"])
    df["Crx Fraction"] = flst
    df["Crx Detected Genomes"] = dlst
    df["Crx Total Genomes"] = tlst

    df["Attributal"] = df["Fraction"] # Attribution is based on specific signal.
    df["Attributal"].where(df["Attributal"] >= args["min_frac"], 0, inplace=True) # Apply limit
    df["Crx Attributal"] = df["Crx Fraction"] 
    df["Crx Attributal"].where(df["Crx Attributal"] >= args["min_frac"], 0, inplace=True) # Apply limit
    
    df["Portion"] = df["Fraction"]
    df["Portion"].where(df["Portion"] >= args["min_frac"], 0, inplace=True) # apply limit
    for source in sources:
        if df.loc[source]["Attributal"] > 0:
            df.loc[source]["Portion"] = df.loc[source]["Attributal"] + df.loc[source]["Crx Attributal"]
        else:
            df.loc[source]["Portion"] = 0
    if args["drop_env"]: # default is false; thus, if the keyword "environmental" doesn't appear in the index (because a custom db is used) all is good
        df.drop("environmental",inplace=True
    df["Portion"] = df["Portion"]/df["Portion"].sum()

    ###  clean up
    if args["aggregate_human"]: # default is false; thus, if the keywords "human" | "wastewater" don't appear in the index (because a custom db is used) all is good
        df.loc["wastewater"]["Portion"] = df.loc["wastewater"]["Portion"] + df.loc["human"]["Attributal"]
        df.loc["wastewater"]["Crx Total Genomes"] = df.loc["wastewater"]["Crx Total Genomes"] + df.loc["human"]["Total Genomes"]
        df.loc["wastewater"]["Crx Detected Genomes"] = df.loc["wastewater"]["Crx Detected Genomes"] + df.loc["human"]["Detected Genomes"]
        dr.drop("human",inplace=True)

    df["Attributal"].where(df["Attributal"] <= 0, 1, inplace=True)
    df["Crx Attributal"].where(df["Crx Attributal"] <= 0, 1, inplace=True)
    
    return df
        
def get_geq(args):
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
        '-o', '--output-dir',
        help='Path to an existing SourceApp output directory',
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
        '--drop-env',
        help='Discard environmental signal from final results. This can significantly impact apportioning results -- you probably want it on if you're apportioning fecal contamination.',
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
    
    if args['sourceapp_database'][-1] == '/': # in the event user provides trailing '/'
        args['sourceapp_database'] = args['sourceapp_database'][:-1]
        
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

    print('Summarizing SourceApp results.', flush=True)
    output_table = summarize(args)
    output_table.set_index("Source", inplace=True)

    raw, app, att, frac = clean_output(output_table, args)

    raw.to_csv(args['output_dir']+'/results.csv', index=True, header=raw.columns)
    att.to_csv(args['output_dir']+'/attributions.csv', index=True, header=["Detection"])
    app.to_csv(args['output_dir']+'/apportions.csv', index=True, header=["Portion"])
    frac.to_csv(args['output_dir']+'/fractions.csv', index=True, header=["Fraction"])

    print('The following results printed to results.csv in output directory:', flush=True)
    print(raw)

    print('Thank you for using SourceApp.', flush=True)

if __name__ == "__main__":
    main()
