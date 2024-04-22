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
def read_trim(args):
    fwd=args['input_files'].split(',')[0]
    rev=args['input_files'].split(',')[1]
    outdir=args['output_dir']
    html=outdir+'/fastp_summary.html'
    json=outdir+'/fastp_summary.json'
    out1=outdir+'/reads.1.fastq'
    out2=outdir+'/reads.2.fastq'
    threads=args['threads']
    skiptrim=args['skip_trimming']
    if skiptrim:
        print("Skipping read trimming...")
        subprocess.run(['cp '+fwd+' '+out1])
        subprocess.run(['cp '+rev+' '+out2])
    else:
        print("Running read trimming...")
        try:
            subprocess.run(["fastp -n 0 -l 70 --detect_adapter_for_pe -i " + fwd + " -I " + rev + " -o " + out1 + " -O " + out2 + " -w " + str(threads) + " -h " + html + " -j " + json], shell=True, check=True)
        except Exception as e:
            print('Error in step 1: read trimming. Exiting . . .')
            sys.exit()

def geq_estimation(args):
    threads=str(args['threads'])
    outdir = args['output_dir']
    fwd=outdir+'/reads.1.fastq'
    rev=outdir+'/reads.2.fastq'
    inputs=fwd + ',' + rev
    outputs=outdir+'/geq.txt'
    usegeq=args['use_geq']
    if usegeq:
        print("Running GEQ estimation...")
        try:
            subprocess.run(["run_microbe_census.py -v -t "+ threads + " -n 10000000 "+ inputs + ' ' + outputs], shell=True, check=True)
        except Exception as e:
            print('Error in step 2: GEQ estimation. Exiting . . .')
            sys.exit()
    else:
        print("Skipping GEQ estimation...")

def read_map(args):
    samthreads='-@'+str(args['threads'])
    threads=args['threads']
    fwd=args['input_files'].split(',')[0]
    rev=args['input_files'].split(',')[1]
    output = args['output_dir'] + '/mappings.bam'
    database = args['sourceapp_database'] + '/database'
    print("Running read mapping...")
    print("This can take a while.")
    try:
        subprocess.run(["bwa mem -t " + str(threads) + ' ' + database + ' ' + fwd + ' ' + rev + " | samtools sort "+ samthreads + " -o " + output + " -"], shell=True, check=True)
    except Exception as e:
        print('Error in step 3: read mapping. Exiting . . .')
        sys.exit()
    # bwa follows optional arguments and then positional required arguments

def read_filter(args):
    bam = args['output_dir'] + '/mappings.bam'
    gdef = args['sourceapp_database'] + '/gdef.txt'
    pid=args['percent_identity']
    qcov=args['query_coverage']
    threads=args['threads']
    trunc=args['limit_threshold']
    output = args['output_dir'] + '/mappings_filtered.txt'
    usegeq=args['use_geq']
    print("Filtering read mapping results...")
    
    if usegeq:
        try:
            subprocess.run(["coverm genome -b "+bam+" --genome-definition "+gdef+" --min-read-percent-identity "+str(pid*100)+
                        " --min-read-aligned-percent "+str(qcov*100)+" --output-format dense -t "+str(threads)+" -m trimmed_mean covered_bases variance "+
                        " -o "+output+" --trim-min "+str(trunc*100)+" --trim-max "+str(100-(trunc*100))], shell=True, check=True)
        except Exception as e:
            print('Error in step 4: filtering of read mappings. Exiting . . .')
            sys.exit()
    else:
        try:
            subprocess.run(["coverm genome -b "+bam+" --genome-definition "+gdef+" --min-read-percent-identity "+str(pid*100)+
                        " --min-read-aligned-percent "+str(qcov*100)+" --output-format dense -t "+str(threads)+" -m relative_abundance "+
                        " -o "+output+" --trim-min "+str(trunc*100)+" --trim-max "+str(100-(trunc*100))], shell=True, check=True)
        except Exception as e:
            print('Error in step 4: filtering of read mappings. Exiting . . .')
            sys.exit()

def summarize(args):
    #produce the final dataframe and make some visuals.
    usegeq=args['use_geq']
    if args["drop_env"]:
        temp_df = pd.read_csv(args['sourceapp_database'] + '/sources.txt',sep="\t",header=0)
        temp_df = temp_df[temp_df["crx"].str.contains("environmental_crx")==False]
        sourcedict = temp_df.iloc[:,0:2].set_index('genome')['source'].to_dict()
    else:
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
                gsum = gsum + (df[df['Genome']==genome].iloc[:,1].sum())/100
            fractions.append([source,gsum,gcount,ncount])
        fractions = pd.DataFrame(fractions, columns=['Source', 'Fraction','Detected Genomes','Total Genomes'])
    return fractions

### Helper functions:
def get_geq(args):
    file = args['output_dir'] + '/geq.txt'
    with open(file) as fh:
        for index, line in enumerate(fh):
            if index == 12:
                censusline = line.split()
    output = float(censusline[1])
    return output
    
def clean_output(table, args):
    if args["drop_env"]:
        table.drop("environmental",inplace=True)
        
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

    plst=[]
    for source in sources:
        if df.loc[source]["Attributal"] > 0:
            plst.append(df.loc[source]["Attributal"] + df.loc[source]["Crx Attributal"])
        else:
            plst.append(0)
            
    df["Total Fraction"] = plst

    if args["aggregate_human"]: # default is false; thus, if the keywords "human" | "wastewater" don't appear in the index (because a custom db is used) all is good
        df.loc["wastewater","Crx Fraction"] = df.loc["wastewater","Crx Fraction"] + df.loc["human","Attributal"]
        df.loc["wastewater","Total Fraction"] = df.loc["wastewater","Total Fraction"] + df.loc["human","Attributal"]
        df.loc["wastewater","Crx Total Genomes"] = df.loc["wastewater","Crx Total Genomes"] + df.loc["human","Total Genomes"]
        df.loc["wastewater","Crx Detected Genomes"] = df.loc["wastewater","Crx Detected Genomes"] + df.loc["human","Detected Genomes"]
        df.drop("human",inplace=True)

    df["Portion"] = df["Total Fraction"]/df["Total Fraction"].sum()
    
    df["Attributal"].where(df["Attributal"] <= 0, 1, inplace=True)
    df["Crx Attributal"].where(df["Crx Attributal"] <= 0, 1, inplace=True)

    df["pid"]=args["percent_identity"]
    df["mask"]=args["limit_threshold"]
    df["qcov"]=args["query_coverage"]
    df["thresh"]=args["min_frac"]
    
    return df
    
### Pipeline:
def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--input-files',
        help='Comma-delimited path to forward and reverse metagenomic reads. Must be in FASTQ format and compressed with gzip',
        metavar='',
        type=str, 
        required=True
        )
    parser.add_argument(
        '-o', '--output-dir',
        help='Path to the desired output directory',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-d', '--sourceapp-database',
        help='Path to directory containing a SourceApp formatted database',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-l', '--limit-threshold',
        help='Sequence breadth needed to consider a genome detected (default=0.1)',
        metavar='',
        type=float,
        default=0.1,
        required=False
        )
    parser.add_argument(
        '-f', '--min-frac',
        help='The minimum read or cell fraction a source must have to be considered detected and therefore apportionable (default=0.0001)',
        metavar='',
        type=float,
        default=0.0001,
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
        '--aggregate-human',
        help='Treat human signal as wastewater signal.',
        action='store_true',
        required=False
        )
    parser.add_argument(
        '--drop-env',
        help='Discard environmental signal from final results. This can significantly impact apportioning results -- you probably want it off if apportioning fecal contamination.',
        action='store_true',
        required=False
        )    
    args=vars(parser.parse_args())

    if args['output_dir'][-1] == '/': # in the event user provides trailing '/'
        args['output_dir'] = args['output_dir'][:-1]

    if args['sourceapp_database'][-1] == '/': # in the event user provides trailing '/'
        args['sourceapp_database'] = args['sourceapp_database'][:-1]
        
    if os.path.isdir(args['output_dir']):
        print('Error: Output directory already exists. Exiting . . .', flush=True)
        sys.exit()
    else:
        os.mkdir(args['output_dir'])

    print('Beginning step 1: read trimming', flush=True)
    read_trim(args)

    if args['use_geq']:
        print('Beginning step 2: GEQ estimation', flush=True)
        geq_estimation(args)
    else:
        print('Skipping step 2: GEQ estimation', flush=True)

    print('Beginning step 3: read mapping', flush=True)
    read_map(args)

    print('Beginning step 4: filtering of read mappings', flush=True)
    read_filter(args)

    print('Beginning step 5: results summarization', flush=True)
    table = summarize(args)
    table.set_index("Source", inplace=True)

    df = clean_output(table, args)

    df.to_csv(args['output_dir']+'/results.csv', index=True, header=df.columns)
    df["Attributal"].to_csv(args['output_dir']+'/attributions.csv', index=True, header=["Detection"]) # ignoring crx signal for attribution
    df["Portion"].to_csv(args['output_dir']+'/apportions.csv', index=True, header=["Portion"]) # portions rely on attribution
    df[["Total Fraction"]].to_csv(args['output_dir']+'/fractions.csv', index=True, header=["Fraction"]) # fractions rely on attribution

    print('The following results printed to results.csv in output directory:', flush=True)
    print(df)

    print('Thank you for using SourceApp.', flush=True)

if __name__ == "__main__":
    main()
