#!/usr/bin/env python3

'''
SourceApp: Python implementation of the Unix-based environmental monitoring tool.

sourceapp.py requires a properly formatted database. This script automates its creation.
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
    threads = args["threads"]
    input_dir = args["input_dir"]
    output_dir = args["output_name"] + "_SourceAppdb"
    try:
        subprocess.run(["checkm2 predict --threads " + str(threads) + " --input " + input_dir + " --output-directory " + output_dir + "/checkm2"], shell=True,check=True)
    except Exception as e:
        print("Error in step 01")
        print(e)
        sys.exit()

def genome_selection(args):
    # toss out bad genomes
    quality = args["genome_quality"]
    output_dir = args["output_name"]
    gdf = pd.read_csv(output_dir + "/checkm2/quality_report.tsv", sep="\t").iloc[:, 0:3]
    gdf.iloc[:, 0] = gdf.iloc[:, 0] + ".fna"
    gdf = gdf[(gdf.iloc[:, 1]/100) - 5 * (gdf.iloc[:, 2]/100) >= quality]
    sdf = pd.read_csv(args["source_associations"], sep="\t")

    for genome in sdf.iloc[:, 0]:  # if a genome is NOT in the good quality list, remove it from source info
        if not gdf.iloc[:, 0].str.contains(genome).any():
            sdf = sdf[~sdf.iloc[:, 0].str.contains(genome)]

    gdf.to_csv(output_dir + "/ginfo.csv", sep=",", header=["genome", "completeness", "contamination"], index=False)
    sdf.to_csv(output_dir + "/sinfo.csv", sep=",", header=None, index=False) # sinfo contains host associations for QC'd genomes


def genome_derep(args):
    input_dir = args["input_dir"]
    output_dir = args["output_name"]
    quality = args["genome_quality"]
    noderep = args["no_dereplication"]
    sdf = pd.read_csv(output_dir + "/sinfo.csv", sep=",",header=None)  # use clean source info to know which genomes passed QC
    sources = sdf.iloc[:, 1].unique()
    ginfo = output_dir + "/ginfo.csv"  # I don't think dRep cares if there are extraneous entries in this file; no need to filter for if/else/for below.
    ani = args["ani"]
    threads = args["threads"]
    subprocess.call(["mkdir " + output_dir + "/final_genomes/"], shell=True)

    if noderep:
        for genome in sdf.iloc[:, 0]:
            subprocess.call(["cp " + genome + " " + output_dir + "/final_genomes/"], shell=True)
        try:
            subprocess.run(["dRep dereplicate " + output_dir + "/drep_all -g "+ output_dir + "/final_genomes/*.fna --S_ani " + str(ani) + " --genomeInfo " + output_dir 
                            + "/ginfo.csv --S_algorithm fastANI -p " + str(threads) + " -comp " + str(quality*100) + " --skip_plots"],shell=True,check=True,stderr=subprocess.DEVNULL)
        except Exception as e:
            print("Error in step 3")
            print(e)
            sys.exit()
        print("no-dereplication has been called, which may result in some genomes from the same source being erroneously flagged as cross-reactive.")
        print("Flagging any cross-reactive genomes for use by sourceapp.py") 
        new_sdf = flag_crx(output_dir)
        new_sdf.to_csv(output_dir + "/sources.txt",index=False,header=None,mode="w")
    else:
        for source in sources:  # if we want to maintain crx genomes, then only dereplicate within sources (thus we'll run dRep N times)
            print("Processing genomes tagged with source: " + source)
            source_df = sdf[sdf.iloc[:, 1] == source].iloc[:, 0]
            source_df.iloc[:] = input_dir + "/" + source_df.iloc[:]
            source_df.to_csv(output_dir + "/" + source + ".glist.txt", index=False, header=None)
            subprocess.run(["echo 'genome,completeness,contamination' > " + output_dir + "/" + source + ".ginfo.csv; while read genome; do sid=$(basename ${genome}); grep ${sid} " 
                            + output_dir + "/ginfo.csv >> " + output_dir + "/" + source + ".ginfo.csv; done < " + output_dir+"/"+source+".glist.txt"],shell=True,check=True)
            try:
                subprocess.run(["dRep dereplicate " + output_dir + "/drep_" + source + " -g " + output_dir + "/" + source + ".glist.txt --S_ani " + str(ani)
                                + " --genomeInfo " + output_dir+"/"+source+".ginfo.csv" + " --S_algorithm fastANI -p " + str(threads) + " -comp " + str(quality*100) 
                                + " --skip_plots"], shell=True,check=True,stderr=subprocess.DEVNULL)
            except Exception as e:
                print("Error in step 3")
                print(e)
                sys.exit()
        print("Collecting non-redundant genome set")
        subprocess.run(["for dir in " + output_dir + "/drep_*; do cp ${dir}/dereplicated_genomes/*.fna " + output_dir + "/final_genomes/; done"],
                        shell=True,check=True)
        print("Flagging any cross-reactive genomes for use by sourceapp.py")
        try:
            subprocess.run(["dRep dereplicate " + output_dir + "/drep_all -g "+ output_dir + "/final_genomes/*.fna --S_ani " + str(ani) + " --genomeInfo " + output_dir 
                            + "/ginfo.csv --S_algorithm fastANI -p " + str(threads) + " -comp " + str(quality*100) + " --skip_plots"],shell=True,check=True,stderr=subprocess.DEVNULL)
        except Exception as e:
            print("Error in step 3")
            print(e)
            sys.exit()
        new_sdf = flag_crx(output_dir)
        new_sdf.to_csv(output_dir + "/sources.txt",index=False,header=None,mode="w")

def build_database(args):
    output_dir = args["output_name"]

    # concatenate genomes
    subprocess.run(["cat " + output_dir + "/final_genomes/*.fna >> " + output_dir + "/database.fna"], shell=True,check=True)

    # build gdef.txt
    subprocess.run(["grep '>' " + output_dir + "/database.fna > " + output_dir + "/contigs.txt"], shell=True,check=True)
    gdef = build_gdef(output_dir)
    gdef.to_csv(output_dir + "/gdef.txt", index=False, header=None, sep="\t")

    # indexing
    file_size = float(os.stat([output_dir + "/database.fna"][0]).st_size)/1073741824
    try: # we are metering our usage of -b based on database size to compromise on memory usage and speed in the case of very large inputs. 
         # i.e., if we want more speed we have to be prepared to provide more memory. we should warn users about memory utilization here
         # they'll need at least RAM >= 2 x input FASTA file size
        if file_size >= 1 and file_size < 30: # 1GB - 30GB
            subprocess.run(["bwa index -b 2500000000 -p " + output_dir + "/database " + output_dir + "/database.fna"],shell=True,check=True,stderr=subprocess.DEVNULL)
        elif file_size >= 30: # greater than 30GB
            subprocess.run(["bwa index -b 7500000000 -p " + output_dir + "/database " + output_dir + "/database.fna"],shell=True,check=True,stderr=subprocess.DEVNULL)
        else: # less than 1GB
            subprocess.run(["bwa index -p " + output_dir + "/database " + output_dir + "/database.fna"],shell=True,check=True,stderr=subprocess.DEVNULL)
    except Exception as e:
        print("Error in step 4")
        print(e)
        sys.exit()

    # clean up
    subprocess.run(["rm " + output_dir + "/final_genome_list.txt " + output_dir + "/*hs.txt " + output_dir + "/sinfo.csv "
         + output_dir + "/*glist.txt " + output_dir + "/contigs.txt " + output_dir + "/*.ginfo.txt"], shell=True,check=True)
    subprocess.run(["rm -r " + output_dir + "/drep* " + output_dir + "/checkm2 " + output_dir + "/final_genomes"],
                    shell=True,check=True)

# helper functions:
def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))

def Fasta_rename_sequences(infile, prefix):
    outfile = prefix + '.rename'
    with open(infile, 'r+') as f, open(outfile, 'w') as o:
        i = 1
        for name, seq in read_fasta(f):
            newName = f'>{prefix}_{i:04}\n{seq}\n'
            o.write(newName)
            i += 1
    _ = subprocess.run(['mv', outfile, infile])

def build_gdef(workdir):
    rhs = pd.read_csv([workdir + "/contigs.txt"][0], header=None).iloc[:,0].str[1:]
    lhs = rhs.str.rsplit("_",n=1).str.get(0)
    lhs = lhs + ".fna"
    gdef = pd.concat([lhs,rhs],axis=1)
    return gdef

def flag_crx(workdir):
    df = pd.read_csv(workdir + "/drep_all/data_tables/Cdb.csv", header=0)
    counts = df['secondary_cluster'].value_counts()
    singletons = counts[counts == 1]
    for cluster in singletons.index:
        df = df[~df.iloc[:,1].str.contains(cluster)]
    crx = df['genome']
    sdf = pd.read_csv(workdir + "/sinfo.csv", sep=",",header=None) 
    for x in range(len(sdf.iloc[:,0])): # all genomes
        genome = sdf.iloc[x,0]
        if crx.str.contains(genome).any(): # genomes found to be concerning
            sdf.iloc[x,1] = sdf.iloc[x,1] + "_crx" # tag those genomes we found concerning
        
    return sdf

### Pipeline:
def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        '-i', '--input-dir',
        help="Path to directory containing input genomes (path/to/dir/*.fna)",
        metavar='',
        type=str,
        required=True
    )
    parser.add_argument(
        '-o', '--output-name',
        help="Path to new directory which SourceApp should build the database in. SourceApp will create the database output directory \
              at the specified path appended with + '_SourceAppdb/'",
        metavar="",
        type=str,
        required=True
    )
    parser.add_argument(
        '-s', '--source-associations',
        help="Text file describing source associations of input genomes",
        metavar="",
        type=str,
        required=True,
    )
    parser.add_argument(
        '-a', '--ani',
        help="ANI threshold for calling genome clusters",
        metavar="",
        type=float,
        default=0.95,
        required=False
    )
    parser.add_argument(
        '-t', '--threads',
        help="Threads available to SourceApp",
        metavar="",
        type=int,
        default=1,
        required=False
    )
    parser.add_argument(
        '-q', '--genome-quality',
        help="Aggregate quality score threshold for accepting input genomes (float, 0.5 default)",
        metavar="",
        type=float,
        default=0.5,
        required=False
    )
    parser.add_argument(
        '--no-dereplication',
        help="Disable genome dereplication. This will create the database using all \
            of the provided genomes which pass quality requirements.",
        action="store_true",
        required=False
    )
    parser.add_argument(
        '-d', '--checkm2_db',
        help="Path to a local installation of the CheckM2 database (.dmnd). If not passed, SourceApp assumes you have let CheckM2 \
             install the database in the default location (~/databases). See 'checkm2 databases -h' for more information.",
        metavar="",
        type=str,
        required=False
    )
    parser.add_argument(
        '-c', '--checkm2_info',
        help="If you've already run CheckM2 yourself, path to the quality_report.tsv output.",
        metavar="",
        type=str,
        required=False
    )
    args = vars(parser.parse_args())

    print("SourceApp_build was called with the following arguments:")
    print(args)

    if args['output_name'][-1] == "/":
        args['output_name'] = args["output_name"][:-1] + "_SourceAppdb"
    else:
        args['output_name'] = args["output_name"] + "_SourceAppdb"
    output_dir = args['output_name']
    subprocess.run(["mkdir " + output_dir], shell=True,check=True)

    if not args["checkm2_info"]:
        if args["checkm2_db"]:
            try:
                subprocess.run(["checkm2 database --setdblocation " + args["checkm2_db"]], shell=True,check=True)
            except Exception as e:
                print("Error in CheckM2 database configuration")
                print(e)
                sys.exit()
        print("Step 01: checking quality of input genomes with CheckM2")
        genome_qc(args)
    else:
        subprocess.run(["mkdir " + output_dir + "/checkm2"], shell=True,check=True)
        subprocess.run(["cp " + args["checkm2_info"] + " " + output_dir + "/checkm2/quality_report.tsv"], shell=True,check=True)

    print("Step 02: finalizing database genome selection")
    genome_selection(args)

    print("Step 03: assessing database redundancy with dRep")
    genome_derep(args)

    # collect representative genomes and rename them
    dir=[output_dir + "/final_genomes"] 
    for genome in os.listdir(dir[0]):
        if genome.endswith(".fna"):
            prefix=os.path.splitext(genome)[0]
            Fasta_rename_sequences([output_dir + "/final_genomes/" + genome][0], prefix)
    
    print("Step 04: indexing database for BWA mem")
    build_database(args)

    print("SourceApp database construction complete")

if __name__ == "__main__":
    main()
