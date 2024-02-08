#!/usr/bin/env python3

'''
SourceApp: Python implementation of the Unix-based environmental monitoring tool.

sourceapp.py requires a properly formatted reference. This script automates creation
of such a database.
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
        subprocess.run(["checkm2 predict --threads " + str(threads) + " --input " + input_dir +
                        " --output-directory " + output_dir + "/checkm2"], shell=True)
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
    sdf.to_csv(output_dir + "/sinfo.csv", sep=",", header=None, index=False)


def genome_derep(args):
    input_dir = args["input_dir"]
    output_dir = args["output_name"]
    noderep = args["no_dereplication"]
    removecrx = args["remove_crx"]
    sdf = pd.read_csv(output_dir + "/sinfo.csv", sep=",",
                      header=None)  # use clean source info to know which genomes passed QC
    sources = sdf.iloc[:, 1].unique()
    ginfo = output_dir + "/ginfo.csv"  # I don't think dRep cares if there are extraneous entries in this file; no need to filter for if/else/for below.
    ani = args["ani"]
    threads = args["threads"]
    subprocess.call(["mkdir " + output_dir + "/final_genomes/"], shell=True)

    if noderep:
        for genome in sdf.iloc[:, 0]:
            subprocess.call(["cp " + genome + " " + output_dir + "/final_genomes/"], shell=True)
    else:
        if removecrx:  # dereplicate ACROSS sources
            sdf.iloc[:, 0] = input_dir + "/" + sdf.iloc[:, 0]
            sdf.iloc[:, 0].to_csv(output_dir + "/glist.txt", index=False, header=None)
            try:
                subprocess.call([
                                    "dRep dereplicate " + output_dir + "/drep -g " + output_dir + "/glist.txt --S_ani " + str(
                                        ani) + " --genomeInfo " + ginfo +
                                    " --S_algorithm fastANI -p " + str(threads)], shell=True)
            except Exception as e:
                print("Error in step 3")
                print(e)
                sys.exit()
            subprocess.call(["cp " + output_dir + "/drep/dereplicated_genomes/*.fna " + output_dir + "/final_genomes/"],
                            shell=True)
        else:  # dereplicate WITHIN sources
            for source in sources:  # if we want to maintain crx genomes, then only dereplicate within sources (thus we'll run dRep N times)
                source_df = sdf[sdf.iloc[:, 1] == source].iloc[:, 0]
                source_df.iloc[:] = input_dir + "/" + source_df.iloc[:]
                source_df.to_csv(output_dir + "/" + source + ".glist.txt", index=False, header=None)
                try:
                    subprocess.call([
                                        "dRep dereplicate " + output_dir + "/drep_" + source + " -g " + output_dir + "/" + source + ".glist.txt" +
                                        " --S_ani " + str(
                                            ani) + " --genomeInfo " + ginfo + " --S_algorithm fastANI -p " + str(
                                            threads)], shell=True)
                except Exception as e:
                    print("Error in step 3")
                    print(e)
                    sys.exit()
            subprocess.call([
                                "for dir in " + output_dir + "/drep_*; do cp ${dir}/dereplicated_genomes/*.fna " + output_dir + "/final_genomes/"],
                            shell=True)

    # collect representative genomes and rename them
    subprocess.call([
                        "for genome in " + output_dir + "/final_genomes/*.fna; do sampleid=$(basename ${genome} | rev | cut -f 2- -d '.' | rev); seqtk rename ${genome} ${sampleid}. > ${genome}.new; done"],
                    shell=True)
    subprocess.call([
                        "for genome in " + output_dir + "/final_genomes/*.new; do sampleid=$(echo ${genome} | rev | cut -f 2- -d '.' | rev); mv ${genome} ${sampleid}; done"],
                    shell=True)


def build_database(args):
    # build source.txt
    output_dir = args["output_name"]
    subprocess.call([
                        "ls " + output_dir + "/final_genomes/*.fna | rev | cut -f 1 -d '/' | rev > " + output_dir + "/final_genome_list.txt"],
                    shell=True)
    subprocess.call(["while read line; do grep ${line} " + args[
        "source_associations"] + " >> rhs.txt; done < " + output_dir + "/final_genome_list.txt"], shell=True)
    subprocess.call(
        ["paste " + output_dir + "/final_genome_list.txt " + output_dir + "/rhs.txt >> " + output_dir + "/sources.txt"],
        shell=True)

    # concatenate genomes
    subprocess.call(["cat " + output_dir + "/final_genomes/*.fna >> " + output_dir + "/database.fna"], shell=True)

    # build gdef.txt
    subprocess.call(["grep '>' " + output_dir + "/database.fna > " + output_dir + "/contigs.txt"], shell=True)
    subprocess.call([
                        "while read line; do genome=$(echo ${line} | rev | cut -f 2- -d '.' | rev); echo ${genome}'.fna' >> " + output_dir + "/lhs.txt; done < " + output_dir + "/contigs.txt"],
                    shell=True)
    subprocess.call(["paste " + output_dir + "/lhs.txt " + output_dir + "/contigs.txt >> " + output_dir + "/gdef.txt"],
                    shell=True)

    # index
    try:
        subprocess.call(["bwa-mem2 index -p database " + output_dir + "/database.fna"],
                        shell=True)
    except Exception as e:
        print("Error in step 4")
        print(e)
        sys.exit()

    # clean up
    subprocess.call(
        ["rm " + output_dir + "/final_genome_list.txt " + output_dir + "/*hs.txt " + output_dir + "/sinfo.csv "
         + output_dir + "/*glist.txt " + output_dir + "/contigs.txt"], shell=True)
    subprocess.call(["rm -r " + output_dir + "/drep* " + output_dir + "/checkm2 " + output_dir + "/final_genomes"],
                    shell=True)


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
        help="Name of the database to be created. SourceApp will create an output directory in the current working \
             directory containing the finished database with the provided string + '_SourceAppdb/'",
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
        '--remove-crx',
        help="Remove genomes found in the same cluster but belonging to different sources.",
        action="store_true",
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
    subprocess.call(["mkdir " + output_dir], shell=True)

    if not args["checkm2_info"]:
        if args["checkm2_db"]:
            try:
                subprocess.call(["checkm2 database --setdblocation " + args["checkm2_db"]], shell=True)
            except Exception as e:
                print("Error in CheckM2 database configuration")
                print(e)
                sys.exit()
        print("Step 01: checking quality of input genomes with CheckM2")
        genome_qc(args)
    else:
        subprocess.call(["mkdir " + output_dir + "/checkm2"], shell=True)
        subprocess.call(["cp " + args["checkm2_info"] + " " + output_dir + "/checkm2/quality_report.tsv"], shell=True)

    print("Step 02: finalizing database genome selection")
    genome_selection(args)

    print("Step 03: assessing database redundancy with dRep")
    genome_derep(args)

    print("Step 04: indexing database for BWA-mem2")
    build_database(args)

    print("SourceApp database construction complete")


if __name__ == "__main__":
    main()
