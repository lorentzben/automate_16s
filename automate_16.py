import subprocess
from subprocess import PIPE
import os
from os import path
import logging
from pathlib import Path
from pathlib import PurePath
import pandas as pd
import csv
import argparse
import qiime2

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
# Logging handler which catches EVERYTHING
file_logger = logging.FileHandler('automate_16s_setup.log')
file_logger.setLevel(logging.DEBUG)
# Logging handler which logs less
console_logger = logging.StreamHandler()
console_logger.setLevel(logging.ERROR)

# Formats the logs so they are pretty
logFormatter = '%(asctime)s- %(name)s - %(lineno)s - %(levelname)s - %(message)s'
formatter = logging.Formatter(logFormatter)
file_logger.setFormatter(formatter)
console_logger.setFormatter(formatter)

#  adds handlers to logger
logger.addHandler(file_logger)
logger.addHandler(console_logger)

def single_or_paired_read(manifest):
    try:
        read_manifest = pd.read_table(manifest, index_col=0, sep='\t')
    except FileNotFoundError:
        logger.critical("that manifest file does not exist")
        exit(1)

    if read_manifest.columns[0] == 'absolute-filepath':
        logger.info("single end analysis")
        return 'single'
    elif read_manifest.columns[0] == 'forward-absolute-filepath':
        logger.info("paired end analsis")
        return 'paired'
    else:
        logger.critical("cannot determine if paired or single end, check manifest file")
        exit(1)
    

def generate_seq_object(manifest, seq_format):
    #manifest = "manifest_trunc.tsv"
    #seq_format = "SingleEndFastqManifestPhred33V2"
    command = "qiime tools import --type 'SampleData[SequencesWithQuality]' --input-path "+manifest+" --output-path single_end_demux.qza --input-format " +seq_format
    result = subprocess.run([command], stderr=PIPE, stdout=PIPE, shell=True)

    logger.info(result.stdout)
    logger.critical(result.stderr)

def qual_control():
    command = "qiime demux summarize --i-data single_end_demux.qza --o-visualization demux_summary.qzv"
    result = subprocess.run([command], stderr=PIPE, stdout=PIPE, shell=True)
    logger.info(result.stdout)
    logger.critical(result.stderr)
    command = "unzip -d inflate demux_summary.qzv"
    result = subprocess.run([command], sterr=PIPE, stdout=PIPE, shell=True)

def calc_qual_cutoff():
    input_file = 'inflate/*/forward-seven-number-summaries.csv'
    summary = pd.read_csv(input_file, index_col=0,)

    mean_qual = summary[4:5]

    average_qual = np.round(mean_qual.mean(axis=1), 0)
    mean_qual_vals = np.array(mean_qual)[0]


    if int(average_qual) < 30:
        print("The Average Quality of these sequences may be a concern would you like to continue?")
        exit(0)


    for i in range(0, len(mean_qual_vals)):
        if mean_qual_vals[i] >= int(average_qual):
            right_cutoff = i+1
            break
    for i in range(0,len(mean_qual_vals)):
        if mean_qual_vals[len(mean_qual_vals)-1-i] >= int(average_qual):
            left_cutoff = len(mean_qual_vals)-i
            break

    logger.info("right cutoff: "+right_cutoff)
    logger.info("left cutoff: " +left_cutoff)

    with open('cutoffs.csv', 'w', newline='') as csvfile:
        fieldnames = ['cutoff', 'value']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        writer.writerow({'cutoff': 'right', 'value': right_cutoff})
        writer.writerow({'cutoff': 'left', 'value': left_cutoff})
        writer.writerow({'cutoff': 'filename', 'value': input_file})

    return(tuple(right_cutoff,left_cutoff))
    



def main(arg):
    single_or_pair = single_or_paired_read(arg.manifest_name)

    if single_or_pair == "single":
        category = "SingleEndFastqManifestPhred33V2"
    elif single_or_pair == "paired":
        category = "PairedEndFastqManifestPhred33V2"
    
    generate_seq_object(arg.manifest_name, category)
    qual_control()
    cutoffs = calc_qual_cutoff()
    right_cutoff = cutoffs[0]
    left_cutoff = cutoffs[1]

    print(right_cutoff, left_cutoff)

if __name__ == "__main__":
    # Build Argument Parser in order to facilitate ease of use for user
    parser = argparse.ArgumentParser(
        description="Checks dependencies are installed, and validates a manifest file for qiime 2")
    parser.add_argument('-n', '--name', action='store', required=True,
                        help="name for the manifest, typically manifest.tsv", dest='manifest_name')
    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s 1.0')

    args = parser.parse_args()
    # print(args)
    main(args)