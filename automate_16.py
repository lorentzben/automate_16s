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
from qiime2.plugins import dada2
import numpy as np
import glob


logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
# Logging handler which catches EVERYTHING
file_logger = logging.FileHandler('automate_16s.log')
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
        logger.critical(
            "cannot determine if paired or single end, check manifest file")
        exit(1)


def generate_seq_object(manifest, seq_format):
    #manifest = "manifest_trunc.tsv"
    #seq_format = "SingleEndFastqManifestPhred33V2"
    command = "qiime tools import --type 'SampleData[SequencesWithQuality]' --input-path " + \
        manifest+" --output-path demux.qza --input-format " + seq_format
    result = subprocess.run([command], stderr=PIPE, stdout=PIPE, shell=True)

    logger.info(result.stdout)
    logger.critical(result.stderr)


def qual_control():
    command = "qiime demux summarize --i-data demux.qza --o-visualization demux_summary.qzv"
    result = subprocess.run([command], stderr=PIPE, stdout=PIPE, shell=True)
    logger.info(result.stdout)
    logger.critical(result.stderr)
    command = "unzip -d inflate demux_summary.qzv"
    result = subprocess.run([command], stderr=PIPE, stdout=PIPE, shell=True)


def calc_qual_cutoff():
    input_file = glob.glob(
        './inflate/*/data/forward-seven-number-summaries.tsv')
    #input_file = 'inflate/*/data/forward-seven-number-summaries.tsv'
    summary = pd.read_table(input_file[0], index_col=0, sep='\t')

    mean_qual = summary[4:5]

    average_qual = np.round(mean_qual.mean(axis=1), 0)
    mean_qual_vals = np.array(mean_qual)[0]

    if int(average_qual) < 30:
        print("The Average Quality of these sequences may be a concern would you like to continue?")
        exit(0)

    for i in range(0, len(mean_qual_vals)):
        if mean_qual_vals[i] >= int(average_qual):
            left_cutoff = i+1
            break
    for i in range(0, len(mean_qual_vals)):
        if mean_qual_vals[len(mean_qual_vals)-1-i] >= int(average_qual):
            right_cutoff = len(mean_qual_vals)-i
            break

    logger.info("right cutoff: "+str(right_cutoff))
    logger.info("left cutoff: " + str(left_cutoff))

    with open('cutoffs.csv', 'w', newline='') as csvfile:
        fieldnames = ['cutoff', 'value']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        writer.writerow({'cutoff': 'right', 'value': right_cutoff})
        writer.writerow({'cutoff': 'left', 'value': left_cutoff})
        writer.writerow({'cutoff': 'filename', 'value': input_file})

    return(right_cutoff, left_cutoff)


def call_denoise(right, left, seq_format):
    if seq_format == 'single':
        command = "qiime dada2 denoise-single --i-demultiplexed-seqs demux.qza --p-trim-left " + str(left)+" --p-trunc-len " + \
            str(right) + " --o-representative-sequences rep-seqs-dada2.qza --o-table table-dada2.qza --o-denoising-stats stats-dada2.qza"
    elif seq_format == 'paired':
        command = "qiime dada2 denoise-paired --i-demultiplexed-seqs demux.qza --p-trim-left " + str(left)+" --p-trunc-len " + \
            str(right) + " --o-representative-sequences rep-seqs-dada2.qza --o-table table-dada2.qza --o-denoising-stats stats-dada2.qza"

    result = subprocess.run([command], stdout=PIPE, stderr=PIPE, shell=True)
    logger.info(result.stdout)
    logger.critical(result.stderr)

    command = "qiime metadata tabulate --m-input-file stats-dada2.qza --o-visualization stats-dada2.qza"

    result = subprocess.run([command], stdout=PIPE, stderr=PIPE, shell=True)
    logger.info(result.stdout)
    logger.critical(result.stderr)


def feature_visualizations(metadata):
    command = "qiime feature-table summarize --i-table table-dada2.qza --o-visualization table.qzv --m-sample-metadata-file " + metadata
    result = subprocess.run([command], stdout=PIPE, stderr=PIPE, shell=True)
    logger.info(result.stdout)
    logger.critical(result.stderr)

    command = "qiime feature-table tabulate-seqs --i-data rep-seqs-dada2.qza --o-visualization rep-seqs.qzv"
    result = subprocess.run([command], stdout=PIPE, stderr=PIPE, shell=True)
    logger.info(result.stdout)
    logger.critical(result.stderr)


def tree_construction():
    command = "qiime phylogeny align-to-tree-mafft-fasttree --i-sequences rep-seqs-dada2.qza --o-alignment aligned-rep-seqs.qza --o-masked-alignment masked-aligned-rep-seqs.qza --o-tree unrooted-tree.qza --o-rooted-tree rooted-tree.qza"
    result = subprocess.run([command], stdout=PIPE, stderr=PIPE, shell=True)
    logger.info(result.stdout)
    logger.critical(result.stderr)


def determine_depth():
    command = "unzip -d inflate table.qzv"
    result = subprocess.run([command], stderr=PIPE, stdout=PIPE, shell=True)

    logger.info(result.stdout)
    logger.critical(result.stderr)

    input_file = glob.glob(
        './inflate/*/data/sample-frequency-detail.csv')

    features = pd.read_csv(input_file[0], index_col=0, header=None)

    total_count = sum(features[1])
    sampling_depth = 0
    perc_features_retain = 0.0

    print("total count: " + str(total_count))

    feature_array = np.array(features)

    for i in range(len(feature_array)-1, -1, -1):
        sampling_depth = feature_array[i][0]
        perc_features_retain = ((sampling_depth * (i+1))/total_count)
        if perc_features_retain > .22:
            sampling_depth = feature_array[i][0]
            print("sampling depth: " + str(sampling_depth) + " % features retained: " +
                  str(round(perc_features_retain, 3)) + " samples retained: " + str(i))
            break

    with open('sampling_depth.csv', 'w', newline='') as csvfile:
        fieldnames = ['stat', 'value']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        writer.writerow({'stat': 'sampling_depth', 'value': sampling_depth})
        writer.writerow({'stat': '%_features_retained',
                         'value': perc_features_retain})
        writer.writerow({'stat': 'filename', 'value': input_file})

    logger.info("sampling_depth: "+str(sampling_depth))
    logger.info("%_features_retained: " + str(perc_features_retain))

    return(sampling_depth)


def diversity_measure(metadata, depth):
    command = "qiime diversity core-metrics-phylogenetic --i-phylogeny rooted-tree.qza --i-table table-dada2.qza --p-sampling-depth " + \
        str(int(depth)) + " --m-metadata-file " + \
        metadata + " --output-dir core-metrics-results"
    result = subprocess.run([command], stdout=PIPE, stderr=PIPE, shell=True)
    logger.info(result.stdout)
    logger.critical(result.stderr)


def alpha_div_calc(metadata):
    command = "qiime diversity alpha-group-significance --i-alpha-diversity core-metrics-results/faith_pd_vector.qza --m-metadata-file " + \
        metadata + "--o-visualization core-metrics-results/faith-pd-group-significance.qzv"
    result = subprocess.run([command], stdout=PIPE, stderr=PIPE, shell=True)
    logger.info(result.stdout)
    logger.critical(result.stderr)

    command = "qiime diversity alpha-group-significance --i-alpha-diversity core-metrics-results/evenness_vector.qza --m-metadata-file " + \
        metadata + "--o-visualization core-metrics-results/evenness-group-significance.qzv"
    result = subprocess.run([command], stdout=PIPE, stderr=PIPE, shell=True)
    logger.info(result.stdout)
    logger.critical(result.stderr)


def beta_div_calc(metadata, item_of_interest):
    command = "qiime diversity beta-group-significance --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza --m-metadata-file " + \
        metadata + " --m-metadata-column "+item_of_interest + \
        " --o-visualization fore-metric-results/unwighted-unifrac-body-site-significance.qzv --p-pairwise"
    result = subprocess.run([command], stdout=PIPE, stderr=PIPE, shell=True)
    logger.info(result.stdout)
    logger.critical(result.stderr)

# TODO put checks in to pickeup from where a failed run left off.


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

    call_denoise(right_cutoff, left_cutoff, single_or_pair)

    feature_visualizations(arg.metadata)

    tree_construction()

    depth = determine_depth()

    diversity_measure(arg.metadata, depth)

    alpha_div_calc(arg.metadata)

    if arg.interest:
        beta_div_calc(arg.metadata, arg.interest)


if __name__ == "__main__":
    # Build Argument Parser in order to facilitate ease of use for user
    parser = argparse.ArgumentParser(
        description="Checks dependencies are installed, and validates a manifest file for qiime 2")
    parser.add_argument('-n', '--name', action='store', required=True,
                        help="name for the manifest, typically manifest.tsv", dest='manifest_name')
    parser.add_argument('-m', '--metadata', action='store', required=True,
                        help="name of the metadata file, usually metadata.tsv", dest='metadata')
    parser.add_argument('-i', "--interest", action='store', required=False,
                        help="item of interest for beta diversity analysis, must match one of the column-names in the metadata file", dest="interest")
    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s 1.0')

    args = parser.parse_args()
    # print(args)
    main(args)
