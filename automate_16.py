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
import numpy as np
import glob
import datetime


logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
# Logging handler which catches EVERYTHING
logfile_name = "automate_"+str(datetime.datetime.now().date()) + \
    '_'+str(datetime.datetime.now().time()).replace(':', '.')+'.log'
file_logger = logging.FileHandler(logfile_name)
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

folder = "inflate_" + str(datetime.datetime.now().date()) + \
    '_'+str(datetime.datetime.now().time()).replace(':', '.')


def unique_folder_name(stub):
    u_folder = stub + "_" + str(datetime.datetime.now().date()) + \
        '_'+str(datetime.datetime.now().time()).replace(':', '.')
    return u_folder


def single_or_paired_read(manifest):
    logger.debug("determining paired or single end design")
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


def generate_seq_object(manifest, seq_format, seq_format2):
    logger.debug("importing fastq files to qiime2 artifact format")
    command = "qiime tools import --type "+seq_format2+" --input-path " + \
        manifest+" --output-path demux.qza --input-format " + seq_format
    result = subprocess.run([command], stderr=PIPE, stdout=PIPE, shell=True)
    if result.returncode == 1:
        logger.critical("There are missing files please review manifest")
        exit(1)
    logger.info(result.stdout)
    logger.error(result.stderr)


def qual_control():
    logger.debug("checking quality of reads")
    command = "qiime demux summarize --i-data demux.qza --o-visualization demux_summary.qzv"
    result = subprocess.run([command], stderr=PIPE, stdout=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)

    command = "unzip -d " + folder + " demux_summary.qzv"
    result = subprocess.run([command], stderr=PIPE, stdout=PIPE, shell=True)

# point of optimization: increase the quality score higher by 3 points? for reverse reads.


def find_cutoffs(dataframe):
    mean_qual = dataframe[4:5]

    average_qual = np.round(mean_qual.mean(axis=1), 0)
    mean_qual_vals = np.array(mean_qual)[0]

    if int(average_qual) < 30:
        print(
            "The Average Quality of these sequences may be a concern would you like to continue?")
        exit(0)

    for i in range(0, len(mean_qual_vals)):
        if mean_qual_vals[i] >= int(average_qual):
            left_cutoff = i+1
            break
    for i in range(0, len(mean_qual_vals)):
        if mean_qual_vals[len(mean_qual_vals)-1-i] >= int(average_qual):
            right_cutoff = len(mean_qual_vals)-i
            break
    return(left_cutoff, right_cutoff)

def find_rev_cutoffs(dataframe):
    mean_qual = dataframe[4:5]

    average_qual = np.round(mean_qual.mean(axis=1), 0)+2
    mean_qual_vals = np.array(mean_qual)[0]

    if int(average_qual) < 30:
        print(
            "The Average Quality of these sequences may be a concern would you like to continue?")
        exit(0)

    for i in range(0, len(mean_qual_vals)):
        if mean_qual_vals[i] >= int(average_qual):
            left_cutoff = i+1
            break
        else:
            left_cutoff = 0
            break
    for i in range(0, len(mean_qual_vals)):
        if mean_qual_vals[len(mean_qual_vals)-1-i] >= int(average_qual):
            right_cutoff = len(mean_qual_vals)-i
            break
        else:
            right_cutoff=(len(mean_qual_vals)-1)
            break
    return(left_cutoff, right_cutoff)

def calc_qual_cutoff(seq_format):
    if seq_format == "single":
        logger.debug("determining left and right cutoffs based on qual score")

        input_file = glob.glob(
            './'+folder+'/*/data/forward-seven-number-summaries.tsv')

        summary = pd.read_table(input_file[0], index_col=0, sep='\t')
        left_cutoff, right_cutoff = find_cutoffs(summary)

        logger.info("right cutoff: "+str(right_cutoff))
        logger.info("left cutoff: " + str(left_cutoff))

        with open('cutoffs.csv', 'w', newline='') as csvfile:
            fieldnames = ['cutoff', 'value']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

            writer.writeheader()
            writer.writerow({'cutoff': 'right', 'value': right_cutoff})
            writer.writerow({'cutoff': 'left', 'value': left_cutoff})
            writer.writerow({'cutoff': 'filename', 'value': input_file})

        return(left_cutoff, right_cutoff)

    elif seq_format == "paired":
        logger.debug(
            "determining forward and revese, left and right cutoffs based on qual score")
        forward_file = glob.glob(
            './'+folder+'/*/data/forward-seven-number-summaries.tsv')
        fr_summary = pd.read_table(forward_file[0], index_col=0, sep='\t')

        forward = find_cutoffs(fr_summary)

        reverse_file = glob.glob(
            './'+folder+'/*/data/reverse-seven-number-summaries.tsv')
        rev_summary = pd.read_table(reverse_file[0], index_col=0, sep='\t')

        reverse = find_rev_cutoffs(rev_summary)

        logger.info("forward cutoffs: "+str(forward))
        logger.info("reverse cutoffs: " + str(reverse))

        with open('cutoffs.csv', 'w', newline='') as csvfile:
            fieldnames = ['cutoff', 'value']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

            writer.writeheader()
            writer.writerow({'cutoff': 'forward left', 'value': forward[0]})
            writer.writerow({'cutoff': 'forward right', 'value': forward[1]})
            writer.writerow({'cutoff': 'reverse left', 'value': reverse[0]})
            writer.writerow({'cutoff': 'reverse right', 'value': reverse[1]})
            writer.writerow({'cutoff': 'filename', 'value': forward_file})
            writer.writerow({'cutoff': 'filename', 'value': reverse_file})

        return(forward, reverse)


def call_denoise(cutoff, seq_format):
    logger.debug("denoising using dada2")
    if seq_format == 'single':
        left = cutoff[0]
        right = cutoff[0]
        command = "qiime dada2 denoise-single --i-demultiplexed-seqs demux.qza --p-trim-left " + str(left)+" --p-trunc-len " + \
            str(right) + " --o-representative-sequences rep-seqs-dada2.qza --o-table table-dada2.qza --o-denoising-stats stats-dada2.qza"
    elif seq_format == 'paired':
        forward_left = cutoff[0][0]
        forward_right = cutoff[0][1]
        rev_left = cutoff[1][0]
        rev_right = cutoff[1][1]
        command = "qiime dada2 denoise-paired --i-demultiplexed-seqs demux.qza --p-trunc-len-f " + str(forward_right)+" --p-trunc-len-r " + \
            str(rev_right) + " --p-trim-left-f " + str(forward_left)+" --p-trim-left-r " + str(rev_left) + \
            " --o-representative-sequences rep-seqs-dada2.qza --o-table table-dada2.qza --o-denoising-stats stats-dada2.qza"

    result = subprocess.run([command], stdout=PIPE, stderr=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)

    command = "qiime metadata tabulate --m-input-file stats-dada2.qza --o-visualization stats-dada2.qzv"

    result = subprocess.run([command], stdout=PIPE, stderr=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)


def feature_visualizations(metadata):
    logger.debug("creating visualization objects")
    command = "qiime feature-table summarize --i-table table-dada2.qza --o-visualization table.qzv --m-sample-metadata-file " + metadata
    result = subprocess.run([command], stdout=PIPE, stderr=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)

    command = "qiime feature-table tabulate-seqs --i-data rep-seqs-dada2.qza --o-visualization rep-seqs.qzv"
    result = subprocess.run([command], stdout=PIPE, stderr=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)


def tree_construction():
    logger.debug("generating phylogenetic tree")
    command = "qiime phylogeny align-to-tree-mafft-fasttree --i-sequences rep-seqs-dada2.qza --o-alignment aligned-rep-seqs.qza --o-masked-alignment masked-aligned-rep-seqs.qza --o-tree unrooted-tree.qza --o-rooted-tree rooted-tree.qza"
    result = subprocess.run([command], stdout=PIPE, stderr=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)


def determine_depth():
    logger.debug("determining the best sampling depth to use ")
    command = "unzip -d "+folder+" table.qzv"
    result = subprocess.run([command], stderr=PIPE, stdout=PIPE, shell=True)

    logger.info(result.stdout)
    logger.error(result.stderr)

    input_file = glob.glob(
        './'+folder+'/*/data/sample-frequency-detail.csv')

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
    logger.debug("writing dept out to file")
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
    logger.debug("calculating general diversity measurments")
    command = "qiime diversity core-metrics-phylogenetic --i-phylogeny rooted-tree.qza --i-table table-dada2.qza --p-sampling-depth " + \
        str(int(depth)) + " --m-metadata-file " + \
        metadata + " --output-dir core-metrics-results --p-n-jobs-or-threads 'auto'"
    result = subprocess.run([command], stdout=PIPE, stderr=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)


def alpha_div_calc(metadata):
    logger.debug('calculating alpha diversity')
    command = "qiime diversity alpha-group-significance --i-alpha-diversity core-metrics-results/faith_pd_vector.qza --m-metadata-file " + \
        metadata + " --o-visualization core-metrics-results/faith-pd-group-significance.qzv"
    result = subprocess.run([command], stdout=PIPE, stderr=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)

    command = "qiime diversity alpha-group-significance --i-alpha-diversity core-metrics-results/evenness_vector.qza --m-metadata-file " + \
        metadata + " --o-visualization core-metrics-results/evenness-group-significance.qzv"
    result = subprocess.run([command], stdout=PIPE, stderr=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)


def beta_div_calc(metadata, item_of_interest):
    logger.debug(
        'calculating beta diversity, only done if column of metadata is provided')
    command = "qiime diversity beta-group-significance --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza --m-metadata-file " + \
        metadata + " --m-metadata-column "+item_of_interest + \
        " --o-visualization core-metrics-results/unweighted-unifrac-" + \
        item_of_interest+"-significance.qzv --p-pairwise"
    result = subprocess.run([command], stdout=PIPE, stderr=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)
    if result.returncode == 1:
        logger.critical(
            "the variable provided does not appear to be a column of the metadata file, please review")
        exit(1)


def generate_result_file(metadata):
    # dada 2 stats extracted from stats-dada2.qzv
    command = "unzip -d "+folder+" stats-dada2.qzv"
    result = subprocess.run([command], stderr=PIPE, stdout=PIPE, shell=True)

    logger.info(result.stdout)
    logger.error(result.stderr)

    command = ('cp ./'+folder+'/*/data/metadata.tsv dada2_stats.tsv')
    result = subprocess.run([command], stderr=PIPE, stdout=PIPE, shell=True)

    logger.info(result.stdout)
    logger.error(result.stderr)

    # Alpha diversity core-metrics-results/evenness-group-significance.qzv AND faith-pd-group-significance.qzv
    # evenness diversity measurement
    command = "cp core-metrics-results/evenness-group-significance.qzv ."
    result = subprocess.run([command], stderr=PIPE, stdout=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)
    
    command = "unzip -d evenness -j -n evenness-group-significance.qzv"
    result = subprocess.run([command], stderr=PIPE, stdout=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)

    # faith diversity measurement
    command = "cp core-metrics-results/faith-pd-group-significance.qzv ."
    result = subprocess.run([command], stderr=PIPE, stdout=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)

    command = "unzip -d faith -j -n faith-pd-group-significance.qzv"
    result = subprocess.run([command], stderr=PIPE, stdout=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)

    # Beta diversity core-metrics-results/unweighted_unifrac_emperor.qzv
    '''
    # bray curtis distance measuremnet
    command = "cp core-metrics-results/bray_curtis_distance_matrix.qza ."
    result = subprocess.run([command], stderr=PIPE, stdout=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)
    
    command = "unzip -d bray -j -n bray_curtis_distance_matrix.qza"
    result = subprocess.run([command], stderr=PIPE, stdout=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)

    
    # jaccard distance measuremnt
    command = "cp core-metrics-results/jaccard_distance_matrix.qza ."
    result = subprocess.run([command], stderr=PIPE, stdout=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)
    
    command = "unzip -d jaccard -j -n jaccard_distance_matrix.qza"
    result = subprocess.run([command], stderr=PIPE, stdout=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)

    
    # unweighted unifrac distance
    command = "cp core-metrics-results/unweighted_unifrac_distance_matrix.qza ."
    result = subprocess.run([command], stderr=PIPE, stdout=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)
    
    command = "unzip -d un-unifrac -j -n unweighted_unifrac_distance_matrix.qza"
    result = subprocess.run([command], stderr=PIPE, stdout=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)

    
    #   weighted unifrac distance
    command = "cp core-metrics-results/weighted_unifrac_distance_matrix.qza ."
    result = subprocess.run([command], stderr=PIPE, stdout=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)
    
    command = "unzip -d unifrac -j -n weighted_unifrac_distance_matrix.qza"
    result = subprocess.run([command], stderr=PIPE, stdout=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)
    '''
    command = "cp "+metadata+" metadata.tsv"
    result = subprocess.run([command], stderr=PIPE, stdout=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)

    # Sequence to render the rnotebook into a html object
    command = 'Rscript -e "rmarkdown::render(\'report.Rmd\', clean=TRUE)"'
    result = subprocess.run([command], stdout=PIPE, stderr=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)
    if result.returncode == 1:
        logger.critical(
            "there was an issue generating the report for this analysis")
        exit(1)


def main(arg):

    single_or_pair = single_or_paired_read(arg.manifest_name)

    if single_or_pair == "single":
        category = "SingleEndFastqManifestPhred33V2"
        cat2 = "SampleData[SequencesWithQuality]"
    elif single_or_pair == "paired":
        category = "PairedEndFastqManifestPhred33V2"
        cat2 = "SampleData[PairedEndSequencesWithQuality]"

    #out: demux.qza
    generate_seq_object(arg.manifest_name, category, cat2)
    qual_control()
    cutoffs = calc_qual_cutoff(single_or_pair)

    print(cutoffs)

    # in: demux.qza out: rep-seqs-dada2.qza table-dada2.qza stats-dada2.qza
    call_denoise(cutoffs, single_or_pair)

    # in: table-dada2.qza out: table.qzv rep-seqs.qzv
    feature_visualizations(arg.metadata)

    # in: rep-seqs-dada2.qza out: aligned-rep-seqs.qza
    tree_construction()

    # in: table.qzv  out: sampling_depth.csv
    depth = determine_depth()

    # in: rooted-tree.qza table-dada2.qza out: core-metrics-results/
    diversity_measure(arg.metadata, depth)

    # in: core-metrics-results/faith_pd_vector.qza out: core-metrics-results/faith-pd-group-significance.qzv
    alpha_div_calc(arg.metadata)

    if arg.interest:
        beta_div_calc(arg.metadata, arg.interest)

    generate_result_file(arg.metadata)

    logger.info('done')


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
