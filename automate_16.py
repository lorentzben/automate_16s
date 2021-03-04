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
    '_'+str(datetime.datetime.now().time()).replace(':', '')+'.log'
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

report_filename = "report_" + str(datetime.datetime.now().date()) + \
    '_'+str(datetime.datetime.now().time()).replace(':', '.')+".html"


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
    command = "qiime tools import \
        --type "+seq_format2+" \
        --input-path " + manifest+" \
        --output-path demux.qza \
        --input-format " + seq_format
    result = subprocess.run([command], stderr=PIPE, stdout=PIPE, shell=True)
    if result.returncode == 1:
        logger.critical("There are missing files please review manifest")
        exit(1)
    logger.info(result.stdout)
    logger.error(result.stderr)


def qual_control():
    logger.debug("checking quality of reads")
    command = "qiime demux summarize \
        --i-data demux.qza \
        --o-visualization demux_summary.qzv"
    result = subprocess.run([command], stderr=PIPE, stdout=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)

    export_demux_command = "qiime tools export \
        --input-path demux_summary.qzv \
        --output-path demux_summary/"
    result = subprocess.run([export_demux_command],
                            stderr=PIPE, stdout=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)

# point of optimization: increase the quality score higher by 3 points? for reverse reads.


def find_cutoffs(dataframe):
    mean_qual = dataframe[4:5]

    average_qual = np.round(mean_qual.mean(axis=1), 0)-1
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

    left_cutoff = 0
    right_cutoff = len(mean_qual_vals)-1

    for i in range(0, len(mean_qual_vals)):
        if mean_qual_vals[i] >= int(average_qual):
            left_cutoff = i+1
            break

    for i in range(0, len(mean_qual_vals)):
        if mean_qual_vals[len(mean_qual_vals)-1-i] >= int(average_qual):
            right_cutoff = len(mean_qual_vals)-i
            break

    return(left_cutoff, right_cutoff)


def calc_qual_cutoff(seq_format):
    if seq_format == "single":
        logger.debug("determining left and right cutoffs based on qual score")

        input_file = "demux_summary/forward-seven-number-summaries.tsv"

        summary = pd.read_table(input_file, index_col=0, sep='\t')
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
        forward_file = "demux_summary/forward-seven-number-summaries.tsv"
        fr_summary = pd.read_table(forward_file, index_col=0, sep='\t')

        forward = find_cutoffs(fr_summary)

        reverse_file = "demux_summary/reverse-seven-number-summaries.tsv"
        rev_summary = pd.read_table(reverse_file, index_col=0, sep='\t')

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
        command = "qiime dada2 denoise-single \
            --i-demultiplexed-seqs demux.qza \
            --p-trim-left " + str(left)+" \
            --p-trunc-len " + str(right) + " \
            --o-representative-sequences rep-seqs-dada2.qza \
            --o-table table-dada2.qza \
            --o-denoising-stats stats-dada2.qza"
    elif seq_format == 'paired':
        forward_left = cutoff[0][0]
        forward_right = cutoff[0][1]
        rev_left = cutoff[1][0]
        rev_right = cutoff[1][1]
        command = "qiime dada2 denoise-paired \
            --i-demultiplexed-seqs demux.qza \
            --p-trunc-len-f " + str(forward_right)+" \
            --p-trunc-len-r " + str(rev_right) + " \
            --p-trim-left-f " + str(forward_left)+" \
            --p-trim-left-r " + str(rev_left) + " \
            --o-representative-sequences rep-seqs-dada2.qza \
            --o-table table-dada2.qza \
            --o-denoising-stats stats-dada2.qza"

    result = subprocess.run([command], stdout=PIPE, stderr=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)

    command = "qiime metadata tabulate \
        --m-input-file stats-dada2.qza \
        --o-visualization stats-dada2.qzv"

    result = subprocess.run([command], stdout=PIPE, stderr=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)


def feature_visualizations(metadata):
    logger.debug("creating visualization objects")
    command = "qiime feature-table summarize \
        --i-table table-dada2.qza \
        --o-visualization table.qzv \
        --m-sample-metadata-file " + metadata
    result = subprocess.run([command], stdout=PIPE, stderr=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)

    command = "qiime feature-table tabulate-seqs \
        --i-data rep-seqs-dada2.qza \
        --o-visualization rep-seqs.qzv"
    result = subprocess.run([command], stdout=PIPE, stderr=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)


def tree_construction():
    logger.debug("generating phylogenetic tree")
    command = "qiime phylogeny align-to-tree-mafft-fasttree \
        --i-sequences rep-seqs-dada2.qza \
        --o-alignment aligned-rep-seqs.qza \
        --o-masked-alignment masked-aligned-rep-seqs.qza \
        --o-tree unrooted-tree.qza \
        --o-rooted-tree rooted-tree.qza"
    result = subprocess.run([command], stdout=PIPE, stderr=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)


def determine_depth():
    logger.debug("determining the best sampling depth to use ")
    export_table_vis_command = "qiime tools export \
        --input-path table.qzv \
        --output-path table_viz"
    result = subprocess.run([export_table_vis_command],
                            stderr=PIPE, stdout=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)

    input_file = "table_viz/sample-frequency-detail.csv"

    features = pd.read_csv(input_file, index_col=0, header=None)

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
    command = "qiime diversity core-metrics-phylogenetic \
        --i-phylogeny rooted-tree.qza \
        --i-table table-dada2.qza \
        --p-sampling-depth " + str(int(depth)) + " \
        --m-metadata-file " + metadata + " \
        --output-dir core-metrics-results \
        --p-n-jobs-or-threads 'auto'"
    result = subprocess.run([command], stdout=PIPE, stderr=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)

    logger.debug("calculating shannon alpha diversity")
    command = "qiime diversity alpha \
    --i-table table-dada2.qza \
    --p-metric shannon \
    --o-alpha-diversity shannon.qza"
    result = subprocess.run([command], stdout=PIPE, stderr=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)

    logger.debug("calculating simpson alpha diversity")
    command = "qiime diversity alpha \
    --i-table table-dada2.qza \
    --p-metric simpson \
    --o-alpha-diversity simpson.qza"
    result = subprocess.run([command], stdout=PIPE, stderr=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)

    logger.debug("calculating chao1 alpha diversity")
    command = "qiime diversity alpha \
    --i-table table-dada2.qza \
    --p-metric chao1 \
    --o-alpha-diversity chao1.qza"
    result = subprocess.run([command], stdout=PIPE, stderr=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)

    logger.debug("calculating ACE alpha diversity")
    command = "qiime diversity alpha \
    --i-table table-dada2.qza \
    --p-metric ace \
    --o-alpha-diversity ace.qza"
    result = subprocess.run([command], stdout=PIPE, stderr=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)

    logger.debug("calculating obs alpha diversity")
    command = "qiime diversity alpha \
    --i-table table-dada2.qza \
    --p-metric observed_otus \
    --o-alpha-diversity obs.qza"
    result = subprocess.run([command], stdout=PIPE, stderr=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)

    logger.debug("calculating phylogenetic diversity")
    command = "qiime diversity alpha-phylogenetic \
    --i-table table-dada2.qza \
    --i-phylogeny rooted-tree.qza \
    --p-metric=faith_pd\
    --o-alpha-diversity faith_pd.qza"
    result = subprocess.run([command], stdout=PIPE, stderr=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)

    logger.debug("calculating euclidian distance")
    calc_euclid_dist_command = "qiime diversity beta \
    --i-table table-dada2.qza \
    --p-metric euclidean \
    --o-distance-matrix core-metrics-results/euclidean_distance_results.qza"
    result = subprocess.run([calc_euclid_dist_command],
                            stdout=PIPE, stderr=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)


def calc_rare_depth():
    # reads sample frequences from previously unzipped archive
    sample_freq = pd.read_csv("table_viz/sample-frequency-detail.csv")
    depth = sample_freq.median()[0]

    return depth


def rarefy_curve_calc(depth, metadata):

    logger.debug("performing alpha rarefying")

    alpha_rare_command = "qiime diversity alpha-rarefaction \
        --i-table table-dada2.qza \
        --i-phylogeny rooted-tree.qza \
        --p-max-depth "+str(int(depth))+" \
        --m-metadata-file "+metadata+" \
        --o-visualization alpha-rarefaction.qzv"
    result = subprocess.run([alpha_rare_command],
                            stdout=PIPE, stderr=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)

    logger.debug("exporting alpha rarefaction so that r-notebook can find it")
    export_rare_command = "qiime tools export \
        --input-path alpha-rarefaction.qzv \
        --output-path alpha-rareplot"
    result = subprocess.run([export_rare_command],
                            stdout=PIPE, stderr=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)


def alpha_div_calc(metadata):
    logger.info("calculating significant features ")
    command = "qiime diversity alpha-group-significance \
    --i-alpha-diversity shannon.qza \
    --m-metadata-file " + metadata + " \
    --o-visualization shannon.qzv"
    result = subprocess.run([command], stdout=PIPE, stderr=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)

    command = "qiime diversity alpha-group-significance \
    --i-alpha-diversity simpson.qza \
    --m-metadata-file " + metadata + " \
    --o-visualization simpson.qzv"
    result = subprocess.run([command], stdout=PIPE, stderr=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)

    command = "qiime diversity alpha-group-significance \
    --i-alpha-diversity chao1.qza \
    --m-metadata-file " + metadata + " \
    --o-visualization chao1.qzv"
    result = subprocess.run([command], stdout=PIPE, stderr=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)

    command = "qiime diversity alpha-group-significance \
    --i-alpha-diversity ace.qza \
    --m-metadata-file " + metadata + " \
    --o-visualization ace.qzv"
    result = subprocess.run([command], stdout=PIPE, stderr=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)

    command = "qiime diversity alpha-group-significance \
    --i-alpha-diversity obs.qza \
    --m-metadata-file " + metadata + " \
    --o-visualization obs.qzv"
    result = subprocess.run([command], stdout=PIPE, stderr=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)

    command = "qiime diversity alpha-group-significance \
    --i-alpha-diversity faith_pd.qza \
    --m-metadata-file " + metadata+" \
    --o-visualization faith_pd.qzv"
    result = subprocess.run([command], stdout=PIPE, stderr=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)

    logger.info("unzipping features into separate dirs")
    command = "qiime tools export \
    --input-path shannon.qzv \
    --output-path shannon"
    result = subprocess.run([command], stdout=PIPE, stderr=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)

    command = "qiime tools export \
    --input-path simpson.qzv \
    --output-path simpson"
    result = subprocess.run([command], stdout=PIPE, stderr=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)

    command = "qiime tools export \
    --input-path chao1.qzv \
    --output-path chao1"
    result = subprocess.run([command], stdout=PIPE, stderr=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)

    command = "qiime tools export \
    --input-path ace.qzv \
    --output-path ace"
    result = subprocess.run([command], stdout=PIPE, stderr=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)

    command = "qiime tools export \
    --input-path obs.qzv \
    --output-path obs"
    result = subprocess.run([command], stdout=PIPE, stderr=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)

    command = "qiime tools export \
    --input-path faith_pd.qzv \
    --output-path faith_pd"
    result = subprocess.run([command], stdout=PIPE, stderr=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)


def beta_div_calc(metadata, item_of_interest):
    logger.debug(
        'calculating beta diversity unweighted, only done if column of metadata is provided')
    command = "qiime diversity beta-group-significance \
        --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
        --m-metadata-file " + metadata + " \
        --m-metadata-column "+item_of_interest + " \
        --o-visualization core-metrics-results/unweighted-unifrac-" + item_of_interest+"-significance.qzv \
        --p-pairwise"
    result = subprocess.run([command], stdout=PIPE, stderr=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)
    if result.returncode == 1:
        logger.critical(
            "the variable provided does not appear to be a column of the metadata file, or is not a categorical item, please review")
        exit(1)

    logger.debug("exporting unweighted unifrac visualization")
    export_unweight_command = "qiime tools export \
        --input-path core-metrics-results/unweighted-unifrac-" + item_of_interest+"-significance.qzv \
        --output-path unweighted-sig/"
    result = subprocess.run([export_unweight_command],
                            stdout=PIPE, stderr=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)

    logger.debug(
        'calculating beta diversity weighted, only done if column of metadata is provided')
    command = "qiime diversity beta-group-significance \
        --i-distance-matrix core-metrics-results/weighted_unifrac_distance_matrix.qza \
        --m-metadata-file " + metadata + " \
        --m-metadata-column "+item_of_interest + " \
        --o-visualization core-metrics-results/weighted-unifrac-" + item_of_interest+"-significance.qzv \
        --p-pairwise"
    result = subprocess.run([command], stdout=PIPE, stderr=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)
    if result.returncode == 1:
        logger.critical(
            "the variable provided does not appear to be a column of the metadata file, or is not a categorical item, please review")
        exit(1)

    logger.debug("exporting weighted unifrac visualization")
    export_weight_command = "qiime tools export \
        --input-path core-metrics-results/weighted-unifrac-" + item_of_interest+"-significance.qzv \
        --output-path weighted-sig/"
    result = subprocess.run([export_weight_command],
                            stdout=PIPE, stderr=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)


def assign_taxonomy():
    classifier = "wget https://data.qiime2.org/2021.2/common/silva-138-99-515-806-nb-classifier.qza -O silva-138-99-515-806-nb-classifier.qza"
    result = subprocess.run([classifier], stdout=PIPE, stderr=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)

    logger.info("taxonomically classifying sequences using sklearn")
    assign_tax = "qiime feature-classifier classify-sklearn \
    --i-classifier 'silva-138-99-515-806-nb-classifier.qza' \
    --i-reads rep-seqs-dada2.qza \
    --p-confidence 0.6 \
    --o-classification taxonomy.qza"
    result = subprocess.run([assign_tax], stdout=PIPE, stderr=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)

    tabulate = "qiime metadata tabulate \
    --m-input-file taxonomy.qza \
    --o-visualization taxonomy.qzv"
    result = subprocess.run([tabulate], stdout=PIPE, stderr=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)


def generate_phylogenetic_trees(metadata, item_interest):
    # pull in metadata file to extract the categories for the item of interest
    metadata_table = pd.read_table(metadata, sep='\t')
    ioi = item_interest

    metadata_table = metadata_table.drop([0, 1])

    ioi_set = set(metadata_table[item_interest])
    # iterates over the items of interest to produce a circular phylogenetic tree per category e.g. CONTROL TREATMENT
    for item in ioi_set:

        # filters/splits the feature table based on the current ioi
        filter_command = "qiime feature-table filter-samples \
            --i-table table-dada2.qza \
            --m-metadata-file metadata.tsv \
            --p-where " + ioi + "=" + item + " \
            --o-filtered-table "+item+"-filtered-table.qza"
        result = subprocess.run(
            [filter_command], shell=True, stdout=PIPE, stderr=PIPE)
        logger.info(result.stdout)
        logger.error(result.stderr)

        # adds taxonomic info needed for plotting
        collapse_command = "qiime taxa collapse \
            --i-table "+item + "-filtered-table.qza \
            --o-collapsed-table collapse-" + item+"-table.qza \
            --p-level 7 \
            --i-taxonomy taxonomy.qza"
        result = subprocess.run(
            [collapse_command], shell=True, stdout=PIPE, stderr=PIPE)
        logger.info(result.stdout)
        logger.error(result.stderr)

        # exports artifact so that the next step can collect it
        export_command = "qiime tools export \
            --input-path collapse-" + item+"-table.qza \
            --output-path collapse-"+item+"-frequency/"
        result = subprocess.run(
            [export_command], shell=True, stdout=PIPE, stderr=PIPE)
        logger.info(result.stdout)
        logger.error(result.stderr)

        # turns feature table into a human-reable format
        biom_command = "biom convert -i collapse-"+item + \
            "-frequency/feature-table.biom -o otu-"+item + \
            "-table.tsv --to-tsv --header-key taxonomy"
        result = subprocess.run(
            [biom_command], shell=True, stdout=PIPE, stderr=PIPE)
        logger.info(result.stdout)
        logger.error(result.stderr)

        # formatting the table so that it is in the correct order
        table = pd.read_table(
            "otu-"+str(item)+"-table.tsv", sep='\t', header=1)
        table = table.drop(columns=['taxonomy'])
        table = table.rename(columns={"#OTU ID": "taxonomy"})
        tax = table.pop("taxonomy")
        insertion_site = len(table.columns)
        table.insert(insertion_site, "taxonomy", tax)
        table.insert(0, "OTU_ID", np.arange(len(table)))
        table.to_csv("otu-"+str(item)+"-mod-table.tsv", sep='\t', index=False)

        # human readable table into compressed computer-readble format
        biom_format_command = "biom convert -i otu-"+str(item) + \
            "-mod-table.tsv -o otu-table-mod.biom --to-hdf5 --table-type=\'OTU table\' --process-obs-metadata taxonomy"
        result = subprocess.run([biom_format_command],
                                shell=True, stdout=PIPE, stderr=PIPE)
        logger.info(result.stdout)
        logger.error(result.stderr)

        # Outputs the current ioi so that it can be annotatted in the graphlan image
        with open("current.txt", "w") as file:
            file.write(item)

        # bash script call to handle the steps within a conda python 2.7.17 envionment
        generate_image_command = "./graph.sh"
        result = subprocess.run(
            [generate_image_command], shell=True, stdout=PIPE, stderr=PIPE)
        logger.info(result.stdout)
        logger.error(result.stderr)

        # renaming otu tables so they have meaning
        rename_table = "cp otu-table-mod.biom otu-table-"+item+"-mod.biom"
        result = subprocess.run(
            [rename_table], shell=True, stdout=PIPE, stderr=PIPE)
        logger.info(result.stdout)
        logger.error(result.stderr)

        # renaming the output of the graping bash script so that it has meaning
        rename_image = "cp image_graph.png image_"+item+"_graph.png"
        result = subprocess.run(
            [rename_image], shell=True, stdout=PIPE, stderr=PIPE)
        logger.info(result.stdout)
        logger.error(result.stderr)

        # renaming the output of the graping bash script so that it has meaning
        rename_image = "cp image_pdf_graph.png image_"+item+"+_pdf_g.png"
        result = subprocess.run(
            [rename_image], shell=True, stdout=PIPE, stderr=PIPE)
        logger.info(result.stdout)
        logger.error(result.stderr)


def lefse_analysis(item_interest):
    # call script to format the qiime data into lefse compatable format
    qiime_to_lefse_command = "Rscript qiime_to_lefse.R " + str(item_interest)
    result = subprocess.run([qiime_to_lefse_command],
                            shell=True, stdout=PIPE, stderr=PIPE)
    logger.info(result.stdout)
    logger.error(result.stderr)

    # call script to run conda env to run lefse in
    lefse_analysis_command = "bash lefse_analysis.sh"
    result = subprocess.run([lefse_analysis_command],
                            shell=True, stdout=PIPE, stderr=PIPE)
    logger.info(result.stdout)
    logger.error(result.stderr)


def generate_result_file(metadata):
    # dada 2 stats extracted from stats-dada2.qzv
    #command = "unzip -d "+folder+" stats-dada2.qzv"
    export_dada2_stats_command = "qiime tools export \
        --input-path stats-dada2.qzv \
        --output-path stats-dada2"
    result = subprocess.run([export_dada2_stats_command],
                            stderr=PIPE, stdout=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)

    command = ('cp stats-dada2/metadata.tsv dada2_stats.tsv')
    result = subprocess.run([command], stderr=PIPE, stdout=PIPE, shell=True)

    logger.info(result.stdout)
    logger.error(result.stderr)

   # Beta diversity core-metrics-results/unweighted_unifrac_emperor.qzv

    command = "cp "+metadata+" metadata.tsv"
    result = subprocess.run([command], stderr=PIPE, stdout=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)

    make_report_command = "./make_report.sh"
    result = subprocess.run([make_report_command],
                            stdout=PIPE, stderr=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)



def main(arg):

    with open('item_of_interest.csv', 'w', newline='') as csvfile:
        fieldnames = ['item name']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        writer.writerow({'item name': arg.interest})

    command = "cp " + arg.metadata+" metadata.tsv"
    result = subprocess.run([command], stdout=PIPE, stderr=PIPE, shell=True)
    logger.info(result.stdout)
    logger.error(result.stderr)

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

    assign_taxonomy()

    rare_depth = calc_rare_depth()

    # in: table-dada2.qza out: dir with rarefaction table
    rarefy_curve_calc(rare_depth, arg.metadata)

    # in: core-metrics-results/faith_pd_vector.qza out: core-metrics-results/faith-pd-group-significance.qzv
    alpha_div_calc(arg.metadata)

    if arg.interest:
        beta_div_calc(arg.metadata, arg.interest)

    generate_phylogenetic_trees(arg.metadata, arg.interest)

    lefse_analysis(arg.interest)

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
