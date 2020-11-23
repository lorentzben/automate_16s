# A script to check the manifest provided, the presence of the fastq files, and that dependencies are all setup
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
import pandas as pd
import datetime


logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
# Logging handler which catches EVERYTHING
logfile_name = "setup_"+str(datetime.datetime.now().date()) + \
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

# TODO add a check for pandas


def check_dependencies():
    # checks to see if qiime is installed, all software versions are saved to log
    try:
        command = subprocess.run(
            ['qiime info'], stdout=PIPE, stderr=PIPE, shell=True)
        logger.info(command.stdout)
    except FileNotFoundError:
        logger.info(
            "It appears qiime2 is possibly not installed or loaded, but this error appears when everything works ok too")
        # exit(1)
    logger.info('all software installed and ready to go')

# TODO determine a way to verify the dir provided is the same as the manifest


def check_provided_dir(manifest, seq_dir):
    try:
        read_manifest = pd.read_table(manifest, index_col=0, sep='\t')
    except FileNotFoundError:
        logger.critical("that manifest file does not exist")
        exit(1)

    path = read_manifest.iloc[[1][0]][0]
    manifest_seq_dir = os.path.split(os.path.split(path)[0])[1]

    if seq_dir != manifest_seq_dir:
        logger.critical(
            "The provided directory is not the same as the directory outlined in the manifest; This check will not function correctly")
        exit(1)


def verify_manifest(manifest, seq_dir):
    try:
        read_manifest = pd.read_table(manifest, index_col=0, sep='\t')
    except FileNotFoundError:
        logger.critical("that manifest file does not exist")
        exit(1)

    # sets current dir and finds the fastq and fastq.gz files in the current directory
    p = Path.cwd()
    list_of_fastq = list(p.glob(seq_dir + '/*.fastq'))
    list_of_gz = list(p.glob(seq_dir+'/*.fastq.gz'))

    fastq_files = []
    gz_files = []
    found = []
    missing = []

    # pulls only the filename and saves to a list
    for item in list_of_fastq:
        filename = os.path.split(item)[1]
        fastq_files.append(filename)

    for item in list_of_gz:
        filename = os.path.split(item)[1]
        gz_files.append(filename)
    if read_manifest.columns[0] == 'forward-absolute-filepath':
        # iterates over the forward reads and then the reverse reads to check to make sure they are all accounted for
        try:
            for item in read_manifest['forward-absolute-filepath']:
                filename = os.path.split(item)[1]
                if filename in fastq_files:
                    logger.info(filename + ' found')
                    found.append(filename)
                else:
                    if filename in gz_files:
                        logger.info(filename + ' found')
                        found.append(filename)
                    else:
                        logger.info(filename + ' missing')
                        missing.append(filename)
        except KeyError:
            logger.info('single read project')

        # try except in the case that the user only has single end reads.
        try:
            for item in read_manifest['reverse-absolute-filepath']:
                filename = os.path.split(item)[1]
                if filename in fastq_files:
                    logger.info(filename + ' found')
                    found.append(filename)
                else:
                    if filename in gz_files:
                        logger.info(filename + ' found')
                        found.append(filename)
                    else:
                        logger.info(filename + ' missing')
                        missing.append(filename)

        except KeyError:
            logger.info("looking for forward only reads")
    else:
        # this case is if there are only single reads and after which we can figure that the manifest file is wrong
        try:
            for item in read_manifest['absolute-filepath']:
                filename = os.path.split(item)[1]
                if filename in fastq_files:
                    logger.info(filename + ' found')
                    found.append(filename)
                else:
                    if filename in gz_files:
                        logger.info(filename + ' found')
                        found.append(filename)
                    else:
                        logger.info(filename + ' missing')
                        missing.append(filename)
        except KeyError:
            logger.critical("headings in the manifest appear to be incorrect")
            exit(1)

    # if missing is not an empty list, i.e. a file listed in the manifest is not detected, it raises an error and
    # creates a list for the user
    if missing != []:

        logger.error("files are missing, please see missing.csv to correct")

        with open('missing.csv', 'w', newline='') as csvfile:
            fieldnames = ['filename']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

            writer.writeheader()

            for filename in missing:
                writer.writerow({'filename': filename})

        exit(0)

    logging.info("the manifest called: " + manifest +
                 " is valid and ready to go")


def main(arg):
    check_dependencies()
    # TODO can come back and make this user-editable, but for right now I will hard-code it.
    check_provided_dir(arg.manifest_name, arg.seq_dir)
    verify_manifest(arg.manifest_name, arg.seq_dir)


if __name__ == "__main__":
    # Build Argument Parser in order to facilitate ease of use for user
    parser = argparse.ArgumentParser(
        description="Checks dependencies are installed, and validates a manifest file for qiime 2")
    parser.add_argument('-n', '--name', action='store', required=True,
                        help="name for the manifest, typically manifest.tsv", dest='manifest_name')
    parser.add_argument('-d', '--dir', action='store', required=True,
                        help="directory where sequences are stored", dest="seq_dir")
    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s 1.0')

    args = parser.parse_args()
    # print(args)
    main(args)
