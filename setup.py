#A script to check the manifest provided, the presence of the fastq files, and that dependencies are all setup
import subprocess
from subprocess import PIPE
import os
from os import path
import logging
from pathlib import Path
from pathlib import PurePath
import pandas as pd 
import csv


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

def check_dependencies():
    #checks to see if nextflow is installed, version is saved to log
    try:
        command = subprocess.run(['nextflow -v'], stdout=PIPE, stderr=PIPE)
        logger.info("nextflow installed: " + command.stdout)
    except FileNotFoundError:
        logger.info("It appears nextflow is possibly not installed or loaded, but this error appears when everything works ok too")
        #exit(1)
    #checks to see if qiime is installed, all software versions are saved to log
    try:
        command = subprocess.run(['qiime info'], stdout=PIPE, stderr=PIPE)
        logger.info(command.stdout)
    except FileNotFoundError:
        logger.info("It appears qiime2 is possibly not installed or loaded, but this error appears when everything works ok too")
        #exit(1)
    logger.info('all software installed and ready to go')

def verify_manifest(manifest):
    read_manifest = pd.read_table(manifest, index_col=0, sep='\t')

    # sets current dir and finds the fastq and fastq.gz files in the current directory
    p = Path.cwd()
    list_of_fastq = list(p.glob('**/*.fastq'))
    list_of_gz = list(p.glob('**/*.fastq.gz'))

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

    # iterates over the forward reads and then the reverse reads to check to make sure they are all accounted for
    for item in read_manifest['forward-absolute-filepath']:
        filename = os.path.split(item)[1]
        if filename in fastq_files:
            found.append(filename)
        if filename in gz_files:
            found.append(filename)
        else:
            missing.append(filename)

    # try except in the case that the user only has single end reads. 
    try:
        for item in read_manifest['reverse-absolute-filepath']:
            filename = os.path.split(item)[1]
            if filename in fastq_files:
                found.append(filename)
            if filename in gz_files:
                found.append(filename)
            else:
                missing.append(filename)

    except KeyError:
        logger.info("looking for forward only reads")

    # if missing is not an empty list, i.e. a file listed in the manifest is not detected, it raises an error and 
    # creates a list for the user
    if missing != []:

        logger.critical("files are missing, please see missing.csv to correct")

        with open('missing.csv', 'w', newline='') as csvfile:
            fieldnames = ['filename']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

            writer.writeheader()

            for filename in missing:
                writer.writerow({'filename': filename})

        exit(0)

    logging.info("the manifest called: " +manifest+ " is valid and ready to go")

def error():
    print("something is wrong and I will try to tell you what the problem is")

def main():
    check_dependencies()
    # TODO can come back and make this user-editable, but for right now I will hard-code it. 
    verify_manifest("manifest.tsv")

main()
