#A script to check the manifest provided, the presence of the fastq files, and that dependencies are all setup
import subprocess
from subprocess import PIPE
import os
from os import path
import logging
from pathlib import Path
from pathlib import PurePath

def setup_logging():
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
        logger.critical("It appears nextflow is not installed or loaded")
        exit(1)
    #checks to see if qiime is installed, all software versions are saved to log
    try:
        command = subprocess.run(['qiime info'], stdout=PIPE, stderr=PIPE)
        logger.info(command.stdout)
    except FileNotFoundError:
        logger.critical("It appears qiime2 is not installed or loaded")
        exit(1)
    logger.info('all software installed and ready to go')

def verify_manifest(mainfest):
    read_manifest = pd.read_csv(maifest, index_col=0,)

    p = Path.cwd()
    list_of_fastq = list(p.glob('**/*.fastq'))
    list_of_gz = list(p.glob('**/*.fastq.gz'))
    found = []
    not_found = []

    #it will make more sense to iterate over the manifest as opposed to the files in the directory
    for item in list_of_fastq:
        print(item)
        found.append(item)
    for item in list_of_gz:
        print(item)
        found.append(item)

    if len(found) not len(read_manifest):
        logger.critical("The files in the manifest were not found in this directory")
    else:
        logging.info("the manifest called: " +manifest+ " is valid and ready to go")

def error():
    print("something is wrong and I will try to tell you what the problem is")