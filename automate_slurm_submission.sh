#!/bin/bash

#SBATCH --partition=batch
#SBATCH --job-name=16s_analysis
#SBATCH --nodes=2
#SBATCH --ntasks=8
#SBATCH --time=24:00:00
#SBATCH --mem=16gb

#Replace this with your UGA email to get notified on completion
#SBATCH --mail-user= ""
#SBATCH --mail-type=BEGIN,END,FAIL

#These values need to be be replaced with your files
MANIFEST="manifest.tsv"
METADATA="metadata_fake.tsv"
SEQS="demultiplexed_seqs"

#this must be a column name in the metadata.tsv file 
VAR_OF_INTEREST="body-site"

module load Python/3.8.2-GCCcore-8.3.0
module load QIIME2/2020.6
module load R/4.0.0-foss-2019b
module load GDAL/3.0.2-foss-2019b-Python-3.7.4
module load PROJ/6.2.1-GCCcore-8.3.0
module load GEOS/3.8.0-GCC-8.3.0-Python-3.7.4

python3 -m pip install --user pandas
CODE= python3 setup.py -n $MANIFEST -d $SEQS

if [[ $CODE -eq 1 ]]
then
   echo "There was an issue validating please review the logfile"
   exit
fi

CODE2= python3 automate_16.py -n $MANIFEST -m $METADATA -i $VAR_OF_INTEREST

if [[ $CODE2 -eq 1 ]]
then
   echo "There was an issue please review the logfile"
   exit
fi