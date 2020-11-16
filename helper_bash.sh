#!/usr/bin/env bash
MANIFEST="manifest.tsv"
METADATA="metadata_fake.tsv"
SEQS="demultiplexed_seqs"

#this must be a column name in the metadata.tsv file 
VAR_OF_INTEREST="treatment"

module load Python/3.8.2-GCCcore-8.3.0
module load Nextflow/20.04.1
module load QIIME2/2020.6

python3 -m pip install --user pandas
python3 setup.py -n $MANIFEST -d $SEQS

python3 automate_16.py -n $MANIFEST -m $METADATA -i $VAR_OF_INTEREST