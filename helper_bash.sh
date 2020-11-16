#!/usr/bin/env bash
$manifest = "manifest.tsv"
$metadata = "metadata.tsv"
$seqs = "seq"

#this must be a column name in the metadata.tsv file 
$var_of_interest = "treatment"

module load Python/3.8.2-GCCcore-8.3.0
module load Nextflow/20.04.1
module load QIIME2/2020.6

python3 -m pip install --user pandas
python3 setup.py -n $manifest -d $seqs

python3 automate_16.py -n $manifest -m $metadata -i $var_of_interest