#!/usr/bin/env bash

module load Python/3.8.2-GCCcore-8.3.0
module load Nextflow/20.04.1

python3 setup.py

nextflow run automate_16s.nf