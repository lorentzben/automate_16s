#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate python2

pip install graphlan
pip install export2graphlan
pip install hclust2
pip install biom-format==2.1.7
pip install h5py
conda install graphlan
pip install numpy
pip install qiime
conda install -c biobakery lefse

LEFSE="$(which lefse.py)"
LEFSE_DIR="${LEFSE::-8}"
cp plot_res.py $LEFSE_DIR
cp plot_cladogram.py $LEFSE_DIR


conda deactivate
