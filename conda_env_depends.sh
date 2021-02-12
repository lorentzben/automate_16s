#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate python2

pip install graphlan
pip install export2graphlan
pip install hclust2
pip install biom-format==2.1.7
pip install h5py
conda install graphlan

conda deactivate
