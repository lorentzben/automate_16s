#A script to determine the --p-sampling-depth parameter for the core metrics phylogenetic method in qiime2  requires pandas
#TODO we need to find the qiime file and unpack it so we can get the summary out of it to do this processing.
import csv
import pandas as pd 
import numpy as np

#the stat we will use to determine cutoff will be aiming for the equation of (Sampling depth * samples)/(sum of all samples) between
#22 and 26%, we will calc this by adding in samples one by one until that range is met. 

#this comes from the compressed archive table.qzv which comes from feature-table summarize
#the dir structure is UUID/data/sample-frequency-detail.csv
input_file = 'sample-frequency-detail.csv'

features = pd.read_csv(input_file, index_col=0,header=None)

total_count = sum(features[1])
sampling_depth = 0 
perc_features_retain = 0.0

print(total_count)
