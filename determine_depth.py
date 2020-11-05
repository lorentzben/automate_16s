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

print("total count" + str(total_count))

feature_array = np.array(features)

for i in range(len(feature_array)-1, -1,-1):
    sampling_depth = feature_array[i][0]
    perc_features_retain = ((sampling_depth * (i+1))/total_count)
    if perc_features_retain > .22:
        sampling_depth = feature_array[i][0]
        print("sampling depth: " + str(sampling_depth) + " % features retained: " + str(round(perc_features_retain,3)) + " samples retained: " + str(i))
        break


with open('sampling_depth.csv', 'w', newline='') as csvfile:
    fieldnames = ['stat', 'value']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

    writer.writeheader()
    writer.writerow({'stat': 'sampling_depth' , 'value': sampling_depth})
    writer.writerow({'stat': '%_features_retained', 'value': perc_features_retain})
    writer.writerow({'stat': 'filename', 'value': input_file})

        
    
    