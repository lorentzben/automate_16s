#A method to use the 5 number summary generated by qiime demux summarize to determine what cutoffs are nessecary, requires pandas
#TODO we need to find the qiime file and unpack it so we can get the summary out of it to do this processing.
import csv
import pandas as pd 
import numpy as np

input_file = 'forward-seven-number-summaries.csv'
#input_file = 'moving_pics.csv'

summary = pd.read_csv(input_file, index_col=0,)

mean_qual = summary[4:5]

average_qual = np.round(mean_qual.mean(axis=1), 0)
mean_qual_vals = np.array(mean_qual)[0]


if int(average_qual) < 30:
    print("The Average Quality of these sequences may be a concern would you like to continue?")
    exit(0)


for i in range(0, len(mean_qual_vals)):
    if mean_qual_vals[i] >= int(average_qual):
        right_cutoff = i+1
        break
for i in range(0,len(mean_qual_vals)):
    if mean_qual_vals[len(mean_qual_vals)-1-i] >= int(average_qual):
        left_cutoff = len(mean_qual_vals)-i
        break

print(right_cutoff)
print(left_cutoff)

with open('cutoffs.csv', 'w', newline='') as csvfile:
    fieldnames = ['cutoff', 'value']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

    writer.writeheader()
    writer.writerow({'cutoff': 'right', 'value': right_cutoff})
    writer.writerow({'cutoff': 'left', 'value': left_cutoff})
    writer.writerow({'cutoff': 'filename', 'value': input_file})