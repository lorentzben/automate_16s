Automate 16s
-------------------------------------------------
This project's aim is focused on automating processes of 16s rRNA analysis that are possible to make the process more accessible. 

## Prerequisities
* Linux based system
* qiime 2
* nextflow v
* python3 

## Install

```shell
$ git clone https://github.com/lorentzben/automate_16s.git
```

After cloning a folder called MiniProject will be created. Inside will be: 


## Running 
A manifest is required for analysis if using demultiplexed reads, i.e. each sample has two fastq files if using paired-end reads. 

The manifest will have three columns and be saved as a .tsv:
```shell
sample-id     forward-absolute-filepath       reverse-absolute-filepath
sample-1      $PWD/some/filepath/sample0_R1.fastq.gz  $PWD/some/filepath/sample1_R2.fastq.gz
sample-2      $PWD/some/filepath/sample2_R1.fastq.gz  $PWD/some/filepath/sample2_R2.fastq.gz
```
for more information see the [qiime manifest tutorial](https://docs.qiime2.org/2020.8/tutorials/importing/)

A design document, which is also nessecary for analysis, will have sample names, and treatment details. It should be in QIIME 2 format which can be found [here](https://docs.qiime2.org/2020.8/tutorials/metadata/) and verified with [Keemei](https://keemei.qiime2.org/). 

```shell
$ python3 
```

## Output
TODO generate results for alpha diversity of a sample dataset

## Version
* Version 1.0

## Author
* Ben Lorentz

## Future Plans
* Perfom analysis from raw data to general figures with minimal user input
* Produce clear logfiles to ensure users can report parameters used
* Offer 'tweaking' of parameters 