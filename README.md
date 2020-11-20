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

## Known R Issue
If you are using a fresh R installation and are on a computer without admin/sudo access, you will have to set up an R profile where your packages are installed. This script does not setup the library, this will need to be done manually. The easiest way is to open the R shell and attemp to install a package where it will prompt you through the library creation. 
```shell
$ R
> install.packages("ggplot2")
Warning in install.packages("ggplot2") :
  'lib = "/apps/eb/R/3.6.2-foss-2019b/lib64/R/library"' is not writable
Would you like to use a personal library instead? (yes/No/cancel) yes
Would you like to create a personal library
‘~/R/x86_64-pc-linux-gnu-library/3.6’
to install packages into? (yes/No/cancel) yes
--- Please select a CRAN mirror for use in this session ---
Secure CRAN mirrors
...
Selection: 1
```
Once this is done, the script should run correctly, if you would like to change the mirror that r gets packages from the first couple line in report.Rmd can be updated. 
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