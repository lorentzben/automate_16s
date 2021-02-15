Automate 16s
-------------------------------------------------
This project's aim is focused on automating processes of 16s rRNA analysis that are possible to make the process more accessible. 

## Prerequisities
* Linux based system
* qiime 2
* python3 

## Install

```shell
$ git clone https://github.com/lorentzben/automate_16s.git
```

After cloning a folder called automate_16s will be created. Inside will be: 
* README.md
* setup.py
* automate_16.py
* automate_single.sh
* automate_slurm_submission.sh
* report.Rmd
* EXAMPLE_MANIFEST.tsv
* EXAMPLE_METADATA.tsv

automate_single.sh and/or automate_slurm_submission.sh must be edited with relevant information including filenames for manifest and metadata, email for notification of completed jobs, and the directory where fastq sequences are located. 

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

There are multiple ways of invoking the program
```shell
$ sbatch automate_slurm_submission.sh 

$ sbatch automate_single.sh

$ python3 automate_16s.py -n MANIFEST -m METADATA -i ITEM_OF_INTEREST

```

## Output

report.nb.html (can be downloaded from cluster to local machine using [winscp](https://winscp.net/eng/docs/task_download) or [scp](https://www.garron.me/en/articles/scp.html) to examine)
setup_YYYY_MM_DD_HH_MM_SSSSSS.log
automate_YY_MM_DD_HH_MM_SSSSSS.log

---

## Creating Manifest

To speed up the manifest creation the command below will list the compressed files in the directory and save this to a txt file that can be opened in a program like Excel to format in the shape of EXAMPLE_MANIFEST.tsv 

```shell
ls $SEQUENCE_DIR/*.fastq.gz > sequences.tsv

```

#### Useful files

core-metrics-results/

table.qza

## Example Sequences
To test installation of the program it can be helpful to complete a diagnostic test run with a subsample of a full dataset. I used the data from the MiSeq SOP for the program Mothur written by Dr. Patrick Schloss. The sequences can be found [here](https://mothur.s3.us-east-2.amazonaws.com/wiki/miseqsopdata.zip), I left the mock samples out for my testing. The example metadata and manifest files included are what I used when running these tests. 

## Version
* Version 1.0

## Author
* Ben Lorentz
* Benjamin.Lorentz@uga.edu

## Future Plans
* Offer 'tweaking' of parameters
* Automated Consolidation of logfiles