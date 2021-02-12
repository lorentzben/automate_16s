#! /usr/bin/env bash

dt=$(date '+%d-%m-%Y_%H.%M.%S');
module load R/4.0.0-foss-2019b
module load GDAL/3.0.2-foss-2019b-Python-3.7.4

Rscript -e "rmarkdown::render('report.Rmd', output_file='report_$dt.html', clean=TRUE)"