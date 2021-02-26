#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate python2
lefse_files=`ls *lefse_formatted.txt`
for eachfile in $lefse_files
do
    format_input.py $eachfile lefse_formatted.in -c 2 -u 1
    run_lefse.py lefse_formatted.in lefse_result.res 
    plot_res.py lefse_result.res "${eachfile::-4}"_res.png --dpi 100 
    plot_cladogram.py lefse_result.res "${eachfile::-4}"_cladogram.png --format png --dpi 200
    cp lefse_result.res "${eachfile::-4}"_result.res
done
conda deactivate