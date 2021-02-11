#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate python2
format_input.py lefse_formatted.txt lefse_formatted.in -c 2 -u 1
run_lefse.py lefse_formatted.in lefse_result.res 
plot_res.py lefse_result.res lefse_analysis.png
plot_cladogram.py lefse_result.res lefse_cladogram.png --format png

conda deactivate