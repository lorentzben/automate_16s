#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate python2
format_input.py lefse_formatted.txt lefse_formatted.in -c 2 -u 1
run_lefse.py lefse_formatted.in lefse_result.res 
plot_res.py lefse_result.res lefse_analysis.png --dpi 150 
plot_cladogram.py lefse_result.res lefse_cladogram.png --format png
#TODO use export2graphlan to make a better cladogram
export2graphlan.py -i lefse_formatted.txt -o lefse_result.res -t tree.nwk -a annot_lefse.txt --external_annotations 2,3,4,5,6 --fname_row 0 --skip_rows 1 
graphlan_annotate.py --annot annot_lefse.txt tree.nwk annotated_tree_lefse.txt
graphlan.py annotated_tree_lefse.txt image_graph_lefse.png --format png --size 16 --dpi 300 --pad 0.2
conda deactivate