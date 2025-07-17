#!/bin/sh
# Calls Docker container for segmenting prostate from T2 volume

path_input=$1 # Path to directory with T2 volume (as a nifti file)
path_output=$2 # Path to directory where contour will be saved 

sudo docker run --ipc="host" -v "$path_input":/data_in -v "$path_output":/data_out localhost/autoseg_prostate
