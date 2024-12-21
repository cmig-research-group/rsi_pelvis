#!/bin/sh
# Calls Docker container for segmenting prostate from T2 volume

path_seg_dir=$1 # Path to directory where segmentation output will be saved

sudo docker run -i --rm -u 0:0 -v "$path_seg_dir":/input:ro -v "$path_seg_dir":/output 696021503801.dkr.ecr.us-east-1.amazonaws.com/prostate-segmentation:dev
