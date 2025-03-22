#!/bin/sh
# Calls Docker container for segmenting prostate from T2 volume

path_seg_dir=$1 # Path to directory where segmentation output will be saved

# Mirror contents from path_seg_dir (main container) to shared-volume (sidecar container)
rsync -avq --delete $path_seg_dir/ /shared-volume/

# TODO change this to an non-user account
echo $GHCR_LOGIN | docker login ghcr.io -u stephen-bluerocksoft --password-stdin

echo docker run -i --rm -u 0:0 -v /shared-volume:/input -v /shared-volume:/output ghcr.io/precision-health/prostate-segmentation:dev
docker run -i --rm -u 0:0 -v /shared-volume:/input -v /shared-volume:/output ghcr.io/precision-health/prostate-segmentation:dev

# After the job finishes, sync back the results from shared-volume to path_seg_dir
rsync -avq --delete /shared-volume/ $path_seg_dir/
