# WARNING: it is recommended to pull the latest version of GCparagon container from cloud.sylabs.io instead of
# building the current code into a container (newer changes might not be reflected in the def file)

# build image using def file
sudo singularity build gcparagon_0.6.9.sif gcparagon-mainBranch.def

# test with home paths only
singularity run gcparagon_0.6.9.sif --bam /home/<INPUT>.bam --temporary-directory /home/tmp --out-dir /home/test_output --preset 1 --threads 4 --reference-genome-build hg38

# test with binding mnt directory and wrong reference genome build:
singularity run -B /mnt:/mnt gcparagon.sif --bam /mnt/<INPUT>.bam --temporary-directory /mnt/tmp --out-dir /mnt/test_output --preset 1 --threads 4 --reference-genome-build hg19