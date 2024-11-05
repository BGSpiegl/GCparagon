# WARNING: it is recommended to pull the latest version of GCparagon container from cloud.sylabs.io instead of
# building the current code into a container (newer changes might not be reflected in the def file)

# build image using def file
sudo singularity build gcparagon_0.6.12.sif gcparagon-mainBranch.def

# test with home paths only
singularity run gcparagon_0.6.12.sif --bam /home/<USER>/<INPUT>.bam --temporary-directory /home/<USER>/tmp --out-dir /home/<USER>/test_output --preset 1 --threads 4 --reference-genome-build hg38

# test with binding mnt directory and wrong reference genome build:
singularity run -B /mnt gcparagon_0.6.12.sif --bam /mnt/<INPUT>.bam --temporary-directory /mnt/tmp --out-dir /mnt/test_output --preset 1 --threads 4 --reference-genome-build hg19

# test BAM outputting
singularity run -B /mnt gcparagon_0.6.12.sif --bam /mnt/<INPUT>.bam --temporary-directory /mnt/tmp --out-dir /mnt/test_output --preset 1 --threads 4 --output-bam
