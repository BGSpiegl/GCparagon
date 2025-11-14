# WARNING: it is recommended to pull the latest version of GCparagon container from cloud.sylabs.io instead of
# building the current code into a container (newer changes might not be reflected in the def file)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# BUILDING AND PUSHING TO REPO
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# to build a sif file, run:
sudo singularity build gcparagon_<VERSION>.sif gcparagon-mainBranch.def

# after building the image, you need ot sign it:
singularity sign gcparagon_<VERSION>.sif

# if you want to push to a repo, you need to make sure that your access token is not expired.
# Best - create a new one at https://cloud.sylabs.io/tokens. For this you need to run:
singularity remote login
# (and then paste the newly created access token)

# finally, push the sif to your repository:
singularity push gcparagon_<VERSION>.sif library://<YOUR_USER>/gcparagon/<SIF_NAME>:latest


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# PULLING AND TESTING
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# pull the latest image from cloud.sylabs.io (recommended; requires installed singularity)
singularity pull library://bgspiegl/gcparagon/gcparagon_<VERSION>
# (you can also specify the long version including the Unique ID which is the file's sha26 checksum)

# after pulling, verify the image:
singularity verify -a gcparagon_<VERSION>.sif

# build image using def file (not recommended)
sudo singularity build gcparagon_<VERSION>.sif gcparagon-mainBranch.def

# test with home paths only (no paths bound)
singularity run gcparagon_<VERSION>.sif --bam /home/<USER>/<INPUT>.bam --temporary-directory /home/<USER>/tmp --out-dir /home/<USER>/test_output --preset 1 --threads 4 --reference-genome-build hg38

# test with binding mnt directory and wrong reference genome build:
singularity run -B /mnt gcparagon_<VERSION>.sif --bam /mnt/<INPUT>.bam --temporary-directory /mnt/tmp --out-dir /mnt/test_output --preset 1 --threads 4 --reference-genome-build hg19

# test BAM outputting
singularity run -B /mnt gcparagon_<VERSION>.sif --bam /mnt/<INPUT>.bam --temporary-directory /mnt/tmp --out-dir /mnt/test_output --preset 1 --threads 4 --output-bam
