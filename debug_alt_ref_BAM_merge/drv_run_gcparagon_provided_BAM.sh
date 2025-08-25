#!/usr/bin/env bash

set -e -o pipefail

# ACTIVATE CONDA ENVIRONMENT
eval "$(conda shell.bash hook)"
conda activate GCparagon

output_dir='/home/servitorbeta/Downloads/VILLANUEVA_DEBUG-gcparagon/DEBUG_output'
input_BAM='/home/servitorbeta/Downloads/VILLANUEVA_DEBUG-gcparagon/CF24_0064.rh.bam'
temp_dir='/home/servitorbeta/temp_dir'
gcparagon_sif_path='/home/servitorbeta/GitHub_cloned_repos/GCparagon/singularity_definition_file/gcparagon_0.6.13-ubuntu-22_04-container_latest.sif'

# PYTHON§ SCRIPT MODE
# run GCparagon normally to verify it works with the BAM file
direct_output="${output_dir}/python3_call"
gcparagon --temporary-directory "${temp_dir}" --bam "${input_BAM}" --out-dir "${direct_output}" --output-bam --reference-genome-build hg38 --threads 3 --preset 1
# SINGULARITY
# use most recent singularity image file
singularity_output="${output_dir}/sif_call"
singularity run "${gcparagon_sif_path}" --temporary-directory "${temp_dir}" --bam "${input_BAM}" --out-dir "${singularity_output}" --output-bam --reference-genome-build hg38 --threads 3 --preset 1
