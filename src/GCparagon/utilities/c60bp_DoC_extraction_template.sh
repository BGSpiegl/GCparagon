#!/bin/bash

PYTHON_PATH="$(which python3)"
INPUT_TAGGED_BAM=''  # <<<---- YOUR BAM PATH HERE!
INPUT_SAMPLE_NAME=''  # <<<---- YOUR SAMPLE NAME HERE!
COVERAGE_REGIONS_BED=''  # <<<---- YOUR REGION DEFINING BED FILE PATH HERE! CAN BE CREATED FOR TSS COORDINATES USING
# THE 'HK.txt' AND 'PAU.txt' FILES IN COMBINATION WITH THE 'hg38_MANE.bed' FROM THE "accessory_files/TSS/" DIRECTORY
OUTPUT_DIRECTORY_ORIGINAL_COVERAGE='TSS_original'  # <<<---- YOUR FOLDER NAME FOR ORIGINAL COVERAGE OUTPUT HERE!
OUTPUT_DIRECTORY_CORRECTED_COVERAGE='TSS_corrected'  # <<<---- YOUR FOLDER NAME FOR CORRECTED COVERAGE OUTPUT HERE!


# -- EXECUTE DoC EXTRACTION -----------------------------------------
"${PYTHON_PATH}" analyse_coverage_pysam_presum_c60.py -bf "${INPUT_TAGGED_BAM}" -ibe "${COVERAGE_REGIONS_BED}" -c 4 -o "${OUTPUT_DIRECTORY_ORIGINAL_COVERAGE}" -on "${INPUT_SAMPLE_NAME}"
"${PYTHON_PATH}" analyse_coverage_pysam_presum_c60.py -bf "${INPUT_TAGGED_BAM}" -ibe "${COVERAGE_REGIONS_BED}" -c 4 -o "${OUTPUT_DIRECTORY_CORRECTED_COVERAGE}" -on "${INPUT_SAMPLE_NAME}" -gc True
