#!/bin/bash

# WARNING: if pip . install is used, both hg19 and hg38 2bit reference files MUST be downloaded first!

# The hg19 analysis set can be downloaded like this:
# ----------------------------------------------------------------------------------------------------------------------
set -e  # exit as soon as a command fails
dwnld_script_dir="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"  # to change into target dir later
two_bit_ref_url="https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit"
# or analysis set which must be converted to 2bit format (WARNING: I had a version exception raised with twobitreader in Pyton3):
#two_bit_ref_url="https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/analysisSet/hg19.p13.plusMT.no_alt_analysis_set.fa.gz"
echo "changing directory to ${dwnld_script_dir}..."
cd "${dwnld_script_dir}" || exit 1
wget --ftps-implicit --auth-no-challenge --no-cookies --no-directories --no-host-directories "${two_bit_ref_url}"
# usually you can also download a chromosome sizes file: wget --ftps-implicit --auth-no-challenge --no-cookies --no-directories --no-host-directories "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes"
# For more information regarding the reference genome build, read the information provided here: https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/analysisSet/
# If you downloaded a FASTA file (not a 2bit version), the downloaded fasta.gz (fa.gz) file must be (unzipped and) converted to 2bit for use with GCparagon:
#gunzip hg19.p13.plusMT.no_alt_analysis_set.fa.gz && faToTwoBit -noMask -long -ignoreDups hg19.p13.plusMT.no_alt_analysis_set.fa hg19.p13.plusMT.no_alt_analysis_set.2bit