#!/bin/bash
# bash script: execute to download the standard analysis set reference file (without the hs38d1 decoys) for hg38
# already in 2bit format
set -e  # exit as soon as a command fails

dwnld_script_dir="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"  # change into target dir
two_bit_ref_url="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.2bit"
echo "changing directory to ${dwnld_script_dir}..."
cd "${dwnld_script_dir}" || exit 1
wget --ftps-implicit --auth-no-challenge --no-cookies --no-directories --no-host-directories "${two_bit_ref_url}"

# to download a reference including the hs38d1 decoy sequences, use the code below:
# ----------------------------------------------------------------------------------------------------------------------
#decoy_ref="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz"
#wget --ftps-implicit --auth-no-challenge --no-cookies --no-directories --no-host-directories "${decoy_ref}"
#gzip -d GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz
#fasta_to_two_bit_path="$( which faToTwoBit )"
#"${fasta_to_two_bit_path}" -noMask -long -ignoreDups GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.2bit
# # 2bit conversion command above ignores lowercase masking; uses 64 bit for index; ignores duplicate sequence names (uses first instance found only), the latter which should not occur

# The hg19 analysis set can be downloaded like this:
# ----------------------------------------------------------------------------------------------------------------------
#set -e  # exit as soon as a command fails
#dwnld_script_dir="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"  # to change into target dir later
# two_bit_ref_url="https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit"
# or analysis set which must be converted to 2bit format (WARNING: I had a version exception raised with twobitreader in Pyton3):
#two_bit_ref_url="https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/analysisSet/hg19.p13.plusMT.no_alt_analysis_set.fa.gz"
#echo "changing directory to ${dwnld_script_dir}..."
#cd "${dwnld_script_dir}" || exit 1
#wget --ftps-implicit --auth-no-challenge --no-cookies --no-directories --no-host-directories "${two_bit_ref_url}"
# usually you can also download a chromosome sizes file: wget --ftps-implicit --auth-no-challenge --no-cookies --no-directories --no-host-directories "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes"
# For more information regarding te reference genome build, read the information provided here: https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/analysisSet/
# If you downloaded a FASTA file (not a 2bit version), the downloaded fasta.gz (fa.gz) file must be (unzipped and) converted to 2bit for use with GCparagon:
#gunzip hg19.p13.plusMT.no_alt_analysis_set.fa.gz && faToTwoBit -noMask -long -ignoreDups hg19.p13.plusMT.no_alt_analysis_set.fa hg19.p13.plusMT.no_alt_analysis_set.2bit