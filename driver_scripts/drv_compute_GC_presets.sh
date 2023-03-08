#!/bin/bash
#SBATCH --job-name=GCparagon                         # Job name
#SBATCH --nodes=1                                    # Run computation on a single NODE (i.e. 'server')
#SBATCH --mem=20G                                    # memory to use
#SBATCH --ntasks=1                                   # Run a single task (you can define tasks and allocate specific resources to them below in the command section)
#SBATCH --cpus-per-task=12                           # Number of CPUs to use for multithread/processing job
#SBATCH --output=GCparagon_dev-%j.log                # Standard output log
#SBATCH --error=GCparagon_dev-%j.err                 # Error log
# ------------------------------ COMMAND SECTION -------------------------------------------
# EXECUTION: either run this script with sbatch on a slurm cluster OR just run in a bash shell
# exit code 1: script exits with 1 if a path definition was incorrect
# exit code 2: script exits with 2 if no input data could be found in 'preset_computation'

echo "ENTERED PRESET COMPUTATION SCRIPT"

# CONDA WARNING - the conda that you use must be the one that is referenced in your .bashrc file!
#                  (and the env must be from that conda)
# for conda environment activation to work, you must source the user's .bashrc file' first:
eval "$(conda shell.bash hook)"
# Still old conda will be used. Mamba install does not seem to work
conda activate GCparagon_py3.10  # NAME OF CONDA ENV GOES HERE

echo "conda environment activated"

# analysis definitions
n_processes=12
python3_path="$(which python3)"
script_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )  # assertion: script in driver_scripts
content_root="$(dirname "${script_dir}")"
profiling_script="benchmark_mprof.py"
profiling_script_path="${content_root}/${profiling_script}"
if [ ! "${profiling_script_path}" ]; then
  echo "ERROR: profiling script not found at: ${profiling_script_path}"
  exit 1
fi
GCparagon_script="${content_root}/GCparagon.py"
if [ ! "${GCparagon_script}" ]; then
  echo "ERROR: GCparagon.py script not found at: ${GCparagon_script}"
  exit 1
fi
TWOBIT_REF_GENOME="${content_root}/2bit_reference/hg38.analysisSet.2bit"
if [ ! "${TWOBIT_REF_GENOME}" ]; then
  echo "ERROR: 2bit version of hg38 reference genome not found at: ${TWOBIT_REF_GENOME}." \
"Please download using the EXECUTE_reference_download.sh script there to obtain it!"
  exit 1
fi
bam_dir="${content_root}/preset_computation"
if [ ! -d "${bam_dir}" ]; then
  echo "ERROR: test BAM directory not found at: ${bam_dir}"
  exit 1
fi
# fetch BAM files & exit if none found
test_bam_wildcard="${bam_dir}/*.bam"
test_bam_files=$(eval ls $test_bam_wildcard)
bam_counter=0
for _b in $test_bam_files
do
  ((bam_counter++))
done
if (( bam_counter == 0 )); then
  echo "ERROR: test BAM directory contained no BAM files"
  # exit 2
else
  echo "info: found ${bam_counter} BAM files for processing"
fi
# create output directory if it does not exist
test_output_dir="${bam_dir}/benchmark_results"
if [ ! -d "${bam_dir}" ]; then
  mkdir -p "${bam_dir}"
fi
# create a temporary directory
tmp_out_dir="$(mktemp -d)"

# TASK: profile each one of the B01, H01, C01, P01 samples from the publication
declare -a presets=(3 )  #1  2 for local test preset 1 suffices



for test_bam in $test_bam_wildcard
do
  sample_id="$(basename -- "${test_bam}" | cut -d '.' -f 1)"
  echo "sample name is: ${sample_id}"
  for preset in "${presets[@]}"
  do
    echo "preset is: ${preset}"
    preset_out_dir="${bam_dir}/${preset}"  # subdir with sample_id will be created by the GCparagon script
    if [ ! -d "${preset_out_dir}" ]; then
      mkdir -p "${preset_out_dir}"
    fi
    "${python3_path}" "${profiling_script_path}" --iter 2 --sampling-frequency 10 --output-path "${test_output_dir}" \
    --script "${GCparagon_script}" --use-parameter-preset "${preset}" --bam "${test_bam}" \
    --two-bit-reference-genome "${TWOBIT_REF_GENOME}" --out-dir "${preset_out_dir}" \
    --temporary-directory "${tmp_out_dir}" --write-chunk-exclusion --threads "${n_processes}"
  done
done
