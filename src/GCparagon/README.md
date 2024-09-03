# GCparagon

```
_____________________________________________________________________________
v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v
                                                                             
 _.--'"`'--._    _.--'"`'--._    _.--'"`'--._    _.--'"`'--._    _.--'"`'--._
`.'|`|"':-.  '-:`.'|`|"':-.  '-:`.'|`|"':-.  '-:`.'|`|"':-.  '-:`.'|`|"':-   
.  | |  | |'.  '.  | |  | |'.  '.  | |  | |'.  '.  | |  | |'.  '.  | |  | |'.
 '.| |  | |  '.  '.| |  | |  '.  '.| |  | |  '.  '.| |  | |  '.  '.| |  | |  
.  `.:_ | :_.' '.  `.:_ | :_.' '.  `.:_ | :_.' '.  `.:_ | :_.' '.  `.:_ | :_.
  `-..,..-'       `-..,..-'       `-..,..-'       `-..,..-'       `-..,..-'  
                                                                             
 ██████╗  ██████╗██████╗  █████╗ ██████╗  █████╗  ██████╗  ██████╗ ███╗   ██╗
██╔════╝ ██╔════╝██╔══██╗██╔══██╗██╔══██╗██╔══██╗██╔════╝ ██╔═══██╗████╗  ██║
██║  ███╗██║     ██████╔╝███████║██████╔╝███████║██║  ███╗██║   ██║██╔██╗ ██║
██║   ██║██║     ██╔═══╝ ██╔══██║██╔══██╗██╔══██║██║   ██║██║   ██║██║╚██╗██║
╚██████╔╝╚██████╗██║     ██║  ██║██║  ██║██║  ██║╚██████╔╝╚██████╔╝██║ ╚████║
 ╚═════╝  ╚═════╝╚═╝     ╚═╝  ╚═╝╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝  ╚═════╝ ╚═╝  ╚═══╝
                                                                             
 _.--'"`'--._    _.--'"`'--._    _.--'"`'--._    _.--'"`'--._    _.--'"`'--._
`.'|`|"':-.  '-:`.'|`|"':-.  '-:`.'|`|"':-.  '-:`.'|`|"':-.  '-:`.'|`|"':-   
.  | |  | |'.  '.  | |  | |'.  '.  | |  | |'.  '.  | |  | |'.  '.  | |  | |'.
 '.| |  | |  '.  '.| |  | |  '.  '.| |  | |  '.  '.| |  | |  '.  '.| |  | |  
.  `.:_ | :_.' '.  `.:_ | :_.' '.  `.:_ | :_.' '.  `.:_ | :_.' '.  `.:_ | :_.
  `-..,..-'       `-..,..-'       `-..,..-'       `-..,..-'       `-..,..-'  
                                                                             
^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^
-----------------------------------------------------------------------------
```

## Contents:

-------------------------------------------------------------------------------------------------------------------
- [Description](#description)
  - [How to use cfDNA fragment weights](#how-to-use-cfdna-fragment-weights)
  - [Software Dependencies](#software-dependencies)
    - [Conda installation](#conda-installation)
    - [Singularity container](#singularity-container)
  - [Hardware Requirements](#hardware-requirements)
  - [Repository Structure](#repository-structure)
  - [Required Files](#required-files)


- [Installation](#installation)
  - [Singularity Installation Procedure](#singularity-installation-procedure)
  - [Conda Installation Procedure](#conda-installation-procedure)

- [Usage](#usage)
  - [Examples](#examples)
  - [Performance](#performance)
  - [Full Commandline Description](#full-commandline-description)
  - [Parameter Presets](#parameter-presets)

- [Result of GC Bias Correction](#result-of-gc-bias-correction)
  - [Output](#output)
  - [Genomic Region Preselection](#genomic-region-preselection)
  - [Optimization of Combining Preselected Regions](#optimization-of-combining-results-from-preselected-regions)
  - [Validation](#validation-v060)

- [Copyright](#copyright)
- [Software License](#software-license)
- [Contributors](#contributors)


## Description

-------------------------------------------------------------------------------------------------------------------

GCparagon is a Python commandline tool for the rapid calculation and correction of fragment length specific GC biases
in WGS cfDNA sequencing datasets for liquid biopsy applications. Code was developed for UNIX machines but can be also used on Windows or MAC using the available Apptainer/SingularityCE image (requires a working [singularity installation][singularity_installation]).

GCparagon takes an aligned BAM file as input and processes the alignments in predefined genomic intervals to estimate 
GC bias in the sample.
Important: the input BAM files must have been created using an alignment algorithm that conforms with the SAM format 
specification (e.g. [BWA MEM][bwa_mem]).
GCparagon uses the observed template length (TLEN column in BAM) to estimate fragment length.
(Preferably TLEN#1 as shown in [SAMv1.pdf][samtools_spec]).

The algorithm assigns weight values to cfDNA fragments based on their length and GC base count. Weights can either
be read as 'GC'-tags from alignments in the output BAM file (enable tagged BAM writing using the `--output-bam` flag),
or from one of the `*_gc_weights_*.txt.gz` output files.
The latter can be loaded in Python using the `numpy.loadtxt()` function or the `load_txt_to_matrix_with_meta()` function 
from [plot_distributions.py](src/GCparagon/utilities/plot_distributions.py).
The tag string can be redefined using the `--tag-name` parameter.

Latest version is v0.6.8


### How to use cfDNA fragment weights
Instead of counting fragment occurrences or their attributes, the user can sum the GC bias correction weights of these 
fragments to obtain an unbiased result for signal extraction. An example could be depth of coverage computation for 
specific groups of transcription start sites as shown for samples B01, C01, P01, and H01 in the section 
"[Result of GC Bias Correction](#result-of-gc-bias-correction)".
To this end, Faruk Kapidzic created a Pysam fork which can directly use the tags.


### Software Dependencies

- UNIX system (server or HPC cluster recommended)

The GCparagon commandline tool was tested on an Ubuntu 20.04.5 LTS and Ubuntu 22.04.4 LTS operating systems,
using a Python3.10 conda environment.

#### Conda installation
This software environment can be installed using the provided [GCparagon_py3.10_env.yml](conda_env/GCparagon_py3.10_env.yml) file.

Using the [GCparagon_py3.10_env.yml](conda_env/GCparagon_py3.10_env.yml) file, the following dependencies will be installed into the conda env named `GCparagon_py3.10`:
  - samtools=1.16
  - bedtools=2.30
  - python=3.10
  - pip=22.3
  - numpy=1.23
  - pysam=0.19
  - natsort=8.2
  - py2bit=0.3
  - cycler=0.11
  - pandas=1.5
  - scipy=1.9
  - ucsc-fatotwobit=377
  - twobitreader=3.1
  - plotly_express=0.4
  - python-kaleido=0.2
  - psutil=5.9
  - requests=2.28
  - memory_profiler
  - pybedtools
  - polars
  - scikit-learn
  - matplotlib
  - singularity

You can create the environment using the following command:

`conda env create -f GCparagon_py3.10_env.yml`

#### Singularity container
As an alternative, a singularity image can be either downloaded or created (using the singularity executable 
from the aforementioned conda environment from OUTSIDE conda) using the 
[gcparagon.def](singularity_definition_file/gcparagon.def) singularity definitions file. 
A tested and pre-built Apptainer/SingularityCE container can be pulled from `cloud.sylabs.io` into the current 
directory and verified using a functioning singularity installation:

`singularity pull library://bgspiegl/gcparagon/gcparagon-ubuntu-22_04-container:latest && 
singularity verify gcparagon-ubuntu-22_04-container_latest.sif`

Note: this requires a working installation of singularity. The easiest way to use singularity is via conda:

`conda install singularity`

If you get an error like the following "*Unable to get library client configuration: remote has no library client*",
it is likely that no remote was specified for your installation of Apptainer/SingularityCE. You might want to 
use the `remote` command to specify cloud.sylabs.io as remote endpoint to be able to pull the latest GCparagon container 
(see [here][singularity_remotes]).

GCparagon can then be run from the singularity container file like this 
(from within the singularity image file's parent directory):

`singularity run gcparagon.sif <YOUR PARAMETERS HERE>` 

NOTE: do not run singularity commands with the activated conda environment that includes singularity!
Rather use the absolute path to the singularity executable from outside conda like this:

`/home/<USERNAME>/miniforge-pypy3/bin/singularity run ...`

As a workaround, you can create a symlink to that executable in your local bin folder:

`sudo ln -s /home/<USERNAME>/miniforge-pypy3/bin/singularity /usr/local/bin`

After that, singularity is available via the expected `singularity` command.

IMPORTANT: programs run from inside Apptainer/SingularityCE images (i.e., GCparagon) do not have access to all 
directories of your file system per default. Examples of exceptions are the user's `$HOME` directory, `\sys`, `\proc`, 
and `\var\tmp`. To make an inaccessible path available to the GCparagon image, follow the steps detailed in the 
[Apptainer documentation][apptainer_files_doc] or the [SingularityCE documentation][singularity_binds_and_mounts_doc] 
(see the `--bind` Examples).

An example of using your local directories under `/mnt` would be:

`singularity run -B /mnt:/mnt gcparagon.sif <YOUR PARAMETERS HERE INCLUDING PATHS UNDER MNT>` 


### Hardware Requirements

- 12 cores are used per default; more cores are better; can be run on a single core
- for 12 cores: \> 8.5 GiB of RAM, 16 GB recommended (max. observed memory usage for the preset 2 computation was 8.5 GiB @ 24 cores. 
A preset 3 computation is expected to use more RAM!)
- SSD scratch drive for `--temporary-directory` with at least twice the input BAM file's size in free space for the 
tagging procedure

Computation time might increase significantly if hardware requirements are not met.
Computation may terminate if drive space and/or RAM size requirements are not met!

### Repository Structure

Instructions for FastA reference genome sequence download can be found 
[here](src/GCparagon/2bit_reference/EXECUTE_hg38_reference_download.sh) [and here for hg19](src/GCparagon/2bit_reference/EXECUTE_hg19_reference_download_hg19.sh).

Code for genomic regions exclusion list creation can be found in 
[GC_correction_exclusion_list_hg19_creation.info](accessory_files/GC_correction_exclusion_list_hg19_creation.info).
Code for genomic interval preselection for a currently unsupported build can be found in folder 
[accessory_files/genomic_interval_preselection-shifted16x_hg19](accessory_files/genomic_interval_preselection-shifted16x_hg19).

Users successfully creating reference files for other genome builds are kindly asked to create a branch + pull request, where 
they also provide these resources (i.e., a URL to the 2bit or FastA reference sequence file,
<GENOME_BUILD>.genome_file.tsv,
<GENOME_BUILD>_minimalExclusionListOverlap_1Mbp_intervals_33pcOverlapLimited.FGCD.bed, 
<GENOME_BUILD>_reference_GC_content_distribution.tsv,
and optionally a different reference fragment length distribution file 
<ISOLATION_PROTOCOL>_<SAMPLE_TYPE>_reference_fragment_length_distribution.tsv in case the plasmaSeq ccfDNA 
reference is inadequate for the analyzed samples).

Note: a separate fragment length distribution file is required for urine samples because of the higher
fragmentation of cfDNA in urine samples (compared to blood plasma cfDNA).

Results of the correction benchmark including the Griffin algorithm can be found in [validation](validation).

Results from the [profile_command.py](src/GCparagon/profile_command.py) script are stored in 
[preset_computation/benchmark_results](preset_computation/benchmark_results).


### Required Files

GCparagon requires a 2bit version of the reference genome sequence which was used to create the aligned, SAM format 
specification conforming input BAM file.
The reference genome used to create the 4 BAM files in plots can be downloaded using the 
[EXECUTE_reference_download_hg38.sh](src/GCparagon/2bit_reference/EXECUTE_hg38_reference_download.sh) bash script.
It downloads the hg38 lowercase-masked standard analysis set reference file in 2bit format from 
[https://hgdownload.soe.ucsc.edu][hg38_std_analysis_set].

Alternatively, you can download a hg38 reference genome file in FastA.gz format which is converted into the 2bit format
containing decoys from NCBI's FTP server at [ftp.ncbi.nlm.nih.gov][hg38_decoy_analysis_set]
(execute 
[EXECUTE_reference_download_hg19.sh](src/GCparagon/2bit_reference/EXECUTE_hg19_reference_download.sh) to download the hg19 2bit reference)

GCparagon uses preselected genomic regions for GC bias computation. These are provided for hg19 and for hg38 
via [BED file](accessory_files/hg38_minimalExclusionListOverlap_1Mbp_intervals_33pcOverlapLimited.FGCD.bed).
Please see [Genomic Region Preselection](#genomic-region-preselection) section for more information.




## Installation

-------------------------------------------------------------------------------------------------------------------

GCparagon can be used out of the box by running `python3 GCparagon.py` using an appropriate Python 3.10+ 
software environment with all [dependencies](#software-dependencies) installed. It is recommended though to install 
software dependencies via the provided [GCparagon_py3.10_env.yml](conda_env/GCparagon_py3.10_env.yml) file in the 
[conda_env](conda_env) subdirectory by following the installation steps described below.

NOTE: an issue was reported for Ubuntu 20 where the conda/mamba solvers could not create the environment by 
directly using the YAML file.
If you experience the same issue, it can possibly be solved by manually installing dependencies one-by-one.

The author recommends to use [mamba/micromamba][mamba install] for environment creation/resolving of dependencies.
In case you are using conda version 22 or above, your installation can be easily upgraded with the [mamba solver][libmamba solver].
Mamba can be added to an older existing [conda installation][conda install].

For a detailed list of dependencies (manual installation, not recommended!) please go to 
[Software Dependencies](#software-dependencies)

### Singularity Installation Procedure
See here: [Singularity container](#singularity-container)

### Conda Installation Procedure
First, move into the directory where you want to store the GCparagon code and clone the [GitHub][github repo] 
repository:

`git clone https://github.com/BGSpiegl/GCparagon`

After making sure that conda is available on your system and up to date, move inside the cloned repository

`cd GCparagon`

and create the GCparagon software environment using the [GCparagon_py3.10_env.yml](conda_env/GCparagon_py3.10_env.yml) 
file:

`mamba env create -f conda_env/GCparagon_py3.10_env.yml` 

OR

`conda env create -f conda_env/GCparagon_py3.10_env.yml`

Activate the new environment:

`conda activate GCparagon`

Download the 2bit reference files - if `pip .` install is used in the next step, both the hg19 and the hg38 2bit reference files MUST be 
already available! Download them using the information provided in `src/GCparagon/2bit_reference/EXECUTE_reference_download_hg38.sh` 
and `src/GCparagon/2bit_reference/EXECUTE_reference_download_hg19.sh`!

Run `pip` to make GCparagon directly executable from the 
console:

`pip install .`

After successful setup, GCparagon should be available via the `gcparagon` command. For a detailed help, type:
`gcparagon --help`

To make the reference build definition work with the gcparagon command, you might need to copy the contents of 
the `src/GCparagon/2bit_reference/` directory and the `accessory_files/` directory to the conda env lib 
site-packages folder of your user like this:

`mkdir ~/miniconda3/envs/GCparagon/lib/python3.10/site-packages/GCparagon/2bit_reference/ && 
cp src/GCparagon/2bit_reference/* ~/miniconda3/envs/GCparagon/lib/python3.10/site-packages/GCparagon/2bit_reference/`

and:

`mkdir ~/miniconda3/envs/GCparagon/lib/python3.10/accessory_files/ && 
cp accessory_files/*_minimalExclusionListOverlap_* accessory_files/*_reference_* ~/miniconda3/envs/GCparagon/lib/python3.10/accessory_files/`

This will be added to the setup in a future version, if possible.

See [Usage](#usage) for a complete explanation of commandline options.
Default output created by GCparagon is described [here](#result-of-gc-bias-correction).
There are several options available to alter plotting behaviour or to keep intermediate data created during
simulation rounds. Note that **per default, the tagged BAM file is _NOT_ output**.
To activate BAM output, use the `--output-bam` flag.

Be mindful of setting the `--temporary-directory` to a reasonable path!
(I.e. high IOPS hardware if available +
sufficient storage space available for tagged BAM etc.)

To gain access to the four BAM files used in the publication, go to EGA, create an account and ask for access to dataset
[EGAS00001006963].

To recreate the **tagged** BAM files and matrix visualisations for the three presets and four samples from [EGAS00001006963],
first download the 2bit version of the reference genome sequence which was used to create the four aligned BAM files.
The hg38 reference genome used to create the 4 BAM files in plots can be downloaded using the 
[EXECUTE_reference_download.sh](src/GCparagon/2bit_reference/EXECUTE_hg38_reference_download.sh) bash script:

`bash src/GCparagon/2bit_reference/EXECUTE_reference_download.sh`

After the download has finished, you can run the [driver script](driver_scripts/drv_compute_GC_presets.sh) from
**within the activated conda environment**:

`bash driver_scripts/drv_compute_GC_presets.sh`

You might want to do this inside a tmux session from which you can detach. Preset 3 computations will take around 
50 minutes for each sample.


## Usage

-------------------------------------------------------------------------------------------------------------------
Run the GCparagon.py script with installed [dependencies](#software-dependencies) using an appropriate Python3.10+ 
interpreter.

### Examples

#### Basic
The most basic call after downoading the 2bit version of the reference genome is as follows:

`python3 src/GCparagon/correct_GC_bias.py --bam <INPUT_BAM>`

OR:

`gcparagon --bam <INPUT_BAM>`
(available only if `pip install .` was run)

This minimalistic setup uses the parent directory of the input BAM fle as output directory.
The `-b`/`--bam` parameter is always required (BAM file path to hg38 aligned cfDNA paired-end sequencing reads).

#### Output tagged BAM file
To output a GC correction weights tagged BAM file, set the `--output-bam` flag:

`gcparagon --bam <INPUT_BAM> --output-bam`

#### Select reference genome build
To use hg19 reference genome-specific files for the bias computation, set the `-rgb`/`--reference-genome-build` parameter:

`gcparagon --bam <INPUT_BAM> --reference-genome-build hg19`

The option sets the following three parameters:
 - `-rtb`/`--two-bit-reference-genome`
 - `-c`/`--intervals-bed`
 - `-rgcd`/`--reference-gc-content-distribution-table`

These parameters can be redefined separately.

#### Temporary directory
It is recommended to set `--temporary-directory` to be located on SSD hardware:

`gcparagon --bam <INPUT_BAM> --output-bam --temporary-directory <PATH_TO_TEMP_DIR>`

If not set, the temporary directory will default to the output of Python's `tempfile.gettempdir()`. All created files 
are saved to the temporary directory first before being moved to the output directory after successful
script execution. 

#### Further options
Rich customization options are available:
To increase the number of logical cores used by GCparagon, use the `-t`/`--threads` flag:

`gcparagon --bam <INPUT_BAM> --output-bam --threads 24`

To get a quick estimate of the GCbias, the user can set a lower preset

`gcparagon --bam <INPUT_BAM> --threads 24 --preset 1`

The `--preset 2` setup is recommended though.


### Performance

A benchmark of GCparagon, preset2 (v0.5.5) against the Griffin algorithm highlighted the superior computation speed of GCparagon:

| Sample |    Griffin table output    | GCparagon table output | GCparagon tagged BAM output | 
|:------:|:--------------------------:|:----------------------:|:---------------------------:|
|        |          hh:mm:ss          |   hh:mm:ss; speedup    |      hh:mm:ss; speedup      |
|  P01   |          11:54:20          | 0:04:58; ***143.8x***  |    0:23:03; ***31.0x***     |
|  H01   |          06:51:43          |    0:05:14;  78.7x     |       0:20:09; 20.4x        |
|  C01   |          05:13:52          |    0:04:32;  69.2x     |       0:17:09; 18.3x        |
|  B01   |          04:51:39          |    0:04:43;  61.8x     |       0:15:32; 18.8x        |


Concerning the time to output of the correction or bias table,
**GCparagon was up to 144x but at least 62 times faster than Griffin**.
When comparing Griffin table output time to duration of GCparagon bias computation and tagged BAM output, 
GCparagon was up to 31x but at least 18x faster than Griffin.
The update to GCparagon v0.6.0 resulted in increased FGCD correction and cDoC to Griffin results
but also slightly increased the computation time compared to the v0.5.5 benchmark. 

The GC bias computation time depends linearly on the portion of the input data which is processed.
The average DoC of the visualized samples is between 10x and 30x.
For preset 1 and preset 2, the duration of bias computation was found to be no longer than 3 and 16 minutes 
respectively.
The plot below shows the computation time in relation to the preset choice.

![linregress_comp_time](https://github.com/BGSpiegl/GCparagon/blob/including_EGAS00001006963_results/preset_computation/benchmark_results/GCparagon_computation_time_presetsAndSamples_linRegress_SSD-ref.png?raw=true)

The variability of the preset 1 calculation time over the actual fragments processed measured for more than 600 samples
is shown in the figure below (from outdated v0.5.4).

![preset1_comp_time](https://github.com/BGSpiegl/GCparagon/blob/including_EGAS00001006963_results/preset_computation/benchmark_results/GCparagon_computation_time_preset1_vs_processed_fragments_cluster_683BAMs.png?raw=true)

The amount of consumed memory is roughly independent of the number of processed fragments
(result of two iterations shown; from outdated v0.5.4).
In general, memory consumption depends on the DoC of a sample and the chosen number of simulation rounds,i.e., the chosen parameter preset.
The user can expect a memory usage between 4 and 9 GiB for default settings (12 cores). Reducing cores reduces memory consumption.

![linregress_comp_time](https://github.com/BGSpiegl/GCparagon/blob/including_EGAS00001006963_results/preset_computation/benchmark_results/GCparagon_memory_consumption_presets_SSD-ref.png?raw=true)

Memory consumption over time can be visualized using the [profile_command.py](src/GCparagon/profile_command.py) script 
(figure: P01, preset 2; v0.6.0):

![memory_consumption_over_time](https://github.com/BGSpiegl/GCparagon/blob/including_EGAS00001006963_results/validation/00_result_plots_preset2_for_readme/memory_usage.png?raw=true)

Writing of tagged BAM files also uses multiprocessing. This step usually takes longer than the bias 
computation itself. A test run using 12 cores and parameter preset 1 for a 30 GB BAM file took 25 minutes 
(computing GC weights + writing tagged BAM file).

Always make sure that there is enough space on the drive(s) containing the temporary directory and the final output 
directory before running GCparagon with `--output-bam`!


### Full Commandline Description

```
The following argument is required: -b/--bam

usage: correct_GC_bias.py [-h] [--version] -b List[File] [List[File] ...]
                          [-rgb {hg38,hg19}] [-rtb File] [-c File] [-rgcd File]
                          [-ec File] [-cw File] [-wm File] [-p Integer] [-to]
                          [-rep Integer] [-mafl Integer] [-mifl Integer] [-t Integer]
                          [-rs Integer] [-sp File] [-nf Integer] [-mf Integer] [-anf]
                          [-ucmaf Float] [-dw DEFAULT_FRAGMENT_WEIGHT]
                          [-reto UNALIGNED_EXTRACTION_TIMEOUT] [-do] [-odm String]
                          [-ods Integer] [-sw] [-sk String] [-si Integer] [-v]
                          [-o File] [-tmp File] [-np] [-os] [-k] [-ob] [-our]
                          [-fp Integer] [-tg String] [-wie] [-nfp] [-sf]

_____________________________________________________________________________
v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v
                                                                             
 _.--'"`'--._    _.--'"`'--._    _.--'"`'--._    _.--'"`'--._    _.--'"`'--._
`.'|`|"':-.  '-:`.'|`|"':-.  '-:`.'|`|"':-.  '-:`.'|`|"':-.  '-:`.'|`|"':-   
.  | |  | |'.  '.  | |  | |'.  '.  | |  | |'.  '.  | |  | |'.  '.  | |  | |'.
 '.| |  | |  '.  '.| |  | |  '.  '.| |  | |  '.  '.| |  | |  '.  '.| |  | |  
.  `.:_ | :_.' '.  `.:_ | :_.' '.  `.:_ | :_.' '.  `.:_ | :_.' '.  `.:_ | :_.
  `-..,..-'       `-..,..-'       `-..,..-'       `-..,..-'       `-..,..-'  
                                                                             
 ██████╗  ██████╗██████╗  █████╗ ██████╗  █████╗  ██████╗  ██████╗ ███╗   ██╗
██╔════╝ ██╔════╝██╔══██╗██╔══██╗██╔══██╗██╔══██╗██╔════╝ ██╔═══██╗████╗  ██║
██║  ███╗██║     ██████╔╝███████║██████╔╝███████║██║  ███╗██║   ██║██╔██╗ ██║
██║   ██║██║     ██╔═══╝ ██╔══██║██╔══██╗██╔══██║██║   ██║██║   ██║██║╚██╗██║
╚██████╔╝╚██████╗██║     ██║  ██║██║  ██║██║  ██║╚██████╔╝╚██████╔╝██║ ╚████║
 ╚═════╝  ╚═════╝╚═╝     ╚═╝  ╚═╝╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝  ╚═════╝ ╚═╝  ╚═══╝
                                                                             
 _.--'"`'--._    _.--'"`'--._    _.--'"`'--._    _.--'"`'--._    _.--'"`'--._
`.'|`|"':-.  '-:`.'|`|"':-.  '-:`.'|`|"':-.  '-:`.'|`|"':-.  '-:`.'|`|"':-   
.  | |  | |'.  '.  | |  | |'.  '.  | |  | |'.  '.  | |  | |'.  '.  | |  | |'.
 '.| |  | |  '.  '.| |  | |  '.  '.| |  | |  '.  '.| |  | |  '.  '.| |  | |  
.  `.:_ | :_.' '.  `.:_ | :_.' '.  `.:_ | :_.' '.  `.:_ | :_.' '.  `.:_ | :_.
  `-..,..-'       `-..,..-'       `-..,..-'       `-..,..-'       `-..,..-'  
                                                                             
-----------------------------------------------------------------------------
                                                                             
             GCparagon (v0.6.7) maintained by @BGSpiegl
                 Copyright (c) 2023 Benjamin Spiegl
            GitHub: https://github.com/BGSpiegl/GCparagon

^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^
_____________________________________________________________________________

options:
  -h, --help            show this help message and exit
  --version             show program's version number and exit

Input (required):
  -b List[File] [List[File] ...], --bam List[File] [List[File] ...]
                        Path to sorted BAM file for which the fragment length-
                        dependent GC-content-based over-representation (= 'GC-bias')
                        should be computed and/or corrected. WARNING: don't use
                        unaligned BAM files (uBAM) or multi-sample/run BAM files! If
                        the BAM's index file is not found on runtime, GCparagon tries
                        to create it. The alignment algorithm used for creating the
                        input BAM file MUST follow the SAM format specifications! The
                        TLEN column is used by GCparagon. [ PARAMETER REQUIRED ]
  -rgb {hg38,hg19}, --reference-genome-build {hg38,hg19}
                        The version of the reference genome build used for creating
                        the input BAM file. Currently only hg19 and hg38 are
                        available. Used to define -rtb/--two-bit-reference-genome,
                        --c/-intervals-bed and -rgcd/--reference-gc-content-
                        distribution-table. The definition is overridden if any of
                        these parameters are specified on the command line. [ DEFAULT:
                        'hg38' ]
  -rtb File, --two-bit-reference-genome File
                        Path to 2bit version of the reference genome FastA file which
                        was used for read alignment of the input BAM file. If the 2bit
                        version is missing, one can create the file using the
                        following command: 'faToTwoBit <PATH_TO_REF_FASTA> -long
                        <PATH_TO_OUT_2BIT>' (see
                        genome.ucsc.edu/goldenPath/help/twoBit.html for more details)
                        [ DEFAULT: '/mnt/NVMeScratch/PycharmProjects/GCparagon_public/
                        src/GCparagon/2bit_reference/hg38.analysisSet.2bit' ]
  -c File, --intervals-bed File
                        Path to BED file containing predefined genomic intervals to
                        process. These should have been selected based on minimal
                        overlap with exclusion-masked regions of the reference genome
                        build used for read alignment earlier (i.e., creation of
                        --bam). Since v0.5.6, the table also contains expected GC
                        content counts which should be computed using the fragment
                        length distribution from
                        'accessory_files/reference_fragment_lenght_distribution.tsv'.
                        The GC content distributions are used to create an optimized
                        consolidated weight matrix using a weighted mean. The weights
                        for each region are selected such that the reference GC
                        content distribution in
                        'accessory_files/hg38_reference_GC_content_distribution.tsv'
                        is approximated better than by naive averaging region GC
                        content distributions. This optimization is might fail in
                        cases where the majority of preselected genomic regions are
                        used to estimate GC bias. In these cases, the naive approach
                        is used as a fallback option. [ DEFAULT: '/mnt/NVMeScratch/Pyc
                        harmProjects/GCparagon_public/accessory_files/hg38_minimalExcl
                        usionListOverlap_1Mbp_intervals_33pcOverlapLimited.FGCD.bed' ]
  -rgcd File, --reference-gc-content-distribution-table File
                        Path to TSV file containing two data columns with header:
                        'gc_percentage', and 'relative_frequency'. This table defines
                        a GC content distribution (0% GC to 100% GC) as relative
                        frequencies of these percentage bins (summing up to 1). If a
                        custom reference genome is used, this file should be created
                        anew from the simulated genome-wide ideal fragment GC content
                        as simulated assuming a fragment length distribution as the
                        one stored in 'accessory_files/plasmaSeq_ccfDNA_reference_frag
                        ment_length_distribution.tsv'! The provided information is
                        used to optimize the combination of correction weight matrices
                        from different genomic intervals to achieve a linear
                        combination of these regions which resembles the reference GC
                        content distribution defined here. This file can be computed
                        following the instructions of gGenomic intervals preselection
                        script 4: 'accessory_files/genomic_interval_preselection-
                        shifted16x_hg38/04-GI-
                        preselection_simulate_genomewide_reference_FGCD_hg38.py'. [
                        DEFAULT: '/mnt/NVMeScratch/PycharmProjects/GCparagon_public/ac
                        cessory_files/hg38_reference_GC_content_distribution.tsv' ]
  -ec File, --exclude-intervals File
                        Path to library file (BED-like) holding DoC-specific
                        definition of bad intervals (intervals must be exact genomic
                        locus match for exclusion, DO NOT expect bedtools intersect-
                        like behavior!). If the bad intervals library is left default,
                        the bad intervals library with the most recent time stamp in
                        the parent directory of the default library/BED file is used.
                        The bad intervals library is intended to speed up the sample
                        processing by excluding intervals with insufficient DoC form
                        the beginning. Excluded intervals were observed to appear most
                        frequently close to centromeres.
  -cw File, --correction-weights File
                        Optional input for --tag-only mode: a matrix file
                        ('*_gc_weights.txt.gz') containing correction weights to be
                        used in tag-only mode ('--tag-only' flag must be set) to
                        create a new, GC-bias-corrected BAM file with weighted
                        fragments (GC-tag).
  -wm File, --weights-mask File
                        Optional path to a weights matrix mask file. These are usually
                        named '<SAMPLE_ID>_gc_bias_computation_mask.txt.gz'. If none
                        is defined (default behaviour), either the currently create
                        mask (when computing GC bias correction weights) or (for
                        --tag-only) a mask file in the same input directory as the
                        correction weights matrix defined via --correction-weights
                        parameter is used based on the naming of the file. If none is
                        specified or found, correction weights are reduced to rows and
                        columns containing non-default values (values other than 1.
                        and 0.).

Output options:
  -v, --verbose         This flag can be set to provide more output, especially
                        information about fragments that fall outside of the defined
                        fragment length window.
  -o File, --out-dir File
                        Path to which output is moved from --temporary-directory after
                        each processing step (GC-bias computation, BAM tagging). The
                        directory will be created if it does not exist. Make sure that
                        it is empty if it exists, otherwise the whole directory will
                        be deleted before writing to it in the GC-bias computation
                        step! If none is provided, a new subdirectory named
                        'GC_bias_correction_GCparagonv0.6.7' will be created in the
                        input BAM's parent directory and used as output directory. The
                        output for each sample will be gathered in a subdirectory of
                        this --out-dir which will be named after the sample. The
                        output directory may be located on slow hardware such as a USB
                        drive or a network storage since everything is stored in
                        --temporary-directory first and moved after completion of all
                        defined phases of the GC bias computation.
  -tmp File, --temporary-directory File
                        Directory to which all files will be written as temporary
                        files in a subdirectory named after the sample (sample id is
                        extracted from the BAM file name) during processing. Directory
                        will be created if non-existent. Subdirectory for the sample
                        will be deleted if it exists initially. If not specified, this
                        directory is identical to the output of Python's tempfile
                        module's gettempdir() function. Permanent non-temporary output
                        files will be MOVED to the --out-dir using shutil's move
                        function from this directory. The temporary directory should
                        be located on a high performance hardware (high IOPs)! [
                        DEFAULT: /tmp ]
  -np, --no-plots       Flag suppresses creation of fragment length distribution plot
                        and heatmaps for observed, expected, correction, and
                        computation mask matrices.Has no effect in --tag-only mode.
  -os, --output-simulations
                        Optional flag for GC-bias computation for plotting individual
                        simulation results (simulated fragments and iteration-specific
                        masks). The simulated fragment attribute distributions and
                        computation masks are plotted for all simulations then.
  -k, --keep-interval-data
                        Optional flag which can be used to save intermediate data per
                        genomic interval.
  -ob, --output-bam     Optional flag to activate writing of the GC-correction-
                        weights-tagged BAM file AFTER COMPUTING GC BIAS (--tag-only
                        flag is not set), either using the statistics computed from
                        the input BAM file or a correction weights matrix specified
                        via --correction-weights. Is implicit if --tag-only is used.
                        WARNING: currently, the output BAM won't contain unaligned
                        reads!
  -our, --output-unaligned-reads
                        Optional flag to activate writing of unaligned reads to a
                        separate BAM file. Per default, unaligned reads are not
                        output. Setting this flag only has an effect if either the
                        --output-bam flag was set or GCparagon was started in the
                        --tag-only mode.
  -fp Integer, --float-precision Integer
                        Optional parameter for GC-bias computation number of digits
                        after the comma for floating point data to be stored in text-
                        based matrices, e.g. for correction weights data. Choose
                        according to expected depth of coverage -> if you would expect
                        10,000, you can go for 5 or even 6 digits. Otherwise this will
                        not have an effect. If you compute signals that are sums over
                        many regions, multiply the expected DoC with how many signals
                        you sum up to get an estimate of which precision you would
                        need to definitively be able to rule out any influence by
                        rounding errors. These should average out though. [ DEFAULT: 6
                        ]
  -tg String, --tag-name String
                        Name of the GC-bias correction weight tag that will be added
                        to alignments in the BAM file. If none is provided, the
                        default tag will be used. Must not be longer than 2
                        characters! [ DEFAULT: GC ]
  -wie, --write-interval-exclusion
                        Optional flag for writing an updated version of the library
                        listing intervals marked for exclusion from the analysis. Per
                        default, genomic intervals are marked for exclusion if drawing
                        fragments of a specific size repeatedly fails (at least 55
                        times (for strict reference N base handling, 33 times
                        otherwise) or 1/3 of number of fragments that need to be
                        drawn, whichever is higher) due to getting only poly-N
                        sequences. In general, the frequency of these exclusion events
                        is dependent on the DoC of the sample, which can be
                        substituted by the number of fragments estimated to be
                        obtained from all predefined intervals in BAM file in a first
                        approximation. WARNING: don't mix exclusion-marked interval
                        libraries computed from different (predefined) interval BED
                        files! If the user places the output BED file library in the
                        default directory, the new library will be used per default
                        for future computations. Genomic intervals will be marked for
                        exclusion depending on a data set's fragment length
                        distribution and sequencing depth.
  -nfp, --no-focused-plots
                        Optional flag to deactivate focusing of matrix plots on non-
                        default values (focus uses a border of up to 10 default
                        values). Only has an effect if --no-plots flag is not set.
  -sf, --show-figures   Optional flag to display plots in an interactive browser
                        window in addition to saving them to a file.

Processing options:
  -p Integer, --preset Integer
                        Optional parameter preset to use for GC bias computation. Must
                        be an integer int the rangeof 0-3 (inclusive). A preset value
                        of 0 leaves parameters at default if not defined differently
                        by the user (unchanged parameters will match preset 1). Other
                        integer values from 1 to 3 define presets with increasing
                        input data usage and required processing time (durations
                        preset 1-3: 1-3 min, 5-10 min, and ~1h depending on file size.
                        Maximum across 4 samples and 2 iterations each computed using
                        12 cores and the profile_command.py script. Maximum memory
                        consumption for any preset should stay below 4 GiB. If preset
                        is not zero, any customized parameters conflicting with the
                        preset will be ignored. A non-zero preset will set the
                        following parameters: number of simulations, the target number
                        of processed fragments, minimum number of fragment attribute
                        combination occurrences, and the options for outlier detection
                        and smoothing. Noise within the resulting correction weights
                        is reduced when selecting a higher preset value. Preset 3 will
                        attempt to process all genomic intervals (target number of
                        fragments set to 100B) within the limits of the maximum
                        allowedexclusion marked regions overlap (per default default
                        ~2.7 Gb of reference are processed). NOTE: the percentage of
                        total GC bias corrected fragments in the dataset for presets 1
                        vs. 3 increases only from 99.837% to 99.938% (average across 4
                        samples). Other fragment weights default to 1.0). The primary
                        advantage of processing more fragments is the reduction of
                        noise in computed weights and a better reconstruction of the
                        reference fragment GC content distribution. It is recommended
                        to use a higher preset for a 'preprocess-once,analyze often'
                        scenario and/or when a high bias is expected/observed (e.g.
                        FastQC average GC percentage). Correction by preset 1, 2, and
                        3 was found to yield 100.74%, 99.98%, and 99,91% of the raw
                        fragment count respectively (average percentage across 4
                        samples). [ DEFAULT: 2 ]
  -to, --tag-only       Optional flag which makes the software switch to tag-only
                        mode. A correction weights matrix must be specified in this
                        case via the '--correction-weights' flag. A valid samtools
                        path must be available via the system path variable or
                        provided using --samtools-path. Be mindful of setting the
                        temporary directory correctly for your system! (e.g. should be
                        set to output of 'echo $TEMP' on HPC clusters)
  -rep Integer, --repetition-of-simulation Integer
                        (PRESET precedence if specified) This value can be left at
                        default if the target number of processed fragments is
                        sufficiently high (e.g. >=5M). The lower the number of target
                        fragments, the stronger is the effect of increasing the number
                        of simulation rounds. Increasing this value increases the
                        computation time almost accordingly (scales linearly). [
                        DEFAULT: 6 ]
  -mafl Integer, --maximum-fragment-length Integer
                        Defines upper length limit for fragments which should be
                        included in computation. This parameter does not impact
                        computation speed. It only increases plotting times for
                        matrices by a few seconds and memory consumption. [ DEFAULT:
                        550bp ]
  -mifl Integer, --minimum-fragment-length Integer
                        Defines lower length limit for fragments which should be
                        included in computation. Must be positive integer. A value
                        below the sequenceable fragment length of the device used to
                        create the dataset is not recommended. [ DEFAULT: 20bp ]
  -t Integer, --threads Integer
                        Total number of threads to be used for BAM processing. If the
                        --single-thread-processes flag was set, this number
                        corresponds to the number of processes spawned for BAM
                        processing. For BAM tagging, multiple threads are used for the
                        sort/merge operations so fewer processes might be used
                        simultaneously. Should be lower than the total number of
                        logical cores available on the hardware. Will be reduced to
                        max. available number of logical cores if is set higher by the
                        user. [ DEFAULT: 12 ]
  -rs Integer, --random-seed Integer
                        Optional random seed to be used for genomic sampling patterns.
                        Warning: the notion that all computed numbers will turn out
                        identical when using the same random seed using different
                        interpreters or different machines should be discarded right
                        away. Might only be useful when repeatedly running the script
                        within the same python interpreter instance! [ DEFAULT: 99
                        (randomly drawn from 0-999) ]
  -sp File, --samtools-path File
                        Optional input: path to specific samtools executable. A valid
                        path is required for creating the tagged BAM output. By
                        default, this path will be used: '/home/benjamin/mambaforge-
                        pypy3/envs/GCparagon/bin/samtools' (empty or None if path is
                        not found). Code tested with samtools version 1.16.1 using
                        htslib 1.16 [ PARAMETER REQUIRED IF DEFAULT VALUE IS
                        EMPTY/NONE ]
  -nf Integer, --target-fragment-number Integer
                        (PRESET precedence if specified) GC-bias computation will stop
                        after surpassing this threshold for processed fragments. Still
                        running subprocesses will be finished and results included so
                        usually this value is overshot by up to several million
                        fragments depending on the amount of processes chosen and the
                        DoC inside defined bins. Increasing this value will reduce the
                        noise in computed weights. Concerning the number of corrected
                        fragments, processing more than 5 million fragments will only
                        increase the number of corresponding computed weights only
                        miniscule. Doubling the target processed fragment amount
                        typically leads only to an increase in corrected fragments by
                        less than one percent. Five million fragments (preset 1)
                        should be enough to correct between 99.5% and 99.9% of all DNA
                        fragments observed in the dataset based on GC content and
                        fragment length. To reach >99.9% of corrected fragments, this
                        parameter should be increased. [ DEFAULT: 5,000,000 ]
  -mf Integer, --minimum-fragment-occurrences Integer
                        (PRESET precedence if specified) This parameter defines the
                        minimum number of fragment occurrences for a specific
                        length/GC-content attribute combination to be regarded in
                        correction weights computation (= mask definition). Higher
                        values result in less extreme weight outliers, especially for
                        the low-GC-content mononucleosomal fragment length range. The
                        absolute lowest supported value is 2 based on visual
                        inspection of resulting weights matrices for different
                        samples. If this value is too low (e.g. 1), strong 'salt-and-
                        pepper'-type noise was observed for rare attribute
                        combinations along with very high weight outliers. A value of
                        10 here means that a particular attribute combination must
                        occur at least once per million fragments in the dataset for a
                        '--number-of-fragments-to-process' value of 10,000,000. As a
                        rule of thumb, one can set this to number of million target
                        fragments (i.e. set to 10 for the target value of 10M
                        processed fragments as in the example above) [ DEFAULT: 3 ]
  -anf, --allow-n-base-fragments
                        Per default, any fragment containing N-bases (as determined
                        from the read alignment positions and the reference genome
                        sequence) is excluded from the analysis. This parameter was
                        not found to cause any problems for Illumina NovaSeq data. If
                        such fragments have to be included, this flag can be set to
                        allow for up to 1/3 N-bases for fragments. Parameter mainly
                        influences the simulation step and how many times random
                        fragment drawing must be repeated for individual genomic
                        intervals. Also can lead to fewer intervals being discarded
                        (and marked as bad genomic interval) if flag is set.
  -ucmaf Float, --unclipped-min-aln-fraction Float
                        This parameter defines the minimum unclipped fraction of an
                        alignment to be counted in the observed fragment attributes
                        matrix O_gc. This might affect how many small fragments are
                        observed and effectively corrected. [ DEFAULT: 0.75 ]
  -dw DEFAULT_FRAGMENT_WEIGHT, --default-weight DEFAULT_FRAGMENT_WEIGHT
                        Parameter redefines the weight which is assigned to fragments
                        with fragment length + GC base count combinations that lie
                        outside of the non-default range of the (computed) weights
                        matrix. Should be 1.0 for GCparagon. Can be e.g. 0.0 for other
                        algorithms like Griffin. Choose according to the source of
                        your weights matrix! [ DEFAULT: 1.0 ]
  -reto UNALIGNED_EXTRACTION_TIMEOUT, --reads-extraction-timeout UNALIGNED_EXTRACTION_TIMEOUT
                        Sets the timout in seconds for unaligned reads extraction.
                        Only has an effect if '--output-bam' or '--tag-only' and '--
                        output-unaligned-reads' is set. [ DEFAULT: 1800 seconds ]

Post-processing options:
  -do, --detect-outliers
                        (PRESET precedence if specified) If this flag is set, extreme
                        outliers will be detected and limited to a threshold value
                        that is computed from the fragment weights. The default method
                        to detect outliers is Q3 + 8x inter-quartile range (IQR).
                        Values above this threshold will be limited to the threshold.
                        It is highly recommended to detect and limit outliers.
  -odm String, --outlier-detection-method String
                        (PRESET precedence if specified) If the --detect-outliers flag
                        is set, the detection method can be set here. Either a method
                        based on the inter-quartile range or a method based on
                        standard deviation can be selected. Must be one of {'IQR',
                        'SD'}. [ DEFAULT: IQR ]
  -ods Integer, --outlier-detection-stringency Integer
                        (PRESET precedence if specified) If the --detect-outliers flag
                        is set, this parameter defines how stringent the outlier
                        detection threshold is set. Must be an integer in the range of
                        1-7 (inclusive). [ DEFAULT: 2 ]
  -sw, --smooth         (PRESET precedence if specified) If this flag is set, computed
                        weights will also be smoothed. An additional matrix is output
                        containing these post-processed values. If plotting is set to
                        true, also a visualisation of the smoothed weights will be
                        created.It is recommended to smooth weights if not the entire
                        dataset is processed (like is done in preset 3).
  -sk String, --smooth-kernel String
                        (PRESET precedence if specified) If the '--smooth' flag is
                        set, the type of kernel used in the 2D convolution operation
                        can be set here. In general, a Gaussian kernel makes more
                        sense because it assigns directly adjacent values a higher
                        weight in computing the smoothed value of the current
                        position. Must be one of {'gauss', 'constant'}. [ DEFAULT:
                        gauss ]
  -si Integer, --smoothing-intensity Integer
                        (PRESET precedence if specified) If the '--smooth' flag is
                        set, the smoothing intensity defines the range of the 2D
                        kernel used in the smoothing operation. Must be an integer in
                        the range of 1-10 (inclusive). [ DEFAULT: 5 ]
```


### Parameter Presets
Parameter presets are defined using `-p`/`--preset`.
The following table shows pre-defined parameters for each preset along with the average computation time across the 4 
samples from [EGAS00001006963].
(Preset 2 is default)

|     Preset | target fragment number | simulation rounds | minimum attribute pair count | outlier detection | weights smoothing | smoothing strength | est. computation time |
|-----------:|-----------------------:|------------------:|-----------------------------:|:-----------------:|:-----------------:|:------------------:|----------------------:|
|          1 |              5,000,000 |                 6 |                            2 |        on         |        on         |         5          |               1-3 min |
| DEFAULT: 2 |             50,000,000 |                 4 |                           10 |        on         |        on         |         2          |              5-10 min |
|          3 |         99,999,999,999 |                 4 |                           20 |        on         |        on         |         2          |              ~50 min* |

*depends on DoC of BAM file


## Result of GC Bias Correction

-------------------------------------------------------------------------------------------------------------------
GC bias correction results using default parameter preset 2 (5–10 min) of two human cfDNA samples (paired-end WGS, 
[EGAS00001006963]) are shown below.
For each sample, the original fragment GC content distribution (= "FGCD", dotted colored lines), the GCparagon-corrected
FGCD (solid colored lines), the correction achieved with the Griffin algorithm (grey dashed lines), and the simulated 
FGCD across the entire analyzable genome (GRCh38, black lines) are displayed.
The GC content of fragments was estimated either:
 - by using the read sequence, if the template length is (shorter or) equal to the read sequence length,
 - or from slices of the reference genome using the leftmost alignment position and the template length otherwise.

![p01_preset2_correction](https://github.com/BGSpiegl/GCparagon/blob/including_EGAS00001006963_results/validation/03_genome-wide_correction_fidelity/FGCD_correction_plots/preset2/B01_GCparagon_GC-content-comparison_GC-bias-correction_SPLINE_cfDNAref.png?raw=true)

(B01, GCparagon preset 2 vs. Griffin correction; fragment GC content in 2 %GC bins, spline interpolated)


![p01_preset2_correction](https://github.com/BGSpiegl/GCparagon/blob/including_EGAS00001006963_results/validation/03_genome-wide_correction_fidelity/FGCD_correction_plots/preset2/P01_GCparagon_GC-content-comparison_GC-bias-correction_SPLINE_cfDNAref.png?raw=true)

(P01, GCparagon preset 2 vs. Griffin correction; fragment GC content in 2 %GC bins, spline interpolated)

The maximum residual bias seems to depend on the initial intensity of GC bias and on the representativeness of
preselected genomic intervals that were used for computing GC bias.
In general, more intervals are used for shallow sequenced samples than for deeper sequenced ones to 
achieve a comparable number of processed fragments.
With increasing preset number, the number of genomic intervals used
for GC bias computation also increases. This also allows for a better representation of the reference FGCD by linear 
combination of weight matrices computed from individual preseleted genomic intervals.

When applied to multiple transcription start sites (TSSs) of genes which are generally expected to be inactive in adult 
humans ([975 genes](accessory_files/TSSs/PAU.txt) as derived from the [protein atlas](https://www.proteinatlas.org/)),
and active genes ([1179 "housekeeping" genes](accessory_files/TSSs/HK.txt)), GC bias manifests as changes in the average
central 61 bp fragment depth of coverage (cDoC) across these 5' -> 3' oriented sites. Active genes are expected 
to show a nucleosome depleted region (unprotected -> decrease in coverage) slightly upstream to the TSS, whereas 
unexpressed or lowly expressed genes should show an almost flat cDoC profile.

Examples of positive (P01, +5.0%) and negative GC bias (B01, -2.2%) on the average DoC for expressed and 
unexpressed genes is shown below (fragment coverage in silico reduced to their central 61 bp).
Original coverage before correction is shown as thin-colored lines.
cDoC of mono-nucleosomal fragments after GC bias correction with GCparagon is shown as thick-colored lines.

For these plots, only mono-nucleosomal fragments with observed template length between 110-210bp were included to create 
a surrogate signal to assess the average nucleosome positioning across combined regions (active genes: expected 
increased positioning of nculeosomes downstream to the TSS; inactive genes: negligible positioning).
The original H01 sample has the lowest deviation of average GC content from the expected 40.4% and shows the weakest
GC bias. Hence, the original and corrected DoC profiles are very similar. We also corrected the same samples using the 
Griffin algorithm and processed fragments accordingly. Correction resulting from Griffin is shown as back lines.

![doc_corr_res_tsss](https://github.com/BGSpiegl/GCparagon/blob/including_EGAS00001006963_results/validation/02_loci_overlay_central_coverages/coverage_overlay_plot_output/DoC_bias_correction_effect_TSSs_preset2.png?raw=true)

The DoC increase/decrease before correction for samples showing a positive/negative GC bias (P01/B01) results from the local sequence GC content.
Similarly, the artificial DoC increase in the following figure downstream to position 0 (= TSS) results from increased GC content of 5'-UTR sequences 
compared to the immediate upstream core promoter.
Promoter sequences can contain the [TATA-box] element, usually located 25bp upstream to position zero (approx. every 3rd promoter).

![doc_corr_res_tfbss](https://github.com/BGSpiegl/GCparagon/blob/including_EGAS00001006963_results/validation/02_loci_overlay_central_coverages/coverage_overlay_plot_output/DoC_bias_correction_effect_TFBSs_preset2.png?raw=true)

Similarly, many transcription factor binding motifs show an increased GC content. Original cDoC is increased/decreased 
for samples showing a positive/negative GC bias (P01/B01) with the most extreme distortion observed for the LYL1 locus 
of P01 which shows the most intense GC bias (+5.0% GC).


### Output

Default outputs are:

- log files (stdout and stderr logged individually)
- fragment length distribution (plot only; can be computed from `*_observed_attributes_matrix.txt.gz`)

(P01, default preset 2 examples shown)

![p01_frag_length_dist](https://github.com/BGSpiegl/GCparagon/blob/including_EGAS00001006963_results/preset_computation/target_output-preset2/P01/P01.fragment_length_distribution.png?raw=true)

- observed fragment attributes (plot, data in `*_observed_attributes_matrix.txt.gz`)

![p01_observed_atts](https://github.com/BGSpiegl/GCparagon/blob/including_EGAS00001006963_results/preset_computation/target_output-preset2/P01/P01.O_gc.heatmap.png?raw=true)

- simulated fragment attributes using reference genome and fragment length distribution (plot, data in 
`*_simulated_attributes_matrix.txt.gz`)

![p01_simmed_atts](https://github.com/BGSpiegl/GCparagon/blob/including_EGAS00001006963_results/preset_computation/target_output-preset2/P01/P01.S_gc.heatmap.png?raw=true)

- weights computation mask (plot; data in `*_gc_bias_computation_mask.txt.gz`)

![p01_comp_mask](https://github.com/BGSpiegl/GCparagon/blob/including_EGAS00001006963_results/preset_computation/target_output-preset2/P01/P01.Mask.heatmap.png?raw=true)

- correction weights matrix (plot; data in `*_gc_weights_*simsMean.txt.gz`)

![p01_w_gc](https://github.com/BGSpiegl/GCparagon/blob/including_EGAS00001006963_results/preset_computation/target_output-preset2/P01/P01.W_gc.heatmap.png?raw=true)

- correction weights matrix, extreme outliers capped at threshold (plot; data in 
`*_gc_weights_*simsMean.*outliersRemoved.txt.gz`)

![p01_w_gc_ol](https://github.com/BGSpiegl/GCparagon/blob/including_EGAS00001006963_results/preset_computation/target_output-preset2/P01/P01.W_gc_outliers_removed.heatmap.png?raw=true)

- correction weights matrix, extreme outliers capped at threshold, local smoothing applied (plot; data in 
`*_gc_weights_*simsMean.*outliersRemoved.*gaussSmoothed.txt.gz`)

![p01_w_gc_ol_sm](https://github.com/BGSpiegl/GCparagon/blob/including_EGAS00001006963_results/preset_computation/target_output-preset2/P01/P01.W_gc_outliers_removed_smoothed.heatmap.png?raw=true)


### Genomic Region Preselection

The code uses up to 2,342 
[preselected 1 Mb genomic intervals](accessory_files/hg38_minimalExclusionListOverlap_1Mbp_intervals_33pcOverlapLimited.FGCD.bed) of 
hg38 reference genome for processing (2,204 regions for hg19).
Preselection was carried on the basis of an 
[exclusion listed regions BED file](accessory_files/hg38_GCcorrection_ExclusionList.merged.sorted.bed) to:
- reduce the overlap of 1 Mb genomic intervals with regions found to be problematic for short read sequencing or sequence alignment
- minimize inclusion of assembly gaps and errors
- slightly optimize the placement of intervals along the genome while ensuring preselection of only non-overlapping regions
- enforce a threshold for maximum percentage of exclusion listed bases per genomic interval (33% max. overlap 
implemented)

Creation of the exclusion listed regions BED file, starting from the 
[ENCODE exclusion listed regions v2 BED file](accessory_files/bad_regions/hg38-ENCODE_blacklist.v2.bed),
is documented [here](accessory_files/GC_correction_exclusion_list_creation.info).
The hg19 version is documented [here](accessory_files/GC_correction_exclusion_list_hg19_creation.info).
The final GRCh38 exclusion list includes 385,022,382 reference positions (=12.68% of GRCh38).

Currently, only a preselection of GRCh38 and hg19 genomic 1 Mb regions is available. To specify the genome 
build for your analysis, set the `-rgb`/`--reference-genome-build` parameter to either `hg38` or `hg19`. Per 
default, `hg38` is used.
Preselection for other genome builds currently has to be carried out by the user using their own exclusion list.
The following files need to be created:
 - *<GENOME_BUILD>.2bit*
 - *<GENOME_BUILD>_minimalExclusionListOverlap_1Mbp_intervals_33pcOverlapLimited.FGCD.bed*
 - *<GENOME_BUILD>.genome_file.tsv*
 - *<GENOME_BUILD>_reference_GC_content_distribution.tsv*
 - optional: *<ISOLATION_PROTOCOL>_<SAMPLE_TYPE>_reference_fragment_length_distribution.tsv*
(if the plasmaSeq ccfDNA reference fragment length distribution is inadequate)

Once created,
some of these files need to be passed to GCparagon via the following commandline parameters to take effect:
 - `--reference-gc-content-distribution-table <path_to/GENOME_BUILD>_reference_GC_content_distribution.tsv`
 - `--intervals-bed <path_to/GENOME_BUILD>_minimalExclusionListOverlap_1Mbp_intervals_33pcOverlapLimited.FGCD.bed`
 - `--two-bit-reference-genome <path_to/GENOME_BUILD>.2bit`

These files can be created using the following code/scripts in ascending order:
 - [reference genome build download instructions](src/GCparagon/2bit_reference/EXECUTE_hg38_reference_download.sh)
 - [Exclusion list creation info](accessory_files/GC_correction_exclusion_list_hg19_creation.info)
 - [01-GI-preselection_test_Mbp_genomic intervals_against_ExclusionList_hg19.py](accessory_files/genomic_interval_preselection-shifted16x_hg19/01-GI-preselection_test_Mbp_intervals_against_ExclusionList_hg19.py) 
(creation of equal-sized genomic intervals, intersected with the exclusion marked regions list)
 - [02-GI-preselection_select_low_scoring_regions_from_overlapping_hg19.py](accessory_files/genomic_interval_preselection-shifted16x_hg19/02-GI-preselection_select_low_scoring_regions_from_overlapping_hg19.py) 
(selection of least-exclusion-list-overlapping regions)
 - [03-GI-preselection_compute_genomic_interval_fragment_GC_content_hg19.py](accessory_files/genomic_interval_preselection-shifted16x_hg19/03-GI-preselection_compute_genomic_interval_fragment_GC_content_hg19.py) 
(computing fragment GC content distributions (FGCDs) for each preselected genomic interval (GI);
also detects bad simulation intervals and writes them to a separate BED file 
'<GENOME_BUILD>_genomic_intervals_failed_fragment_drawing.bed')
 - [04-GI-preselection_simulate_genomewide_reference_FGCD_hg19.py](accessory_files/genomic_interval_preselection-shifted16x_hg19/04-GI-preselection_simulate_genomewide_reference_FGCD_hg19.py) 
(computing the genome-wide reference fragment GC content distribution = "reference FGCD" based on the reference 
fragment length distribution).

 - Optionally,
the diversity of fragment GC content distributions among preselected genomic intervals can be visualized using 
[05-GI-preselection_visualize_preselected_intervals_GC_content_distributions_hg19.py](accessory_files/genomic_interval_preselection-shifted16x_hg19/05-GI-preselection_visualize_preselected_intervals_GC_content_distributions_hg19.py)


![hg19_GI_preselection](https://github.com/BGSpiegl/GCparagon/blob/dev-hg19ExclusionList/accessory_files/genomic_interval_preselection-shifted16x_hg19/hg19_FGCD_intervals_vs_reference/GC_content_reference_genome_vs_preselected_intervals.51bins.png?raw=true)

It is recommended to use at least the ENCODE exclusion list to restrict genomic interval preselection.
The size of preselected genomic intervals should be uniform and fit the application (i.e., larger genomic intervals for 
shallow sequenced WGS data). A tradeoff was made for GRCh38 by choosing 1 Mb genomic intervals.

Not including assembly gaps might cause a significant performance hit because GCparagon would then try and repeatedly 
fail to draw fragments from gapped regions containing only Ns unnecessarily.
This can be prevented by including, for example, the UCSC gap track regions in the region exclusion.
The gap track must match the genome build of interest!
Gap tracks can be downloaded using the UCSC table browser [UCSC Table Browser][table browser].

An IGV screenshot visualizing the distribution of 
[hg38 preselected 1Mb genomic intervals](accessory_files/hg38_minimalExclusionListOverlap_1Mbp_intervals_33pcOverlapLimited.FGCD.bed)
across the whole genome, and additionally for two chromosomes, is provided 
[here](accessory_files/genomic_interval_preselection-shifted16x_hg19/IGV_composite_preselected_intervals_hg38.png).


### Optimization of Combining Results from Preselected Regions
From v0.6.0 on, fragment attribute counts from preselected genomic intervals are combined using a linear combination 
such that their simulated cfDNA fragment GC content distribution approximates the FGCD of the reference genome without 
dropping information gathered from individual regions. The range of weights from which the linear combination of 
genomic interval FGCDs is created is: [0.1, 10].

A plot showing the relative frequency of selected weights within this range is created by default:

![p01_int_w](https://github.com/BGSpiegl/GCparagon/blob/including_EGAS00001006963_results/validation/00_result_plots_preset2_for_readme/P01.preselectedRegionWeights_refGCreconstruction.png?raw=true)

The weights should not accumulate at the lower end of the scale because this would mean an over-fitting to a small 
subset of genomic intervals that themselves resemble the reference GC content best in a specific linear combination. 

Another plot shows the improvement of the reference GC content approximation by the computed linear combination over 
a naive combination of GC correction weights:

![p01_recons_GC](https://github.com/BGSpiegl/GCparagon/blob/including_EGAS00001006963_results/validation/00_result_plots_preset2_for_readme/P01.GCcontent_RefVsReconstructed.png?raw=true)


### Validation v0.6.0+
Access the validation code via the [including_EGAS00001006963_results][validation branch]!

To recreate the presented GC correction results, run [this driver script](driver_scripts/drv_compute_GC_presets.sh) 
after setting up the conda env and downloading the 2bit reference genome file and the [EGAS00001006963] BAM files.
The BAMs used in plots can be requested for download from EGA via the accession [EGAS00001006963].
If required, a new EGA account can be created for free.

To compute the validation results for the depth of coverage signals, one can either use the npz exported, normalized and
averaged DoC values for each loci group, preset and sample (mean of 10th to 90th percentile of normalized cDoC values, 
per base position) in the subdirectories of the 
[coverage signals validation directory](validation/02_loci_overlay_central_coverages) and run the two plotting code 
scripts for [TSS loci groups](validation/02-2_create_TSS_subplot_figure_from_extracted_cDoC.py) and 
[TFBS loci groups](validation/02-3_create_TFBS_subplot_figure_from_extracted_cDoC.py)
OR create these coverage *.npz files themselves using their own custom code 
(extract from 110bp-210 bp fragments the central 61 bp fragment DoC signal for each locus and compute central 80% DoC mean 
for each position in the window).
To process Griffin output (bias matrices), one must run the scripts in the [validation directory](validation) in numerical 
order (i.e. from [01-1_transform_Griffin_GC_bias_weights.py](validation/01-1_transform_Griffin_GC_bias_weights.py) to 
[03-3_plot_per_fragment_GC_content_correction_fidelity.py](validation/03-3_plot_per_fragment_GC_content_correction_fidelity.py)).

The full validation requires the GCparagon preset benchmark results to be available.
To create the genome-wide correction fidelity plots, the transformed Griffin bias matrices and Griffin weights-tagged 
BAM files must be available (compute matrices with 
[01_transform_Griffin_GC_bias_weights.py](validation/01-1_transform_Griffin_GC_bias_weights.py)).
All results created by these scrips are provided along with the code, except the Griffin correction weights-tagged BAM files.
These can be created with the [02-1_tag_BAMs_with_Griffin_result.py](validation/02-1_tag_BAMs_with_Griffin_result.py) script.
Some of the validation template scripts must be completed with the appropriate directory and/or code paths before being 
called (see "TODO" comments at the beginning of scripts).


## Copyright

-------------------------------------------------------------------------------------------------------------------
- Original work on GCparagon.py and accessory code Copyright (c) 2023 Benjamin Spiegl
- Original work on profile_command.py Copyright (c) 2023 Marharyta Papakina and Benjamin Spiegl

## Software license

-------------------------------------------------------------------------------------------------------------------
[MIT License](src/GCparagon/LICENSE)

Intended for research use only.

## Contributors

-------------------------------------------------------------------------------------------------------------------
- Benjamin Spiegl ([BGSpiegl][github user])
- Marharyta Papakina
- Sebastian Röner
- Isaac Lazzeri

-------------------------------------------------------------------------------------------------------------------

GCparagon developed at the [D&F Research Center of Molecular Biomedicine][molbiomed graz], 
[Institute of Human Genetics][humgen graz], [Medical University of Graz, Austria][mug]

GCparagon uses the [ENCODE Blacklist (The Boyle Lab, GitHub)][encode_exclusion_listed_regions_url]:
Amemiya, H.M., Kundaje, A. & Boyle, A.P. The ENCODE Blacklist: Identification of Problematic Regions of the Genome. Sci 
Rep 9, 9354 (2019). https://doi.org/10.1038/s41598-019-45839-z

GCparagon uses resources from the [UCSC Genome browser][genome browser]

[github user]: https://github.com/BGSpiegl
[github repo]: https://github.com/BGSpiegl/GCparagon
[TATA-box]: https://en.wikipedia.org/wiki/TATA_box
[hardware_requirements]: https://github.com/BGSpiegl/GCparagon#hardware-requirements
[software_dependencies]: https://github.com/BGSpiegl/GCparagon#software-dependencies
[usage]: https://github.com/BGSpiegl/GCparagon#usage
[presets]: https://github.com/BGSpiegl/GCparagon#presets
[repository_structure]: https://github.com/BGSpiegl/GCparagon#repository-structure
[molbiomed graz]: https://www.medunigraz.at/en/research-centers-and-institutes/diagnostic-and-research-center-for-molecular-biomedicine
[humgen graz]: https://humangenetik.medunigraz.at/en/
[mug]: https://www.medunigraz.at/en/
[mamba install]: https://mamba.readthedocs.io/en/latest/installation.html
[conda install]: https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html
[hg38_std_analysis_set]: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/
[hg38_decoy_analysis_set]: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/
[bwa_mem]: https://github.com/lh3/bwa
[samtools_spec]: https://samtools.github.io/hts-specs/SAMv1.pdf
[EGAS00001006963]: https://ega-archive.org/studies/EGAS00001006963
[genome browser]: https://genome.ucsc.edu/
[encode_exclusion_listed_regions_url]: https://github.com/Boyle-Lab/Blacklist/
[table browser]: https://genome.ucsc.edu/cgi-bin/hgTables
[libmamba solver]: https://conda.github.io/conda-libmamba-solver/user-guide/
[validation branch]: https://github.com/BGSpiegl/GCparagon/tree/including_EGAS00001006963_results
[singularity_installation]: https://docs.sylabs.io/guides/3.0/user-guide/installation.html
[apptainer_files_doc]: https://apptainer.org/docs/user/latest/quick_start.html#working-with-files
[singularity_binds_and_mounts_doc]: https://docs.sylabs.io/guides/3.9/user-guide/bind_paths_and_mounts.html
[singularity_remotes]: https://docs.sylabs.io/guides/3.5/user-guide/endpoint.html