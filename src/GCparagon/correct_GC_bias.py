#!/usr/bin/env python3

# full imports
import os
import gc
import sys
import time
import math
import shutil
import random
import logging
import tempfile
import datetime
import traceback
import contextlib
import multiprocessing

# import aliases
import numpy as np
import pandas as pd
import subprocess as sp
import multiprocessing as mp
import multiprocessing.connection as mp_connection

# partial imports
from pathlib import Path
from collections import deque
from natsort import humansorted
from pysam import AlignmentFile  # coordinates in pysam are always 0-based (following python convention)
from scipy.optimize import nnls, minimize
from scipy.stats import trim_mean
from twobitreader import TwoBitFile, TwoBitSequence
from typing import Union, Dict, List, Tuple, Optional, Any
from argparse import ArgumentParser, RawDescriptionHelpFormatter
OneOf = Union

# TODO: add hg19 compatibility and cmdline flag specifying hg19 instead of hg38;
#                 increment MINOR_RELEASE number and reset patch number -> v0.6.0;
#                 create new release and upload to zenodo

# version
MAJOR_RELEASE = 0
MINOR_RELEASE = 6
PATCH_NUMBER = 0
VERSION_STRING = f'v{MAJOR_RELEASE}.{MINOR_RELEASE}.{PATCH_NUMBER}'

# GitHub link
github_url = 'https://github.com/BGSpiegl/GCparagon'

# dependencies (see "GCparagon_py3.10_env.yml" file for conda env creation):
#   - samtools=1.16
#   - bedtools=2.30
#   - python=3.10
#   - pip=22.3
#   - numpy=1.23
#   - pysam=0.19
#   - natsort=8.2
#   - py2bit=0.3
#   - cycler=0.11
#   - pandas=1.5
#   - scipy=1.9
#   - ucsc-fatotwobit=377
#   - twobitreader=3.1
#   - plotly_express=0.4
#   - python-kaleido=0.2
#   - psutil=5.9
#   - requests=2.28
#   - matplotlib=3.6
#   - memory_profiler
#   - pybedtools
#   - polars

# TIMEOUT DEFINITIONS
READS_EXTRACTION_TIMEOUT = 1800  # seconds; wait a maximum of 30 minutes for single pass through entire file,
# extracting mates where at least one dir not align
# default definitions for analysis
MAX_FRAGMENT_LENGTH = 1200  # maximum allowed fragment length; THIS MUST BE CHANGED TO ENABLE USING LONGER FRAGMENTS!
DEFAULT_MIN_FRAGMENT_LENGTH = 20  # do not set to 0!
DEFAULT_MAX_FRAGMENT_LENGTH = 550
DEFAULT_TAG_NAME = 'GC'
DEFAULT_WEIGHT = 1.0
DEFAULT_MIN_OCCURRENCES = 3
ABSOLUTE_MIN_OCCURRENCES = 2
TAGGING_INTERVAL_SIZE = 50 * 10 ** 6
# estimated HW capacities
max_logical_cores = multiprocessing.cpu_count()
max_physical_cores = max_logical_cores // 2 if max_logical_cores > 1 else 1
# PARAMETERS DEFINING RUNTIME OF GC-BIAS COMPUTATION:
# ----------------------------------------------------------------------------------------------------------------------
DEFAULT_NUMBER_PROCESSES = min(12, max_logical_cores)  # limit default value to meaningful amount
DEFAULT_TARGET_NUMBER_FRAGMENTS_PROCESSED = 5 * 10 ** 6  # 5 million
DEFAULT_PRESET = 2  # 1 is fastest; 2 gives best original count reconstruction
DEFAULT_SIMULATION_REPETITIONS = 6
# ----------------------------------------------------------------------------------------------------------------------
DEFAULT_FLOAT_PRECISION = 6
DEFAULT_FRAGMENT_N_CONTENT_THRESHOLD = 0.3
DEFAULT_MAX_INTERVAL_PERCENTAGE_EXCLUSIONLIST_OVERLAP = 1 / 3 * 100.  # of exclusion-listed regions for 1 Mbp intervals
DEFAULT_MIN_UNCLIPPED_ALN_FRACTION = 0.75  # alignment-clipping limitation: regard alignment faulty
# if more than 25% are clipped and ignore mates of fragment
# POSTPROCESSING DEFAULTS:
DEFAULT_SMOOTHING_INTENSITY = 5
DEFAULT_SMOOTHING_KERNEL = 'gauss'
DEFAULT_OUTLIER_DETECTION_STRINGENCY = 2
DEFAULT_OUTLIER_DETECTION_METHOD = 'IQR'
# random numbers:
RANDOM_SEED = random.randint(1, 999)
# PATH DEFINITIONS:
SOURCE_CODE_ROOT_PATH = Path(__file__).parent
sys.path.append(str(SOURCE_CODE_ROOT_PATH))  # to enable relative imports
SOURCE_CODE_ROOT_DIR = str(SOURCE_CODE_ROOT_PATH)
DEFAULT_SAMTOOLS_PATH = shutil.which('samtools')
DEFAULT_TEMPORARY_DIRECTORY = tempfile.gettempdir()
EXPECTED_TWO_BIT_REFERENCE_GENOME_PATH = SOURCE_CODE_ROOT_PATH / '2bit_reference/hg38.analysisSet.2bit'
PREDEFINED_1MBP_INTERVALS_TO_PROCESS = SOURCE_CODE_ROOT_PATH.parent.parent / \
    'accessory_files/hg38_minimalExclusionListOverlap_1Mbp_intervals_33pcOverlapLimited.FGCD.bed'
DEFAULT_REFERENCE_GENOME_TARGET_GC_CONTENT_DISTRIBUTION = SOURCE_CODE_ROOT_PATH.parent.parent / \
    'accessory_files/hg38_reference_GC_content_distribution.tsv'
TIMESTAMP_FORMAT = '%Y-%m-%d_%H-%M-%S'
BAD_INTERVALS_FILE_PATTERN_AS_PATH = Path(f'bad_intervals_{TIMESTAMP_FORMAT}.bed')

# define custom types
BadIntervalsDict = Dict[Tuple[str, int, int], List[int]]
GenomicIntervalList = List[Tuple[str, int, int, int]]

# module imports:
from utilities.plot_GC_matrices import plot_statistic_matrices, limit_extreme_outliers, smooth_2d_gc_weights
from utilities.plot_distributions import plot_fragment_length_dists, load_txt_to_matrix_with_meta, visualize_weights, \
    visualize_reconstruction_result
from utilities.secure_file_handling import AtomicOpen
from utilities.gc_logging import set_up_logging, set_new_log_paths, log, gib_cmd_logger

LOGGER = gib_cmd_logger()  # basically just to stop linter to assume it is None (was set to None in earlier version)


# definitions of specific exceptions
class OutOfGenomicBoundsError(Exception):
    """
    To be raised in conditions where a genomic locus (= chromosome/scaffold/contig + coordinate) does not exist in a
    given reference genome build, whichever that may be.
    """
    pass


# commandline interface
def get_cmdline_args():
    random_seed = random.randint(1, 1000)
    commandline_parser = ArgumentParser(
        description="""
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
                                                                             
"""
        f'             GCparagon ({VERSION_STRING}) maintained by @BGSpiegl\n'
        f'                 Copyright (c) 2023 Benjamin Spiegl\n'
        f'            GitHub: {github_url}\n\n'
        '^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^v^'
        'v^v^v^v^v^v^\n____________________________________________________'
        '_________________________',
        formatter_class=RawDescriptionHelpFormatter,
        epilog='GCparagon developed at the D&R Center of Molecular Biomedicine '
               '(https://www.medunigraz.at/en/research-centers-and-institutes/'
               'diagnostic-and-research-center-for-molecular-biomedicine), Institute '
               'of Human Genetics (https://humangenetik.medunigraz.at/en), Medical '
               'University of Graz, Austria (https://www.medunigraz.at/en)\n'
               'ASCII art created using textkool.com. DNA adapted from dariusz '
               'szenfeld (https://ascii.co.uk/art/dna)')
    # define argument groups
    input_args = commandline_parser.add_argument_group('Input (required)')
    output_args = commandline_parser.add_argument_group('Output options')
    processing_args = commandline_parser.add_argument_group('Processing options')
    postprocessing_args = commandline_parser.add_argument_group('Post-processing options')
    # input options
    input_args.add_argument('-b', '--bam', dest='input_bams', nargs='+', required=True, metavar='List[File]',
                            help='Path to sorted BAM file for which the fragment length-dependent GC-content-based '
                                 "over-representation (= 'GC-bias') should be computed and/or corrected. WARNING: "
                                 "don't use unaligned BAM files (uBAM) or multi-sample/run BAM files! If the BAM's "
                                 "index file is not found on runtime, GCparagon tries to create it. The alignment "
                                 "algorithm used for creating the input BAM file MUST follow the SAM format "
                                 "specifications! The TLEN column is used by GCparagon. [ PARAMETER REQUIRED ]")
    input_args.add_argument('-rtb', '--two-bit-reference-genome', dest='two_bit_reference_file',
                            default=EXPECTED_TWO_BIT_REFERENCE_GENOME_PATH,
                            help='Path to 2bit version of the reference genome FastA file which was used for read '
                                 'alignment of the input BAM file. If the 2bit version is missing, one can create the '
                                 'file using the following command: '
                                 "'faToTwoBit <PATH_TO_REF_FASTA> -long <PATH_TO_OUT_2BIT>' "
                                 "(see genome.ucsc.edu/goldenPath/help/twoBit.html for more details)", metavar='File')
    input_args.add_argument('-c', '--intervals-bed', dest='genomic_intervals_bed_file',
                            default=PREDEFINED_1MBP_INTERVALS_TO_PROCESS,
                            help='Path to BED file containing predefined genomic intervals to process. These should '
                                 'have been selected based on minimal overlap with exclusion-masked regions of the '
                                 'reference genome build used for read alignment earlier (i.e., creation of --bam). '
                                 'Since v0.5.6, the table also contains expected GC content counts which should be '
                                 "computed using the fragment length distribution from 'accessory_files/"
                                 "reference_fragment_lenght_distribution.tsv'. The GC content distributions are used "
                                 "to create an optimized consolidated weight matrix using a weighted mean. The weights "
                                 "for each region are selected such that the reference GC content distribution in "
                                 f"[ DEFAULT: '{PREDEFINED_1MBP_INTERVALS_TO_PROCESS}' ]", metavar='File')
    input_args.add_argument('-rgcd', '--reference-gc-content-distribution-table', dest='ref_gc_dist_path', 
                            default=DEFAULT_REFERENCE_GENOME_TARGET_GC_CONTENT_DISTRIBUTION,
                            help="Path to TSV file containing two data columns with header: 'gc_percentage', and "
                                 "'relative_frequency'. This table defines a GC content distribution (0%% GC to 100%% "
                                 "GC) as relative frequencies of these percentage bins which (summing up to 1). If a "
                                 "custom reference genome is used, this file should be created anew from the simulated "
                                 "genome-wide ideal fragment GC content as simulated assuming a fragment length "
                                 "distribution as the one stored in 'accessory_files/"
                                 "reference_fragment_lenght_distribution.tsv'! The provided information is used to "
                                 "optimize the combination of correction weight matrices from different genomic "
                                 "intervals to achieve a linear combination of these regions which resembles the "
                                 "reference GC content distribution defined here.", metavar='File')
    input_args.add_argument('-ec', '--exclude-intervals', dest='exclude_genomic_intervals_bed_file',
                            help='Path to library file (BED-like) holding DoC-specific definition of bad intervals '
                                 '(intervals must be exact genomic locus match for exclusion, DO NOT expect bedtools '
                                 'intersect-like behavior!). If the bad intervals library is left default, the bad '
                                 'intervals library with the most recent time stamp in the parent directory of the '
                                 'default library/BED file is used. The bad intervals library is intended to speed up '
                                 'the sample processing by excluding intervals with insufficient DoC form the '
                                 'beginning. Excluded intervals were observed to appear most frequently close to '
                                 'centromeres.', metavar='File')
    input_args.add_argument('-cw', '--correction-weights', dest='correction_weights', metavar='File',
                            help="Optional input for --tag-only mode: a matrix file ('*_gc_weights.txt.gz') containing "
                                 "correction weights to be used in tag-only mode ('--tag-only' flag must be set) to "
                                 'create a new, GC-bias-corrected BAM file with weighted fragments (GC-tag).')
    input_args.add_argument('-wm', '--weights-mask', dest='weights_mask', metavar='File',
                            help="Optional path to a weights matrix mask file. These are usually named '<SAMPLE_ID>_gc"
                                 "_bias_computation_mask.txt.gz'. If none is defined (default behaviour), either the "
                                 'currently create mask (when computing GC bias correction weights) or '
                                 '(for --tag-only) a mask file in the same input directory as the correction weights '
                                 'matrix defined via --correction-weights parameter is used based on the naming of the '
                                 'file. If none is specified or found, correction weights are reduced to rows and '
                                 'columns containing non-default values (values other than 1. and 0.).')
    # processing options
    # PRESET DEFINITIONS
    processing_args.add_argument('-p', '--preset', type=int, dest='parameter_preset_number',
                                 choices=range(0, 4, 1), default=DEFAULT_PRESET, metavar='Integer',
                                 help='Optional parameter preset to use for GC bias computation. Must be an integer '
                                      'int the rangeof 0-3 (inclusive). A preset value of 0 leaves parameters at '
                                      'default if not defined differently by the user (unchanged parameters will match '
                                      'preset 1). Other integer values from 1 to 3 define presets with increasing '
                                      'input data usage and required processing time (durations preset 1-3: 1-3 min, '
                                      '5-10 min, and ~1h depending on file size. Maximum across 4 samples and 2 '
                                      'iterations each computed using 12 cores and the profile_command.py script. '
                                      'Maximum memory consumption for any preset should stay below 4 GiB. If preset is '
                                      'not zero, any customized parameters conflicting with the preset will be '
                                      'ignored. A non-zero preset will set the following parameters: number of '
                                      'simulations, the target number of processed fragments, minimum number of '
                                      'fragment attribute combination occurrences, and the options for outlier '
                                      'detection and smoothing. Noise within the resulting correction weights is '
                                      'reduced when selecting a higher '
                                      'preset value. Preset 3 will attempt to process all genomic intervals (target '
                                      'number of fragments set to 100B) within the limits of the maximum allowed'
                                      'exclusion marked regions overlap (per default default ~1.7 Gb of reference are '
                                      'processed). NOTE: the percentage of total GC bias corrected fragments in the '
                                      'dataset for presets 1 vs. 3 increases only from 99.837%% to 99.938%% (average '
                                      'across 4 samples). Other fragment weights default to 1.0). The primary '
                                      'advantage of processing more fragments is the reduction of noise in computed '
                                      "weights and a better reconstruction of the reference fragment GC content "
                                      "distribution. It is recommended to use a higher preset for a 'preprocess-once,"
                                      "analyze often' scenario and/or when a high bias is expected/observed (e.g. "
                                      'FastQC average GC percentage). Correction by preset 1, 2, and 3 was found to '
                                      'yield 100.74%%, 99.98%%, and 99,91%% of the raw fragment count respectively '
                                      f'(average percentage across 4 samples). [ DEFAULT: {DEFAULT_PRESET} ]')
    # individual processing options
    processing_args.add_argument('-to', '--tag-only', dest='only_tag_bam', action='store_true',
                                 help='Optional flag which makes the software switch to tag-only mode. A correction '
                                      "weights matrix must be specified in this case via the '--correction-weights' "
                                      'flag. A valid samtools path must be available via the system path variable or '
                                      'provided using --samtools-path. Be mindful of setting the temporary directory '
                                      "correctly for your system! (e.g. should be set to output of 'echo $TEMP' on HPC "
                                      "clusters)")
    processing_args.add_argument('-rep', '--repetition-of-simulation', type=int, default=DEFAULT_SIMULATION_REPETITIONS,
                                 dest='n_simulations', metavar='Integer',
                                 help='(PRESET precedence if specified) This value can be left at default if the '
                                      'target number of processed fragments is sufficiently high (e.g. >=5M). The '
                                      'lower the number of target fragments, the stronger is the effect of increasing '
                                      'the number of simulation rounds. Increasing this value increases the '
                                      'computation time almost accordingly (scales linearly). [ DEFAULT: '
                                      f'{DEFAULT_SIMULATION_REPETITIONS} ]')
    processing_args.add_argument('-mafl', '--maximum-fragment-length', dest='upper_limit_fragment_length', type=int,
                                 default=DEFAULT_MAX_FRAGMENT_LENGTH, metavar='Integer',
                                 help=f'Defines upper length limit for fragments which should be included in '
                                      f'computation. This parameter does not impact computation speed. It only '
                                      f'increases plotting times for matrices by a few seconds and memory consumption. '
                                      f'[ DEFAULT: {DEFAULT_MAX_FRAGMENT_LENGTH}bp ]')
    processing_args.add_argument('-mifl', '--minimum-fragment-length', dest='lower_limit_fragment_length', type=int,
                                 default=DEFAULT_MIN_FRAGMENT_LENGTH, metavar='Integer',
                                 help=f'Defines lower length limit for fragments which should be included in '
                                      f'computation. Must be positive integer. A value below the sequenceable '
                                      f'fragment length of the device used to create the dataset is not recommended.'
                                      f' [ DEFAULT: {DEFAULT_MIN_FRAGMENT_LENGTH}bp ]')
    processing_args.add_argument('-t', '--threads', dest='total_number_threads', type=int,
                                 metavar='Integer',
                                 help=f'Total number of threads to be used for BAM processing. If the '
                                      '--single-thread-processes flag was set, this number corresponds to the number '
                                      'of processes spawned for BAM processing. For BAM tagging, multiple threads are '
                                      'used for the sort/merge operations so fewer processes might be used '
                                      'simultaneously. Should be lower than the total number of logical cores '
                                      'available on the hardware. Will be reduced to max. available number of logical '
                                      f'cores if is set higher by the user. [ DEFAULT: {DEFAULT_NUMBER_PROCESSES} ]')
    processing_args.add_argument('-rs', '--random-seed', default=random_seed, type=int, dest='random_seed',
                                 help='Optional random seed to be used for genomic sampling patterns. '
                                      'Warning: the notion that all computed numbers will turn out identical when '
                                      'using the same random seed using different interpreters or different machines '
                                      'should be discarded right away. Might only be useful when repeatedly running '
                                      'the script within the same python interpreter instance! [ DEFAULT: '
                                      f'{random_seed} (randomly drawn from 0-999) ]', metavar='Integer')
    processing_args.add_argument('-sp', '--samtools-path', dest='samtools_path', default=DEFAULT_SAMTOOLS_PATH,
                                 help='Optional input: path to specific samtools executable. A valid path is required '
                                      'for creating the tagged BAM output. By default, this path will be used: '
                                      f"'{DEFAULT_SAMTOOLS_PATH}' (empty or None if path is not found). "
                                      "Code tested with samtools version 1.16.1 using htslib 1.16 "
                                      "[ PARAMETER REQUIRED IF DEFAULT VALUE IS EMPTY/NONE ]",
                                 metavar='File')
    processing_args.add_argument('-nf', '--target-fragment-number', dest='process_n_fragments',
                                 default=DEFAULT_TARGET_NUMBER_FRAGMENTS_PROCESSED, type=int, metavar='Integer',
                                 help='(PRESET precedence if specified) GC-bias computation will stop after surpassing '
                                      'this threshold for processed fragments. Still running subprocesses will be '
                                      'finished and results included so usually this value is overshot by up to '
                                      'several million fragments depending on the amount of processes chosen and the '
                                      'DoC inside defined bins. Increasing this value will reduce the noise in '
                                      'computed weights. Concerning the number of corrected fragments, processing more '
                                      'than 5 million fragments will only increase the number of corresponding '
                                      'computed weights only miniscule. Doubling the target processed fragment amount '
                                      'typically leads only to an increase in corrected fragments by less than one '
                                      'percent. Five million fragments (preset 1) should be enough to correct between '
                                      '99.5%% and 99.9%% of all DNA fragments observed in the dataset based on GC '
                                      'content and fragment length. To reach >99.9%% of corrected fragments, this '
                                      'parameter should be increased. [ DEFAULT: '
                                      f'{DEFAULT_TARGET_NUMBER_FRAGMENTS_PROCESSED:,} ]')
    processing_args.add_argument('-mf', '--minimum-fragment-occurrences', dest='min_frag_occurs',
                                 default=DEFAULT_MIN_OCCURRENCES, type=int, metavar='Integer',
                                 help='(PRESET precedence if specified) This parameter defines the minimum number of '
                                      'fragment occurrences for a specific length/GC-content attribute combination to '
                                      'be regarded in correction weights computation (= mask definition). Higher '
                                      'values result in less extreme weight outliers, especially for the '
                                      'low-GC-content mononucleosomal fragment length range. The absolute lowest '
                                      f'supported value is {ABSOLUTE_MIN_OCCURRENCES} based on visual inspection of '
                                      'resulting weights matrices for different samples. If this value is too low '
                                      "(e.g. 1), strong 'salt-and-pepper'-type noise was observed for rare attribute "
                                      'combinations along with very high weight outliers. A value of 10 here means '
                                      'that a particular attribute combination must occur at least once per million '
                                      "fragments in the dataset for a '--number-of-fragments-to-process' value of "
                                      '10,000,000. As a rule of thumb, one can set this to number of million target '
                                      'fragments (i.e. set to 10 for the target value of 10M processed fragments as '
                                      f'in the example above) [ DEFAULT: {DEFAULT_MIN_OCCURRENCES} ]')
    processing_args.add_argument('-anf', '--allow-n-base-fragments', action='store_true', dest='allow_n_base_fragments',
                                 help="Per default, any fragment containing N-bases (as determined from the read "
                                      "alignment positions and the reference genome sequence) is excluded from the "
                                      "analysis. This parameter was not found to cause any problems for Illumina "
                                      "NovaSeq data. If such fragments have to be included, this flag can be set to "
                                      "allow for up to 1/3 N-bases for fragments. Parameter mainly influences the "
                                      "simulation step and how many times random fragment drawing must be repeated for "
                                      "individual genomic intervals. Also can lead to fewer intervals being discarded "
                                      "(and marked as bad genomic interval) if flag is set.")
    processing_args.add_argument('-ucmaf', '--unclipped-min-aln-fraction', dest='min_unclipped_aln_fraction',
                                 default=DEFAULT_MIN_UNCLIPPED_ALN_FRACTION, type=float, metavar='Float',
                                 help='This parameter defines the minimum unclipped fraction of an alignment to be '
                                      'counted in the observed fragment attributes matrix O_gc. This might affect how '
                                      'many small fragments are observed and effectively corrected. [ DEFAULT: '
                                      f'{DEFAULT_MIN_UNCLIPPED_ALN_FRACTION} ]')
    processing_args.add_argument('-dw', '--default-weight', default=DEFAULT_WEIGHT, type=float,
                                 dest='default_fragment_weight',
                                 help='Parameter redefines the weight which is assigned to fragments with fragment '
                                      'length + GC base count combinations that lie outside of the non-default range '
                                      'of the (computed) weights matrix. Should be 1.0 for GCparagon. '
                                      'Can be e.g. 0.0 for other algorithms like Griffin. Choose according to the '
                                      f'source of your weights matrix! [ DEFAULT: {DEFAULT_WEIGHT} ]')
    processing_args.add_argument('-reto', '--reads-extraction-timeout', dest='unaligned_extraction_timeout',
                                 default=READS_EXTRACTION_TIMEOUT,
                                 help='Sets the timout in seconds for unaligned reads extraction. '
                                      "Only has an effect if '--output-bam' or '--tag-only' and "
                                      "'--output-unaligned-reads' is set. "
                                      f"[ DEFAULT: {READS_EXTRACTION_TIMEOUT} seconds ]")
    # post-processing options
    postprocessing_args.add_argument('-do', '--detect-outliers', action='store_true', dest='detect_outliers',
                                     help='(PRESET precedence if specified) If this flag is set, extreme outliers will '
                                          'be detected and limited to a threshold value that is computed from the '
                                          'fragment weights. The default method to detect outliers is '
                                          'Q3 + 8x inter-quartile range (IQR). Values above this threshold will be '
                                          'limited to the threshold. It is highly recommended to detect and limit '
                                          'outliers.')
    postprocessing_args.add_argument('-odm', '--outlier-detection-method', dest='outlier_method',
                                     default=DEFAULT_OUTLIER_DETECTION_METHOD, choices=['IQR', 'SD'],
                                     help='(PRESET precedence if specified) If the --detect-outliers flag is set, the '
                                          'detection method can be set here. Either a method based on the '
                                          'inter-quartile range or a method based on standard deviation can be '
                                          "selected. Must be one of {'IQR', 'SD'}. "
                                          f'[ DEFAULT: {DEFAULT_OUTLIER_DETECTION_METHOD} ]',
                                     metavar='String')
    postprocessing_args.add_argument('-ods', '--outlier-detection-stringency', choices=range(0, 8, 1), type=int,
                                     dest='outlier_stringency', default=DEFAULT_OUTLIER_DETECTION_STRINGENCY,
                                     help='(PRESET precedence if specified) If the --detect-outliers flag is set, this '
                                          'parameter defines how stringent the outlier detection threshold is set. '
                                          'Must be an integer in the range of 1-7 (inclusive). [ DEFAULT: '
                                          f'{DEFAULT_OUTLIER_DETECTION_STRINGENCY} ]', metavar='Integer')
    postprocessing_args.add_argument('-sw', '--smooth', action='store_true', dest='smooth_weights',
                                     help='(PRESET precedence if specified) If this flag is set, computed weights will '
                                          'also be smoothed. An additional matrix is output containing these '
                                          'post-processed values. If plotting is set to true, also a visualisation of '
                                          'the smoothed weights will be created.It is recommended to smooth weights if '
                                          'not the entire dataset is processed (like is done in preset 3).')
    postprocessing_args.add_argument('-sk', '--smooth-kernel', choices=['gauss', 'constant'],
                                     default=DEFAULT_SMOOTHING_KERNEL, dest='smoothing_kernel',
                                     help="(PRESET precedence if specified) If the '--smooth' flag is set, the type of "
                                          'kernel used in the 2D convolution operation can be set here. In general, a '
                                          'Gaussian kernel makes more sense because it assigns directly adjacent '
                                          'values a higher weight in computing the smoothed value of the current '
                                          "position. Must be one of {'gauss', 'constant'}. "
                                          f'[ DEFAULT: {DEFAULT_SMOOTHING_KERNEL} ]', metavar='String')
    postprocessing_args.add_argument('-si', '--smoothing-intensity', dest='smoothing_intensity', type=int,
                                     choices=range(1, 11, 1), default=DEFAULT_SMOOTHING_INTENSITY,
                                     help="(PRESET precedence if specified) If the '--smooth' flag is set, the "
                                          'smoothing intensity defines the range of the 2D kernel used in the '
                                          'smoothing operation. Must be an integer in the range of 1-10 (inclusive). '
                                          f'[ DEFAULT: {DEFAULT_SMOOTHING_INTENSITY} ]', metavar='Integer')
    # output options
    output_args.add_argument('-v', '--verbose', action='store_true', dest='verbose',
                             help='This flag can be set to provide more output, especially information about fragments '
                                  'that fall outside of the defined fragment length window.')
    output_args.add_argument('-o', '--out-dir', dest='output_directory', metavar='File',
                             help='Path to which output is moved from --temporary-directory after each processing step '
                                  '(GC-bias computation, BAM tagging). The directory will be created if it does '
                                  'not exist. Make sure that it is empty if it exists, otherwise the whole directory '
                                  'will be deleted before writing to it in the GC-bias computation step! '
                                  'If none is provided, a new subdirectory named '
                                  f"'GC_bias_correction_GCparagon{VERSION_STRING}' will be created in the input BAM's "
                                  'parent directory and used as output directory. The output for each sample will be '
                                  'gathered in a subdirectory of this --out-dir which will be named after the sample. '
                                  'The output directory may be located on slow hardware such as a USB drive or a '
                                  'network storage since everything is stored in --temporary-directory first and moved '
                                  'after completion of all defined phases of the GC bias computation.')
    output_args.add_argument('-tmp', '--temporary-directory', dest='temporary_directory', metavar='File',
                             default=DEFAULT_TEMPORARY_DIRECTORY,
                             help='Directory to which all files will be written as temporary files in a subdirectory '
                                  'named after the sample (sample id is extracted from the BAM file name) during '
                                  'processing. Directory will be created if non-existent. Subdirectory for the sample '
                                  'will be deleted if it exists initially. If not specified, this directory is '
                                  "identical to the output of Python's tempfile module's gettempdir() function. "
                                  "Permanent non-temporary output files will be MOVED to the --out-dir using shutil's "
                                  'move function from this directory. The temporary directory should be located on a '
                                  f'high performance hardware (high IOPs)! [ DEFAULT: {DEFAULT_TEMPORARY_DIRECTORY} ]')
    output_args.add_argument('-np', '--no-plots', action='store_false', dest='plot_result',
                             help='Flag suppresses creation of fragment length distribution plot and heatmaps for '
                                  'observed, expected, correction, and computation mask matrices.')
    output_args.add_argument('-os', '--output-simulations', action='store_true', dest='output_simulation_results',
                             help='Optional flag for GC-bias computation for plotting individual simulation results '
                                  '(simulated fragments and iteration-specific masks). The simulated fragment '
                                  'attribute distributions and computation masks are plotted for all simulations then.')
    output_args.add_argument('-k', '--keep-interval-data', dest='keep_interval_data', action='store_true',
                             help='Optional flag which can be used to save intermediate data per genomic interval.')
    output_args.add_argument('-ob', '--output-bam', dest='output_corrected_bam', action='store_true',
                             help='Optional flag to activate writing of the GC-correction-weights-tagged BAM file '
                                  'AFTER COMPUTING GC BIAS (--tag-only flag is not set), either using the statistics '
                                  'computed from the input BAM file or a correction weights matrix specified via '
                                  "--correction-weights. WARNING: currently, the output BAM won't contain "
                                  'unaligned reads!')
    output_args.add_argument('-our', '--output-unaligned-reads', dest='output_unaligned_reads', action='store_true',
                             help='Optional flag to activate writing of unaligned reads to a separate BAM file. '
                                  'Per default, unaligned reads are not output. Setting this flag '
                                  'only has an effect if either the --output-bam flag was set or GCparagon was '
                                  'started in the --tag-only mode.')
    output_args.add_argument('-fp', '--float-precision', dest='floating_point_precision',
                             default=DEFAULT_FLOAT_PRECISION, type=int, metavar='Integer',
                             help='Optional parameter for GC-bias computation number of digits after the comma for '
                                  'floating point data to be stored in text-based matrices, e.g. for correction '
                                  'weights data. Choose according to expected depth of coverage -> if you would '
                                  'expect 10,000, you can go for 5 or even 6 digits. Otherwise this will not have an '
                                  'effect. If you compute signals that are sums over many regions, multiply the '
                                  'expected DoC with how many signals you sum up to get an estimate of which precision '
                                  'you would need to definitively be able to rule out any influence by rounding '
                                  f'errors. These should average out though. [ DEFAULT: {DEFAULT_FLOAT_PRECISION} ]')
    output_args.add_argument('-tg', '--tag-name', dest='gc_tag_name', default=DEFAULT_TAG_NAME, metavar='String',
                             help='Name of the GC-bias correction weight tag that will be added to alignments in the '
                                  'BAM file. If none is provided, the default tag will be used. Must not be longer '
                                  f'than 2 characters! [ DEFAULT: {DEFAULT_TAG_NAME} ]')
    output_args.add_argument('-wie', '--write-interval-exclusion', dest='write_updated_bad_intervals_library',
                             action='store_true',
                             help='Optional flag for writing an updated version of the library listing intervals '
                                  'marked for exclusion from the analysis. Per default, genomic intervals are marked '
                                  'for exclusion if drawing fragments of a specific size repeatedly fails (at least 55 '
                                  'times (for strict reference N base handling, 33 times otherwise) or 1/3 of number '
                                  'of fragments that need to be drawn, whichever is higher) due '
                                  'to getting only poly-N sequences. In general, the frequency of these exclusion '
                                  'events is dependent on the DoC of the sample, which can be substituted by the '
                                  'number of fragments estimated to be obtained from all predefined intervals in BAM '
                                  "file in a first approximation. WARNING: don't mix exclusion-marked interval "
                                  'libraries computed from different (predefined) interval BED files! If the user '
                                  'places the output BED file library in the default directory, the new library will '
                                  'be used per default for future computations. Genomic intervals will be marked for '
                                  "exclusion depending on a data set's fragment length distribution and sequencing "
                                  "depth.")
    output_args.add_argument('-nfp', '--no-focused-plots', action='store_true', dest='dont_focus_plots',
                             help='Optional flag to deactivate focusing of matrix plots on non-default values (focus '
                                  'uses a border of up to 10 default values). Only has an effect if --no-plots flag is '
                                  'not set.')
    output_args.add_argument('-sf', '--show-figures', action='store_true', dest='show_plots',
                             help='Optional flag to display plots in an interactive browser window in addition to '
                                  'saving them to a file.')
    return commandline_parser.parse_args()


def create_exception_stack_trace(e):
    exception_list = traceback.format_stack()
    exception_list = exception_list[:-2]
    exception_list.extend(traceback.format_tb(sys.exc_info()[2]))
    exception_list.extend(traceback.format_exception_only(sys.exc_info()[0], sys.exc_info()[1]))
    exception_str = "Traceback (most recent call last):\n"
    exception_str += "".join(exception_list)
    exception_str = exception_str[:-1]  # remove last newline
    return exception_str


@contextlib.contextmanager
def silently_open_alignment_file(*args, **kwargs):
    error_sink = os.open('/dev/null', os.O_WRONLY)
    f_aln = error_sink  # linter silencing
    try:
        if LOGGER:
            f_aln = AlignmentFile(*args, **kwargs)
            yield f_aln  # yield because f_aln can be iterated over -> so a generator makes sense
        else:  # handle the default streams
            std_out_saved = os.dup(1)  # save current stdout stream target
            os.dup2(error_sink, 1)  # redirect stdout to /dev/null (NULL-pointer)
            std_err_saved = os.dup(2)  # save current stderr stream target
            os.dup2(error_sink, 2)  # redirect stderr to /dev/null (NULL-pointer)
            # redirect both to /dev/null (NULL-pointer)
            f_aln = AlignmentFile(*args, **kwargs)
            # redirect back to previous files/streams (does not seem to work)
            os.dup2(std_out_saved, 1)  # load former stdout stream target
            os.dup2(std_err_saved, 2)  # load former stderr stream target
            yield f_aln  # yield because f_aln can be iterated over -> so a generator makes sense
    finally:
        f_aln.close()


# GC-bias computation functions:
# ----------------------------------------------------------------------------------------------------------------------
def read_bed_file(bed_path: str, header=False) -> List[Tuple[str, int, int, Any]]:
    """

    :param bed_path:
    :return:
    """
    with AtomicOpen(bed_path, 'rt') as f_bed:
        if header:
            _hdr = f_bed.readline()
        return_data = [(chrom, int(region_start), int(region_stop),
                        meta_info)  # list of ints -> exclusion_marked_bases followed by fragment GC content counts
                       for chrom, region_start, region_stop, *meta_info
                       in filter(lambda x: x not in ('', '\n', None),
                                 [bed_line.strip().split() for bed_line in f_bed.readlines()])]
        return return_data


def read_bad_genomic_intervals_bed_file(bed_path: str) -> BadIntervalsDict:
    bad_intervals = {}
    interval_lengths = set()
    try:
        for chrm, strt, stp, meta_inf_list in read_bed_file(bed_path=bed_path):
            interval_lengths.update({stp - strt})
            # meta_inf_list: yields of samples (total number of reads) as comma-delimited list of integers;
            # added for cutoff -> some samples with low yield are more likely to run into the problem of randomly not
            # drawing enough N-bases-free fragments than higher yield samples
            # -> improved pre-rejection of genomic intervals:
            # reject all intervals that failed in the fragment drawing process for samples with HIGHER yield!
            # I don't think this has been implemented yet
            if bad_intervals.get((chrm, strt, stp)) is None:
                bad_intervals.update({(chrm, strt, stp): list(map(lambda v: int(v), meta_inf_list[0].split(',')))})
            else:  # interval locus exists already (= multiple entries! Not expected)
                bad_intervals[(chrm, strt, stp)].extend(list(map(lambda v: int(v), meta_inf_list[0].split(','))))
    except ValueError:  # not enough values to unpack (expected 4, got X)
        log(message=f"Could not load bad regions from BED file '{bed_path}' (requires column 4 and column 5 to contain "
                    "values that can be cast to int!). Returning no bad regions instead.",
            log_level=logging.WARNING, logger_name=LOGGER)
    if len(interval_lengths) > 1:
        log(message=f"Genomic intervals of different length encountered in file '{bed_path}'. "
                    f"A BED file containing mixed bad intervals was provided!",
            log_level=logging.WARNING, logger_name=LOGGER)
    return bad_intervals


def create_region_label(chrm: str ,start: int, end: int):
    return f"{chrm}_{start:,}-{end:,}"


def read_scored_regions_bed_file(bed_path: str):
    scored_regions = []
    region_gc_content_distributions = {}
    try:
        for chrm, strt, stp, meta_info in read_bed_file(bed_path=bed_path, header=True):
            scored_regions.append((chrm, strt, stp, int(meta_info[0])))
            # process meta_info
            assert len(meta_info[1:]) == 101
            assert region_gc_content_distributions.get(create_region_label(chrm=chrm, start=strt, end=stp)) is None
            component_counts = np.array(list(map(lambda s: int(s), meta_info[1:])))
            region_gc_content_distributions[create_region_label(chrm=chrm, start=strt, end=stp)] = \
                component_counts / component_counts.sum()  # NORMALIZE TO 1!!!
        return scored_regions, region_gc_content_distributions
    except ValueError:  # not enough values to unpack (expected 4, got X)
        log(message=f"Could not load scored regions from BED file '{bed_path}' (requires column 4 to contain values "
                    "that can be cast to int!). Terminating ..", log_level=logging.ERROR, logger_name=LOGGER)
        sys.exit(1)


def save_matrix_to_txt(matrix: OneOf[np.array, np.matrix], filename: OneOf[str, Path], output_dir: str,
                       gzipped=True, float_data_precision=6, min_frag_length=DEFAULT_MIN_FRAGMENT_LENGTH,
                       max_frag_length=DEFAULT_MAX_FRAGMENT_LENGTH, verbose=False, report_saved_path=False) \
        -> Optional[Path]:
    """

    :param report_saved_path:
    :param verbose:
    :param matrix:
    :param filename:
    :param output_dir:
    :param gzipped:
    :param float_data_precision:
    :param min_frag_length:
    :param max_frag_length:
    :return:
    """
    if gzipped and filename.lower()[-3:] != '.gz':
        filename += '.gz'
    output_path = Path(output_dir) / filename
    # check for data type and select formatter accordingly
    if matrix.dtype in (np.float64, np.float32, np.float16, float):
        formatter_str = f'%1.{float_data_precision}f'
    else:
        if matrix.dtype not in (bool, np.uint64, np.int64, np.int32, np.uint32, np.int16, np.uint16, int, np.uint):
            log(message="A numpy dtype of numpy.float, numpy.integer or bool was expected. Don't know how to "
                        f"handle dtype {matrix.dtype}.", log_level=logging.CRITICAL, close_handlers=True,
                logger_name=LOGGER)
            sys.exit(2)
        formatter_str = '%2.1u'
    if verbose:
        log(message=f"Saving statistic matrix '{filename}' to output directory",
            log_level=logging.INFO, logger_name=LOGGER)
    np.savetxt(output_path, matrix, delimiter=' |', fmt=formatter_str,
               header=f"rows representing fragment lengths ({min_frag_length} bp to {max_frag_length} bp) and columns "
                      f"representing GC base count from zero to max. fragment length included")
    if report_saved_path:
        return output_path


def rhui(num: OneOf[int, float]):
    """

    :param num:
    :return:
    """
    return int(math.floor(num + 0.5))


def check_region_within_genomic_bounds(chrom: str, start: int, stop: int, chromosome_sizes: Dict[str, int],
                                       raise_error=False) -> bool:
    outside = False
    if start not in range(0, chromosome_sizes[chrom]) or stop not in range(1, chromosome_sizes[chrom] + 1):
        outside = True

    if outside:
        if raise_error:
            log(message=f"Region '{chrom}:{start}-{stop}' not within genomic scaffold "
                        f"{chrom}:0-{chromosome_sizes[chrom]}", log_level=logging.CRITICAL, close_handlers=True,
                logger_name=LOGGER)
            raise OutOfGenomicBoundsError()
        return False
    return True


def compute_observed_attributes_matrix(bam_file: str, two_bit_reference_path: OneOf[str, Path], chromosome: str,
                                       tmp_dir: str, sample_id: str, start_coord: int, stop_coord: int,
                                       save_individual_matrices=False, strict_n_ref_bases_handling=True,
                                       frag_n_cont_thresh=DEFAULT_FRAGMENT_N_CONTENT_THRESHOLD,
                                       min_frag_len=DEFAULT_MIN_FRAGMENT_LENGTH, float_precision=6,
                                       max_frag_len=DEFAULT_MAX_FRAGMENT_LENGTH,
                                       min_unclipped_aln_fraction=DEFAULT_MIN_UNCLIPPED_ALN_FRACTION) \
        -> Tuple[OneOf[np.array, np.matrix], int, int]:
    """

    :param two_bit_reference_path:
    :param strict_n_ref_bases_handling:
    :param sample_id:
    :param frag_n_cont_thresh:
    :param tmp_dir:
    :param bam_file:
    :param chromosome:
    :param start_coord:
    :param stop_coord:
    :param save_individual_matrices:
    :param float_precision:
    :param min_frag_len:
    :param max_frag_len:
    :param min_unclipped_aln_fraction:
    :return:
    """
    with silently_open_alignment_file(bam_file, mode='rb') as f_aln:
        # Probably the old intervals loading strategy was better (slicing 500k chars is easier than permanent IO)
        observed_attributes_matrix = np.zeros((max_frag_len - min_frag_len + 1, max_frag_len + 1), dtype=float)
        ref_genome_handle = TwoBitFile(str(two_bit_reference_path))
        chromosome_handle = ref_genome_handle[chromosome]
        try:  # assert: template has between 1 and 2 segments (single-end or paired-end data; SAM format supports more)
            # The leftmost segment has a plus sign and the rightmost has a minus sign. It is set as 0 for
            # single-segment template or when the information is unavailable; template length must be positive and
            # within defined range
            # LEGACY FILTER:
            # filtered_alignments = filter(lambda p: (min_frag_len <= p.template_length <= max_frag_len) and
            #                              p.is_proper_pair and not p.is_supplementary and not p.is_secondary,
            #                              (read for read in f_aln.fetch(chromosome, start_coord, stop_coord,
            #                                                            multiple_iterators=True)))

            # binary filter:
            exclude_flags = np.uint32(3852)  # = 256 + 2048 + 512 + 1024 + 4 + 8
            exclude_flags_binary = bin(exclude_flags)
            # paired_flag = np.uint32(1)  # only paired!
            # paired_flag_binary = bin(paired_flag)  # not used -> replaced by simply checking a.is_paired
            # -> not necessary to check for "mapped" attribute (filtered out if "read unmapped")
            # complete alignment filter:
            # --------------------------
            # EXCLUDE if:
            # read unmapped = 4               '0b100'
            # mate unmapped = 8               '0b1000'
            # not primary (secondary) = 256   '0b100000000'
            # vendor/QC fail = 512            '0b1000000000'
            # PCR or optical duplicate = 1024 '0b10000000000'
            # supplementary = 2048            '0b100000000000'
            # = 3852
            # REQUIRE THAT:
            # alignment is paired = 1 '0b1'
            # mates map to different strands
            #    a.is_forward != a.mate_is_forward
            # TLEN column is (positive and) between defined fragment length limits (inclusive)
            #    min_frag_len <= a.template_length <= max_frag_len
            filtered_alignments = filter(lambda a:
                                         # bin(np.uint32(a.flag) & paired_flag) == paired_flag_binary and
                                         a.is_paired and
                                         bin(~np.uint32(a.flag) & exclude_flags) == exclude_flags_binary and
                                         a.is_forward != a.mate_is_forward and
                                         (min_frag_len <= a.template_length <= max_frag_len),
                                         f_aln.fetch(chromosome, start_coord, stop_coord, multiple_iterators=True))
            # binary filter UPGRADE:
            #    (circumvent fragment length distribution based non-proper-pair-flagging of alignments)
            # EXCHANGE PROPER PAIRED: instead of the aln_seg.is_proper_pair (2), use fragment length (already there)
            #                         use exclusion if (read reverse strand + mate reverse strand) = 48 or none
            #                         of both -> use inequality of these!
            # # binary filter:
            #  exclude_flags = np.uint32(3852)  # = 256 + 2048 + 512 + 1024 + 4 + 8
            #  exclude_flags_binary = bin(exclude_flags)
            #  paired_flag = np.uint32(1)  # only paired!
            #  paired_flag_binary = bin(prop_paired_r1_flags)
            #
            # # complete alignment filter:
            # filtered_alignments = filter(
            #                 lambda a: bin(np.uint32(a.flag) & paired_flag) == paired_flag_binary and
            #                 bin(~np.uint32(p.flag) & exclude_flags) == exclude_flags_binary and
            #                 a.is_forward != a.mate_is_forward and
            #                 (min_frag_len <= a.template_length <= max_frag_len),
            #                 f_alns.fetch(chromosome, start_coord, stop_coord, multiple_iterators=True))
        except OSError as e:  # maybe truncated file
            log(message=f"Encountered a problem while trying to parse file '{bam_file}':  {e}. File is likely either "
                        f"truncated or these coordinates were out of bounds: {chromosome}:{start_coord}-{stop_coord}"
                        f". Terminating..", log_level=logging.ERROR, close_handlers=True, logger_name=LOGGER)
            sys.exit(2)
        n_fragments_processed = 0
        ignored_fragments = 0
        # file must still be open, otherwise we will get an OSError; implicitly increment via enumeration
        if strict_n_ref_bases_handling:  # all fragments containing N-bases are rejected; extensive looping here!
            for n_fragments_processed, aln_segment in enumerate(filtered_alignments, start=1):
                frag_start = min(aln_segment.pos, aln_segment.pnext)
                frag_length = aln_segment.template_length  # pysam: "the observed query template length"
                # -> filtered for inside flength range if PE reads; SE reads not supported at the moment
                try:  # below yields str, takes fragment length, computed GC bases and increments counts matrix
                    if aln_segment.reference_length < aln_segment.query_length * min_unclipped_aln_fraction:
                        # the reference_length (= alen) comparison (instead of flength which reduces the exclusion
                        # solely to fragments shorter than "min_unclipped_aln_fraction times the read length") gets rid
                        # of artifactual alignment observations in the low GC base count space that account for some
                        # noise in the bias matrix
                        raise IndexError  # min. 3/4 must be identical to ref seq per default or fragment is discarded
                    fragment_sequence = chromosome_handle[frag_start:frag_start + frag_length].upper()
                    gc_count = gc_count_rejecting_n_containing(f_seq=fragment_sequence)  # might return high number
                    observed_attributes_matrix[frag_length - min_frag_len, gc_count] += 1
                except IndexError:  # out of bounds for extreme fragment lengths OR N-containing fragments above thrsh.
                    ignored_fragments += 1
                    continue
        else:  # allow N bases until a certain threshold; extensive looping here!
            for n_fragments_processed, aln_segment in enumerate(filtered_alignments, start=1):
                frag_start = min(aln_segment.pos, aln_segment.pnext)
                frag_length = aln_segment.template_length  # pysam: "the observed query template length"
                # -> filtered for positive
                try:  # below yields str, takes fragment length, computed GC bases and increments counts matrix
                    if frag_length < int(aln_segment.query_length * min_unclipped_aln_fraction):
                        raise IndexError  # min. 3/4 must be identical to ref seq per default or fragment is discarded
                    fragment_sequence = chromosome_handle[frag_start:frag_start + frag_length].upper()
                    gc_count = safe_gc_base_count_inference_thresh(f_seq=fragment_sequence, f_len=frag_length,
                                                                   threshold=frag_n_cont_thresh)
                    observed_attributes_matrix[frag_length - min_frag_len, gc_count] += 1
                except IndexError:  # N-containing fragments above upper threshold
                    ignored_fragments += 1
                    continue
    # give feedback and store the observed attributes matrix
    if ignored_fragments:
        log(message=f"Done processing {n_fragments_processed:,} fragments for region '{chromosome}:"
                    f"{start_coord:,}-{stop_coord:,}'. Of these, {ignored_fragments:,} fragments either did not fall "
                    f"into the selected fragment length interval of [{min_frag_len}, {max_frag_len}] bp, or were "
                    f"(soft-)clipped beyond {min_unclipped_aln_fraction:.1%}.",
            log_level=logging.DEBUG, logger_name=LOGGER)
    if save_individual_matrices:
        target_path = Path(tmp_dir)
        target_path.mkdir(parents=True, exist_ok=True)
        interval_str = f"{chromosome}_{start_coord}-{stop_coord}"
        save_matrix_to_txt(matrix=observed_attributes_matrix, output_dir=str(target_path), max_frag_length=max_frag_len,
                           filename=f'{sample_id}_observed_attributes_matrix_{interval_str}.txt.gz',
                           float_data_precision=float_precision, gzipped=True, min_frag_length=min_frag_len)
    return observed_attributes_matrix, n_fragments_processed, ignored_fragments


def gc_count_rejecting_n_containing(f_seq):
    n_count = f_seq.count('N')
    if n_count:  # too many N-bases compared to fragment length
        return 99999999  # will lead to an IndexError
    return f_seq.count('G') + f_seq.count('C')


def safe_gc_base_count_inference_thresh(f_seq, f_len, threshold=0.3):
    """
    Execution takes around 30-60 micro seconds
    :param threshold:
    :param f_seq:
    :param f_len:
    :return:
    :raises: nothing but may lead to IndexError
    """
    n_count = f_seq.count('N')
    if n_count / f_len >= threshold:  # too many N-bases compared to fragment length
        return 99999999  # will lead to an IndexError
    try:
        return rhui(f_len / (f_len - n_count) * (f_seq.count('G') + f_seq.count('C')))
    except ZeroDivisionError:  # if all bases returned are 'N's
        return 99999999  # will lead to an IndexError


def simulate_fragment_attributes(two_bit_reference_path: OneOf[str, Path], tmp_dir: str, chromosome: str,
                                 start_coord: int, sample_id: str, stop_coord: int, expected_yield: int,
                                 statistic_matrix: OneOf[str, np.ndarray, np.matrix],
                                 float_precision=6, save_individual_matrices=False, strict_n_ref_bases_handling=True,
                                 random_seed=RANDOM_SEED, simulation_repetitions=DEFAULT_SIMULATION_REPETITIONS,
                                 min_frag_len=DEFAULT_MIN_FRAGMENT_LENGTH, max_frag_len=DEFAULT_MAX_FRAGMENT_LENGTH,
                                 frag_n_cont_thresh=DEFAULT_FRAGMENT_N_CONTENT_THRESHOLD) \
        -> OneOf[Tuple[np.ndarray, Tuple[np.ndarray]], Tuple[str, int, int, int]]:
    """
    May add intervals to bad intervals library
    :param max_frag_len:
    :param strict_n_ref_bases_handling:
    :param sample_id:
    :param frag_n_cont_thresh:
    :param expected_yield:
    :param tmp_dir:
    :param min_frag_len:
    :param save_individual_matrices:
    :param float_precision:
    :param two_bit_reference_path:
    :param chromosome:
    :param start_coord: inclusive
    :param stop_coord: exclusive
    :param statistic_matrix:
    :param simulation_repetitions:
    :param random_seed:
    :return:
    """
    ref_genome_handle = TwoBitFile(str(two_bit_reference_path))
    # check if interval coordinates are valid
    try:
        if start_coord < 0 or stop_coord > ref_genome_handle.sequence_sizes()[chromosome]:
            log(message=f"Check your coordinates! Your provided {chromosome}:{start_coord}-{stop_coord}",
                log_level=logging.CRITICAL, close_handlers=True, logger_name=LOGGER)
            raise OutOfGenomicBoundsError()
    except KeyError:
        log(message=f"The following contig/scaffold was not found in the 2bit reference genome file: {chromosome}",
            log_level=logging.CRITICAL, close_handlers=True, logger_name=LOGGER)
        sys.exit(2)
    chromosome_handle = ref_genome_handle[chromosome]
    ref_seq_interval_slice = chromosome_handle[start_coord:stop_coord].upper()
    if isinstance(statistic_matrix, np.ndarray) or isinstance(statistic_matrix, np.matrix):
        observed_attributes_matrix = statistic_matrix
    elif isinstance(statistic_matrix, str) and Path(statistic_matrix).is_file():
        observed_attributes_matrix, frag_length_range = load_txt_to_matrix_with_meta(statistic_matrix,
                                                                                     loading_logger=LOGGER)
        if frag_length_range.start != min_frag_len:
            log(message="Mismatching fragment length range starts: the observed attributes matrix stated a different "
                        f"fragment length start ({frag_length_range.start}bp) than the internal program logic "
                        f"({min_frag_len} bp). Will use the loaded value.",
                log_level=logging.WARNING, logger_name=LOGGER)
            min_frag_len = frag_length_range.start
    else:
        log(message="Cannot use provided attribute 'statistic_matrix'. Must be either a path to a saves matrix file or "
                    "the matrix as np.array itself.",
            log_level=logging.CRITICAL, close_handlers=True, logger_name=LOGGER)
        sys.exit(2)
    bad_interval = False
    n_fragment_length_rows = observed_attributes_matrix.shape[0]
    s_gc_content_columns = observed_attributes_matrix.shape[1]
    fragments = np.sum(observed_attributes_matrix, axis=1, dtype=int).flatten()
    # create averaged matrix for in-place manipulation
    simulated_attributes_matrix = np.zeros((n_fragment_length_rows, s_gc_content_columns), dtype=float)
    # create raw matrices for in-place manipulation
    raw_simulated_matrices = []
    for m_idx in range(simulation_repetitions):
        raw_simulated_matrices.append(np.zeros((n_fragment_length_rows, s_gc_content_columns), dtype=float))
    random_number_generator = np.random.default_rng(seed=random_seed)  # use random seed for reproducibility here!
    if save_individual_matrices:  # check and create before iterating
        Path(tmp_dir).mkdir(parents=True, exist_ok=True)
    interval_str = f"{chromosome}_{start_coord}-{stop_coord}"
    for sim_iter_idx in range(simulation_repetitions):
        for length_index, amount_fragments in enumerate(fragments):  # simulate each fragment length separately
            if not amount_fragments:  # no fragments observed for current frag length (observed_attributes[f_len] = 0)
                continue
            actual_fragment_length = length_index + min_frag_len
            # create 2 random values for each fragment
            rand_ints = random_number_generator.integers(low=0,
                                                         high=stop_coord - start_coord - actual_fragment_length + 1,
                                                         size=amount_fragments)
            rand_ints.sort()  # for sweeping-like file pointer placement (if that even matters for twobitreader)
            sampling_failure_threshold = max(int(amount_fragments / 3), 55 if strict_n_ref_bases_handling else 33)
            # drawing fragment must fail at least 55 times (default is strict handling; 33 times otherwise) or
            # one third of all required draws, whichever is higher
            unsorted_randoms_index = 0  # index of new fall-back random positions for each fragment length
            unsorted_randoms = random_number_generator.integers(
                low=0, high=stop_coord - start_coord - actual_fragment_length + 1,
                size=sampling_failure_threshold + 1)  # could be rare fragments -> threshold +1 as minimum
            # 1/3 of target fragment number - whichever is higher - for interval to be marked for exclusion.
            # Approach if not strict_n_ref_bases_handling: linear extrapolation
            # subtract the N bases from the fragment length, compute GC content of the analyzable portion and
            # extrapolate the number of GC-bases according to the fragment's actual length.
            # Still not perfect but the best we can do for a fast algorithm
            # (implemented in 'safe_gc_base_count_inference_thresh()')
            if strict_n_ref_bases_handling:  # faster
                gc_count_iterator = map(lambda q: gc_count_rejecting_n_containing(f_seq=q),
                                        map(lambda s: ref_seq_interval_slice[s:s + actual_fragment_length].upper(),
                                            rand_ints))
            else:
                gc_count_iterator = map(lambda q: safe_gc_base_count_inference_thresh(f_seq=q[0], f_len=q[1],
                                                                                      threshold=frag_n_cont_thresh),
                                        map(lambda s: (ref_seq_interval_slice[s:s + actual_fragment_length].upper(),
                                                       actual_fragment_length), rand_ints))
            for gc_count in gc_count_iterator:
                try:
                    simulated_attributes_matrix[length_index, gc_count] += 1  # gc_count is 99999999 if any/all bases=N
                    raw_simulated_matrices[sim_iter_idx][length_index, gc_count] += 1
                    # since we do not expect 100M bp fragments; Code fails for such samples -> designed for short reads
                except IndexError:  # all-Ns-fragment: slow routine here -> get another random position
                    backup_seq = 'N' * actual_fragment_length  # initialize for while loop
                    while backup_seq.count('N') if strict_n_ref_bases_handling else \
                            (backup_seq.count('N') / actual_fragment_length >= frag_n_cont_thresh):  # keep drawing
                        try:
                            cur_start = unsorted_randoms[unsorted_randoms_index]  # get new random fragment start
                        except IndexError:  # all random integers consumed -> bad interval! (should not occur)
                            log(message=f"Too many attempts of drawing random fragments ({sampling_failure_threshold} "
                                        f"attempts) for interval '{chromosome}:{start_coord}-{stop_coord}' were in vain. "
                                        f"Triggered for fragments of {actual_fragment_length}bp length. "
                                        f"Discarding interval ..", log_level=logging.WARNING, logger_name=LOGGER)
                            bad_interval = True  # stop redrawing once all fall-back random integers have been used up
                            break
                        unsorted_randoms_index += 1  # can max. be 33, then clause below should trigger interval marking
                        if unsorted_randoms_index >= sampling_failure_threshold:  # should always be the reason why a
                            # interval is marked as bad
                            log(message=f"Too many attempts of drawing random fragments "
                                        f"for interval '{chromosome}:{start_coord}-{stop_coord}' were in vain "
                                        f"(threshold was {sampling_failure_threshold} tries). Discarding interval ..",
                                log_level=logging.WARNING, logger_name=LOGGER)
                            bad_interval = True  # stop redrawing once the threshold has been reached
                            break
                        backup_seq = ref_seq_interval_slice[cur_start:cur_start + actual_fragment_length].upper()
                    if bad_interval:
                        break  # for loop -> bad interval encountered!
                    if strict_n_ref_bases_handling:
                        gc_count = gc_count_rejecting_n_containing(f_seq=backup_seq)
                    else:
                        gc_count = safe_gc_base_count_inference_thresh(f_seq=backup_seq, f_len=actual_fragment_length,
                                                                       threshold=frag_n_cont_thresh)
                    raw_simulated_matrices[sim_iter_idx][length_index, gc_count] += 1
            if bad_interval:
                return chromosome, start_coord, stop_coord, expected_yield
        if save_individual_matrices:  # save snapshot of cumulative simulated_attributes_matrix
            simulated_file = f"{sample_id}_cumulative_simulated_attributes_matrix_{sim_iter_idx + 1}_{interval_str}.txt.gz"
            save_matrix_to_txt(matrix=simulated_attributes_matrix, filename=simulated_file, output_dir=tmp_dir,
                               float_data_precision=float_precision, gzipped=True, max_frag_length=max_frag_len,
                               min_frag_length=min_frag_len)
    if save_individual_matrices:  # store individual raw matrices
        for raw_matrix_idx, simulated_raw_matrix in enumerate(raw_simulated_matrices):
            simulated_file = f"{sample_id}_simulated_raw_attributes_matrix_{raw_matrix_idx + 1}_{interval_str}.txt.gz"
            save_matrix_to_txt(matrix=simulated_raw_matrix, filename=simulated_file, output_dir=tmp_dir,
                               float_data_precision=float_precision, gzipped=True, max_frag_length=max_frag_len,
                               min_frag_length=min_frag_len)
    return simulated_attributes_matrix, tuple(raw_simulated_matrices)


def consolidate_results(observed_attributes_matrices_sum: np.array, simulated_attributes_matrices_sum: np.array,
                        simulated_attributes_raw_matrix_sums: List[np.array], n_ogc_summed: int, n_sims: int,
                        n_sgc_summed: int, bam_file: str, tmp_dir: Optional[str], min_frag_len: int, max_frag_len: int,
                        min_frag_occurs: int, ignored_fragments: int, sample_id: str,
                        focus_nondefault_values: Optional[int], precision=6, plot_result=False, output_all=False,
                        show_plots=False) -> Tuple[Tuple[Path, np.array], Tuple[Path, np.array], Tuple[range, range]]:
    """

    :param focus_nondefault_values: plots will focus on these values using the integer value of the parameter as border
    :param ignored_fragments:
    :param output_all:
    :param simulated_attributes_raw_matrix_sums: list of n simulated attribute matrices; used for computing correction
                                                 weights
    :param tmp_dir:
    :param min_frag_occurs:
    :param observed_attributes_matrices_sum:
    :param simulated_attributes_matrices_sum: used to store a "consolidated" version of the n simulation matrices
                                              is not used for correction weights! simulated_attributes_raw_matrix_sums
                                              is used instead to create n correction weight matrices which are then
                                              averaged
    :param n_ogc_summed:
    :param n_sims:
    :param n_sgc_summed:
    :param bam_file:
    :param min_frag_len:
    :param max_frag_len:
    :param sample_id:
    :param precision:
    :param plot_result:
    :return:
    """
    # compute average matrices
    simulated_attributes_matrix_downscaled = simulated_attributes_matrices_sum / n_sims  # genome-wide averaged counts
    save_matrix_to_txt(matrix=simulated_attributes_matrix_downscaled, float_data_precision=precision, verbose=True,
                       filename=f'{sample_id}_simulated_attributes_matrix.txt.gz', output_dir=tmp_dir, gzipped=True,
                       max_frag_length=max_frag_len, min_frag_length=min_frag_len)
    observed_attributes_matrix_sum_path = save_matrix_to_txt(
        matrix=observed_attributes_matrices_sum, output_dir=tmp_dir, gzipped=True, verbose=True,
        filename=f'{sample_id}_observed_attributes_matrix.txt.gz', float_data_precision=precision,
        report_saved_path=True, max_frag_length=max_frag_len, min_frag_length=min_frag_len)  # there is 1 row too many!
    if output_all:
        simulated_attributes_matrix_interval_average = simulated_attributes_matrices_sum / n_sgc_summed  # avg per interval
        save_matrix_to_txt(matrix=simulated_attributes_matrix_interval_average, float_data_precision=precision,
                           verbose=True,
                           filename=f'{sample_id}_simulated_attributes_matrix_interval-average.txt.gz', gzipped=True,
                           output_dir=tmp_dir, max_frag_length=max_frag_len, min_frag_length=min_frag_len)
        for sim_idx, raw_sim_mat in enumerate(simulated_attributes_raw_matrix_sums):
            save_matrix_to_txt(matrix=raw_sim_mat, float_data_precision=precision, verbose=True, output_dir=tmp_dir,
                               filename=f'{sample_id}_simulated_attributes_raw_matrix_{sim_idx}.txt.gz', gzipped=True,
                               max_frag_length=max_frag_len, min_frag_length=min_frag_len)
        observed_attributes_matrix_interval_average = observed_attributes_matrices_sum / n_ogc_summed  # average per interval
        save_matrix_to_txt(
            matrix=observed_attributes_matrix_interval_average, output_dir=tmp_dir, gzipped=True, verbose=True,
            float_data_precision=precision, filename=f'{sample_id}_observed_attributes_matrix_interval-average.txt.gz',
            max_frag_length=max_frag_len, min_frag_length=min_frag_len)
    # compute correction matrix
    if min_frag_occurs > 1:  # compute mask based on combination occurrences (f_length and GC base count)
        log(message=f"Mask computation: used a minimum threshold of {min_frag_occurs:,} for total "
                    "counts of individual GC-base count-fragment length combinations.",
            log_level=logging.INFO, logger_name=LOGGER)
    else:
        log(message=f"Masking GC-base-count/fragment-length combinations that were not observed..",
            log_level=logging.INFO, logger_name=LOGGER)
    masks_for_sims = []
    for raw_sim_mat in simulated_attributes_raw_matrix_sums:
        masks_for_sims.append((observed_attributes_matrices_sum >= min_frag_occurs) *  # min_frag_occurs
                              (raw_sim_mat >= min_frag_occurs))  # min_frag_occurs
    if output_all:
        for sim_idx, raw_mask in enumerate(masks_for_sims):
            save_matrix_to_txt(matrix=raw_mask, output_dir=tmp_dir, float_data_precision=precision, gzipped=True,
                               filename=f'{sample_id}_gc_bias_computation_mask_simulation{sim_idx}.txt.gz',
                               max_frag_length=max_frag_len, min_frag_length=min_frag_len)
    # compute correction weights for each simulation iteration separately
    # (high values in high frequency attribute combinations in S_gc due to stochasticity should be reduced)
    correction_weights_matrices = []
    for sim_idx in range(n_sims):
        correction_weights_matrices.append(np.ones(observed_attributes_matrices_sum.shape, dtype=np.float64))
        correction_weights_matrices[sim_idx][masks_for_sims[sim_idx]] = np.divide(
            simulated_attributes_raw_matrix_sums[sim_idx], observed_attributes_matrices_sum,
            where=masks_for_sims[sim_idx])[masks_for_sims[sim_idx]]
    # compute average weights matrix over simulation repetitions
    correction_weights_matrix_average = sum(correction_weights_matrices) / n_sims
    # mask nonsense positions in W_gc matrix with 0:
    for f_len in range(min_frag_len, max_frag_len + 1, 1):
        correction_weights_matrix_average[f_len - min_frag_len, f_len + 1:] = 0.
    # ALWAYS SAVE PRIMARY RESULT!
    correction_weights_matrix_path = save_matrix_to_txt(matrix=correction_weights_matrix_average, output_dir=tmp_dir,
                                                        float_data_precision=precision, max_frag_length=max_frag_len,
                                                        filename=f'{sample_id}_gc_weights_{n_sims}simsMean.txt.gz',
                                                        gzipped=True, min_frag_length=min_frag_len,
                                                        report_saved_path=True)
    # compute corrected values mask
    mask_from_average = (correction_weights_matrix_average != 0.) * (correction_weights_matrix_average != 1.)
    if output_all:
        save_matrix_to_txt(matrix=mask_from_average, float_data_precision=precision, gzipped=True,
                           filename=f'{sample_id}_gc_bias_computation_from_average_mask.txt.gz',
                           output_dir=tmp_dir, max_frag_length=max_frag_len, min_frag_length=min_frag_len)
    # compute total mask from OR-linked simulation masks
    complete_mask = np.zeros(mask_from_average.shape).astype(bool)
    for sim_mask in masks_for_sims:
        complete_mask = complete_mask | sim_mask  # or-link all non-default values
    complete_mask_path = save_matrix_to_txt(matrix=complete_mask, output_dir=tmp_dir, float_data_precision=precision,
                                            filename=f'{sample_id}_gc_bias_computation_mask.txt.gz',
                                            gzipped=True, max_frag_length=max_frag_len, min_frag_length=min_frag_len,
                                            report_saved_path=True)
    # plot matrices if specified by user
    deleted_rows, deleted_columns = range(0), range(0)
    if plot_result:
        # plot fragment length distribution
        plot_fragment_length_dists(matrix_data_frame=None, matrix_file_list=[observed_attributes_matrix_sum_path],
                                   out_dir_path=Path(tmp_dir), normalize_to_dataset_size=True, show_figure=False,
                                   strip_xaxis_end_zeros=True, parent_logger=LOGGER, sample_id=sample_id)
        if focus_nondefault_values is not None:  # create focused plots
            complete_mask_focused, (deleted_rows, deleted_columns) = reduce_matrix(
                matrix_to_trim=complete_mask, trim_dimensions_exclusively_containing=[False],
                border_elements=focus_nondefault_values)
            correction_weights_matrix_average_focused = trim_2d_matrix(matrix=correction_weights_matrix_average,
                                                                       rows=(deleted_rows.start, deleted_rows.stop),
                                                                       columns=(deleted_columns.start,
                                                                                deleted_columns.stop))
            observed_attributes_matrices_sum_focused = trim_2d_matrix(matrix=observed_attributes_matrices_sum,
                                                                      rows=(deleted_rows.start, deleted_rows.stop),
                                                                      columns=(deleted_columns.start,
                                                                               deleted_columns.stop))
            simulated_attributes_matrix_downscaled_focused = trim_2d_matrix(
                matrix=simulated_attributes_matrix_downscaled, rows=(deleted_rows.start, deleted_rows.stop),
                columns=(deleted_columns.start, deleted_columns.stop))
            frq_data = {'W_gc': pd.DataFrame(correction_weights_matrix_average_focused),
                        'O_gc': pd.DataFrame(observed_attributes_matrices_sum_focused),
                        'S_gc': pd.DataFrame(simulated_attributes_matrix_downscaled_focused),
                        'Mask': pd.DataFrame(complete_mask_focused)}
        else:
            frq_data = {'W_gc': pd.DataFrame(correction_weights_matrix_average),
                        'O_gc': pd.DataFrame(observed_attributes_matrices_sum),
                        'S_gc': pd.DataFrame(simulated_attributes_matrix_downscaled),
                        'Mask': pd.DataFrame(complete_mask)}
        if output_all:
            for s_idx, sim_mat in enumerate(simulated_attributes_raw_matrix_sums):
                if focus_nondefault_values is not None:
                    sim_mat_focused = trim_2d_matrix(matrix=sim_mat, rows=(deleted_rows.start, deleted_rows.stop),
                                                     columns=(deleted_columns.start, deleted_columns.stop))
                    frq_data.update({f'S_gc_simulation_{s_idx}': pd.DataFrame(sim_mat_focused)})
                    mask_sim_focused = trim_2d_matrix(matrix=masks_for_sims[s_idx],
                                                      rows=(deleted_rows.start, deleted_rows.stop),
                                                      columns=(deleted_columns.start, deleted_columns.stop))
                    frq_data.update({f'Mask_simulation_{s_idx}': pd.DataFrame(mask_sim_focused)})
                    # plot also difference between observed and expected
                    diff_mat_focused = trim_2d_matrix(matrix=sim_mat - observed_attributes_matrices_sum,
                                                      rows=(deleted_rows.start, deleted_rows.stop),
                                                      columns=(deleted_columns.start, deleted_columns.stop))
                    frq_data.update({f'D_gc_simulation_{s_idx}': pd.DataFrame(diff_mat_focused)})
                else:
                    frq_data.update({f'S_gc_simulation_{s_idx}': pd.DataFrame(sim_mat)})
                    frq_data.update({f'Mask_simulation_{s_idx}': pd.DataFrame(masks_for_sims[s_idx])})
                    # plot also absolute error between observed and expected
                    frq_data.update(
                        {f'D_gc_simulation_{s_idx}': pd.DataFrame(sim_mat - observed_attributes_matrices_sum)})
        for data_category in frq_data.keys():
            plot_statistic_matrices(frq_data=frq_data, data_id_to_show=data_category,
                                    y_tick_label_offset=deleted_rows.start + min_frag_len,
                                    x_tick_label_offset=deleted_columns.start, show_figure=show_plots,
                                    in_file=bam_file, output_dir=tmp_dir, sample_id=sample_id, fig_width=1800,
                                    fig_height=2000, fig_fontsize=50, parent_logger=LOGGER)
    try:  # give feedback about estimated percentage of corrected fragments:
        included_fragments = int(observed_attributes_matrices_sum.sum())  # without ignored fragments!
        r_fragments_corrected = observed_attributes_matrices_sum[complete_mask].sum() / included_fragments
        log(message=f"Estimated percentage of weighted fragments in dataset: {r_fragments_corrected:.4%} "
                    f"(based on {included_fragments:,} processed fragments included in statistics, which is "
                    f"{included_fragments / (included_fragments + ignored_fragments):.2%} of all processed fragments)",
            log_level=logging.INFO, logger_name=LOGGER)
    except ZeroDivisionError:  # this should never happen
        log(message=f"Unable to estimate weighted dataset fraction!",
            log_level=logging.WARNING, logger_name=LOGGER)
    return (correction_weights_matrix_path, correction_weights_matrix_average), \
        (complete_mask_path, complete_mask), \
        (deleted_rows, deleted_columns)


def sort_intervals_by_exclusion_ist_overlap(all_intervals_with_score: GenomicIntervalList,
                                            expected_dataset_fragments: int,
                                            bad_intervals: Optional[BadIntervalsDict], remove_bad_intervals=True,
                                            max_overlap_percentage=DEFAULT_MAX_INTERVAL_PERCENTAGE_EXCLUSIONLIST_OVERLAP) \
        -> List[Tuple[str, int, int]]:
    # sort ascending overlapping bases % (normalized to interval length accounting for possible size differences)
    intervals_passing_filters = sorted(all_intervals_with_score,
                                       key=lambda c: c[3] / (c[2] - c[1]),  # sort: ascending overlap -> exclusion list
                                       reverse=False)
    intervals_passing_filters = list(map(
        lambda t: t[:3],  # discard the exclusion list overlap for further processing
        filter(lambda c: max_overlap_percentage >= c[3] / (c[2] - c[1]) * 100.,
               intervals_passing_filters)))
    log(message=f"{len(all_intervals_with_score) - len(intervals_passing_filters):,} intervals were excluded from "
                f"further analysis based on the {max_overlap_percentage:.1f}% interval overlap with exclusion marked "
                "regions threshold.", log_level=logging.INFO, logger_name=LOGGER)
    if remove_bad_intervals and bad_intervals is not None:  # remove bad intervals from library if any are defined
        # TODO: UNTESTED CLAUSE!
        # -> select next higher or equal target fragment value
        bad_intervals_for_check = []
        pre_bc_removal_intervals = len(intervals_passing_filters)
        for (chrom, start, stop), sample_frag_yield_list in bad_intervals.items():
            if max(sample_frag_yield_list) >= expected_dataset_fragments:
                # ignore bad intervals which were recorded for samples with HIGHER DoC than the current one!
                bad_intervals_for_check.append((chrom, start, stop))
        intervals_passing_filters = list(filter(lambda c: c not in bad_intervals_for_check,
                                                intervals_passing_filters))
        intervals_after_library_exclusion = len(intervals_passing_filters)
        if pre_bc_removal_intervals != intervals_after_library_exclusion:
            log(message=f"{pre_bc_removal_intervals - intervals_after_library_exclusion:,} intervals were excluded "
                        "from further analysis based on the bad intervals library file. Bad intervals are "
                        "automatically replaced by other intervals.", log_level=logging.INFO, logger_name=LOGGER)
    return intervals_passing_filters


def gc_bias_worker(weighted_intervals_to_process: List[Tuple[float, Tuple[str, int, int]]], n_sims: int,
                   sender: mp_connection.Connection, min_frag_len: int, max_frag_len: int,
                   two_bit_genome_file: OneOf[str, Path],
                   chromosome_sizes: Dict[str, int], input_bam: str, sample_id: str, expected_yield: int, precision=6,
                   min_unclipped_aln_fraction=DEFAULT_MIN_UNCLIPPED_ALN_FRACTION, random_seed=RANDOM_SEED, tmp_dir=None,
                   keep_interval_data=False, strict_n_base_exclusion=True):
    """

    :param strict_n_base_exclusion:
    :param sample_id:
    :param expected_yield:
    :param tmp_dir:
    :param two_bit_genome_file:
    :param chromosome_sizes:
    :param input_bam:
    :param weighted_intervals_to_process:
    :param n_sims:
    :param sender:
    :param min_frag_len:
    :param max_frag_len:
    :param precision:
    :param random_seed:
    :param keep_interval_data:
    :param min_unclipped_aln_fraction:
    :return:
    """
    n_discarded_intervals = 0
    n_observed_gc_matrices = 0
    n_expected_gc_matrices = 0
    observed_attributes_cumulated_matrix = np.zeros((max_frag_len - min_frag_len + 1, max_frag_len + 1),
                                                    dtype=float)  # interval is both sides inclusive
    simulated_attributes_cumulated_matrix = np.zeros((max_frag_len - min_frag_len + 1, max_frag_len + 1),
                                                     dtype=float)  # interval is both sides inclusive
    simulated_attributes_cumulated_raw_matrices = []
    for mat_idx in range(n_sims):
        simulated_attributes_cumulated_raw_matrices.append(np.zeros((max_frag_len - min_frag_len + 1, max_frag_len + 1),
                                                                    dtype=float))
    n_fragments_processed = 0
    n_fragments_ignored = 0
    individual_matrices_tmp_dir = str(Path(tmp_dir) / 'data_per_interval') if keep_interval_data else None
    bad_intervals_list = []
    # process all received intervals (combine using weight!)
    for region_weight, (interval_chrom, interval_start, interval_end) in weighted_intervals_to_process:
        check_region_within_genomic_bounds(chrom=interval_chrom, start=interval_start, stop=interval_end,
                                           chromosome_sizes=chromosome_sizes, raise_error=True)
        # COMPUTE O_GC
        gc_matrix_observed, interval_fragments_processed, frag_ignored = compute_observed_attributes_matrix(
            two_bit_reference_path=two_bit_genome_file, tmp_dir=individual_matrices_tmp_dir, bam_file=input_bam,
            start_coord=interval_start, stop_coord=interval_end, save_individual_matrices=keep_interval_data,
            sample_id=sample_id, float_precision=precision, min_frag_len=min_frag_len, chromosome=interval_chrom,
            max_frag_len=max_frag_len, strict_n_ref_bases_handling=strict_n_base_exclusion,
            min_unclipped_aln_fraction=min_unclipped_aln_fraction)  # returns unweighted matrices
        # COMPUTE S_GC; in-place manipulate the expected matrix
        simulated_attributes_matrix, raw_simulated_matrices = simulate_fragment_attributes(
            two_bit_reference_path=two_bit_genome_file, tmp_dir=individual_matrices_tmp_dir, min_frag_len=min_frag_len,
            statistic_matrix=gc_matrix_observed, chromosome=interval_chrom, sample_id=sample_id,
            start_coord=interval_start,
            stop_coord=interval_end, simulation_repetitions=n_sims, save_individual_matrices=keep_interval_data,
            random_seed=random_seed, float_precision=precision, expected_yield=expected_yield,
            strict_n_ref_bases_handling=strict_n_base_exclusion, max_frag_len=max_frag_len)
        # check if we've got a bad interval
        if isinstance(simulated_attributes_matrix, tuple) and len(simulated_attributes_matrix) == 4:
            # if so, discard entire interval info!
            # (reduces the accuracy of the reconstruction but what can you do. It is what it is.)
            bad_interval = (region_weight, simulated_attributes_matrix)  # in this case, simulated_attributes_matrix is
            # a tuple: (chromosome, start_coord, stop_coord, expected_yield)
            bad_intervals_list.append(bad_interval)
            n_discarded_intervals += 1
            continue
        observed_attributes_cumulated_matrix += (gc_matrix_observed * region_weight)
        n_fragments_processed += interval_fragments_processed - frag_ignored  # use actual fragment count?
        n_fragments_ignored += frag_ignored  # use actual fragment count?
        n_observed_gc_matrices += 1
        simulated_attributes_cumulated_matrix += (simulated_attributes_matrix * region_weight)
        for mat_idx in range(n_sims):
            simulated_attributes_cumulated_raw_matrices[mat_idx] += (raw_simulated_matrices[mat_idx] * region_weight)
        n_expected_gc_matrices += n_sims
    sender.send((observed_attributes_cumulated_matrix, simulated_attributes_cumulated_matrix,
                 simulated_attributes_cumulated_raw_matrices, n_observed_gc_matrices,
                 n_expected_gc_matrices, n_fragments_processed, n_fragments_ignored, n_discarded_intervals,
                 bad_intervals_list))
    sender.close()


def compute_lin_comb_reconstructed_gc_distribution_error(
        bad_region_labels: List[str], processed_regions_weighted: List[Tuple[float, Tuple[str, int, int]]],
        reference_distribution: np.array, all_fgcds: Dict[str, np.array]) \
        -> Tuple[np.array, np.array, np.array, float]:
    # insanity check
    assert 0.999 < reference_distribution.sum() < 1.001
    # used_weights = [reg_weight if create_region_label(chrm=hrm, start=strt, end=nd) not in bad_region_labels else 0
    #                 for reg_weight, (hrm, strt, nd) in processed_regions_weighted]
    weighted_fgcds = np.array([(all_fgcds[create_region_label(chrm=hrm, start=strt, end=nd)] * reg_weight)
                               if create_region_label(chrm=hrm, start=strt, end=nd) not in bad_region_labels else
                               np.zeros(len(all_fgcds[create_region_label(chrm=hrm, start=strt, end=nd)]))
                               for reg_weight, (hrm, strt, nd) in processed_regions_weighted], dtype=float)
    original_fgcd_dist = np.array([(all_fgcds[create_region_label(chrm=hrm, start=strt, end=nd)])
                               if create_region_label(chrm=hrm, start=strt, end=nd) not in bad_region_labels else
                               np.zeros(len(all_fgcds[create_region_label(chrm=hrm, start=strt, end=nd)]))
                               for reg_weight, (hrm, strt, nd) in processed_regions_weighted], dtype=float)
    original_fgcd_dist = original_fgcd_dist.sum(axis=0) / len(processed_regions_weighted)
    # assert 0.999 < original_fgcd_dist.sum() < 1.001
    reconstructed_ref_fgcd = weighted_fgcds.sum(axis=0)  # column-wise sum over all rows
    # assert 0.999 < reconstructed_ref_fgcd.sum() < 1.001
    # assert len(reconstructed_ref_fgcd) == len(reference_distribution)
    residual_distribution = reconstructed_ref_fgcd - reference_distribution
    # assert np.abs(residual_distribution).sum() < 1.0
    absolute_reconstruction_error_pc = np.abs(residual_distribution).sum()
    # from aes(): np.abs((components.T * weights).sum(axis=1) - target).sum()
    return original_fgcd_dist, reconstructed_ref_fgcd, residual_distribution, absolute_reconstruction_error_pc


def compute_gc_bias_parallel(weighted_intervals_to_process: List[Tuple[float, Tuple[str, int, int]]], threads: int,
                             simulation_count: int,
                             min_flen: int, max_flen: int, out_dir_sample: str, in_bam: str, sample_name: str,
                             two_bit_reference_file: OneOf[str, Path], chrom_sizes: Dict[str, int],
                             min_frag_occurs: int,
                             target_fragments_processed: int, expected_yield: int, tmp_dir_sample: str,
                             bad_intervals_library_file: Optional[str], plot_focus_border: Optional[int],
                             reference_gc_content_distribution: np.array,
                             interval_gc_content_distributions: Dict[str, np.array],
                             float_precision=6, strict_n_base_exclusion=True, keep_interval_data=False,
                             visualize_matrices=False, output_all=False, write_updated_bad_intervals_library=True,
                             detect_outliers=True, focus_custom_values=True,
                             outlier_detection_method='IQR', outlier_detection_stringency=2, smooth_weights=True,
                             smoothing_kernel='gauss', smoothing_intensity=2, show_plots=False,
                             min_unclipped_aln_fraction=DEFAULT_MIN_UNCLIPPED_ALN_FRACTION, random_seed=RANDOM_SEED) \
        -> Tuple[np.array, np.array]:
    """
    :param interval_gc_content_distributions:
    :param reference_gc_content_distribution:
    :param show_plots:
    :param focus_custom_values:
    :param plot_focus_border:
    :param strict_n_base_exclusion:
    :param smooth_weights:
    :param smoothing_kernel:
    :param smoothing_intensity:
    :param outlier_detection_stringency:
    :param outlier_detection_method:
    :param detect_outliers:
    :param output_all:
    :param bad_intervals_library_file:
    :param write_updated_bad_intervals_library:
    :param expected_yield:
    :param tmp_dir_sample: primary output directory used by this function. Results stored here will be moved later.
    :param min_frag_occurs:
    :param target_fragments_processed:
    :param chrom_sizes:
    :param two_bit_reference_file:
    :param weighted_intervals_to_process:
    :param threads:
    :param simulation_count:
    :param min_flen:
    :param max_flen:
    :param out_dir_sample:
    :param in_bam:
    :param sample_name:
    :param float_precision:
    :param keep_interval_data:
    :param random_seed:
    :param visualize_matrices:
    :param min_unclipped_aln_fraction: minimum fraction of unclipped alignment length for fragment to be counted in O_gc
    :return:
    """
    n_processes = int(threads)
    # -> multiprocessing gives better performance than threading
    # split up sorted intervals in lists for workers
    intervals_for_workers = [[] for _i in range(n_processes)]
    # sort intervals according to chromosome
    sorted_weighted_intervals_to_process = humansorted(weighted_intervals_to_process,
                                                       key=lambda t: (t[1][0], t[1][1]))  # sort by chrom., then start
    n_intervals = len(sorted_weighted_intervals_to_process)
    # create consecutive intervals on as few chromosomes as much as possible for all workers
    # (speedup due to BAM file pointer requiring less jumping?)
    intervals_per_process = math.ceil(len(sorted_weighted_intervals_to_process) / n_processes)
    for split_idx in range(n_processes):
        if split_idx < n_processes - 1:  # not last split
            intervals_for_workers[split_idx] = sorted_weighted_intervals_to_process[
                                               intervals_per_process * split_idx:
                                               intervals_per_process * (split_idx + 1)]
        else:
            intervals_for_workers[split_idx] = sorted_weighted_intervals_to_process[
                                               intervals_per_process * split_idx:]
    # create connections, workers, and processes
    receivers = []
    senders = []
    communication_lines = [mp.Pipe(duplex=False) for _i in range(n_processes)]
    for p_idx in range(n_processes):
        receivers.append(communication_lines[p_idx][0])
        senders.append(communication_lines[p_idx][1])
    # check sample tmp dir
    tmp_sample_output_path = Path(tmp_dir_sample)
    if tmp_sample_output_path.exists():
        log(message=f"Temporary output directory for sample '{sample_name}' exists and will be deleted: Deleting "
                    f"'{tmp_dir_sample}' ..", log_level=logging.WARNING, logger_name=LOGGER)
        shutil.rmtree(tmp_sample_output_path)
    tmp_sample_output_path.mkdir(parents=True, exist_ok=True)
    # create shared value and lock for counting processed fragments
    all_worker_kwargs = [{'weighted_intervals_to_process': intervals_for_workers[w_idx], 'n_sims': simulation_count,
                          'sender': senders[w_idx], 'min_frag_len': min_flen, 'max_frag_len': max_flen,
                          'keep_interval_data': keep_interval_data, 'random_seed': random_seed, 'input_bam': in_bam,
                          'precision': float_precision, 'two_bit_genome_file': two_bit_reference_file,
                          'chromosome_sizes': chrom_sizes,
                          'tmp_dir': tmp_dir_sample, 'expected_yield': expected_yield, 'sample_id': sample_name,
                          'strict_n_base_exclusion': strict_n_base_exclusion,
                          'min_unclipped_aln_fraction': min_unclipped_aln_fraction}
                         # write to TMP first, then move files which should be kept (is much faster on cluster!)
                         for w_idx in range(n_processes)]
    gc_bias_workers = [mp.Process(target=gc_bias_worker, kwargs=worker_kwargs)
                       for worker_kwargs in all_worker_kwargs]
    # start all workers, receive results, wait until finished, and close them
    observed_attributes_matrices_sum = np.zeros((max_flen - min_flen + 1, max_flen + 1), dtype=float)
    simulated_attributes_matrices_sum = np.zeros((max_flen - min_flen + 1, max_flen + 1), dtype=float)
    simulated_attributes_raw_matrix_sums = []
    for _s_idx in range(simulation_count):
        simulated_attributes_raw_matrix_sums.append(np.zeros((max_flen - min_flen + 1, max_flen + 1), dtype=float))
    n_summed_observed_attributes_matrices = 0
    n_summed_s_gc_matrices = 0
    total_processed_fragments = 0
    total_ignored_fragments = 0
    total_discarded_intervals = 0
    all_bad_intervals = {}
    bad_interval_weights = []
    deque(map(lambda w: w.start(), gc_bias_workers), maxlen=0)
    # if memory profiler is applied to code and processes are spawned rather than forked, we get a PicklingError her:
    # "Can't pickle <function gc_bias_worker at 0xsomething>: attribute lookup gc_bias_worker on __main__ failed"
    received_data = [incoming_result.recv() for incoming_result in receivers]  # collect self-terminating worker returns
    for o_gc_sum_mat, s_gc_sum_mat, s_gc_sum_raw_mats, n_observed, n_expected, n_fragments_processed, \
            n_fragments_ignored, n_discarded_intervals, list_of_bad_intervals in received_data:
        # received:
        # (observed_attributes_cumulated_matrix, simulated_attributes_cumulated_matrix,
        #  simulated_attributes_cumulated_raw_matrices, n_observed_gc_matrices, n_expected_gc_matrices,
        #  n_fragments_processed, n_fragments_ignored, n_discarded_intervals, bad_intervals_list)
        observed_attributes_matrices_sum += o_gc_sum_mat  # inplace np.array manipulation
        simulated_attributes_matrices_sum += s_gc_sum_mat  # inplace np.array manipulation; interval weights already
        # included. Note for computing the proper average, weights in discarded_intervals (List[Tuple[float, interval]])
        # must be subtracted from the denominator before the division operation in averaging!
        # The matrices can simply be summed below because their intervals have been scaled already by
        # their corresponding weights.
        for mat_idx, s_gc_sum_raw_mat in enumerate(s_gc_sum_raw_mats):
            simulated_attributes_raw_matrix_sums[mat_idx] += s_gc_sum_raw_mat
        # UPDATE COUNT STATS
        n_summed_observed_attributes_matrices += n_observed
        n_summed_s_gc_matrices += n_expected
        total_processed_fragments += n_fragments_processed
        total_ignored_fragments += n_fragments_ignored
        total_discarded_intervals += n_discarded_intervals
        # process bad regions information
        for bad_region_weight, (chrm, strt, stp, rst) in list_of_bad_intervals:
            bad_interval_weights.append(bad_region_weight)
            if all_bad_intervals.get((chrm, strt, stp)) is None:
                all_bad_intervals.update({(chrm, strt, stp): [int(rst[1])]})
            else:  # insanity check: interval locus already exists (= multiple entries! Not expected)
                raise KeyError(f"the interval {(chrm, strt, stp)} has already been processed and should ot be present "
                               f"during accumulation of consolidated results!")
                # all_bad_intervals[(chrm, strt, stp)].append(int(rst[1]))
    actually_processed_intervals = n_intervals - len(list_of_bad_intervals)
    # divide by applied region weights to finalize weighted mean
    bad_interval_labels = list(map(lambda t: create_region_label(chrm=t[0], start=t[1], end=t[2]),
                                   all_bad_intervals.keys()))
    # used_weights_sum = sum([0 if create_region_label(chrm=hrm, start=strt, end=nd) in bad_interval_labels else reg_wght
    #                         for reg_wght, (hrm, strt, nd) in sorted_weighted_intervals_to_process])
    observed_attributes_matrices_sum *= actually_processed_intervals
    simulated_attributes_matrices_sum *= actually_processed_intervals
    for mat_idx in range(len(simulated_attributes_raw_matrix_sums)):
        simulated_attributes_raw_matrix_sums[mat_idx] *= actually_processed_intervals
    # create user feedback
    log(message="---------------------------------------------------------------------------------",
        log_level=logging.INFO, logger_name=LOGGER)
    if total_processed_fragments - total_ignored_fragments < target_fragments_processed:
        if total_processed_fragments == 0:
            log(message="No fragments were processed. Check your environment for samtools and pysam! Terminating ..",
                log_level=logging.CRITICAL, close_handlers=True, logger_name=LOGGER)
            sys.exit(3)
        log(message="Could not reach target count of processed fragments: "
                    f"{total_processed_fragments-total_ignored_fragments:,}/{target_fragments_processed:,} fragments "
                    "included in statistics." + (" (check W_gc plot for safety)" if visualize_matrices else ''),
            log_level=logging.WARNING, logger_name=LOGGER)
    else:
        log(message=f"Target number of fragments to process reached: "
                    f"{total_processed_fragments - total_ignored_fragments:,}/{target_fragments_processed:,} "
                    f"(incl. ignored fragments).", log_level=logging.INFO, logger_name=LOGGER)
    if total_discarded_intervals:
        log(message=f"Discarded some intervals that contained too many fragments with N-bases: "
                    f"{total_discarded_intervals:,} (NOTE: this reduces the accuracy of the approximation of the "
                    f"reference genome fragment GC content by a linear combination of the processed intervals)",
            log_level=logging.INFO, logger_name=LOGGER)
    log(message=f"In total, {n_summed_observed_attributes_matrices:,} "
                "genomic intervals were included in GC-bias computation.", log_level=logging.INFO, logger_name=LOGGER)
    if total_ignored_fragments:
        log(message=f"Total number of all ignored fragments due to length out of bounds, extensive clipping, or "
                    f"excessive reference base N-content: "
                    f"{total_ignored_fragments:,} (= {total_ignored_fragments / total_processed_fragments:.2%} "
                    f"processed alignments ignored)", log_level=logging.INFO, logger_name=LOGGER)
        log(message=f"Fragments included in statistics: "
                    f"{total_processed_fragments - total_ignored_fragments:,}", log_level=logging.INFO,
            logger_name=LOGGER)
        if total_ignored_fragments / total_processed_fragments > 0.05:
            log(message="More than 5% of all processed fragments with length out of bounds, were extensively clipped, "
                        "or had extensive reference base N-content! You may want to re-run this analysis with a higher "
                        "maximum fragment length!", log_level=logging.WARNING, logger_name=LOGGER)
    else:
        log(message=f"No alignments were rejected.", log_level=logging.INFO, logger_name=LOGGER)
    # give feedback about the reconstruction accuracy of the distribution of the reference genome fragment GC content
    original_fgcd, reconstructed_fgcd, residual_distribution, absolute_error_pc = \
        compute_lin_comb_reconstructed_gc_distribution_error(
            bad_region_labels=bad_interval_labels, processed_regions_weighted=sorted_weighted_intervals_to_process,
            reference_distribution=reference_gc_content_distribution, all_fgcds=interval_gc_content_distributions)
    visualize_reconstruction_result(out_dir=tmp_dir_sample, show=show_plots, original_dist=original_fgcd,
                                    target_dist=reference_gc_content_distribution, abs_err=absolute_error_pc,
                                    reconstructed_dist=reconstructed_fgcd, residual_dists=residual_distribution,
                                    sample_id=sample_name, reduced_bins=True)
    log(message="Actual cumulative reconstruction error of the fragment GC content distribution of the reference "
                f"genome: {absolute_error_pc:.2%}", log_level=logging.INFO, logger_name=LOGGER)
    _ = deque(map(lambda w: w.join(), gc_bias_workers), maxlen=0)  # wait until all processes have finished
    _ = deque(map(lambda w: w.close(), gc_bias_workers), maxlen=0)
    # update bad intervals library
    if write_updated_bad_intervals_library:
        if not all_bad_intervals:
            log(message="No bad intervals encountered. No new bad intervals library will be output.",
                log_level=logging.INFO, logger_name=LOGGER)
        else:
            # re-read bad intervals (might have changed in the meantime by other processes;
            # include also bad intervals form lower-expected fragment count datasets in new library!)
            if bad_intervals_library_file is not None:
                curren_bad_intervals = read_bad_genomic_intervals_bed_file(bed_path=bad_intervals_library_file)
                for k, yield_list in curren_bad_intervals.items():
                    if all_bad_intervals.get(k) is None:
                        all_bad_intervals.update({k: yield_list})
                    else:
                        all_bad_intervals[k].extend(yield_list)  # append received sample yield list to existing list
            newline = '\n'
            # create lines for (possibly) combined bad intervals
            new_lib_buffer = [f"{newline if line_idx else ''}{chrm}\t{strt}\t{stp}\t{max(yld_lst)}"
                              for line_idx, ((chrm, strt, stp), yld_lst) in
                              enumerate(humansorted(all_bad_intervals.items(),
                                                    key=lambda t: (t[0][0], t[0][1]),
                                                    reverse=False))]
            new_lib_path = Path(out_dir_sample) / \
                f'bad_intervals_{time.strftime(TIMESTAMP_FORMAT, time.localtime())}.bed'
            with AtomicOpen(new_lib_path, 'wt') as f_new_bclib:
                f_new_bclib.writelines(new_lib_buffer)
            log(message=f"Updated bad intervals library BED file written to: '{new_lib_path}'",
                log_level=logging.INFO, logger_name=LOGGER)
    # combine results across simulations using specific binary masks for each iteration; output consensus S_gc and mask
    (use_correction_matrix_path, correction_matrix), \
        (weights_mask_path, weights_mask), \
        trimmed_dimensions = consolidate_results(
        observed_attributes_matrices_sum=observed_attributes_matrices_sum, min_frag_len=min_flen, bam_file=in_bam,
        simulated_attributes_matrices_sum=simulated_attributes_matrices_sum, max_frag_len=max_flen,
        simulated_attributes_raw_matrix_sums=simulated_attributes_raw_matrix_sums, plot_result=visualize_matrices,
        tmp_dir=tmp_dir_sample, n_sgc_summed=n_summed_s_gc_matrices, n_ogc_summed=n_summed_observed_attributes_matrices,
        n_sims=simulation_count, sample_id=sample_name, min_frag_occurs=min_frag_occurs, precision=float_precision,
        ignored_fragments=total_ignored_fragments, output_all=output_all, focus_nondefault_values=plot_focus_border)
    use_correction_matrix = correction_matrix
    # detect and set outliers if defined:
    postprocessing_data_id = 'W_gc'
    if detect_outliers:
        postprocessing_data_id = f'{postprocessing_data_id}_outliers_removed'
        use_correction_matrix = limit_extreme_outliers(outliers_matrix=use_correction_matrix,
                                                       outliers_factor=10 - outlier_detection_stringency,
                                                       detection_method=outlier_detection_method, parent_logger=LOGGER)
        use_correction_matrix_path = save_matrix_to_txt(  # always save non-focused matrices
            matrix=use_correction_matrix, output_dir=tmp_dir_sample,
            filename=str(Path('.'.join(use_correction_matrix_path.name.split('.')[:-2] +
                                       [f'{outlier_detection_stringency}{outlier_detection_method.upper()}'
                                        f'outliersRemoved', 'txt', 'gz']))),
            gzipped=True, float_data_precision=float_precision, report_saved_path=True, max_frag_length=max_flen,
            min_frag_length=min_flen)
        if visualize_matrices:
            visualize_correction_matrix = use_correction_matrix
            if focus_custom_values:
                visualize_correction_matrix = trim_2d_matrix(matrix=use_correction_matrix, rows=trimmed_dimensions[0],
                                                             columns=trimmed_dimensions[1])
            plot_statistic_matrices(frq_data={postprocessing_data_id: pd.DataFrame(visualize_correction_matrix)},
                                    data_id_to_show=postprocessing_data_id, in_file=in_bam, output_dir=tmp_dir_sample,
                                    y_tick_label_offset=min_flen + trimmed_dimensions[0].start, sample_id=sample_name,
                                    x_tick_label_offset=trimmed_dimensions[1].start, fig_width=1800, fig_height=2000,
                                    fig_fontsize=50, parent_logger=LOGGER, show_figure=show_plots)
    # create smoothed version of the matrix:
    if smooth_weights:
        postprocessing_data_id = f'{postprocessing_data_id}_smoothed'
        use_correction_matrix = smooth_2d_gc_weights(smooth_matrix=use_correction_matrix, parent_logger=LOGGER,
                                                     smoothing_intensity=smoothing_intensity, default_matrix_value=1.,
                                                     smoothing_kernel=smoothing_kernel, min_flen=min_flen)
        use_correction_matrix_path = save_matrix_to_txt(
            matrix=use_correction_matrix, output_dir=tmp_dir_sample,
            filename=str(Path('.'.join(use_correction_matrix_path.name.split('.')[:-2] +
                                       [f'{smoothing_intensity}I{smoothing_kernel}Smoothed', 'txt', 'gz']))),
            gzipped=True, float_data_precision=float_precision, report_saved_path=True, max_frag_length=max_flen,
            min_frag_length=min_flen)
        if visualize_matrices:
            visualize_correction_matrix = use_correction_matrix
            if focus_custom_values is not None:
                visualize_correction_matrix = trim_2d_matrix(matrix=use_correction_matrix, rows=trimmed_dimensions[0],
                                                             columns=trimmed_dimensions[1])
            plot_statistic_matrices(frq_data={postprocessing_data_id: pd.DataFrame(visualize_correction_matrix)},
                                    data_id_to_show=postprocessing_data_id, in_file=in_bam, output_dir=tmp_dir_sample,
                                    y_tick_label_offset=min_flen + trimmed_dimensions[0].start, sample_id=sample_name,
                                    x_tick_label_offset=trimmed_dimensions[1].start, fig_width=1800, fig_height=2000,
                                    fig_fontsize=50, parent_logger=LOGGER, show_figure=show_plots)
    # move all output from temporary sample dir into output dir after checking that the latter does not exist
    target_path = Path(out_dir_sample)
    target_path.mkdir(parents=True, exist_ok=True)  # ensure parent path exists to be able to move sample dir
    log(message=f"Copying GC-bias computation output from '{tmp_dir_sample}' to '{out_dir_sample}' ..",
        log_level=logging.INFO, logger_name=LOGGER)
    source_path = Path(tmp_dir_sample)
    for f in source_path.iterdir():  # handle all files/dirs separately
        # -> don't just batch delete existing target dir! Will exist probably due to GC-bias computation step
        if f.is_dir():  # subdir -> can be completely moved & target completely deleted
            # (will only concern the bam parts if they were kept in a previous analysis)
            target_dir_path = target_path / f.name
            if target_dir_path.exists():
                log(message=f"Target path for directory '{f.name}' exists and will be completely deleted! "
                            f"Deleting '{target_dir_path}' ..", log_level=logging.WARNING, logger_name=LOGGER)
                shutil.rmtree(target_dir_path)
            _ = shutil.move(f, target_path)
        elif f.is_file():
            target_file_path = target_path / f.name
            if target_file_path.exists():
                os.remove(target_file_path)
            _ = shutil.move(f, target_path)  # just move inside target dir
    if source_path.exists() and target_path.exists() and \
            (len([None for _2 in target_path.iterdir()]) >= len([None for _1 in source_path.iterdir()]) or
             len([None for _1 in source_path.iterdir()]) == 0):  # delete if same or more content OR source is empty
        shutil.rmtree(source_path)  # delete temp dir and all content that still exists
        # (should be gone after successfully moving entire dir)
    moved_correction_matrix_path = target_path / use_correction_matrix_path.name
    if not moved_correction_matrix_path.is_file():
        log(message=f"Moved correction matrix not found at '{moved_correction_matrix_path}'",
            log_level=logging.WARNING, logger_name=LOGGER)
    moved_weights_mask_matrix_path = target_path / weights_mask_path.name
    if not moved_weights_mask_matrix_path.is_file():
        log(message=f"Moved correction matrix not found at '{moved_weights_mask_matrix_path}'",
            log_level=logging.WARNING, logger_name=LOGGER)
    return (moved_correction_matrix_path, use_correction_matrix), (moved_weights_mask_matrix_path, weights_mask)


def load_until_leftmost_not_poly_n(loaded_sequences: List[Optional[str]], loaded_intervals: List[bool], chrom_seq: str,
                                   interval_loading_size=5000000, max_loaded=10, trigger_garbage_collection=False) \
        -> Tuple[List[Optional[str]], List[bool]]:
    if max_loaded < 2:
        log(message=f"Invalid choice for maximum number of loaded genomic intervals! Must be at least 2. Setting it to "
                    f"default of 10 ..",
            log_level=logging.WARNING, close_handlers=True, logger_name=LOGGER)
        max_loaded = 10
    first_loaded_sequence_index = loaded_intervals.index(True)
    while loaded_sequences[first_loaded_sequence_index].count('N') == interval_loading_size:
        loaded_sequences[first_loaded_sequence_index] = None  # unload leftmost poly-N interval
        loaded_intervals[first_loaded_sequence_index] = False
        if first_loaded_sequence_index >= len(loaded_sequences) - 2:  # add rightmost reference interval(s) if < 2 loaded
            for i in range(first_loaded_sequence_index - len(loaded_sequences) + 3):  # load >=2 non-poly-N intervals
                loaded_sequences.append(chrom_seq[interval_loading_size * len(loaded_intervals):
                                                  interval_loading_size * (len(loaded_intervals) + 1)])
                loaded_intervals.append(True)
        first_loaded_sequence_index = loaded_intervals.index(True)
    # check for too many loaded intervals
    if sum(loaded_intervals) > max_loaded:
        remove_n = sum(loaded_intervals) - max_loaded
        # unload interval with most neighboring not loaded entries; if tied, unload leftmost of these intervals
        for _i in range(remove_n):
            unloaded_neighbors = [0] * len(loaded_intervals)
            inverted_loaded_intervals = [not entr for entr in loaded_intervals]
            for pos in range(len(loaded_intervals)):
                unloaded_neighbors[pos] = sum(inverted_loaded_intervals[pos + 1:]) + sum(inverted_loaded_intervals[:pos])
            idx_interval_to_remove = unloaded_neighbors.index(max(unloaded_neighbors))
            loaded_sequences[idx_interval_to_remove] = None
            loaded_intervals[idx_interval_to_remove] = False
    if trigger_garbage_collection:
        gc.collect()
    return loaded_sequences, loaded_intervals


def extend_intervals_for_index(target_index: int, scaffold_length: int, scaffold_name: str,
                            loaded_sequences: List[Optional[str]], loaded_intervals: List[bool],
                            chrom_handle: TwoBitSequence, parent_logger: str, max_loaded=10,
                            interval_loading_size=5000000, trigger_garbage_collection=False) \
        -> Tuple[List[Optional[str]], List[bool]]:
    """

    :param target_index:
    :param scaffold_length:
    :param scaffold_name:
    :param loaded_sequences:
    :param loaded_intervals:
    :param chrom_handle:
    :param parent_logger:
    :param max_loaded:
    :param interval_loading_size:
    :param trigger_garbage_collection:
    :return:
    :raises: AttributeError if scaffold length is smaller than starting coordinate of interval which should be loaded
    """
    last_idx = target_index - len(loaded_intervals)
    for l_idx in range(target_index - len(loaded_intervals) + 1):
        if l_idx == last_idx:  # actually load the sequence
            hypothetical_end_coordinate = interval_loading_size * (len(loaded_intervals) + 1)
            if hypothetical_end_coordinate > scaffold_length:
                if scaffold_length < interval_loading_size * len(loaded_intervals):  # should never occur - illogical error
                    log(message=f"Critical error encountered loading interval  {target_index} for scaffold {scaffold_name}"
                                f". Length of scaffold {scaffold_length:,}bp was smaller than starting coordinate "
                                f"{interval_loading_size * len(loaded_intervals)}bp. Terminating..",
                        log_level=logging.CRITICAL, flush=True, logger_name=parent_logger)
                    raise AttributeError
                loaded_sequences.append(chrom_handle[interval_loading_size * len(loaded_intervals):scaffold_length].upper())
            else:
                loaded_sequences.append(chrom_handle[interval_loading_size * len(loaded_intervals):
                                                     interval_loading_size * (len(loaded_intervals) + 1)].upper())
            loaded_intervals.append(True)
        else:  # just elongate the list without loading anything
            loaded_sequences.append(None)
            loaded_intervals.append(False)
    # check for too many loaded intervals
    if sum(loaded_intervals) > max_loaded:
        remove_n = sum(loaded_intervals) - max_loaded
        # unload interval with most neighboring not loaded entries; if tied, unload leftmost of these intervals
        for _i in range(remove_n):
            unloaded_neighbors = [0] * len(loaded_intervals)
            inverted_loaded_intervals = [not entry for entry in loaded_intervals]
            for pos in range(len(loaded_intervals)):
                unloaded_neighbors[pos] = sum(inverted_loaded_intervals[pos + 1:]) + sum(inverted_loaded_intervals[:pos])
            idx_interval_to_remove = unloaded_neighbors.index(max(unloaded_neighbors))
            # unload leftmost interval with the most not-loaded neighbors
            loaded_sequences[idx_interval_to_remove] = None
            loaded_intervals[idx_interval_to_remove] = False
    if trigger_garbage_collection:
        gc.collect()
    return loaded_sequences, loaded_intervals


def load_specific_interval(target_index: int, scaffold_length: int, scaffold_name: str,
                        loaded_sequences: List[Optional[str]], loaded_intervals: List[bool],
                        chrom_handle: TwoBitSequence, parent_logger: str, max_loaded=10,
                        interval_loading_size=5000000, trigger_garbage_collection=False) \
        -> Tuple[List[Optional[str]], List[bool]]:
    """
    Use only if target_index lies within len(laoded_seqeunces) - 1!
    :param parent_logger:
    :param scaffold_name:
    :param target_index:
    :param scaffold_length:
    :param loaded_sequences:
    :param loaded_intervals:
    :param chrom_handle:
    :param max_loaded:
    :param interval_loading_size:
    :param trigger_garbage_collection:
    :return:
    """
    if target_index >= len(loaded_sequences):
        log(message=f"Cannot load interval {target_index} for scaffold {scaffold_name}. Terminating..",
            log_level=logging.CRITICAL, flush=True, close_handlers=True, logger_name=parent_logger)
        raise AttributeError
    if interval_loading_size * (target_index + 1) > scaffold_length:
        loaded_sequences[target_index] = chrom_handle[interval_loading_size * target_index:scaffold_length].upper()
    else:
        loaded_sequences[target_index] = chrom_handle[interval_loading_size * target_index:
                                                      interval_loading_size * (target_index + 1)].upper()
    loaded_intervals[target_index] = True
    # check for too many loaded intervals
    if sum(loaded_intervals) > max_loaded:
        remove_n = sum(loaded_intervals) - max_loaded
        # unload interval with most neighboring not-loaded-entries; if tied, unload leftmost of these intervals
        for _i in range(remove_n):
            unloaded_neighbors = [0] * len(loaded_intervals)
            inverted_loaded_intervals = [not entr for entr in loaded_intervals]
            for pos in range(len(loaded_intervals)):
                unloaded_neighbors[pos] = sum(inverted_loaded_intervals[pos + 1:]) + sum(inverted_loaded_intervals[:pos])
            idx_interval_to_remove = unloaded_neighbors.index(max(unloaded_neighbors))
            if target_index == idx_interval_to_remove:  # must not be target index!
                dummy = unloaded_neighbors[:]
                dummy[idx_interval_to_remove] = -1  # set to lowest
                idx_interval_to_remove = dummy.index(max(dummy))  # find next maximum
            loaded_sequences[idx_interval_to_remove] = None
            loaded_intervals[idx_interval_to_remove] = False
    if trigger_garbage_collection:
        gc.collect()
    return loaded_sequences, loaded_intervals


def unaligned_bam_worker(bam_path: OneOf[str, Path], output_path: OneOf[str, Path], tag_name: str):
    # binary filter:
    unaligned_read = np.uint32(4)
    unaligned_read_binary = bin(unaligned_read)
    unaligned_mate = np.uint32(8)
    unaligned_mate_binary = bin(unaligned_mate)
    unaligned_flags = np.uint32(unaligned_read + unaligned_mate)
    unaligned_flags_binary = bin(unaligned_flags)
    unaligned_options = (unaligned_read_binary, unaligned_mate_binary, unaligned_flags_binary)
    # -> fetch if "read/mate unmapped"
    # unaligned filter:
    # -----------------
    # REQUIRE THAT:
    # read unmapped = 4               '0b100'
    # mate unmapped = 8               '0b1000'
    # = 12 (are handled in unaligned reads extraction)
    with AlignmentFile(bam_path, mode='rb') as f_in:
        unaligned_pair_segments = filter(lambda a: bin(np.uint32(a.flag) & unaligned_flags) in unaligned_options,
                                         f_in.fetch(until_eof=True, multiple_iterators=True))
        # todo: add functionality to extract aligned but broken pairs! ("extract singular reads" or similar)
        with AlignmentFile(output_path, header=f_in.header, mode='wb') as f_unaligned_tagged:
            for aln in unaligned_pair_segments:
                # if not aln.is_mapped or not aln.mate_is_mapped:  # LEGACY
                aln.set_tag(tag=tag_name, value=0.)  # unaligned reads get a GC correction weight of 0 because the
                # fragment sequence cannot be safely inferred from the read halve
                f_unaligned_tagged.write(aln)


def bam_tagging_worker_single_interval(bam_path: str, correction_weights: np.array, temp_dir: str, gc_bases_offset: int,
                                       fragment_length_range: range, two_bit_reference_path: OneOf[str, Path],
                                       tagging_intervals_list: List[str], reference_lengths: Dict[str, int],
                                       sender_connection: mp_connection.Connection, parent_logger: str,
                                       default_weight: float = 1.0, tag_name=DEFAULT_TAG_NAME,
                                       ref_interval_loading_size: int = 500000, annotation_str=None):
    """

    :param default_weight:
    :param tagging_intervals_list: (chromosome, interval_start, interval_end, interval_size)
    :param gc_bases_offset: this is the offset for the weights matrix column required for weights retrieval
    :param parent_logger:
    :param two_bit_reference_path:
    :param annotation_str:
    :param bam_path:
    :param correction_weights:
    :param temp_dir:
    :param fragment_length_range:
    :param reference_lengths: length of reference chromosomes and other scaffolds
    :param sender_connection:
    :param tag_name:
    :param ref_interval_loading_size:
    :return:
    """
    # get chromosome sequence handle
    reference_handle = TwoBitFile(str(two_bit_reference_path))
    # (interval is sliced again subsequently for faster access)
    min_frag_len = fragment_length_range.start
    # compute corrected alignments -> use buffer of 200000 entries
    tagged_bam_files = []
    with silently_open_alignment_file(bam_path, mode='rb') as input_bam_file:
        for c_idx, (chromosome, start_coord, stop_coord, _ch_len) in enumerate(tagging_intervals_list):
            scaffold_length = reference_lengths[chromosome]
            chromosome_handle = reference_handle[chromosome]  # is large; just slice for interval sequence retrieval
            tagged_bam_file = str(Path(temp_dir) /
                                  '.'.join(Path(bam_path).name.split('.')[:-2] +
                                           [f"{Path(bam_path).name.split('.')[-2]}"  # anno str will be None
                                            f"+{chromosome}-{start_coord}-{stop_coord}.GCcorr"] +
                                           ([annotation_str, 'bam'] if annotation_str else ['bam'])))
            with silently_open_alignment_file(tagged_bam_file, mode='wb', template=input_bam_file) as f_gc_tagged:
                # reference scaffold handle management:
                # ______________________________________________________________________________________________________
                # preload first 5 Mbp of reference sequence; use 2 lists, one stores the sequences, the other stores
                # loading status of each interval as boolean value (intervals are loaded consecutively without gaps)
                if ref_interval_loading_size > scaffold_length:
                    ref_interval_sequences = [chromosome_handle[0:scaffold_length].upper()]
                    loaded_ref_intervals = [True]
                elif scaffold_length < 2 * ref_interval_loading_size:
                    ref_interval_sequences = [chromosome_handle[0:ref_interval_loading_size].upper(),
                                           chromosome_handle[ref_interval_loading_size:scaffold_length].upper()]
                    loaded_ref_intervals = [True, True]
                else:
                    ref_interval_sequences = [chromosome_handle[0:ref_interval_loading_size].upper(),
                                           chromosome_handle[ref_interval_loading_size:2 * ref_interval_loading_size].upper()]
                    loaded_ref_intervals = [True, True]
                # fastforward until no N-contigs are in ref_contig_intervals-deque any more
                try:
                    ref_interval_sequences, loaded_ref_intervals = load_until_leftmost_not_poly_n(
                        loaded_sequences=ref_interval_sequences, loaded_intervals=loaded_ref_intervals, max_loaded=10,
                        chrom_seq=chromosome_handle, interval_loading_size=ref_interval_loading_size)
                except AttributeError:
                    log(message=f"An error occurred when trying to load initial intervals for tagging. Cannot continue "
                                f"processing scaffold '{chromosome}'. Exiting BAM tagging worker..",
                        log_level=logging.ERROR, logger_name=parent_logger)
                    sender_connection.send(-1)
                    sender_connection.close()
                    return -1
                # ______________________________________________________________________________________________________
                # binary filter:
                exclude_flags = np.uint32(12)
                exclude_flags_binary = bin(exclude_flags)
                # -> not necessary to check for "mapped" attribute (filtered out if "read unmapped")
                # complete alignment filter: skip alns outside of interval borders and unaligned reads
                # ---------------------------------------------------------------------------------
                # EXCLUDE if:
                # read unmapped = 4    '0b100'
                # mate unmapped = 8    '0b1000'
                #          SUM  = 12 (excluded alns are handled in unaligned reads extraction)
                # REQUIRE THAT:
                # alignment start lies inside interval -> otherwise processed by previous interval!
                alignments_mapped_to_interval = filter(
                    lambda a: bin(~np.uint32(a.flag) & exclude_flags) == exclude_flags_binary and
                              start_coord <= a.pos < stop_coord,
                    input_bam_file.fetch(chromosome, start_coord, stop_coord, multiple_iterators=True))
                # iterate over alignments from specified reference scaffolds
                aln_buffer = []  # (re-) set alignment buffer
                for aln_seg in alignments_mapped_to_interval:
                    try:
                        frag_size = abs(aln_seg.template_length)
                        if frag_size < min_frag_len:
                            raise IndexError
                        correction_weights_row = correction_weights[frag_size - min_frag_len]
                    except TypeError:  # this should not occur at all! -> unaligned are filtered and are only extracted
                        # if the --output-unaligned flag was set by another process!
                        # in case of unaligned segment or mate: '<' not supported between instances of
                        # 'NoneType' and 'int'
                        aln_seg.set_tag(tag_name, value=0.,  # give unaligned reads a GC weight of 0.
                                        value_type="f", replace=True)  # should not occur when fetching intervals
                        aln_buffer.append(aln_seg)  # no need to check aln buffer; will be written in non-default case
                        continue
                    except IndexError:  # fragment length not in reduced weight matrix -> use default value
                        aln_seg.set_tag(tag_name, value=default_weight,
                                        value_type="f", replace=True)
                        aln_buffer.append(aln_seg)  # no need to check aln buffer; will be written in non-default case
                        continue
                    if frag_size <= abs(aln_seg.query_alignment_length):  # get sequence from aligned read portion
                        frag_seq = aln_seg.query_alignment_sequence.upper()
                        # ^-- This is now a substring of query_sequence, excluding flanking bases that were soft clipped
                    else:  # get fragment sequence from reference genome
                        frag_start_scaffold = min(aln_seg.reference_start, aln_seg.next_reference_start)
                        target_interval_index = frag_start_scaffold // ref_interval_loading_size
                        if target_interval_index >= len(loaded_ref_intervals):  # required interval not loaded -> load!
                            try:
                                ref_interval_sequences, loaded_ref_intervals = extend_intervals_for_index(
                                    target_index=target_interval_index, loaded_sequences=ref_interval_sequences,
                                    parent_logger=parent_logger, loaded_intervals=loaded_ref_intervals,
                                    chrom_handle=chromosome_handle, scaffold_name=chromosome, max_loaded=10,
                                    interval_loading_size=ref_interval_loading_size, scaffold_length=scaffold_length)
                            except AttributeError:
                                log(message=f"Error occurred when trying to load a downstream interval. Cannot continue "
                                            f"processing scaffold '{chromosome}'. Exiting BAM tagging worker..",
                                    log_level=logging.ERROR, logger_name=parent_logger)
                                sender_connection.send(-1)
                                sender_connection.close()
                                return -1
                        elif not loaded_ref_intervals[target_interval_index]:
                            try:  # interval in range but was unloaded/not loaded (shouldn't be necessary)
                                ref_interval_sequences, loaded_ref_intervals = load_specific_interval(
                                    target_index=target_interval_index, loaded_sequences=ref_interval_sequences,
                                    loaded_intervals=loaded_ref_intervals, chrom_handle=chromosome_handle,
                                    max_loaded=10, scaffold_name=chromosome, parent_logger=parent_logger,
                                    interval_loading_size=ref_interval_loading_size, scaffold_length=scaffold_length)
                            except AttributeError:
                                log(message="Error occurred when trying to load an unloaded/intermediate interval. "
                                            f"Cannotcontinue processing scaffold '{chromosome}'. "
                                            "Exiting BAM tagging worker..",
                                    log_level=logging.ERROR, logger_name=parent_logger)
                                sender_connection.send(-1)
                                sender_connection.close()
                                return -1
                        frag_start_interval = frag_start_scaffold % ref_interval_loading_size
                        frag_seq = ref_interval_sequences[target_interval_index][
                                   frag_start_interval:frag_start_interval + frag_size].upper()
                    try:  # to retrieve correction weight
                        gc_column_index = frag_seq.count('G') + frag_seq.count('C') \
                                          - gc_bases_offset  # shift back by trimmed columns!
                        if gc_column_index < 0:
                            raise IndexError
                        corr_factor = correction_weights_row[gc_column_index]
                    except IndexError:  # also: for not-in-proper-pair alns (across breakpoints) or alns mapping..
                        corr_factor = default_weight  # ..far apart due to deletion; also for fragment lengths + GC ..
                        # ..bases pointing to columns/rows that exclusively contain the default corr-factor and were
                        # trimmed (weights matrix is trimmed -> index errors)
                    aln_seg.set_tag(tag_name, value=corr_factor, value_type="f", replace=True)
                    aln_buffer.append(aln_seg)
                    if len(aln_buffer) > 10000:
                        _ = deque(map(lambda a: f_gc_tagged.write(a), aln_buffer), maxlen=0)  # overwrites existing tags
                        aln_buffer = []
                # write remaining buffer elements
                if len(aln_buffer):
                    _ = deque(map(lambda a: f_gc_tagged.write(a), aln_buffer), maxlen=0)  # overwrites existing tags
            # append created BAM file to return list
            tagged_bam_files.append(tagged_bam_file)
    # report path to new BAM file
    sender_connection.send(tagged_bam_files)
    sender_connection.close()


def try_clear_temp_dir_and_exit(tmp_dir: str, exit_code: int, message=None):
    if message is not None:
        log(message=message, log_level=logging.ERROR, logger_name=LOGGER)
    if Path(tmp_dir).is_dir():
        shutil.rmtree(tmp_dir, ignore_errors=True)  # remove full directory with content -> use shutil rm_tree
    if LOGGER:
        for hdlr in logging.getLogger(LOGGER).handlers:
            hdlr.flush()
            hdlr.close()
    sys.exit(exit_code)


def trim_2d_matrix(matrix: np.array, rows=(0, 0), columns=(0, 0)) -> np.array:
    if isinstance(rows, range):
        rows = (rows.start, rows.stop)
    if isinstance(columns, range):
        columns = (columns.start, columns.stop)
    # trim rows:
    if rows[1]:
        matrix = matrix[rows[0]:-rows[1], :]
    else:
        matrix = matrix[rows[0]:, :]
    # trim rows
    if columns[1]:
        matrix = matrix[:, columns[0]:-columns[1]]
    else:
        matrix = matrix[:, columns[0]:]
    return matrix


def reduce_matrix(matrix_to_trim: np.array, trim_dimensions_exclusively_containing: List[float],
                  border_elements: Optional[int]) -> Tuple[np.array, Tuple[range, range]]:
    """
    Remove non-informative matrix entries so that (stored) numpy matrix is smaller (faster/smaller on disc)
    :param border_elements:
    :param matrix_to_trim:
    :param trim_dimensions_exclusively_containing:
    :return:
    """
    if border_elements is None:  # should not occur but casting according to logic for safety
        border_elements = 0
    # trim leading lines
    delete_initial_rows = 0
    for line_idx in range(matrix_to_trim.shape[0]):
        if all(map(lambda v: v in trim_dimensions_exclusively_containing, matrix_to_trim[line_idx, :])):
            delete_initial_rows += 1
        else:
            break
    if delete_initial_rows > border_elements:  # don't remove lines otherwise
        delete_initial_rows -= border_elements
        matrix_to_trim = matrix_to_trim[delete_initial_rows:, :]
    else:
        delete_initial_rows = 0
    # trim tailing lines
    delete_tailing_rows = 0
    for inverse_line_idx in range(1, matrix_to_trim.shape[0] + 1):
        if all(map(lambda v: v in trim_dimensions_exclusively_containing, matrix_to_trim[-inverse_line_idx, :])):
            delete_tailing_rows += 1
        else:
            break
    if delete_tailing_rows > border_elements:  # don't remove lines otherwise
        delete_tailing_rows -= border_elements
        matrix_to_trim = matrix_to_trim[:-delete_tailing_rows, :]
    else:
        delete_tailing_rows = 0
    # trim leading columns
    delete_initial_columns = 0
    for column_idx in range(matrix_to_trim.shape[1]):
        if all(map(lambda v: v in trim_dimensions_exclusively_containing, matrix_to_trim[:, column_idx])):
            delete_initial_columns += 1
        else:
            break
    if delete_initial_columns > border_elements:  # don't remove columns otherwise
        delete_initial_columns -= border_elements
        matrix_to_trim = matrix_to_trim[:, delete_initial_columns:]
    else:
        delete_initial_columns = 0
    # trim tailing columns
    delete_tailing_columns = 0
    for column_idx in range(1, matrix_to_trim.shape[1] + 1):
        if all(map(lambda v: v in trim_dimensions_exclusively_containing, matrix_to_trim[:, -column_idx])):
            delete_tailing_columns += 1
        else:
            break
    if delete_tailing_columns > border_elements:  # don't remove columns otherwise
        delete_tailing_columns -= border_elements
        matrix_to_trim = matrix_to_trim[:, :-delete_tailing_columns]
    else:
        delete_tailing_columns = 0
    return matrix_to_trim, (range(delete_initial_rows, delete_tailing_rows),
                            range(delete_initial_columns, delete_tailing_columns))


def get_unaligned_reads(bam_path: OneOf[str, Path], output_dir: OneOf[str, Path], tag_name: str) \
        -> Tuple[Optional[Path], mp.Process]:
    sample_id = Path(bam_path).stem
    output_bam_unaligned = Path(output_dir) / f'{sample_id}.unaligned.bam'
    output_dir.mkdir(parents=True, exist_ok=True)  # ensure output directory exists
    unaligned_process_handle = mp.Process(target=unaligned_bam_worker,
                                          kwargs={'bam_path': bam_path, 'output_path': output_bam_unaligned,
                                                  'tag_name': tag_name})
    unaligned_process_handle.start()
    return output_bam_unaligned, unaligned_process_handle


def samtools_cat_bams(list_of_bams: List[str], samtools_path: OneOf[str, Path],
                      tmp_dir: OneOf[str, Path], output_bam: Path, keep_input=False):
    concatenation_command = [str(samtools_path), 'cat', '-o', output_bam, '--no-PG', '--threads', '4'] + list_of_bams
    called_concatenation_command = sp.run(concatenation_command)
    try:
        called_concatenation_command.check_returncode()
    except sp.CalledProcessError:
        exit_message = f"Subprocess for concatenating scaffold BAM files failed. Terminating main.."
        try_clear_temp_dir_and_exit(tmp_dir=tmp_dir, exit_code=1, message=exit_message)
    if not Path(output_bam).is_file():
        exit_message = f"Concatenated BAM file '{output_bam}' not found. Terminating main.."
        try_clear_temp_dir_and_exit(tmp_dir=tmp_dir, exit_code=1, message=exit_message)
    # index final BAM file
    create_bam_index(bam_path=output_bam, samtools_path=samtools_path, check_success=True)
    # remove temp dir
    if not keep_input:
        temp_bam_parent = Path(list_of_bams[0]).parent
        for d_bm in list_of_bams:
            if os.path.isfile(d_bm):
                os.remove(d_bm)
            d_idx_putative = str(d_bm)[:-1] + 'i'
            if os.path.isfile(d_idx_putative):
                os.remove(d_idx_putative)
            d_2idx_putative = f'{d_bm}.bai'
            if os.path.isfile(d_2idx_putative):
                os.remove(d_2idx_putative)
        if temp_bam_parent.is_dir():
            shutil.rmtree(temp_bam_parent)


def order_bams(bam_list: List[str]) -> List[str]:
    # extract ref scaffold order form header:
    with AlignmentFile(bam_list[0], mode='rb') as f_scaff_order:
        scaffold_order = f_scaff_order.references
    # if the provided value is not in the header
    bam_paths_with_locus = []
    for bam_nm in bam_list:
        chrm, strt, stp = Path(bam_nm).name.split('.GCcorr')[0].split('+')[-1].split('-')  # chrom-start-stop
        bam_paths_with_locus.append(((chrm, int(strt)), str(bam_nm)))
    return [bam_path for _, bam_path in
            sorted(bam_paths_with_locus, key=lambda s: (scaffold_order.index(s[0][0]), s[0][1]))]


def get_genomic_intervals_for_tagging(bam_for_tagging: OneOf[str, Path], interval_size=TAGGING_INTERVAL_SIZE,
                                      offset=0) -> list:
    ref_lengths = get_reference_tuples(bam=bam_for_tagging)
    whole_genome_regions = [(chrm, offset, r_len) for chrm, r_len in ref_lengths]
    genomic_intervals = []
    for chrom, _strt, stop in whole_genome_regions:
        n_splits = (stop - offset) // interval_size
        cum_size = offset
        for _split_idx in range(0, n_splits, 1):
            genomic_intervals.append((chrom, cum_size, cum_size + interval_size, interval_size))
            cum_size += interval_size
        resudual_bases = stop - n_splits * interval_size
        if resudual_bases * 10 <= interval_size and n_splits:  # remaining part of scaffold is <= 10% of interval size
            # -> just add to last! BUT: we need to have at least one full interval. Otherwise, just add the entire sequence
            #    because it is smaller than one interval
            last_chrom, last_start, last_end, _chk_size = genomic_intervals[-1]
            genomic_intervals[-1] = (last_chrom, last_start, stop, interval_size + resudual_bases)
        else:  # just add as separate interval otherwise
            genomic_intervals.append((chrom, cum_size, stop, stop - cum_size))
    return genomic_intervals


def tag_bam_with_correction_weights_parallel(sample_output_dir: str, two_bit_genome_file: OneOf[str, Path],
                                             threads: int, correction_matrix: np.array, frag_len_range: range,
                                             bam_path: str, ref_lengths: Dict[str, int],
                                             temporary_directory_sample: str, gc_base_limits: range,
                                             output_unaligned=False, default_fragment_weight: float = 1.,
                                             tag_name=DEFAULT_TAG_NAME, samtools_path=DEFAULT_SAMTOOLS_PATH):
    """
    Size increase of BAM file: 6.8 Gb to 6.9 Gb ~= 1.5%
    Test on the 22/11/2022: duration of BAM file tagging was 0:11:30 (h:mm:ss)
    :param default_fragment_weight:
    :param gc_base_limits:
    :param frag_len_range:
    :param threads:
    :param sample_output_dir:
    :param bam_path:
    :param correction_matrix:
    :param ref_lengths:
    :param two_bit_genome_file:
    :param temporary_directory_sample:
    :param tag_name:
    :param samtools_path:
    :param output_unaligned:
    :return:
    """
    gc_start = gc_base_limits.start  # this is the offset for the weights matrix column required for weights retrieval
    n_processes = threads - 1  # 1 process for extracting unaligned reads
    # prepare parallel processes for alignment correction
    receivers = []
    senders = []
    all_worker_kwargs = []
    intervals_per_proc = []
    interval_lengths_per_proc = []
    proc_list_done = []
    for _idx in range(n_processes):
        all_worker_kwargs.append([])
        recv, sender = mp.Pipe(duplex=False)
        receivers.append(recv)
        senders.append(sender)
        intervals_per_proc.append([])
        interval_lengths_per_proc.append([])
        proc_list_done.append(False)
    # reference_contigs -> ref_lengths! compute total number of bases;
    # add contigs until over threshold, then record base+/-; regard in next iteration (add fewer/more bases)
    # -> as soon as the next added contig would overrun the threshold,
    #    look for the most fitting contig and add that one instead
    genomic_intervals = get_genomic_intervals_for_tagging(bam_for_tagging=bam_path,
                                                          interval_size=TAGGING_INTERVAL_SIZE, offset=0)
    target_base_sum_per_process = 1 + sum(map(lambda c: c[3], genomic_intervals)) // n_processes
    for scaff_idx, (scaff, strt, stp, scaff_length) in enumerate(genomic_intervals):
        target_proc_idx = scaff_idx % n_processes
        # check if adding to list is possible; if so, just add to current (base_offset does not change)
        if sum(interval_lengths_per_proc[target_proc_idx]) + scaff_length <= target_base_sum_per_process:
            interval_lengths_per_proc[target_proc_idx].append(scaff_length)
            intervals_per_proc[target_proc_idx].append((scaff, strt, stp, scaff_length))
        else:  # the current contig should be added somewhere else and choose a better suited list
            bp_overshot = [sum(interval_lengths_per_proc[list_search_index]) + scaff_length - target_base_sum_per_process
                           for list_search_index in range(n_processes)]
            final_index = bp_overshot.index(min(bp_overshot))
            interval_lengths_per_proc[final_index].append(scaff_length)
            intervals_per_proc[final_index].append((scaff, strt, stp, scaff_length))
    # feedback to user
    log(message='\n'.join([f"process ID {cid} processes contigs: {', '.join(map(lambda e: e[0], chrms))}\n"
                           f"This corresponds to {sum(interval_lengths_per_proc[cid]):,} bp (bp compared to target of "
                           f"{target_base_sum_per_process:,}bp = "
                           f"{sum(interval_lengths_per_proc[cid]) / target_base_sum_per_process:.1%})"
                           for cid, chrms in enumerate(intervals_per_proc)]), log_level=logging.DEBUG, logger_name=LOGGER)
    # start unaligned reads extraction
    worker_output_path = Path(temporary_directory_sample) / 'scaffold_BAMs_pre-merging'
    worker_output_path.mkdir(parents=True, exist_ok=True)
    unaligned_bam = None
    unaligned_extraction_handle = None
    if output_unaligned:
        unaligned_bam, unaligned_extraction_handle = get_unaligned_reads(
            bam_path=bam_path, output_dir=worker_output_path, tag_name=tag_name)  # uses 1 process
    # create worker kwargs
    worker_output_dir = str(worker_output_path)
    for w_idx in range(n_processes):
        all_worker_kwargs[w_idx] = {'bam_path': bam_path, 'correction_weights': correction_matrix,
                                    'tag_name': tag_name, 'temp_dir': worker_output_dir,
                                    'fragment_length_range': frag_len_range,
                                    'two_bit_reference_path': two_bit_genome_file,
                                    'tagging_intervals_list': intervals_per_proc[w_idx],
                                    'sender_connection': senders[w_idx],
                                    'reference_lengths': ref_lengths,
                                    'parent_logger': LOGGER,
                                    'gc_bases_offset': gc_start,
                                    'default_weight': default_fragment_weight}
    # create worker processes
    bam_tagging_workers = [mp.Process(target=bam_tagging_worker_single_interval, kwargs=worker_kwargs)
                           for worker_kwargs in all_worker_kwargs]
    # start all workers, receive results, wait until finished, and close them
    _ = deque(map(lambda w: w.start(), bam_tagging_workers), maxlen=0)  # ~1GB RAM usage using 12 dual-thread processes
    tagged_scaffold_bam_files = []
    try:
        _ = [tagged_scaffold_bam_files.extend(incoming_result.recv())
             for incoming_result in receivers]  # get returns from workers
    except TypeError:  # " 'int' object is not iterable" occurs if one of the bam tagging workers returns with error
        log(message="Error in BAM tagging worker detected. Cannot proceed. Terminating..",
            log_level=logging.ERROR, close_handlers=True, logger_name=LOGGER)
        sys.exit(3)
    _ = deque(map(lambda w: w.join(), bam_tagging_workers), maxlen=0)
    _ = deque(map(lambda w: w.close(), bam_tagging_workers), maxlen=0)
    # sanity check: do we have any scaffold/contig BAM files for merging?
    if len(tagged_scaffold_bam_files) == 0:
        log(message="Did not get any scaffold/contig-wise BAM files for merging. Cannot proceed. Terminating..",
            log_level=logging.ERROR, close_handlers=True, logger_name=LOGGER)
        sys.exit(3)
    # define final output path
    tagged_bam_file_path = Path(temporary_directory_sample) / \
        '.'.join(Path(bam_path).name.split('.')[:-2] +
                 [f"{Path(bam_path).name.split('.')[-2]}", "GCtagged", "bam"])
    tagged_bam_file = tagged_bam_file_path
    # concatenate BAM files and index
    tagged_scaffold_bam_files_in_order = order_bams(bam_list=tagged_scaffold_bam_files)
    if output_unaligned and unaligned_extraction_handle is not None:
        unaligned_extraction_handle.join(timeout=READS_EXTRACTION_TIMEOUT)  # wait for max. 30 minutes
        if unaligned_extraction_handle.is_alive():
            log(message=f"Unaligned reads extraction terminated due to timeout "
                        f"({READS_EXTRACTION_TIMEOUT} seconds). You can increase this via the "
                        f"'--reads-extraction-timeout' commandline parameter. Continuing ..",
                log_level=logging.WARNING, logger_name=LOGGER)
            unaligned_extraction_handle.terminate()
        else:
            unaligned_extraction_handle.close()
    # move from temporary sample dir into output dir after checking that the latter does not exist
    target_path = Path(sample_output_dir)  # this is the sample output dir
    target_path.mkdir(parents=True, exist_ok=True)  # ensure parent path exists to be able to move sample dir
    log(message=f"Combining tagged scaffold BAMs '{temporary_directory_sample}' to '{sample_output_dir}' ..",
        log_level=logging.INFO, logger_name=LOGGER)
    if unaligned_bam is not None and unaligned_bam.is_file():
        create_bam_index(bam_path=unaligned_bam, samtools_path=samtools_path, check_success=True)
        # check number of total reads in unaligned BAM file
        with silently_open_alignment_file(unaligned_bam, mode='rb') as f_ubam:
            index_statistics = f_ubam.get_index_statistics()
        if sum([total for _cont, _mp, _ump, total in index_statistics]) == 0:
            log(message=f"There were no unaligned reads detected in the input file. Nothing will be output.",
                log_level=logging.INFO, logger_name=LOGGER)
            unaligned_bam.unlink()  # delete empty uBAM
            Path(f'{unaligned_bam}.bai').unlink()  # delete empty index
        else:  # there were unaligned reads
            if output_unaligned:
                _ = shutil.move(unaligned_bam, target_path)  # copy uBAM to target dir
                _ = shutil.move(f'{unaligned_bam}.bai', target_path)  # move index for index-free samtools cat
            else:
                unaligned_bam.unlink()  # delete empty uBAM
                Path(f'{unaligned_bam}.bai').unlink()  # delete empty index
    # combine individual scaffold BAM files to
    log(message=f"Merging tagged scaffold BAMs into temporary output file '{tagged_bam_file}' ..",
        log_level=logging.INFO, logger_name=LOGGER)
    samtools_cat_bams(list_of_bams=tagged_scaffold_bam_files_in_order, samtools_path=samtools_path,
                      keep_input=False, tmp_dir=temporary_directory_sample, output_bam=tagged_bam_file)  # also indexes
    # move output files from temporary directory to target directory
    log(message=f"Moving GC-bias computation output from '{temporary_directory_sample}' to '{sample_output_dir}' ..",
        log_level=logging.INFO, logger_name=LOGGER)
    for f in Path(temporary_directory_sample).iterdir():  # handle all files/dirs separately
        # -> don't just batch delete existing target dir! Will exist probably due to GC-bias computation step
        if f.is_dir():  # subdir -> can be completely moved & target completely deleted
            # (will only concern the bam parts if kept)
            target_dir_path = target_path / f.name
            if target_dir_path.exists():
                log(message=f"Target path for directory '{f.name}' exists and will be completely deleted! "
                            f"Deleting '{target_dir_path}' ..", log_level=logging.WARNING, logger_name=LOGGER)
                shutil.rmtree(target_dir_path)
            _ = shutil.move(f, target_path)
        elif f.is_file():
            target_file_path = target_path / f.name
            if target_file_path.exists():
                os.remove(target_file_path)
            _ = shutil.move(f, target_path)  # just move inside target dir
    source_path = Path(temporary_directory_sample)
    if source_path.exists() and target_path.exists() and \
            (len([None for _2 in target_path.iterdir()]) == len([None for _1 in source_path.iterdir()]) or
             len([None for _1 in source_path.iterdir()]) == 0):  # same content OR source is empty
        shutil.rmtree(source_path)  # delete temp dir and all content that still exists
        # (should be gone after successfully moving entire dir)


def reduce_weights_for_tagging(weights_path: Path, mask_path: Optional[Path], sample_id: str,
                               mask_matrix: Optional[np.array], correction_matrix: Optional[np.array],
                               weights_flen_range=None, mask_flen_range=None, default_weight: float = 1.) \
        -> Tuple[np.array, range, range]:
    """
    :param weights_path:
    :param mask_path: str or Path to mask matrix file created during tagging procedure. If available, the
    matrix will be used to reduce the weights matrix to non-default values. If not provided, the weights matrix parent
    directory will be searched for such a matrix file (assertion: naming of files not changed). The weights matrix
    will be reduced in size based on it's non-default values instead, if no mask matrix was provided/found made
    available. Downside of this fallback behaviour: potentially non-default weights merely resulting from blurring
    weights will be included in the correction procedure.
    :param sample_id:
    :param mask_matrix:
    :param correction_matrix:
    :param weights_flen_range:
    :param mask_flen_range:
    :param default_weight: weight assigned to attribute combinations which lie outside reduced weights.
                           1. for GCparagon [DEFAULT]; may be 0. for other algorithms like Griffin
    :return:
    """
    reduced_mask = None
    deleted_rows = range(0)
    deleted_columns = range(0)
    if correction_matrix is None:  # try to load correction weights from file if no matrix given
        if not weights_path.is_file():
            raise AttributeError(f"correction weights matrix file {weights_path} does not exist")
        correction_matrix, weights_flen_range = load_txt_to_matrix_with_meta(filename=weights_path, to_dtype=np.float64,
                                                                             loading_logger=LOGGER)
    elif weights_flen_range is None:
        raise AttributeError(f"correction_matrix was provided but weights_flen_range was not provided!")
    if mask_matrix is None:  # try to load from mask file if no mask matrix is given
        if mask_path is None or not mask_path.is_file():
            try:  # to find mask based on correction_matrix name
                candidate_matrix = sorted(list(Path(weights_path.parent).glob(
                    f"{sample_id}*_gc_bias_computation_mask.txt*")), key=lambda c: len(c.name), reverse=True)[0]
                mask_matrix, mask_flen_range = load_txt_to_matrix_with_meta(filename=candidate_matrix, to_dtype=bool,
                                                                            loading_logger=LOGGER)
            except (IndexError, AttributeError):  # IndexError out of range; 'NoneType' object no attribute 'parent'
                pass  # there is no weights file -> only matrix; behaviour not optimal!
                # Will potentially use correction weights that are non-default just because of blurring
        else:
            mask_matrix, mask_flen_range = load_txt_to_matrix_with_meta(filename=mask_path,
                                                                        to_dtype=bool,
                                                                        loading_logger=LOGGER)
    if mask_flen_range is not None and \
            mask_flen_range == weights_flen_range:  # sanity check dimensions - discard mask if different
        reduced_mask, (deleted_rows, deleted_columns) = reduce_matrix(matrix_to_trim=mask_matrix, border_elements=0,
                                                                      trim_dimensions_exclusively_containing=[False])
    if reduced_mask is None:  # workaround: reduce based on correction weights themselves
        reduced_weights, (deleted_rows, deleted_columns) = reduce_matrix(
            matrix_to_trim=correction_matrix, trim_dimensions_exclusively_containing=[default_weight, 0.],
            border_elements=0)
        resulting_flen_range = range(weights_flen_range.start + deleted_rows.start,
                                     weights_flen_range.stop - deleted_rows.stop)  # shift by removed rows
    else:  # default behavior
        resulting_flen_range = range(mask_flen_range.start + deleted_rows.start,
                                     mask_flen_range.stop - deleted_rows.stop)  # shift by removed rows
        reduced_weights = trim_2d_matrix(matrix=correction_matrix, rows=deleted_rows, columns=deleted_columns)
    resulting_gc_bases_range = range(deleted_columns.start,
                                     correction_matrix.shape[1] - deleted_columns.stop - 1)  # -1 because zero start idx
    # the RANGE ALWAYS NAMES THE ACTUAL GC base counts of the column! last column: GC base count 550 but 551 elements!
    # removed  - 1!
    # resulting_gc_bases_range are the actual GC base counts that define the first and last column
    # everything outside that automatically receives a default weight! Default needs to be set appropriate to algorithm!
    # assert resulting_gc_bases_range.stop - resulting_gc_bases_range.start + 1 == reduced_weights.shape[1]  # valid
    return reduced_weights, resulting_flen_range, resulting_gc_bases_range


def get_reference_tuples(bam: OneOf[str, Path]) -> Tuple[Tuple[Any, Any]]:
    with silently_open_alignment_file(bam, mode='rb', threads=1) as input_bam_file:
        # allow for any error output resulting from BAM access to be printed here
        return tuple(zip(input_bam_file.references, input_bam_file.lengths))


def get_reference_contig_lengths(bam: OneOf[str, Path]):
    reference_tuples = get_reference_tuples(bam=bam)
    # create chromosome/scaffold length dictionary
    reference_lengths = {}.fromkeys([r_n for r_n, _r_l in reference_tuples])
    for ref_name, ref_len in reference_tuples:
        reference_lengths[ref_name] = ref_len
    return reference_lengths


def sufficient_aligned_reads_available(bam_path: str, target_fragments_processed: int,
                                       intervals_for_processing_with_score: GenomicIntervalList) -> Tuple[bool, int]:
    """
    Only standard-non-sex-scaffolds are counted. Number of mappable GRCh38 bases (without exclusion list): 2,848,552,403
    :param bam_path:
    :param target_fragments_processed:
    :param intervals_for_processing_with_score:
    :return: bool if enough fragments can be expected from te computation
    """
    std_chroms = [f'chr{idx}' for idx in range(1, 23, 1)]
    with silently_open_alignment_file(bam_path, mode='rb') as f_in_bam:
        index_statistics = f_in_bam.get_index_statistics()
    mapped_std_reads = sum(filter(lambda x: x is not None,
                                  [idx_stat.mapped if idx_stat.contig.lower() in std_chroms else None
                                   for idx_stat in index_statistics]))
    # get number of bases in all intervals = max. portion of reference genome we can compute over
    bases_in_all_std_intervals = sum(filter(
        lambda x: x is not None, [stop - start if chrom in std_chroms else None
                                  for chrom, start, stop, _score in intervals_for_processing_with_score]))
    fraction_intervaled_genome = bases_in_all_std_intervals / 2848552403  # second number is "usable" bases of hg38
    expected_number_fragments_in_intervals = int(round(mapped_std_reads / 2 * fraction_intervaled_genome, ndigits=0))
    # paired-end data -> divide by 2
    # expected_number_fragments_in_intervals = 0  # coding instruction: if exception, catch it here and set value to 0
    return expected_number_fragments_in_intervals >= target_fragments_processed, expected_number_fragments_in_intervals


def create_bam_index(bam_path: OneOf[Path, str], samtools_path: OneOf[Path, str], check_success=True):
    bam_path = Path(bam_path)
    indexing_command = [str(samtools_path), 'index', str(bam_path)]
    ran_indexing_subp = sp.run(indexing_command)
    try:
        ran_indexing_subp.check_returncode()
    except sp.CalledProcessError:
        log(message=f"BAM indexing process for file '{Path(bam_path).name}' returned with error. Cannot proceed."
                    f"Terminating ..",
            log_level=logging.ERROR, close_handlers=True, logger_name=LOGGER)
        sys.exit(2)
    if check_success:
        # check for existence of putative BAM index:
        if not any([Path(potential_index).is_file()
                    for potential_index in (f'{bam_path}.bai',
                                            Path(bam_path.parent) / f'{bam_path.stem}.bai')]):
            log(message=f"No BAM index found for file '{bam_path}'. Cannot proceed. Terminating ..",
                log_level=logging.ERROR, close_handlers=True, logger_name=LOGGER)
            sys.exit(2)


def fix_bam_index(bam_path: OneOf[Path, str], samtools_path: str, silent=True):
    with AlignmentFile(bam_path, 'rb') as f_aln:
        try:
            f_aln.get_index_statistics()  # check_index(self) could also be used
        except ValueError:  # ValueError: mapping information not recorded in index or index not available
            # create missing index
            if not silent:
                log(message=f"GCparagon requires index statistics to be present in BAM index file. Creating missing "
                            f"index file for bam '{bam_path}' ..", log_level=logging.INFO, logger_name=LOGGER)
            create_bam_index(bam_path=bam_path, samtools_path=samtools_path)
        except AttributeError:  # if htsfile is SAM formatted and thus has no index
            log(message=f"input BAM file is actually a SAM file. Code requires a BAM file. Terminating ..",
                log_level=logging.ERROR, close_handlers=True, logger_name=LOGGER)
            sys.exit(2)


def manage_bad_intervals(bad_genomic_intervals_bed: Optional[str]) \
        -> Tuple[Optional[Dict[Tuple[str, int, int],
                               List[int]]],
                 Optional[str]]:
    if bad_genomic_intervals_bed is not None:  # figure out most recent one!
        candidate_bad_library_files = Path(bad_genomic_intervals_bed).parent.glob(
            Path(bad_genomic_intervals_bed).name)  # yields generator object -> no iteration if empty/nothing found
        lib_dates = []
        for bc_library in candidate_bad_library_files:  # may contain timestamp
            cur_lib = Path(bc_library)
            try:
                lib_dates.append(time.strptime('_'.join(cur_lib.stem.split('_')[-2:]), TIMESTAMP_FORMAT))
            except ValueError:  # is raised if the expected timestamp is not present
                pass
        # find most recent time stamp, assuming default output naming -> just use input as-is otherwise
        try:
            most_recent_lib_date = max(lib_dates)
            new_bad_genomic_intervals_bed_path = Path(
                time.strftime(str(Path(bad_genomic_intervals_bed).parent / BAD_INTERVALS_FILE_PATTERN_AS_PATH),
                              most_recent_lib_date))
            if new_bad_genomic_intervals_bed_path.is_file():
                log(message=f"Using bad intervals library file '{new_bad_genomic_intervals_bed_path}'",
                    log_level=logging.INFO, logger_name=LOGGER)
                bad_genomic_intervals_bed = new_bad_genomic_intervals_bed_path
        except ValueError:  # max() arg is an empty sequence
            pass  # leave bad_genomic_intervals_bed as was input
    if bad_genomic_intervals_bed is None or not Path(bad_genomic_intervals_bed).exists():  # if none was found/path is not valid
        return None, None
    bad_intervals = read_bad_genomic_intervals_bed_file(bed_path=bad_genomic_intervals_bed)  # requires score in col 4 and integer field in col 5
    return bad_intervals, bad_genomic_intervals_bed


def read_gc_distribution(ref_table_path: OneOf[str, Path]):
    gc_dist = []
    with open(ref_table_path, 'rt') as f_gc:
        hdr_content = f_gc.readline().strip().split('\t')  # gc_percentage	relative_frequency
        if hdr_content != ['gc_percentage', 'relative_frequency']:
            raise ValueError(f"expected header content like: "
                             f"['gc_percentage', 'relative_frequency'] but got {hdr_content}.")
        gc_pc_range = tuple(range(0, 101, 1))
        for line_idx, data_line in enumerate(f_gc.readlines()):
            gc_content, relative_freq = data_line.strip().split('\t')
            assert int(gc_content) == gc_pc_range[line_idx]  # relative frequencies in ascending GC content order
            gc_dist.append(float(relative_freq))
    return np.array(gc_dist)


def infer_intervals_for_n_fragments(intervals: List[Tuple[str, int, int]], bam_path: OneOf[str, Path],
                                    target_fragment_count: int, flength_range: range,
                                    estimate_from_n_intervals: int = 8, repetitions: int = 3):
    log(message=f"Inferring required number of genomic intervals to reach {target_fragment_count:,} processed "
                f"fragments ..", log_level=logging.INFO, logger_name=LOGGER)
    start_time = time.localtime()
    max_num_intervals = len(intervals)
    rep_results = []
    n_frag_results = []
    for rep in range(repetitions):
        # create random intervals list
        n_frags = 0
        random_intervals = []
        for i in range(estimate_from_n_intervals):
            not_drawn = True
            while not_drawn:
                random_interval = random.choice(intervals)
                if random_interval in random_intervals:
                    continue  # draw again
                random_intervals.append(random_interval)
                not_drawn = False
        assert len(set(random_intervals)) == estimate_from_n_intervals
        with AlignmentFile(bam_path, 'rb') as f_aln:
            for chrom, start, end, *rest in random_intervals:
                exclude_flags = np.uint32(3852)  # = 256 + 2048 + 512 + 1024 + 4 + 8
                exclude_flags_binary = bin(exclude_flags)
                # -> not necessary to check for "mapped" attribute (filtered out if "read unmapped")
                # complete alignment filter:
                # --------------------------
                # EXCLUDE if:
                # read unmapped = 4               '0b100'
                # mate unmapped = 8               '0b1000'
                # not primary (secondary) = 256   '0b100000000'
                # vendor/QC fail = 512            '0b1000000000'
                # PCR or optical duplicate = 1024 '0b10000000000'
                # supplementary = 2048            '0b100000000000'
                # = 3852
                # REQUIRE THAT:
                # alignment is paired = 1 '0b1'
                # mates map to different strands
                #    a.is_forward != a.mate_is_forward
                # TLEN column is (positive and) between defined fragment length limits (inclusive)
                #    min_frag_len <= a.template_length <= max_frag_len
                filtered_alignments = filter(lambda a:
                                             # bin(np.uint32(a.flag) & paired_flag) == paired_flag_binary and
                                             a.is_paired and
                                             bin(~np.uint32(a.flag) & exclude_flags) == exclude_flags_binary and
                                             a.is_forward != a.mate_is_forward and
                                             (flength_range.start <= a.template_length <= flength_range.stop),
                                             f_aln.fetch(chrom, start, end, multiple_iterators=True))
                for _aln in filtered_alignments:
                    n_frags += 1
        # estimated number of required GI to reach target_fragment_count:
        use_n_intervals = math.ceil(target_fragment_count / n_frags * estimate_from_n_intervals) + 1
        rep_results.append(use_n_intervals)
        n_frag_results.append(n_frags)
    # consolidate:
    actually_use_n_intervals = round(np.mean(rep_results))
    actually_required_fragments = round(np.mean(n_frag_results))
    if actually_use_n_intervals * 1.1 < 1.0:
        actually_use_n_intervals += 1  # require one additional interval
        actually_required_fragments += actually_required_fragments / (actually_use_n_intervals - 1)  # estimate
    else:
        actually_use_n_intervals *= 1.1  # add 10%
        actually_required_fragments *= 1.1
    # cast to closest int
    actually_use_n_intervals = round(actually_use_n_intervals)
    actually_required_fragments = round(actually_required_fragments)
    elapsed_time = datetime.timedelta(seconds=time.mktime(time.localtime()) - time.mktime(start_time))
    log(message=f"Estimated in {elapsed_time} (h:mm:ss) from {estimate_from_n_intervals:,} genomic intervals that "
                f"{actually_use_n_intervals:,} intervals will be needed to reach the target number of "
                f"{target_fragment_count:,} processed fragments (est. number: "
                f"{round(actually_required_fragments / estimate_from_n_intervals * actually_use_n_intervals):,}; "
                f"average over {repetitions} repetitions.)",
        log_level=logging.INFO, logger_name=LOGGER)
    # limit insanity
    if actually_use_n_intervals < 10:
        if max_num_intervals < 10:
            log(message=f"the number of estimated genomic intervals was below 10. Will use all {max_num_intervals:,} "
                        "pre-defined intervals instead.",
                log_level=logging.WARNING, logger_name=LOGGER)
            actually_use_n_intervals = max_num_intervals
        else:
            log(message=f"the number of estimated genomic intervals was below 10. Will use 10 intervals instead!",
                log_level=logging.WARNING, logger_name=LOGGER)
            actually_use_n_intervals = 10
    elif actually_use_n_intervals > max_num_intervals:
        log(message=f"the number of estimated genomic intervals was higher than the number of available genomic "
                    f"intervals. Will use all {max_num_intervals:,} pre-defined intervals.",
            log_level=logging.WARNING, logger_name=LOGGER)
        actually_use_n_intervals = max_num_intervals
    return actually_use_n_intervals


def aes(target: np.array, components: np.array, weights: np.array) -> float:
    """
    Asserts target sums up to 1 and each row of components to sum up to 1;
    The product of components and weights must also result in 1.0 if summed.
    :param target: target distribution array relative to which the erorr is computed
    :param components: the components that will be weighted by weights to approximate the target distribution
    :param weights: np.array of weights with its dimension matching the number of components
    :return: sum of absolute errors
    """
    return np.abs((components.T * weights).sum(axis=1) - target).sum()


def objective_function_mse(scale_factor: float, to_scale: np.array, target_func: np.array):
    # Ensure positive and non-zero scaling factor
    if scale_factor <= 0:
        return np.inf  # Return a large value to indicate infeasibility
    scaled = scale_factor * to_scale
    mean_squared_diff = np.mean((scaled - target_func) ** 2)
    return np.sum(mean_squared_diff)


def reconstruct_distribution_nnls(target_distribution: np.array,
                                  components: Dict[str, np.array],
                                  component_order: List[str], verbose=True) \
        -> Tuple[Dict[str, float], float]:
    start_time = time.perf_counter_ns()
    # force all components to individually sum up to one:
    for cmp_lab, cmp in components.items():  # TypeError: 'cell' object is not callable
        if not (0.999999 <= cmp.sum() <= 1.00000001):
            components[cmp_lab] = cmp / cmp.sum()
    if isinstance(target_distribution, (list, tuple)):  # not a np.array
        target_distribution = np.array(target_distribution)
    # create normalized target distribution:
    normalized_target_distribution = target_distribution / target_distribution.sum()
    ordered_components = np.array([components[cmp] for cmp in component_order])
    # make sure most weights are not zero
    n_components = len(component_order)
    if ordered_components.shape != (len(component_order), len(target_distribution)):  # check dimensions
        if ordered_components.T.shape != (len(component_order), len(target_distribution)):
            raise ValueError(f"the expected dimensions of the components matrix ({ordered_components.shape}) seems to "
                             "be broken (not transposed, I checked that).")
        ordered_components = ordered_components.T  # try transposed version...
    # compute initial error
    initial_re = aes(target=normalized_target_distribution, components=ordered_components,
                     weights=np.array([1 / n_components] * n_components))
    component_weight_boundaries = (0.1, 10)
    # compute component fit:
    residuals_from_component_fit = []
    for cmp_idx, o_cmp in enumerate(ordered_components):  # scale each component linearly
        # Initial guess for the scaling factor:
        initial_guess = normalized_target_distribution.sum() / o_cmp.sum()
        # Use scipy.optimize.minimize to find the optimal scaling factor
        # Define constraint to ensure that the scaling factor is positive
        positive_constraint = {'type': 'ineq',  # "[..] inequality means that it is to be non-negative."
                               'fun': lambda x: x}
        result = minimize(fun=objective_function_mse, x0=initial_guess, args=(o_cmp, normalized_target_distribution),
                          method='SLSQP', constraints=positive_constraint)
        optimal_scale_factor = result.x[0]
        scaled_component_residual = aes(target=normalized_target_distribution,
                                        components=np.array([o_cmp]),
                                        weights=np.array([optimal_scale_factor]))
        residuals_from_component_fit.append((scaled_component_residual, component_order[cmp_idx]))
    # such that sum of multiplied components is one
    min_weight, max_weight = component_weight_boundaries
    residual_values = [res[0] for res in residuals_from_component_fit]
    _5pc, upper_bound = np.percentile(residual_values, (5, 95))  # to make more resilient against outliers
    lower_bound = min(residual_values)  # try lowering lower bound
    # transform residuals to weights
    weights_per_component = {}.fromkeys(component_order)
    for cmp_idx, o_cmp in enumerate(component_order):  # scale each component linearly
        res, cmp_lab = residuals_from_component_fit[cmp_idx]
        assert o_cmp == cmp_lab
        if res <= lower_bound:
            weights_per_component[o_cmp] = max_weight  # low deviation -> gets highest weight!
        elif upper_bound <= res:
            weights_per_component[o_cmp] = min_weight  # large deviation -> gets lowest weight
        else:  # scaled somewhere in the middle
            weights_per_component[o_cmp] = ((max_weight - min_weight) *
                                            (1 - (res - lower_bound) / (upper_bound - lower_bound))
                                            + min_weight)  # res == upper_bound -> min_weight AND
            # res == lower_bound -> max_weight
    ordered_weights = np.array([weights_per_component[o_cmp] for o_cmp in component_order])
    total_weights_sum = ordered_weights.sum()
    final_ordered_weights = ordered_weights / total_weights_sum
    # compute final weights which can directly be used for the computation of a weighted mean across GIs
    for o_cmp in component_order:
        weights_per_component[o_cmp] /= total_weights_sum
    # compute improvement:
    residual_error = aes(target=normalized_target_distribution,
                         components=ordered_components,
                         weights=final_ordered_weights)

    elapsed_time_ns = time.perf_counter_ns() - start_time
    if residual_error > initial_re:
        log(message=f"Could not achieve a reduction of reconstruction error (took {elapsed_time_ns / 10**6:,.2f} ms). "
                    f"Will use default weights to combine results from genomic intervals.",
            log_level=logging.INFO, logger_name=LOGGER)
        fallback_weights_per_component = {}.fromkeys(weights_per_component)
        for comp in fallback_weights_per_component.keys():
            fallback_weights_per_component[comp] = 1 / n_components
        return fallback_weights_per_component, initial_re
    log(message=f"Successfully reduced the AES reconstruction error of the genomic GC content from {initial_re:.3f} to "
                f"{residual_error:.3f} (took {elapsed_time_ns / 10**6:,.2f} ms).", log_level=logging.INFO,
        logger_name=LOGGER)
    if verbose:
        # default case: reduced residual error between weighted mean of preselected GI FGCDs and the reference FGCD
        log(message=f"Range of weights used for minimizing MSE: [{min_weight}, {max_weight}]).\n"
                    f"(Range of residuals used for the computation of component weights: "
                    f"[{lower_bound:.6f}, {upper_bound:.6f}])\n"
                    f"Range of weights (before weighted mean): [{min(ordered_weights):.3f}, "
                    f"{max(ordered_weights):.3f}].\n"
                    f"Range of actual weights is: [{min(final_ordered_weights):.6f}, "
                    f"{max(final_ordered_weights):.6f}].",
            log_level=logging.INFO, logger_name=LOGGER)
    return weights_per_component, residual_error


def preselect_genomic_intervals(genomic_intervals_sorted: List[Tuple[str, int, int]], reference_fgcd: np.array,
                                interval_fgcds, bam_file: OneOf[str, Path], target_fragment_number: int,
                                fragment_length_range: range, output_path: OneOf[str, Path], sample_name: str,
                                show_figures: bool = False) -> Tuple[List[Tuple[float, Tuple[str, int, int]]], float]:
    # infer number of required intervals based on # fragments in 10 intervals
    inferred_number_of_required_intervals = infer_intervals_for_n_fragments(
        intervals=genomic_intervals_sorted, bam_path=bam_file, target_fragment_count=target_fragment_number,
        flength_range=fragment_length_range, repetitions=3, estimate_from_n_intervals=8)
    log(message=f"Will use {inferred_number_of_required_intervals:,} genomic intervals for GC bias computation.",
        log_level=logging.INFO, logger_name=LOGGER)
    # select genomic intervals
    preselected_intervals = genomic_intervals_sorted[:inferred_number_of_required_intervals]
    # create interval FGCD subset
    preselected_interval_fgcds = {}
    interval_order = []
    ordered_intervals = []
    for pi_c, pi_s, pi_e, *rest in preselected_intervals:  # rest should be empty...
        region_id = create_region_label(chrm=pi_c, start=pi_s, end=pi_e)
        interval_order.append(region_id)
        if not rest:  # expect this to be an empty list: []
            ordered_intervals.append((pi_c, pi_s, pi_e))
        else:
            ordered_intervals.append((pi_c, pi_s, pi_e, *rest))
        try:
            preselected_interval_fgcds[region_id] = interval_fgcds[region_id]
        except KeyError:
            log(message=f"There was no precomputed FGCD available for region '{region_id}'. Recompute!",
                log_level=logging.ERROR, logger_name=LOGGER, close_handlers=True)
            sys.exit(2)  # TODO: compute the interval instead of terminating !!
    # create linear combination for best reference FGCD representation
    interval_weights_per_component, reconstruction_residual = reconstruct_distribution_nnls(
        target_distribution=reference_fgcd, components=preselected_interval_fgcds, component_order=interval_order,
        verbose=False)  # supress additional debugging info output
    visualize_weights(region_weights=np.array(list(interval_weights_per_component.values()), dtype=float),
                      sample_label=sample_name, out_dir=output_path, compute_skew=True, compute_curtosis=True,
                      show_figure=show_figures)
    return list(zip([interval_weights_per_component[create_region_label(chrm=hrm, start=strt, end=nd)]
                     for hrm, strt, nd in ordered_intervals], ordered_intervals)), reconstruction_residual


def main() -> int:
    global LOGGER  # commandline logging only

    cmd_args = get_cmdline_args()
    # input options
    input_bams = cmd_args.input_bams
    two_bit_reference_file = cmd_args.two_bit_reference_file
    genomic_intervals_bed_file = cmd_args.genomic_intervals_bed_file
    exclude_genomic_intervals_bed_file = cmd_args.exclude_genomic_intervals_bed_file
    correction_weights_matrix_path = cmd_args.correction_weights
    mask_path = cmd_args.weights_mask
    ref_gc_dist_path = cmd_args.ref_gc_dist_path
    # processing options
    preset_number = cmd_args.parameter_preset_number
    n_simulations = cmd_args.n_simulations
    upper_limit_fragment_length = cmd_args.upper_limit_fragment_length
    lower_limit_fragment_length = cmd_args.lower_limit_fragment_length
    total_number_threads = cmd_args.total_number_threads
    random_seed = cmd_args.random_seed
    samtools_path = cmd_args.samtools_path
    process_n_fragments = cmd_args.process_n_fragments
    min_frag_occurs = cmd_args.min_frag_occurs
    only_tag_bam = cmd_args.only_tag_bam  # deactivates bias computation; a weights matrix must be supplied!
    strict_n_base_exclusion = not cmd_args.allow_n_base_fragments
    min_unclipped_aln_fraction = cmd_args.min_unclipped_aln_fraction
    default_fragment_weight = cmd_args.default_fragment_weight
    # post-processing options
    detect_outliers = cmd_args.detect_outliers
    outliers_method = cmd_args.outlier_method
    outlier_stringency = cmd_args.outlier_stringency
    smooth_weights = cmd_args.smooth_weights
    smoothing_kernel = cmd_args.smoothing_kernel
    smoothing_intensity = cmd_args.smoothing_intensity
    # output options
    verbose = cmd_args.verbose
    output_directory = cmd_args.output_directory
    temporary_directory = cmd_args.temporary_directory
    plot_result = cmd_args.plot_result
    output_simulation_results = cmd_args.output_simulation_results
    keep_interval_data = cmd_args.keep_interval_data
    output_corrected_bam = cmd_args.output_corrected_bam
    floating_point_precision = cmd_args.floating_point_precision
    gc_tag_name = cmd_args.gc_tag_name
    write_updated_bad_intervals_library = cmd_args.write_updated_bad_intervals_library
    focus_plots = not cmd_args.dont_focus_plots
    show_plots = cmd_args.show_plots
    output_unaligned_reads = cmd_args.output_unaligned_reads
    # processing settings
    compute_bias = not only_tag_bam
    if correction_weights_matrix_path is not None:
        correction_weights_matrix_path = Path(correction_weights_matrix_path)
    if mask_path is not None:
        mask_path = Path(mask_path)
    np.seterr(all='raise')
    # CHECK INPUT:
    # check and fix correctable parameters
    exit_after_warnings = 0
    print_warnings = []
    if two_bit_reference_file is None:
        if not EXPECTED_TWO_BIT_REFERENCE_GENOME_PATH.is_file():
            print("cannot proceed - no two-boit reference file defined and default expected file not present under "
                  f"{EXPECTED_TWO_BIT_REFERENCE_GENOME_PATH}. Terminating ..")
            sys.exit(3)
        two_bit_reference_file = EXPECTED_TWO_BIT_REFERENCE_GENOME_PATH
    if floating_point_precision < 2:
        print_warnings.append(f"Floating pint precision was set to {floating_point_precision} but needs to "
                              f"be at least 3! Setting to default of {DEFAULT_FLOAT_PRECISION} instead. Continuing ..")
        floating_point_precision = DEFAULT_FLOAT_PRECISION
    if lower_limit_fragment_length < 1:
        print_warnings.append(f"Lower limit for fragment lengths was set to {lower_limit_fragment_length} but must be "
                              f"at least 1. Setting to default value of {DEFAULT_MIN_FRAGMENT_LENGTH} instead ..")
        lower_limit_fragment_length = DEFAULT_MIN_FRAGMENT_LENGTH
    if upper_limit_fragment_length <= lower_limit_fragment_length:
        print_warnings.append(f"Upper limit for fragment lengths was set to {upper_limit_fragment_length} but must be "
                              f"at least 1 higher than lower limit of {lower_limit_fragment_length}. Setting to "
                              f"default value of {DEFAULT_MAX_FRAGMENT_LENGTH} instead ..")
        upper_limit_fragment_length = DEFAULT_MAX_FRAGMENT_LENGTH
    if lower_limit_fragment_length > MAX_FRAGMENT_LENGTH or upper_limit_fragment_length > MAX_FRAGMENT_LENGTH:
        raise ValueError(f"Maximum allowed fragment length (= {MAX_FRAGMENT_LENGTH:,} bp) violation: one of the "
                         "specified minimum or maximum fragment length specified exceeded defined limits: "
                         f"[{lower_limit_fragment_length:,} bp, {upper_limit_fragment_length:,} bp].")
    if upper_limit_fragment_length < lower_limit_fragment_length:
        raise ValueError(f"Maximum fragment length was smaller than minimum fragment length: "
                         f"[{lower_limit_fragment_length:,} bp, {upper_limit_fragment_length:,} bp].")
    if not (0. <= min_unclipped_aln_fraction <= 1.):
        print_warnings.append(f"Minimum unclipped alignment fraction was set to {min_unclipped_aln_fraction} but must "
                              "be a floating point value between 0 and 1. Setting to default value of "
                              f"{DEFAULT_MIN_UNCLIPPED_ALN_FRACTION} instead ..")
        min_unclipped_aln_fraction = DEFAULT_MIN_UNCLIPPED_ALN_FRACTION
    # TODO: make preset not taking precedence over defined parameters -> make preset settings customizable!
    #       If user customizes, check if there was also a preset != 0 and add "<preset-str>-CUSTOMIZED" to output dir!
    # set preset parameters if defined
    if preset_number:  # 1, 2, or 3; 0 means no changes relative to default parameters
        match preset_number:
            case 1:
                min_frag_occurs = 2
                process_n_fragments = 5 * 10 ** 6
                n_simulations = 6
                smooth_weights = True
                smoothing_kernel = 'gauss'
                smoothing_intensity = 5
                detect_outliers = True
                outliers_method = 'IQR'
                outlier_stringency = 2
            case 2:
                min_frag_occurs = 10
                process_n_fragments = 50 * 10 ** 6
                n_simulations = 4
                smooth_weights = True
                smoothing_kernel = 'gauss'
                smoothing_intensity = 2
                detect_outliers = True
                outliers_method = 'IQR'
                outlier_stringency = 2
            case 3:
                min_frag_occurs = 20
                process_n_fragments = 99999999999  # usually for 30x WGS and default bins, ~170M frags can be processed
                n_simulations = 4
                smooth_weights = True
                smoothing_kernel = 'gauss'
                smoothing_intensity = 2
                detect_outliers = True
                outliers_method = 'IQR'
                outlier_stringency = 2
    elif min_frag_occurs < ABSOLUTE_MIN_OCCURRENCES:  # ABSOLUTE_MIN_OCCURRENCES:
        print_warnings.append(f"Lower limit for fragment length/GC-bases combinations was set to {min_frag_occurs} "
                              f"but must be at least {ABSOLUTE_MIN_OCCURRENCES}. Setting to suggested minimal value "
                              f"of {ABSOLUTE_MIN_OCCURRENCES} instead ..")
        min_frag_occurs = ABSOLUTE_MIN_OCCURRENCES
    if n_simulations < 1:
        print_warnings.append(f"Number of simulations per interval was set to {n_simulations} but must be at least 1. "
                              f"Setting to default of {DEFAULT_SIMULATION_REPETITIONS} instead ..")
        n_simulations = DEFAULT_SIMULATION_REPETITIONS
    # find cpu count boundaries (asserts hyper-threading architecture)
    available_logical_cores = len(os.sched_getaffinity(0))
    max_efficiently_usable_physical_cores = available_logical_cores // 2 if available_logical_cores > 1 else 1
    max_efficiently_usable_threads = max_efficiently_usable_physical_cores
    if total_number_threads is None:  # only optimize number of threads if no number was specified
        total_number_threads = DEFAULT_NUMBER_PROCESSES
        if max_efficiently_usable_physical_cores * 2 < 4:  # general HW check
            print_warnings.append('GCparagon requires at least 4 logical cores for being able to run. Only '
                                  f'{max_efficiently_usable_physical_cores * 2} were estimated to be available. '
                                  f'Exiting..')
            exit_after_warnings = 1
        if total_number_threads > max_efficiently_usable_threads:
            print_warnings.append(f'CPUs to use in multiprocessing operations was set to {total_number_threads} but '
                                  f'number of available efficiently usable logical cores estimated to be only '
                                  f'{max_efficiently_usable_threads}. Setting to {max_efficiently_usable_threads}.')
            total_number_threads = max_efficiently_usable_threads
    elif total_number_threads > available_logical_cores:  # limit to 75% of max. available cores
        print_warnings.append(f'CPUs to use in multiprocessing operations was set to {total_number_threads} but number '
                              f'of available logical cores is only {available_logical_cores}. Setting to '
                              f'{int(available_logical_cores/4*3)} (= 75%).')
        total_number_threads = int(available_logical_cores/4*3)
    # check unfixable parameters
    two_bit_reference_file_path = Path(two_bit_reference_file)
    if not two_bit_reference_file_path.is_file():
        raise AttributeError(f"2bit reference genome file '{two_bit_reference_file}' does not exist!")
    if not samtools_path or not Path(samtools_path).exists() or not Path(samtools_path).is_file():
        raise AttributeError("path to samtools executable either not found or not accessible. Please provide a valid "
                             "and accessible path using '-sp' or '--samtools-path'.")
    if not ref_gc_dist_path.is_file():
        raise FileNotFoundError(f"the reference GC content distribution file could not be found/accessed!")

    if compute_bias:  # load information that needs to be loaded only once
        # choose most recent bad intervals library version if multiple are present in parent directory
        bad_intervals, exclude_genomic_intervals_bed_file = manage_bad_intervals(
            bad_genomic_intervals_bed=exclude_genomic_intervals_bed_file)
        ref_gc_dist = read_gc_distribution(ref_table_path=ref_gc_dist_path)
        # read interval list (should be >=150Mbp in total)
        generally_processable_intervals_with_score, interval_gc_content_distributions = read_scored_regions_bed_file(
            bed_path=genomic_intervals_bed_file)

    # process all input BAM files sequentially
    for input_bam in input_bams:
        log_file, error_log = None, None
        try:
            if not os.path.isfile(input_bam):
                raise AttributeError(f"input BAM file '{input_bam}' does not exist!")
            # manage imago parameters
            input_bam_path = Path(input_bam)
            input_bam_parent_path = input_bam_path.parent
            sample_id = input_bam_path.stem
            if only_tag_bam and not Path(correction_weights_matrix_path).is_file():
                print_warnings.append('input argument --correction-weights missing. '
                                      'Tag-only-mode not possible. Exiting..')
                exit_after_warnings = 1
            if not output_directory or output_directory == str(input_bam_parent_path):
                print_warnings.append('Output directory is either input BAM parent directory or was None. Setting '
                                      "it to subdirectory of input BAM parent directory: 'GC_correction_output'")
                output_directory = str(input_bam_parent_path / f'GC_bias_correction_GCparagon{VERSION_STRING}')
            # set up target output directory and logfile
            start_time = time.localtime()
            sample_out_dir_path = Path(output_directory) / sample_id
            if sample_out_dir_path.exists() and compute_bias:  # do NOT delete if tag only mode is active!
                print_warnings.append(f"Output path for GC bias computation exists. Deleting completely: "
                                      f"'{sample_out_dir_path}'")
                shutil.rmtree(sample_out_dir_path)  # will fail if read-only files are present!
            try:
                sample_out_dir_path.mkdir(parents=True, exist_ok=True)  # ensure output directory for sample exists
            except FileExistsError:  # path is a file -> delete it!
                print_warnings.append(f"Output directory for sample {sample_id} path is a file. Deleting file ..")
                os.remove(sample_out_dir_path)
                sample_out_dir_path.mkdir(parents=True, exist_ok=True)
            sample_out_dir = str(sample_out_dir_path)
            # set up logging (cmdline handler + file handler are created)
            log_file = sample_out_dir_path / (f"{sample_id}_GCbiasCorrection_"
                                              f"{time.strftime('%d-%m-%Y_%H-%M', start_time)}.log")
            error_log = log_file.parent / f"{log_file.stem}.err"
            if LOGGER != 'GCparagon':  # entered when specifying logfiles for the first time
                LOGGER = set_up_logging(logfile_path=log_file, logger_name='GCparagon', verbose=verbose)

            else:  # logger exists - delete handlers of logger and create new ones (and add to logger)
                assert logging.getLogger('GCparagon')
                set_new_log_paths(logfile_path=log_file,
                                  logger_name='GCparagon',  # a logger instance with this name MUST exist!
                                  verbose=verbose)
            if print_warnings:
                for warning_message in print_warnings:
                    log(message=warning_message, log_level=logging.WARNING, logger_name=LOGGER)
            if exit_after_warnings:
                return exit_after_warnings
            if temporary_directory:
                if temporary_directory == input_bam_parent_path:
                    log(message="Temporary directory is identical to input BAM parent directory. Setting it to "
                                "subdirectory of input BAM parent directory: 'GC_correction_tmp'",
                        log_level=logging.WARNING, logger_name=LOGGER)
                    sample_temp_dir_path = input_bam_parent_path / 'GC_correction_tmp'
                else:
                    sample_temp_dir_path = Path(temporary_directory) / sample_id
            else:  # define a temporary directory
                sample_temp_dir_path = Path(sample_out_dir) / 'GC_correction_tmp'  # required for BAM merging (samtools)
            sample_temp_dir = str(sample_temp_dir_path)
            # check if index is there. if not, create it!
            fix_bam_index(bam_path=input_bam, samtools_path=samtools_path, silent=False)
            # get reference scaffolds and lengths
            reference_contig_lengths = get_reference_contig_lengths(bam=input_bam)  # stderr enabled for AlignmentFile
            # TASK 1: compute the GC bias present in the BAM file of the current sample
            #         WARNING: don't use multi-sample/run BAM files!)
            log(message=f"Number of available physical cores: {max_efficiently_usable_physical_cores:,}. "
                        f"Will use {total_number_threads:,} threads.",
                log_level=logging.INFO, logger_name=LOGGER)
            correction_weights_matrix = None
            weights_mask = None
            if compute_bias:
                # Integrate information about all encountered bad intervals so far using bad regions library file which
                # contains bad regions that were recorded during the analyiss of previous samples. Remove regions based
                # on "target fragments processed" which acts as a lower boundary for how many fragments are expected to
                # be available for processing by GCparagon in a dataset.
                expect_sufficient_fragment_count, n_expected_fragments = sufficient_aligned_reads_available(
                    bam_path=input_bam, target_fragments_processed=process_n_fragments,
                    intervals_for_processing_with_score=generally_processable_intervals_with_score)
                log(message='GCparagon (GC-bias computation) Started.\n' +
                            f"|---------------------------------------------------------------------------------\n"
                            f"|   Configuration for processing sample {sample_id} was:\n"
                            f"|   ++++++++++++++++++++++++++++++++++++++++{'+' * len(sample_id)}\n" +
                            (f'|   Using parameter preset {preset_number}\n'
                             if preset_number in (1, 2, 3) else '') +
                            f"|   Minimum fragment length: {lower_limit_fragment_length:,}bp\n"
                            f"|   Maximum fragment length: {upper_limit_fragment_length:,}bp\n"
                            f"|   Minimum number of specific fragment attribute combination occurrences: "
                            f"{min_frag_occurs:,}\n"
                            f"|   Target number of fragments to process: " +
                            ('all' if process_n_fragments == 9999999999 else f'{process_n_fragments:,}') +
                            f" ({'' if expect_sufficient_fragment_count else 'not '}expected to be reached)\n" +
                            ('' if expect_sufficient_fragment_count or not n_expected_fragments else
                             f"|   Number of fragments estimated from BAM index statistics that will be processed: "
                             f"{n_expected_fragments:,}\n") +
                            f"|   Repetitions for simulation of expected GC-content per interval: {n_simulations:,}\n"
                            f"|   Random seed was: {random_seed}\n"
                            f"|   Temporary data will be written to: {sample_temp_dir}\n"
                            f"|   Final results will be moved from temporary path to directory: {sample_out_dir}\n"
                            f"|---------------------------------------------------------------------------------",
                    log_level=logging.INFO, logger_name=LOGGER)
                if not expect_sufficient_fragment_count:
                    if n_expected_fragments == 0:
                        log(message="Number of expected aligned fragments from all intervals could not be determined "
                                    "from the BMA index. The BAM index might not contain information about mapped "
                                    "reads per contig/scaffold or you likely have provided an unaligned BAM file "
                                    "(uBAM).", log_level=logging.WARNING, logger_name=LOGGER)
                    else:
                        log(message=f"Number of aligned fragments expected from all intervals "
                                    f"(= {n_expected_fragments:,}) over the genome was lower than the target fragment "
                                    f"count!",
                            log_level=logging.WARNING, logger_name=LOGGER)
                # sort intervals based on exclusion list overlap
                # check all predefined intervals against bad intervals library (remove exactly matching intervals)
                sorted_eligible_genomic_intervals = sort_intervals_by_exclusion_ist_overlap(
                    all_intervals_with_score=generally_processable_intervals_with_score, bad_intervals=bad_intervals,
                    expected_dataset_fragments=n_expected_fragments,
                    max_overlap_percentage=DEFAULT_MAX_INTERVAL_PERCENTAGE_EXCLUSIONLIST_OVERLAP)
                log(message=f"Predefined genomic intervals loaded: {len(sorted_eligible_genomic_intervals):,} genomic "
                            "intervals available for processing.", log_level=logging.INFO, logger_name=LOGGER)
                if len(sorted_eligible_genomic_intervals) < 10:
                    log(message="The number of predefined genomic intervals is low! You might consider decreasing "
                                "the size of your intervals and increase their number.",
                        log_level=logging.WARNING, logger_name=LOGGER)
                target_genomic_intervals_with_weights, _reconstruction_error = preselect_genomic_intervals(
                    genomic_intervals_sorted=sorted_eligible_genomic_intervals,  # = (chrom, start, end)
                    reference_fgcd=ref_gc_dist, bam_file=input_bam, target_fragment_number=process_n_fragments,
                    interval_fgcds=interval_gc_content_distributions, output_path=sample_out_dir_path,
                    fragment_length_range=range(lower_limit_fragment_length, upper_limit_fragment_length),
                    show_figures=show_plots, sample_name=sample_id)
                # get intervals with weights from estimated reference genome fragment GC content reconstruction!
                # -> linear combination of selected intervals GC content that best approximates the reference GC content
                # following a typical cfDNA fragment length distribution (provided for blood plasma cfDNA, plasmaSeq
                # protocol in file: 'reference_fragment_length_distribution.tsv')
                # NOW: compute GC bias!
                ((correction_weights_matrix_path, correction_weights_matrix),
                 (mask_path, weights_mask)) = compute_gc_bias_parallel(
                    visualize_matrices=plot_result, output_all=output_simulation_results, in_bam=input_bam,
                    out_dir_sample=sample_out_dir, min_unclipped_aln_fraction=min_unclipped_aln_fraction,
                    max_flen=upper_limit_fragment_length,
                    weighted_intervals_to_process=target_genomic_intervals_with_weights,
                    min_flen=lower_limit_fragment_length, reference_gc_content_distribution=ref_gc_dist,
                    simulation_count=n_simulations, threads=total_number_threads, tmp_dir_sample=sample_temp_dir,
                    keep_interval_data=keep_interval_data, random_seed=random_seed, sample_name=sample_id,
                    chrom_sizes=reference_contig_lengths, float_precision=floating_point_precision,
                    two_bit_reference_file=two_bit_reference_file_path, min_frag_occurs=min_frag_occurs,
                    target_fragments_processed=process_n_fragments, expected_yield=n_expected_fragments,
                    write_updated_bad_intervals_library=write_updated_bad_intervals_library,
                    bad_intervals_library_file=exclude_genomic_intervals_bed_file,
                    strict_n_base_exclusion=strict_n_base_exclusion, plot_focus_border=10 if focus_plots else None,
                    interval_gc_content_distributions=interval_gc_content_distributions,
                    detect_outliers=detect_outliers, outlier_detection_method=outliers_method,
                    outlier_detection_stringency=outlier_stringency, smooth_weights=smooth_weights,
                    smoothing_kernel=smoothing_kernel, smoothing_intensity=smoothing_intensity, show_plots=show_plots)
                # compute end time and give feedback
                log(message=f"Correction weights matrix averaged from {n_simulations:,} simulations written to file: "
                            f"'{correction_weights_matrix_path}'", log_level=logging.INFO, logger_name=LOGGER)
                computation_end_time = time.localtime()
                elapsed_time = datetime.timedelta(seconds=time.mktime(computation_end_time) - time.mktime(start_time))
                log(message='GCparagon (GC-bias computation) Finished Successfully. '
                            f"Elapsed time: {elapsed_time} (h:mm:ss)", log_level=logging.INFO, logger_name=LOGGER)
            # TASK 2: create a GC-weights-tagged version of the input BAM file
            if output_corrected_bam or only_tag_bam:
                start_time = time.localtime()
                if compute_bias:  # add separation line in this case for prettier output
                    log(message="---------------------------------------------------------------------------------",
                        log_level=logging.INFO, logger_name=LOGGER)
                log(message='GCparagon (adding weights to BAM file) Tagging Started.',
                    log_level=logging.INFO, logger_name=LOGGER)
                default_flen_range = range(lower_limit_fragment_length, upper_limit_fragment_length)
                weights_matrix_for_tagging, flen_range, gc_bas_range = reduce_weights_for_tagging(
                    weights_path=correction_weights_matrix_path, mask_path=mask_path, mask_matrix=weights_mask,
                    sample_id=sample_id, correction_matrix=correction_weights_matrix,
                    default_weight=default_fragment_weight, weights_flen_range=default_flen_range,
                    mask_flen_range=default_flen_range)
                tag_bam_with_correction_weights_parallel(
                    bam_path=input_bam, tag_name=gc_tag_name, samtools_path=samtools_path,
                    correction_matrix=weights_matrix_for_tagging, gc_base_limits=gc_bas_range,
                    output_unaligned=output_unaligned_reads,
                    two_bit_genome_file=two_bit_reference_file_path,
                    temporary_directory_sample=sample_temp_dir, threads=total_number_threads,
                    default_fragment_weight=default_fragment_weight,
                    ref_lengths=reference_contig_lengths, frag_len_range=flen_range, sample_output_dir=sample_out_dir)
                # compute duration of execution and message user
                computation_end_time = time.localtime()
                elapsed_time = datetime.timedelta(seconds=time.mktime(computation_end_time) - time.mktime(start_time))
                log(message='GCparagon (adding weights to BAM file) Finished Successfully. '
                            f"Elapsed time: {elapsed_time} (h:mm:ss)", log_level=logging.INFO, logger_name=LOGGER)
            # after BAM has been processed successfully, close logger handlers
            current_logger = logging.getLogger(LOGGER)
            for hdlr in current_logger.handlers:
                hdlr.flush()
                if hdlr in current_logger.handlers:
                    current_logger.removeHandler(hdlr)
                hdlr.close()
        except Exception as e:
            log(message=f"The following exception occurred when trying to process BAM file {input_bam}:\n"
                        f"{create_exception_stack_trace(e)}\n"
                        f"Continuing to processing next BAM file..", log_level=logging.ERROR, logger_name=LOGGER)
        finally:
            # close logging handlers
            final_logger = logging.getLogger(LOGGER)
            for hdlr in final_logger.handlers:
                hdlr.flush()
                if hdlr in final_logger.handlers:
                    final_logger.removeHandler(hdlr)  # local variable 'hdlr' referenced before assignment
                hdlr.close()
            # delete empty log files
            for lg_fl in (log_file, error_log):
                if lg_fl is not None and lg_fl.stat().st_size == 0:
                    lg_fl.unlink()
            continue
    logging.shutdown()
    return 0


if __name__ == "__main__":
    sys.exit(main())
