#!/usr/bin/env python3
import multiprocessing
import time
import datetime
import os
import sys
import math
import shutil
import tempfile
import random
import gc
import logging
import contextlib
import numpy as np
import pandas as pd
import subprocess as sp
from pathlib import Path
import multiprocessing as mp
import multiprocessing.connection as mp_connection
from pysam import AlignmentFile  # coordinates in pysam are always 0-based (following python convention)
from typing import Union, Dict, List, Tuple, Optional, Any
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from collections import deque
from twobitreader import TwoBitFile, TwoBitSequence
from natsort import humansorted

# version
MAJOR_RELEASE = 0
MINOR_RELEASE = 5
PATCH_NUMBER = 4
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

# default definitions for analysis
DEFAULT_MIN_FRAGMENT_LENGTH = 20  # do not set to 0!
DEFAULT_MAX_FRAGMENT_LENGTH = 550
DEFAULT_TAG_NAME = 'GC'
DEFAULT_MIN_OCCURRENCES = 3
ABSOLUTE_MIN_OCCURRENCES = 2
TAGGING_CHUNK_SIZE = 50*10**6
# estimated HW capacities
max_logical_cores = multiprocessing.cpu_count()
max_physical_cores = max_logical_cores // 2 if max_logical_cores > 1 else 1
# PARAMETERS DEFINING RUNTIME OF GC-BIAS COMPUTATION:
# ----------------------------------------------------------------------------------------------------------------------
DEFAULT_NUMBER_PROCESSES = min(12, max_logical_cores)  # limit default value to meaningful amount
DEFAULT_TARGET_NUMBER_FRAGMENTS_PROCESSED = 5*10**6  # 5 million
DEFAULT_PRESET = 1  # showed best results in GC curve comparison against randomly drawn 150 bp sequences from reference
DEFAULT_SIMULATION_REPETITIONS = 6
# ----------------------------------------------------------------------------------------------------------------------
DEFAULT_FLOAT_PRECISION = 6
DEFAULT_FRAGMENT_N_CONTENT_THRESHOLD = 0.3
DEFAULT_MAX_CHUNK_PERCENTAGE_BLACKLIST_OVERLAP = 1/3*100.  # of blacklisted regions for default 1 Mbp chunks processing
DEFAULT_MIN_UNCLIPPED_ALN_FRACTION = 0.75
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
PREDEFINED_1MBP_CHUNKS_TO_PROCESS = SOURCE_CODE_ROOT_PATH.parent.parent / \
                                    'accessory_files/hg38_minimalBlacklistOverlap_1Mbp_chunks_33pcOverlapLimited.bed'
TIMESTAMP_FORMAT = '%Y-%m-%d_%H-%M-%S'
BAD_CHUNKS_FILE_PATTERN_AS_PATH = Path(f'bad_chunks_{TIMESTAMP_FORMAT}.bed')

# define custom types
BadChunksDict = Dict[Tuple[str, int, int], List[int]]
ChunksList = List[Tuple[str, int, int, int]]

# module imports:
from utilities.plot_GC_matrices import plot_statistic_matrices, limit_extreme_outliers, smooth_2d_gc_weights
from utilities.plot_distributions import plot_fragment_length_dists, load_txt_to_matrix_with_meta
from utilities.secure_file_handling import AtomicOpen
from utilities.gc_logging import set_up_logging, log, gib_cmd_logger

LOGGER = gib_cmd_logger()  # basically just to stop linter to assume it is None


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
    input_args.add_argument('-b', '--bam', dest='input_bam', required=True, metavar='File',
                            help='Path to sorted BAM file for which the fragment length-dependent GC-content-based '
                                 "over-representation (= 'GC-bias') should be computed and/or corrected. WARNING: "
                                 "don't use unaligned BAM files (uBAM) or multi-sample/run BAM files! If the BAM's "
                                 'index file is not found on runtime, GCparagon tries to create it. '
                                 '[ PARAMETER REQUIRED ]')
    input_args.add_argument('-rtb', '--two-bit-reference-genome', dest='two_bit_reference_file', required=True,
                            help='Path to 2bit version of the reference genome FastA file which was used for read '
                                 'alignment of the input BAM file. If the 2bit version is missing, one can create the '
                                 'file using the following command: '
                                 "'faToTwoBit <PATH_TO_REF_FASTA> -long <PATH_TO_OUT_2BIT>' "
                                 "(see genome.ucsc.edu/goldenPath/help/twoBit.html for more details)"
                                 "[ PARAMETER REQUIRED ]", metavar='File')
    input_args.add_argument('-c', '--chunks-bed', dest='chunks_bed_file', default=PREDEFINED_1MBP_CHUNKS_TO_PROCESS,
                            help='Path to BED file containing chunks to process. Should be selected based on minimal '
                                 'overlap with bad regions of reference genome build used in creation of --bam. '
                                 f"[ DEFAULT: '{PREDEFINED_1MBP_CHUNKS_TO_PROCESS}' ]", metavar='File')
    input_args.add_argument('-ec', '--exclude-chunks', dest='exclude_chunks_bed_file',
                            help='Path to library file (BED-like) holding DoC-specific definition of bad chunks '
                                 '(chunks must be exact genomic locus match for exclusion, DO NOT expect bedtools '
                                 'intersect-like behavior!). If the bad chunks library is left default, the bad chunks '
                                 'library with the most recent time stamp in the parent directory of the default '
                                 'library/BED file is used. The bad chunks library is intended to speed up the sample '
                                 'processing by excluding chunks with insufficient DoC form the beginning. Excluded '
                                 'chunks were observed to appear most frequently close to centromeres.', metavar='File')
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
    processing_args.add_argument('-up', '--use-parameter-preset', type=int, dest='parameter_preset_number',
                                 choices=range(0, 4, 1), default=DEFAULT_PRESET, metavar='Integer',
                                 help='Optional parameter preset to use for GC bias computation. Must be an integer '
                                      'int the range of 0-3 (inclusive). A preset value of 0 leaves parameters at '
                                      'default if not defined differently by the user (unchanged parameters will match '
                                      'preset 1). Other integer values from 1 to 3 define presets with increasing '
                                      'input data usage and required processing time (expected computation times '
                                      'preset 1-3: 2:40, 15:15, and 50:40 (mm:ss)). Computation time of preset 3 '
                                      'depends on the average DoC of the sample. Average across 4 samples and 3 '
                                      'iterations each computed using 12 cores and the benchmark_mprof.py script. '
                                      'Memory consumption preset 1-3: 340 MiB, 290 MiB, and 300 MiB respectively.'
                                      'If preset is not zero, any customized parameters conflicting with the preset '
                                      'will be ignored. A non-zero preset will set the following parameters: number of '
                                      'simulations, the target number of processed fragments, minimum number of '
                                      'fragment attribute combination occurrences, and the options for outlier '
                                      'detection and smoothing. Noise within the resulting correction weights is '
                                      'reduced when selecting a higher preset value. Preset 3 will attempt to process '
                                      'all genomic chunks (target number of fragments set to 100B) within the limits '
                                      'of the maximum allowed blacklisted regions overlap (per default default ~1.7 Gb '
                                      'of reference are processed).\nNOTE: the percentage of total GC bias corrected '
                                      'fragments in the dataset for presets 1 vs. 3 increases only from 99.837%% to '
                                      '99.938%% (average across 4 samples). Other fragment weights default to 1.0). '
                                      'The primary advantage of processing more fragments is the reduction of noise in '
                                      'computed weights. It is recommended to use a higher preset for a '
                                      "'preprocess-once, analyze often' scenario and/or when a high bias is "
                                      'expected/observed (e.g. FastQC average GC percentage). Correction by preset 1, '
                                      '2, and 3 was found to yield 100.39%%, 99.98%%, and 99,94%% of the raw fragment '
                                      'count respectively (average percentage across 4 samples). '
                                      f'[ DEFAULT: {DEFAULT_PRESET} ]')
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
    processing_args.add_argument('-uf', '--upper-fragment-length', dest='upper_limit_fragment_length', type=int,
                                 default=DEFAULT_MAX_FRAGMENT_LENGTH, metavar='Integer',
                                 help=f'Defines upper length limit for fragments which should be included in '
                                      f'computation. This parameter does not impact computation speed. It only '
                                      f'increases plotting times for matrices by a few seconds. '
                                      f'[ DEFAULT: {DEFAULT_MAX_FRAGMENT_LENGTH}bp ]')
    processing_args.add_argument('-lf', '--lower-fragment-length', dest='lower_limit_fragment_length', type=int,
                                 default=DEFAULT_MIN_FRAGMENT_LENGTH, metavar='Integer',
                                 help=f'Defines lower length limit for fragments which should be included in '
                                      f'computation. Must be positive integer. A value below the sequenceable '
                                      f'fragment length of the device used to create the dataset is not recommended.'
                                      f' [ DEFAULT: {DEFAULT_MIN_FRAGMENT_LENGTH}bp ]')
    processing_args.add_argument('-t', '--threads', dest='total_number_threads', type=int,
                                 default=DEFAULT_NUMBER_PROCESSES, metavar='Integer',
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
                                      "alignment positions and the reference genome) is excluded from the analysis. "
                                      "This parameter was not found to cause any problems for Illumina HiSeq/NovaSeq "
                                      "data. If such fragments have to be included, this flag can be set to allow for "
                                      "up to 1/3 N-bases for fragments. Parameter mainly influences the simulation "
                                      "step and how many times random fragment drawing must be repeated for individual "
                                      "chunks. Also can lead to fewer chunks being discarded (and marked as bad chunk) "
                                      "if flag is set.")
    processing_args.add_argument('-mtb', '--multi-thread-bam', dest='multithread_bam_access',
                                 action='store_true',
                                 help='Optional flag to use 2 threads for BAM read/write access. If activated, the '
                                      'parameter --threads corresponds to half the number of parallel processes that '
                                      'will be spawned for BAM processing. Was not found to have a tremendous impact '
                                      'on I/O performance.')
    processing_args.add_argument('-ucmaf', '--unclipped-min-aln-fraction', dest='min_unclipped_aln_fracton',
                                 default=DEFAULT_MIN_UNCLIPPED_ALN_FRACTION, type=float, metavar='Float',
                                 help='This parameter defines the minimum unclipped fraction of an alignment to be '
                                      'counted in the observed fragment attributes matrix O_gc. This might affect how '
                                      'many small fragments are observed and efectively corrected. [ DEFAULT: '
                                      f'{DEFAULT_MIN_UNCLIPPED_ALN_FRACTION} ]')
    # post-processing options
    postprocessing_args.add_argument('-do', '--detect-outliers', action='store_true', dest='detect_outliers',
                                     help='(PRESET precedence if specified) If this flag is set, extreme outliers will '
                                          'be detected and limited to the extreme outliers threshold value computed '
                                          'from weighted attribute combinations. Default method to detect outliers is '
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
                                     metavar='OutlierDetectionBasis')
    postprocessing_args.add_argument('-ods', '--outlier-detection-stringency', choices=range(0, 8, 1), type=int,
                                     dest='outlier_stringency', default=DEFAULT_OUTLIER_DETECTION_STRINGENCY,
                                     help='(PRESET precedence if specified) If the --detect-outliers flag is set, this '
                                          'parameter defines how stringent the outlier detection threshold is set. '
                                          'Must be an integer in the range of 1-7 (inclusive). [ DEFAULT: '
                                          f'{DEFAULT_OUTLIER_DETECTION_STRINGENCY} ]',  metavar='Integer')
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
                                          f'[ DEFAULT: {DEFAULT_SMOOTHING_KERNEL} ]', metavar='KernelType')
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
                                  'If none is provided, a subdirectory of the input BAM file will be used as '
                                  'output directory. The output for each sample will be gathered in a subdirectory of '
                                  '--out-dir which will be named after the sample. Output directory can be located on '
                                  'slow hardware such as a USB drive or a network storage since everything is '
                                  'stored in --temporary-directory first and moved after completion of computation.')
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
    output_args.add_argument('-k', '--keep-chunk-data', dest='keep_chunk_data', action='store_true',
                             help='Optional flag which can be used to save intermediate data per chunk.')
    output_args.add_argument('-ob', '--output-bam', dest='output_corrected_bam', action='store_true',
                             help='Optional flag to activate writing of the GC-correction-weights-tagged BAM file '
                                  'AFTER COMPUTING GC BIAS (--tag-only flag is not set), either using the statistics '
                                  'computed from the input BAM file or a correction weights matrix specified via '
                                  "--correction-weights. WARNING: currently, the output BAM won't contain "
                                  'unaligned reads!')
    output_args.add_argument('-our', '--output-unaligned-reads', dest='output_unaligned_reads', action='store_true',
                             help='Optional flag to activate writing of unaligned reads to a separate BAM file. '
                                  'Unaligned reads are not included in the tagged BAM output file.  This only has an '
                                  'effect if either the --output-bam flag was also set or the --tag-only mode was '
                                  'started.')
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
    output_args.add_argument('-wce', '--write-chunk-exclusion', dest='write_updated_bad_chunks_library',
                             action='store_true',
                             help='Optional flag for writing an updated version of the library listing chunks marked '
                                  'for exclusion from the analysis. Per default, genomic chunks are marked for '
                                  'exclusion if drawing fragments of a specific size repeatedly fails (at least 33 '
                                  'times or 1/3 of number of fragments that need to be drawn, whichever is higher) due '
                                  'to getting only poly-N sequences. In general, the frequency of these exclusion '
                                  'events is dependent on the DoC of the sample, which can be substituted by the '
                                  'number of fragments estimated to be obtained from all predefined chunks in BAM file '
                                  "in a first approximation. WARNING: don't mix exclusion-marked chunk libraries "
                                  'computed from different (predefined) chunk BED files! If the user places the output '
                                  'BED file library in the default directory, the new library will be used per default '
                                  'for future computations. Chunks will be marked for exclusion depending on a data '
                                  "set's fragment length distribution and sequencing depth.")
    output_args.add_argument('-nfp', '--no-focused-plots', action='store_true', dest='dont_focus_plots',
                             help='Optional flag to deactivate focusing of matrix plots on non-default values (focus '
                                  'uses a border of up to 10 default values). Only has an effect if --no-plots flag is '
                                  'not set.')
    return commandline_parser.parse_args()


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
def read_bed_file(bed_path: str) -> List[Tuple[str, int, int, Any]]:
    """

    :param bed_path:
    :return:
    """
    with AtomicOpen(bed_path, 'rt') as f_bed:
        return [(chrom, int(region_start), int(region_stop), meta_info)
                for chrom, region_start, region_stop, *meta_info
                in filter(lambda x: x not in ('', '\n', None),
                          [bed_line.strip().split() for bed_line in f_bed.readlines()])]


def read_bad_chunks_bed_file(bed_path: str) -> BadChunksDict:
    bad_chunks = {}
    chunk_lengths = set()
    try:
        for chrm, strt, stp, rst in read_bed_file(bed_path=bed_path):
            chunk_lengths.update({stp-strt})
            # values: sample yield(s) as comma-deimited list of integers;
            if bad_chunks.get((chrm, strt, stp)) is None:
                bad_chunks.update({(chrm, strt, stp): list(map(lambda v: int(v), rst[0].split(',')))})
            else:  # chunk locus exists already (= multiple entries! Not expected)
                bad_chunks[(chrm, strt, stp)].extend(list(map(lambda v: int(v), rst[0].split(','))))
    except ValueError:   # not enough values to unpack (expected 4, got X)
        log(message=f"Could not load bad regions from BED file '{bed_path}' (requires column 4 and column 5 to contain "
                    "values that can be cast to int!). Returning no bad regions instead.",
            log_level=logging.WARNING, i_log_with=LOGGER)
    if len(chunk_lengths) > 1:
        log(message=f"Chunks of different length encountered in file '{bed_path}'. "
                    f"A BED file containing mixed bad chunks was provided!",
            log_level=logging.WARNING, i_log_with=LOGGER)
    return bad_chunks


def read_scored_regions_bed_file(bed_path: str):
    try:
        scored_regions = [(chrm, strt, stp, int(rst[0]))
                          for chrm, strt, stp, rst in read_bed_file(bed_path=bed_path)]
        return scored_regions
    except ValueError:  # not enough values to unpack (expected 4, got X)
        log(message=f"Could not load scored regions from BED file '{bed_path}' (requires column 4 to contain values "
                    "that can be cast to int!). Terminating ..", log_level=logging.ERROR, i_log_with=LOGGER)
        sys.exit(1)


def save_matrix_to_txt(matrix: Union[np.array, np.matrix], filename: Union[str, Path], output_dir: str,
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
                i_log_with=LOGGER)
            sys.exit(2)
        formatter_str = '%2.1u'
    if verbose:
        log(message=f"Saving statistic matrix '{filename}' to output directory",
            log_level=logging.INFO, i_log_with=LOGGER)
    np.savetxt(output_path, matrix, delimiter=' |', fmt=formatter_str,
               header=f"rows representing fragment lengths ({min_frag_length} bp to {max_frag_length} bp) and columns "
                      f"representing GC content in rounded percentages (0% to 100%)")
    if report_saved_path:
        return output_path


def round_half_up_int(num: Union[int, float]):
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
                i_log_with=LOGGER)
            raise OutOfGenomicBoundsError()
        return False
    return True


def compute_observed_attributes_matrix(bam_file: str, two_bit_reference_path: str, chromosome: str, tmp_dir: str,
                                       sample_id: str, start_coord: int, stop_coord: int, paired_end_data=True,
                                       save_individual_matrices=False, strict_n_ref_bases_handling=True,
                                       multithread_access=True, float_precision=6,
                                       frag_n_cont_thresh=DEFAULT_FRAGMENT_N_CONTENT_THRESHOLD,
                                       min_frag_len=DEFAULT_MIN_FRAGMENT_LENGTH,
                                       max_frag_len=DEFAULT_MAX_FRAGMENT_LENGTH,
                                       min_unclipped_aln_fracton=DEFAULT_MIN_UNCLIPPED_ALN_FRACTION) \
        -> Tuple[Union[np.array, np.matrix], int, int]:
    """

    :param multithread_access:
    :param two_bit_reference_path:
    :param strict_n_ref_bases_handling:
    :param sample_id:
    :param frag_n_cont_thresh:
    :param tmp_dir:
    :param bam_file:
    :param chromosome:
    :param start_coord:
    :param stop_coord:
    :param paired_end_data:
    :param save_individual_matrices:
    :param float_precision:
    :param min_frag_len:
    :param max_frag_len:
    :param min_unclipped_aln_fracton:
    :return:
    """
    with silently_open_alignment_file(bam_file, mode='rb', threads=2 if multithread_access else 1) as f_aln:
        # Probably the old chunks loading strategy was better (slicing 500k chars is easier than permanent IO)
        observed_attributes_matrix = np.zeros((max_frag_len - min_frag_len + 1, max_frag_len + 1), dtype=np.uint64)
        ref_genome_handle = TwoBitFile(two_bit_reference_path)
        chromosome_handle = ref_genome_handle[chromosome]
        try:  # assert: template has between 1 and 2 segments (single-end or paired-end data; SAM format supports more)
            if paired_end_data:
                # The leftmost segment has a plus sign and the rightmost has a minus sign. It is set as 0 for
                # single-segment template or when the information is unavailable; template length must be positive and
                # within defined range
                filtered_alignments = filter(lambda p: (min_frag_len <= p.template_length <= max_frag_len) and
                                             p.is_proper_pair and not p.is_supplementary and not p.is_secondary,
                                             (read for read in f_aln.fetch(chromosome, start_coord, stop_coord,
                                                                           multiple_iterators=True)))
            else:  # NOT SUPPORTED!
                filtered_alignments = filter(lambda p: (min_frag_len <= p.template_length <= max_frag_len) and
                                             not p.is_supplementary and not p.is_secondary,
                                             (read for read in f_aln.fetch(chromosome, start_coord, stop_coord,
                                                                           multiple_iterators=True)))
        except OSError as e:  # maybe truncated file
            log(message=f"Encountered a problem while trying to parse file '{bam_file}':  {e}. File is likely either "
                        f"truncated or these coordinates were out of bounds: {chromosome}:{start_coord}-{stop_coord}"
                        f". Terminating..", log_level=logging.ERROR, close_handlers=True, i_log_with=LOGGER)
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
                    if frag_length < int(aln_segment.query_length * min_unclipped_aln_fracton):
                        raise IndexError  # min. 3/4 must be identical to ref seq per default or fragment is discarded
                    fragment_sequence = chromosome_handle[frag_start:frag_start + frag_length].upper()
                    gc_count = gc_count_rejecting_n_containing(f_seq=fragment_sequence)  # might return high number
                    observed_attributes_matrix[frag_length-min_frag_len, gc_count] += 1
                except IndexError:  # out of bounds for extreme fragment lengths OR N-containing fragments above thrsh.
                    ignored_fragments += 1
                    continue
        else:  # allow N bases until a certain threshold; extensive looping here!
            for n_fragments_processed, aln_segment in enumerate(filtered_alignments, start=1):
                frag_start = min(aln_segment.pos, aln_segment.pnext)
                frag_length = aln_segment.template_length  # pysam: "the observed query template length"
                # -> filtered for positive
                try:  # below yields str, takes fragment length, computed GC bases and increments counts matrix
                    if frag_length < int(aln_segment.query_length * min_unclipped_aln_fracton):
                        raise IndexError  # min. 3/4 must be identical to ref seq per default or fragment is discarded
                    fragment_sequence = chromosome_handle[frag_start:frag_start + frag_length].upper()
                    gc_count = safe_gc_base_count_inference_thresh(f_seq=fragment_sequence, f_len=frag_length,
                                                                   threshold=frag_n_cont_thresh)
                    observed_attributes_matrix[frag_length-min_frag_len, gc_count] += 1
                except IndexError:  # N-containing fragments above upper threshold
                    ignored_fragments += 1
                    continue
    # give feedback and store the observed attributes matrix
    if ignored_fragments:
        log(message=f"I am done processing {n_fragments_processed:,} fragments for region '{chromosome}:"
                    f"{start_coord:,}-{stop_coord:,}'. Of these, {ignored_fragments:,} fragments did not fall into "
                    f"the selected fragment length interval of [{min_frag_len}, {max_frag_len}] bp.",
            log_level=logging.DEBUG, i_log_with=LOGGER)
    if save_individual_matrices:
        target_path = Path(tmp_dir)
        target_path.mkdir(parents=True, exist_ok=True)
        chunk_str = f"{chromosome}_{start_coord}-{stop_coord}"
        save_matrix_to_txt(matrix=observed_attributes_matrix, output_dir=str(target_path), max_frag_length=max_frag_len,
                           filename=f'{sample_id}_observed_attributes_matrix_{chunk_str}.txt.gz',
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
        return round_half_up_int(f_len / (f_len - n_count) * (f_seq.count('G') + f_seq.count('C')))
    except ZeroDivisionError:  # if all bases returned are 'N's
        return 99999999  # will lead to an IndexError


def simulate_fragment_attributes(two_bit_reference_path: str, tmp_dir: str, chromosome: str, start_coord: int,
                                 sample_id: str, stop_coord: int, expected_yield: int,
                                 statistic_matrix: Union[str, np.ndarray, np.matrix],
                                 float_precision=6, save_individual_matrices=False, strict_n_ref_bases_handling=True,
                                 random_seed=RANDOM_SEED, simulation_repetitions=DEFAULT_SIMULATION_REPETITIONS,
                                 min_frag_len=DEFAULT_MIN_FRAGMENT_LENGTH, max_frag_len=DEFAULT_MAX_FRAGMENT_LENGTH,
                                 frag_n_cont_thresh=DEFAULT_FRAGMENT_N_CONTENT_THRESHOLD) \
        -> Union[Tuple[np.ndarray, Tuple[np.ndarray]], Tuple[str, int, int, int]]:
    """
    May add chunks to bad chunks library
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
    ref_genome_handle = TwoBitFile(two_bit_reference_path)
    # check if chunk coordinates are valid
    try:
        if start_coord < 0 or stop_coord > ref_genome_handle.sequence_sizes()[chromosome]:
            log(message=f"Check your coordinates! Your provided {chromosome}:{start_coord}-{stop_coord}",
                log_level=logging.CRITICAL, close_handlers=True, i_log_with=LOGGER)
            raise OutOfGenomicBoundsError()
    except KeyError:
        log(message=f"The following contig/scaffold was not found in the BAM file: {chromosome}",
            log_level=logging.CRITICAL, close_handlers=True, i_log_with=LOGGER)
        sys.exit(2)
    chromosome_handle = ref_genome_handle[chromosome]
    ref_seq_chunk_slice = chromosome_handle[start_coord:stop_coord].upper()
    if isinstance(statistic_matrix, np.ndarray) or isinstance(statistic_matrix, np.matrix):
        observed_attributes_matrix = statistic_matrix
    elif isinstance(statistic_matrix, str) and Path(statistic_matrix).is_file():
        observed_attributes_matrix, frag_length_range = load_txt_to_matrix_with_meta(statistic_matrix,
                                                                                     loading_logger=LOGGER)
        if frag_length_range.start != min_frag_len:
            log(message="Mismatching fragment length range starts: the observed attributes matrix stated a different "
                        f"fragment length start ({frag_length_range.start}bp) than the internal program logic "
                        f"({min_frag_len} bp). Will use the loaded value.",
                log_level=logging.WARNING, i_log_with=LOGGER)
            min_frag_len = frag_length_range.start
    else:
        log(message="Cannot use provided attribute 'statistic_matrix'. Must be either a path to a saves matrix file or "
                    "the matrix as np.array itself.",
            log_level=logging.CRITICAL, close_handlers=True, i_log_with=LOGGER)
        sys.exit(2)
    bad_chunk = False
    n_fragment_length_rows = observed_attributes_matrix.shape[0]
    s_gc_content_columns = observed_attributes_matrix.shape[1]
    fragments = np.sum(observed_attributes_matrix, axis=1, dtype=np.uint64).flatten()
    # create averaged matrix for in-place manipulation
    simulated_attributes_matrix = np.zeros((n_fragment_length_rows, s_gc_content_columns), dtype=np.uint64)
    # create raw matrices for in-place manipulation
    raw_simulated_matrices = []
    for m_idx in range(simulation_repetitions):
        raw_simulated_matrices.append(np.zeros((n_fragment_length_rows, s_gc_content_columns), dtype=np.uint64))
    random_number_generator = np.random.default_rng(seed=random_seed)  # use random seed for reproducibility here!
    if save_individual_matrices:  # check and create before iterating
        Path(tmp_dir).mkdir(parents=True, exist_ok=True)
    chunk_str = f"{chromosome}_{start_coord}-{stop_coord}"
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
            sampling_failure_threshold = max(int(amount_fragments/3), 55 if strict_n_ref_bases_handling else 33)
            # drawing fragment must fail at least 33 times (55 if strict handling) or one third of all required draws
            # whichever is higher
            unsorted_randoms_index = 0  # index of new fall-back random positions for each fragment length
            unsorted_randoms = random_number_generator.integers(
                low=0, high=stop_coord - start_coord - actual_fragment_length + 1,
                size=sampling_failure_threshold + 1)  # could be rare fragments -> threshold +1 as minimum
            # 1/3 of target fragment number - whichever is higher - for chunk to be marked for exclusion.
            # Approach if not strict_n_ref_bases_handling: linear extrapolation
            # subtract the N bases from the fragment length, compute GC content of the analyzable portion and
            # extrapolate the number of GC-bases according to the fragment's actual length.
            # Still not perfect but the best we can do for a fast algorithm
            # (implemented in 'safe_gc_base_count_inference_thresh()')
            if strict_n_ref_bases_handling:  # faster
                gc_count_iterator = map(lambda q: gc_count_rejecting_n_containing(f_seq=q),
                                        map(lambda s: ref_seq_chunk_slice[s:s + actual_fragment_length].upper(),
                                            rand_ints))
            else:
                gc_count_iterator = map(lambda q: safe_gc_base_count_inference_thresh(f_seq=q[0], f_len=q[1],
                                                                                      threshold=frag_n_cont_thresh),
                                        map(lambda s: (ref_seq_chunk_slice[s:s + actual_fragment_length].upper(),
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
                        except IndexError:  # all random integers consumed -> bad chunk! (should not occur)
                            log(message=f"Too many attempts of drawing random fragments ({sampling_failure_threshold} "
                                        f"attempts) for chunk '{chromosome}:{start_coord}-{stop_coord}' were in vain. "
                                        f"Triggered for fragments of {actual_fragment_length}bp length. "
                                        f"Discarding chunk ..", log_level=logging.WARNING, i_log_with=LOGGER)
                            bad_chunk = True  # stop redrawing once all fall-back random integers have been used up
                            break
                        unsorted_randoms_index += 1  # can max. be 33, then clause below should trigger chunk marking
                        if unsorted_randoms_index >= sampling_failure_threshold:  # should always be the reason why a
                            # chunk is marked as bad
                            log(message=f"Too many attempts of drawing random fragments "
                                        f"for chunk '{chromosome}:{start_coord}-{stop_coord}' were in vain "
                                        f"(threshold was {sampling_failure_threshold} tries). Discarding chunk ..",
                                log_level=logging.WARNING, i_log_with=LOGGER)
                            bad_chunk = True  # stop redrawing once the threshold has been reached
                            break
                        backup_seq = ref_seq_chunk_slice[cur_start:cur_start + actual_fragment_length].upper()
                    if bad_chunk:
                        break  # for loop -> bad chunk encountered!
                    if strict_n_ref_bases_handling:
                        gc_count = gc_count_rejecting_n_containing(f_seq=backup_seq)
                    else:
                        gc_count = safe_gc_base_count_inference_thresh(f_seq=backup_seq, f_len=actual_fragment_length,
                                                                       threshold=frag_n_cont_thresh)
                    raw_simulated_matrices[sim_iter_idx][length_index, gc_count] += 1
            if bad_chunk:
                return chromosome, start_coord, stop_coord, expected_yield
        if save_individual_matrices:  # save snapshot of cumulative simulated_attributes_matrix
            simulated_file = f"{sample_id}_cumulative_simulated_attributes_matrix_{sim_iter_idx + 1}_{chunk_str}.txt.gz"
            save_matrix_to_txt(matrix=simulated_attributes_matrix, filename=simulated_file, output_dir=tmp_dir,
                               float_data_precision=float_precision, gzipped=True, max_frag_length=max_frag_len,
                               min_frag_length=min_frag_len)
    if save_individual_matrices:  # store individual raw matrices
        for raw_matrix_idx, simulated_raw_matrix in enumerate(raw_simulated_matrices):
            simulated_file = f"{sample_id}_simulated_raw_attributes_matrix_{raw_matrix_idx + 1}_{chunk_str}.txt.gz"
            save_matrix_to_txt(matrix=simulated_raw_matrix, filename=simulated_file, output_dir=tmp_dir,
                               float_data_precision=float_precision, gzipped=True, max_frag_length=max_frag_len,
                               min_frag_length=min_frag_len)
    return simulated_attributes_matrix, tuple(raw_simulated_matrices)


def consolidate_results(observed_attributes_matrices_sum: np.array, simulated_attributes_matrices_sum: np.array,
                        simulated_attributes_raw_matrix_sums: List[np.array], n_ogc_summed: int, n_sims: int,
                        n_sgc_summed: int, bam_file: str, tmp_dir: Optional[str], min_frag_len: int, max_frag_len: int,
                        min_frag_occurs: int, ignored_fragments: int, sample_id: str,
                        focus_nondefault_values: Optional[int], precision=6, plot_result=False, output_all=False) \
        -> Tuple[Tuple[Path, np.array],
                 Tuple[Path, np.array],
                 Tuple[range, range]]:
    """

    :param focus_nondefault_values: plots will focus on these values using the integer value of the parameter as border
    :param ignored_fragments:
    :param output_all:
    :param simulated_attributes_raw_matrix_sums:
    :param tmp_dir:
    :param min_frag_occurs:
    :param observed_attributes_matrices_sum:
    :param simulated_attributes_matrices_sum:
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
        report_saved_path=True, max_frag_length=max_frag_len, min_frag_length=min_frag_len)
    if output_all:
        simulated_attributes_matrix_chunk_average = simulated_attributes_matrices_sum / n_sgc_summed  # avg per chunk
        save_matrix_to_txt(matrix=simulated_attributes_matrix_chunk_average, float_data_precision=precision,
                           verbose=True,
                           filename=f'{sample_id}_simulated_attributes_matrix_chunk-average.txt.gz', gzipped=True,
                           output_dir=tmp_dir, max_frag_length=max_frag_len, min_frag_length=min_frag_len)
        for sim_idx, raw_sim_mat in enumerate(simulated_attributes_raw_matrix_sums):
            save_matrix_to_txt(matrix=raw_sim_mat, float_data_precision=precision, verbose=True, output_dir=tmp_dir,
                               filename=f'{sample_id}_simulated_attributes_raw_matrix_{sim_idx}.txt.gz', gzipped=True,
                               max_frag_length=max_frag_len, min_frag_length=min_frag_len)
        observed_attributes_matrix_chunk_average = observed_attributes_matrices_sum / n_ogc_summed  # average per chunk
        save_matrix_to_txt(
            matrix=observed_attributes_matrix_chunk_average, output_dir=tmp_dir, gzipped=True, verbose=True,
            float_data_precision=precision, filename=f'{sample_id}_observed_attributes_matrix_chunk-average.txt.gz',
            max_frag_length=max_frag_len, min_frag_length=min_frag_len)
    # compute correction matrix
    if min_frag_occurs > 1:  # compute mask based on combination occurrences (f_length and GC base count)
        log(message=f"Mask computation: used a minimum threshold of {min_frag_occurs:,} for total "
                    "counts of individual GC-base count-fragment length combinations.",
            log_level=logging.INFO, i_log_with=LOGGER)
    else:
        log(message=f"Masking GC-base-count/fragment-length combinations that were not observed..",
            log_level=logging.INFO, i_log_with=LOGGER)
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
    for f_len in range(min_frag_len, max_frag_len, 1):
        correction_weights_matrix_average[f_len - min_frag_len, f_len+1:] = 0.
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
                                   out_dir_path=Path(tmp_dir), normalize_to_dataset_size=True,
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
                        {f'D_gc_simulation_{s_idx}': pd.DataFrame(sim_mat-observed_attributes_matrices_sum)})
        for data_category in frq_data.keys():
            plot_statistic_matrices(frq_data=frq_data, data_id_to_show=data_category,
                                    y_tick_label_offset=deleted_rows.start + min_frag_len,
                                    x_tick_label_offset=deleted_columns.start,
                                    in_file=bam_file, output_dir=tmp_dir, sample_id=sample_id, fig_width=1800,
                                    fig_height=2000, fig_fontsize=32, parent_logger=LOGGER)
    try:  # give feedback about estimated percentage of corrected fragments:
        included_fragments = int(observed_attributes_matrices_sum.sum())  # without ignored fragments!
        r_fragments_corrected = observed_attributes_matrices_sum[complete_mask].sum() / included_fragments
        log(message=f"Estimated percentage of weighted fragments in dataset: {r_fragments_corrected:.4%} "
                    f"(based on {included_fragments:,} processed fragments included in statistics, which is "
                    f"{included_fragments / (included_fragments + ignored_fragments):.2%} of all processed fragments)",
            log_level=logging.INFO, i_log_with=LOGGER)
    except ZeroDivisionError:  # this should never happen
        log(message=f"Unable to estimate weighted dataset fraction!",
            log_level=logging.WARNING, i_log_with=LOGGER)
    return (correction_weights_matrix_path, correction_weights_matrix_average),\
           (complete_mask_path, complete_mask),\
           (deleted_rows, deleted_columns)


def sort_chunks_by_blacklist_overlap(all_chunks: ChunksList, expected_dataset_fragments: int,
                                     bad_chunks: Optional[BadChunksDict], remove_bad_chunks=True,
                                     max_overlap_percentage=DEFAULT_MAX_CHUNK_PERCENTAGE_BLACKLIST_OVERLAP) \
        -> List[Tuple[str, int, int]]:
    # sort ascending overlapping bases % (normalized to chunk length to account for possible chunk size differences)
    chunks_passing_filters = sorted(all_chunks, key=lambda c: c[3] / (c[2] - c[1]), reverse=False)
    chunks_passing_filters = list(map(lambda t: t[:3],  # discard the blacklist overlap for further processing
                                      filter(lambda c: max_overlap_percentage >= c[3] / (c[2]-c[1]) * 100.,
                                             chunks_passing_filters)))
    log(message=f"{len(all_chunks) - len(chunks_passing_filters):,} chunks were excluded from further "
                f"analysis based on the {max_overlap_percentage:.1f}% chunk overlap with blacklisted regions "
                "threshold.", log_level=logging.INFO, i_log_with=LOGGER)
    if remove_bad_chunks and bad_chunks is not None:  # remove bad chunks from library if any are defined
        # -> select next higher or equal target fragment value
        bad_chunks_for_check = []
        pre_bc_removal_chunks = len(chunks_passing_filters)
        for (chrom, start, stop), sample_frag_yield_list in bad_chunks.items():
            if max(sample_frag_yield_list) >= expected_dataset_fragments:  # ignore bad chunks from higher DoC samples
                bad_chunks_for_check.append((chrom, start, stop))
        chunks_passing_filters = list(filter(lambda c: c not in bad_chunks_for_check,
                                             chunks_passing_filters))
        chunks_after_library_exclusion = len(chunks_passing_filters)
        if pre_bc_removal_chunks != chunks_after_library_exclusion:
            log(message=f"{pre_bc_removal_chunks - chunks_after_library_exclusion:,} chunks were excluded from "
                        "further analysis based on the bad chunks library file. Bad chunks are automatically "
                        "replaced by other chunks.", log_level=logging.INFO, i_log_with=LOGGER)
    return chunks_passing_filters


def gc_bias_worker(chunks_to_process: List[Tuple[str, int, int]], n_sims: int, sender: mp_connection.Connection,
                   min_frag_len: int, max_frag_len: int, two_bit_genome_file: str, chromosome_sizes: Dict[str, int],
                   target_fragment_count: int, mproc_lock: mp.Lock, shared_counter: mp.Value, input_bam: str,
                   sample_id: str, expected_yield: int, precision=6,
                   min_unclipped_aln_fracton=DEFAULT_MIN_UNCLIPPED_ALN_FRACTION, random_seed=RANDOM_SEED, tmp_dir=None,
                   keep_chunk_data=False, strict_n_base_exclusion=True, use_multithreading=True):
    """

    :param use_multithreading:
    :param strict_n_base_exclusion:
    :param sample_id:
    :param expected_yield:
    :param tmp_dir:
    :param target_fragment_count:
    :param two_bit_genome_file:
    :param chromosome_sizes:
    :param input_bam:
    :param mproc_lock:
    :param shared_counter:
    :param chunks_to_process:
    :param n_sims:
    :param sender:
    :param min_frag_len:
    :param max_frag_len:
    :param precision:
    :param random_seed:
    :param keep_chunk_data:
    :param min_unclipped_aln_fracton:
    :return:
    """
    discarded_chunks = 0
    n_observed_gc_matrices = 0
    n_expected_gc_matrices = 0
    observed_attributes_cumulated_matrix = np.zeros((max_frag_len - min_frag_len + 1, max_frag_len + 1),
                                                    dtype=np.uint64)
    simulated_attributes_cumulated_matrix = np.zeros((max_frag_len - min_frag_len + 1, max_frag_len + 1),
                                                     dtype=np.uint64)
    simulated_attributes_cumulated_raw_matrices = []
    for mat_idx in range(n_sims):
        simulated_attributes_cumulated_raw_matrices.append(np.zeros((max_frag_len - min_frag_len + 1, max_frag_len + 1),
                                                                    dtype=np.uint64))
    fragments_processed = 0
    fragments_ignored = 0
    individual_matrices_tmp_dir = str(Path(tmp_dir) / 'data_per_chunk') if keep_chunk_data else None
    bad_chunks_list = []
    # process chunks either until global number of processed fragments suffices or we processed all chunks received here
    for chunk_chrom, chunk_start, chunk_end in chunks_to_process:
        check_region_within_genomic_bounds(chrom=chunk_chrom, start=chunk_start, stop=chunk_end,
                                           chromosome_sizes=chromosome_sizes, raise_error=True)
        gc_matrix_observed, chunk_fragments_processed, frag_ignored = compute_observed_attributes_matrix(
            two_bit_reference_path=two_bit_genome_file, tmp_dir=individual_matrices_tmp_dir, bam_file=input_bam,
            start_coord=chunk_start, stop_coord=chunk_end, save_individual_matrices=keep_chunk_data,
            sample_id=sample_id, float_precision=precision, min_frag_len=min_frag_len, chromosome=chunk_chrom,
            max_frag_len=max_frag_len, strict_n_ref_bases_handling=strict_n_base_exclusion,
            multithread_access=use_multithreading, min_unclipped_aln_fracton=min_unclipped_aln_fracton)
        # in-place manipulate the expected matrix
        simulated_attributes_matrix, raw_simulated_matrices = simulate_fragment_attributes(
            two_bit_reference_path=two_bit_genome_file, tmp_dir=individual_matrices_tmp_dir, min_frag_len=min_frag_len,
            statistic_matrix=gc_matrix_observed, chromosome=chunk_chrom, sample_id=sample_id, start_coord=chunk_start,
            stop_coord=chunk_end, simulation_repetitions=n_sims, save_individual_matrices=keep_chunk_data,
            random_seed=random_seed, float_precision=precision, expected_yield=expected_yield,
            strict_n_ref_bases_handling=strict_n_base_exclusion, max_frag_len=max_frag_len)
        if isinstance(simulated_attributes_matrix, tuple) and len(simulated_attributes_matrix) == 4:  # got a bad chunk
            # -> discard entire chunk info
            bad_chunk = simulated_attributes_matrix
            bad_chunks_list.append(bad_chunk)
            discarded_chunks += 1
            continue
        observed_attributes_cumulated_matrix += gc_matrix_observed
        fragments_processed += chunk_fragments_processed - frag_ignored
        fragments_ignored += frag_ignored
        n_observed_gc_matrices += 1
        simulated_attributes_cumulated_matrix += simulated_attributes_matrix
        for mat_idx in range(n_sims):
            simulated_attributes_cumulated_raw_matrices[mat_idx] += raw_simulated_matrices[mat_idx]
        n_expected_gc_matrices += n_sims
        with mproc_lock:
            shared_counter.value += chunk_fragments_processed
            if shared_counter.value > target_fragment_count:
                break  # end iteration and send results
    sender.send((observed_attributes_cumulated_matrix, simulated_attributes_cumulated_matrix,
                 simulated_attributes_cumulated_raw_matrices, n_observed_gc_matrices,
                 n_expected_gc_matrices, fragments_processed, fragments_ignored, discarded_chunks, bad_chunks_list))
    sender.close()


def compute_gc_bias_parallel(chunks_to_process: List[Tuple[str, int, int]], threads: int, simulation_count: int,
                             min_flen: int, max_flen: int, out_dir_sample: str, in_bam: str, sample_name: str,
                             two_bit_reference_file: str, chrom_sizes: Dict[str, int], min_frag_occurs: int,
                             target_fragments_processed: int, expected_yield: int, tmp_dir_sample: str,
                             bad_chunks_library_file: Optional[str], plot_focus_border: Optional[int],
                             float_precision=6, strict_n_base_exclusion=True, keep_chunk_data=False,
                             visualize_matrices=False, output_all=False, write_updated_bad_chunks_library=True,
                             use_multithreading=True, detect_outliers=True, focus_custom_values=True,
                             outlier_detection_method='IQR', outlier_detection_stringency=2, smooth_weights=True,
                             smoothing_kernel='gauss', smoothing_intensity=2,
                             min_unclipped_aln_fracton=DEFAULT_MIN_UNCLIPPED_ALN_FRACTION, random_seed=RANDOM_SEED) \
        -> Tuple[np.array, np.array]:
    """
    :param focus_custom_values:
    :param plot_focus_border:
    :param use_multithreading:
    :param strict_n_base_exclusion:
    :param smooth_weights:
    :param smoothing_kernel:
    :param smoothing_intensity:
    :param outlier_detection_stringency:
    :param outlier_detection_method:
    :param detect_outliers:
    :param output_all:
    :param bad_chunks_library_file:
    :param write_updated_bad_chunks_library:
    :param expected_yield:
    :param tmp_dir_sample:
    :param min_frag_occurs:
    :param target_fragments_processed:
    :param chrom_sizes:
    :param two_bit_reference_file:
    :param chunks_to_process:
    :param threads:
    :param simulation_count:
    :param min_flen:
    :param max_flen:
    :param out_dir_sample:
    :param in_bam:
    :param sample_name:
    :param float_precision:
    :param keep_chunk_data:
    :param random_seed:
    :param visualize_matrices:
    :param min_unclipped_aln_fracton: minimum fraction of unclipped alignment length for fragment to be counted in O_gc
    :return:
    """
    n_processes = int(threads // 1.5 if use_multithreading else threads)  # observed mostly 1 thread active per process
    # -> multiprocessing gives better performance than threading
    # split up sorted chunks in lists for workers
    chunks_for_workers = [[] for _i in range(n_processes)]
    for c_idx, cur_chunk in enumerate(chunks_to_process):
        chunks_for_workers[c_idx % n_processes].append(cur_chunk)
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
                    f"'{tmp_dir_sample}' ..", log_level=logging.WARNING, i_log_with=LOGGER)
        shutil.rmtree(tmp_sample_output_path)
    tmp_sample_output_path.mkdir(parents=True, exist_ok=True)
    # create shared value and lock for counting processed fragments
    shared_fragment_counter = mp.Value('i', 0)  # integer value, initialized as zero
    multiproc_lock = mp.Lock()
    all_worker_kwargs = [{'chunks_to_process': chunks_for_workers[w_idx], 'n_sims': simulation_count,
                          'sender': senders[w_idx], 'min_frag_len': min_flen, 'max_frag_len': max_flen,
                          'keep_chunk_data': keep_chunk_data, 'random_seed': random_seed,
                          'precision': float_precision, 'two_bit_genome_file': two_bit_reference_file,
                          'chromosome_sizes': chrom_sizes, 'target_fragment_count': target_fragments_processed,
                          'mproc_lock': multiproc_lock, 'shared_counter': shared_fragment_counter, 'input_bam': in_bam,
                          'tmp_dir': tmp_dir_sample, 'expected_yield': expected_yield, 'sample_id': sample_name,
                          'strict_n_base_exclusion': strict_n_base_exclusion, 'use_multithreading': use_multithreading,
                          'min_unclipped_aln_fracton': min_unclipped_aln_fracton}
                         # write to TMP first, then move files which should be kept (is much faster on cluster!)
                         for w_idx in range(n_processes)]
    gc_bias_workers = [mp.Process(target=gc_bias_worker, kwargs=worker_kwargs)
                       for worker_kwargs in all_worker_kwargs]
    # start all workers, receive results, wait until finished, and close them
    observed_attributes_matrices_sum = np.zeros((max_flen - min_flen + 1, max_flen + 1), dtype=np.uint64)
    simulated_attributes_matrices_sum = np.zeros((max_flen - min_flen + 1, max_flen + 1), dtype=np.uint64)
    simulated_attributes_raw_matrix_sums = []
    for _s_idx in range(simulation_count):
        simulated_attributes_raw_matrix_sums.append(np.zeros((max_flen - min_flen + 1, max_flen + 1), dtype=np.uint64))
    n_summed_observed_attributes_matrices = 0
    n_summed_s_gc_matrices = 0
    total_ignored_fragments = 0
    total_discarded_chunks = 0
    all_bad_chunks = {}
    deque(map(lambda w: w.start(), gc_bias_workers), maxlen=0)
    # if memory profiler is applied to code and processes are spawned rather than forked, we get a PicklingError her:
    # "Can't pickle <function gc_bias_worker at 0xsomething>: attribute lookup gc_bias_worker on __main__ failed"
    received_data = [incoming_result.recv() for incoming_result in receivers]  # collect self-terminating worker returns
    for o_gc_sum_mat, s_gc_sum_mat, s_gc_sum_raw_mats, n_observed, n_expected, _fragments_processed, \
            ignored_fragments, discarded_chunks, list_of_bad_chunks in received_data:
        observed_attributes_matrices_sum += o_gc_sum_mat  # inplace np.array manipulation
        simulated_attributes_matrices_sum += s_gc_sum_mat  # inplace np.array manipulation
        for mat_idx, s_gc_sum_raw_mat in enumerate(s_gc_sum_raw_mats):
            simulated_attributes_raw_matrix_sums[mat_idx] += s_gc_sum_raw_mat
        n_summed_observed_attributes_matrices += n_observed
        n_summed_s_gc_matrices += n_expected
        total_ignored_fragments += ignored_fragments
        total_discarded_chunks += discarded_chunks
        for chrm, strt, stp, rst in list_of_bad_chunks:
            if all_bad_chunks.get((chrm, strt, stp)) is None:
                all_bad_chunks.update({(chrm, strt, stp): [int(rst[1])]})
            else:  # chunk locus already exists (= multiple entries! Not expected)
                all_bad_chunks[(chrm, strt, stp)].append(int(rst[1]))
    # create user feedback
    log(message="---------------------------------------------------------------------------------",
        log_level=logging.INFO, i_log_with=LOGGER)
    if shared_fragment_counter.value < target_fragments_processed:
        if shared_fragment_counter.value == 0:
            log(message="No fragments were processed. Check your environment for samtools and pysam! Terminating ..",
                log_level=logging.CRITICAL, close_handlers=True, i_log_with=LOGGER)
            sys.exit(3)
        log(message=f"Could not reach target count of processed fragments: {shared_fragment_counter.value:,}/"
                    f"{target_fragments_processed:,} fragments included in statistics. GC-bias weights matrix "
                    "might be subject to stochastic noise!" + (" (see plot if created)" if visualize_matrices else ''),
            log_level=logging.WARNING, i_log_with=LOGGER)
    else:
        log(message=f"Target number of fragments to process reached: {shared_fragment_counter.value:,}/"
                    f"{target_fragments_processed:,} (incl. ignored fragments).",
            log_level=logging.INFO, i_log_with=LOGGER)
    if total_discarded_chunks:
        log(message=f"Number of discarded chunks that contained too many fragments with N-bases: "
                    f"{total_discarded_chunks:,}", log_level=logging.INFO, i_log_with=LOGGER)
    log(message=f"In total, {n_summed_observed_attributes_matrices:,} "
                "genomic chunks were included in GC-bias computation.", log_level=logging.INFO, i_log_with=LOGGER)
    if total_ignored_fragments:
        log(message=f"Total number of fragments with length out of bounds: {total_ignored_fragments:,}\n"
                    f"Effective number of fragments included in statistics: "
                    f"{shared_fragment_counter.value - shared_fragment_counter.value:,}",
            log_level=logging.INFO, i_log_with=LOGGER)
        if total_ignored_fragments / 50000 > 1.:
            log(message=f"More than 50k fragments with length out of bounds detected! You may want to re-run this "
                        f"analysis with a higher maximum fragment length!",
                log_level=logging.WARNING, i_log_with=LOGGER)
    _ = deque(map(lambda w: w.join(), gc_bias_workers), maxlen=0)  # wait until all processes have finished
    _ = deque(map(lambda w: w.close(), gc_bias_workers), maxlen=0)
    # update bad chunks library
    if write_updated_bad_chunks_library:
        if not all_bad_chunks:
            log(message="No bad chunks encountered. No new bad chunks library will be output.",
                log_level=logging.INFO, i_log_with=LOGGER)
        else:
            # re-read bad chunks (might have changed in the meantime by other processes;
            # include also bad chunks form lower-expected fragment count datasets in new library!)
            if bad_chunks_library_file is not None:
                curren_bad_chunks = read_bad_chunks_bed_file(bed_path=bad_chunks_library_file)
                for k, yield_list in curren_bad_chunks.items():
                    if all_bad_chunks.get(k) is None:
                        all_bad_chunks.update({k: yield_list})
                    else:
                        all_bad_chunks[k].extend(yield_list)  # append received sample yield list to existing list
            newline = '\n'
            # create lines for (possibly) combined bad chunks
            new_lib_buffer = [f"{newline if line_idx else ''}{chrm}\t{strt}\t{stp}\t{max(yld_lst)}"
                              for line_idx, ((chrm, strt, stp), yld_lst) in
                              enumerate(humansorted(all_bad_chunks.items(),
                                                    key=lambda t: (t[0][0], t[0][1]),
                                                    reverse=False))]
            new_lib_path = Path(out_dir_sample) / \
                f'bad_chunks_{time.strftime(TIMESTAMP_FORMAT, time.localtime())}.bed'
            with AtomicOpen(new_lib_path, 'wt') as f_new_bclib:
                f_new_bclib.writelines(new_lib_buffer)
            log(message=f"Updated bad chunks library BED file written to: '{new_lib_path}'",
                log_level=logging.INFO, i_log_with=LOGGER)
    # combine results
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
                                                       outliers_factor=10-outlier_detection_stringency,
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
                                    fig_fontsize=32, parent_logger=LOGGER)
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
                                    fig_fontsize=32, parent_logger=LOGGER)
    # move all output from temporary sample dir into output dir after checking that the latter does not exist
    target_path = Path(out_dir_sample)
    target_path.mkdir(parents=True, exist_ok=True)  # ensure parent path exists to be able to move sample dir
    log(message=f"Copying GC-bias computation output from '{tmp_dir_sample}' to '{out_dir_sample}' ..",
        log_level=logging.INFO, i_log_with=LOGGER)
    source_path = Path(tmp_dir_sample)
    for f in source_path.iterdir():  # handle all files/dirs separately
        # -> don't just batch delete existing target dir! Will exist probably due to GC-bias computation step
        if f.is_dir():  # subdir -> can be completely moved & target completely deleted
            # (will only concern the bam parts if they were kept in a previous analysis)
            target_dir_path = target_path / f.name
            if target_dir_path.exists():
                log(message=f"Target path for directory '{f.name}' exists and will be completely deleted! "
                            f"Deleting '{target_dir_path}' ..", log_level=logging.WARNING, i_log_with=LOGGER)
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
            log_level=logging.WARNING, i_log_with=LOGGER)
    moved_weights_mask_matrix_path = target_path / weights_mask_path.name
    if not moved_weights_mask_matrix_path.is_file():
        log(message=f"Moved correction matrix not found at '{moved_weights_mask_matrix_path}'",
            log_level=logging.WARNING, i_log_with=LOGGER)
    return (moved_correction_matrix_path, use_correction_matrix), (moved_weights_mask_matrix_path, weights_mask)


def load_until_leftmost_not_poly_n(loaded_sequences: List[Optional[str]], loaded_chunks: List[bool], chrom_seq: str,
                                   chunk_loading_size=5000000, max_loaded=10, trigger_garbage_collection=False) \
        -> Tuple[List[Optional[str]], List[bool]]:
    if max_loaded < 2:
        log(message=f"Invalid choice for maximum number of loaded genomic chunks! Must be at least 2. Setting it to "
                    f"default of 10 ..",
            log_level=logging.WARNING, close_handlers=True, i_log_with=LOGGER)
        max_loaded = 10
    first_loaded_sequence_index = loaded_chunks.index(True)
    while loaded_sequences[first_loaded_sequence_index].count('N') == chunk_loading_size:
        loaded_sequences[first_loaded_sequence_index] = None  # unload leftmost poly-N chunk
        loaded_chunks[first_loaded_sequence_index] = False
        if first_loaded_sequence_index >= len(loaded_sequences) - 2:  # add rightmost reference chunk(s) if < 2 loaded
            for i in range(first_loaded_sequence_index - len(loaded_sequences) + 3):  # load >=2 non-poly-N chunks
                loaded_sequences.append(chrom_seq[chunk_loading_size * len(loaded_chunks):
                                                  chunk_loading_size * (len(loaded_chunks) + 1)])
                loaded_chunks.append(True)
        first_loaded_sequence_index = loaded_chunks.index(True)
    # check for too many loaded chunks
    if sum(loaded_chunks) > max_loaded:
        remove_n = sum(loaded_chunks) - max_loaded
        # unload chunk with most neighboring not loaded entries; if tied, unload leftmost of these chunks
        for _i in range(remove_n):
            unloaded_neighbors = [0] * len(loaded_chunks)
            inverted_loaded_chunks = [not entr for entr in loaded_chunks]
            for pos in range(len(loaded_chunks)):
                unloaded_neighbors[pos] = sum(inverted_loaded_chunks[pos+1:]) + sum(inverted_loaded_chunks[:pos])
            idx_chunk_to_remove = unloaded_neighbors.index(max(unloaded_neighbors))
            loaded_sequences[idx_chunk_to_remove] = None
            loaded_chunks[idx_chunk_to_remove] = False
    if trigger_garbage_collection:
        gc.collect()
    return loaded_sequences, loaded_chunks


def extend_chunks_for_index(target_index: int, scaffold_length: int, scaffold_name: str,
                            loaded_sequences: List[Optional[str]], loaded_chunks: List[bool],
                            chrom_handle: TwoBitSequence, parent_logger: logging.Logger, max_loaded=10,
                            chunk_loading_size=5000000, trigger_garbage_collection=False) \
        -> Tuple[List[Optional[str]], List[bool]]:
    """

    :param target_index:
    :param scaffold_length:
    :param scaffold_name:
    :param loaded_sequences:
    :param loaded_chunks:
    :param chrom_handle:
    :param parent_logger:
    :param max_loaded:
    :param chunk_loading_size:
    :param trigger_garbage_collection:
    :return:
    :raises: AttributeError if scaffold length is smaller than starting coordinate of chunk which should be loaded
    """
    last_idx = target_index - len(loaded_chunks)
    for l_idx in range(target_index - len(loaded_chunks) + 1):
        if l_idx == last_idx:  # actually load the sequence
            hypothetical_end_coordinate = chunk_loading_size * (len(loaded_chunks) + 1)
            if hypothetical_end_coordinate > scaffold_length:
                if scaffold_length < chunk_loading_size * len(loaded_chunks):  # should never occur - illogical error
                    log(message=f"Critical error encountered loading chunk  {target_index} for scaffold {scaffold_name}"
                                f". Length of scaffold {scaffold_length:,}bp was smaller than starting coordinate "
                                f"{chunk_loading_size * len(loaded_chunks)}bp. Terminating..",
                        log_level=logging.CRITICAL, flush=True, i_log_with=parent_logger)
                    raise AttributeError
                loaded_sequences.append(chrom_handle[chunk_loading_size * len(loaded_chunks):scaffold_length].upper())
            else:
                loaded_sequences.append(chrom_handle[chunk_loading_size * len(loaded_chunks):
                                                     chunk_loading_size * (len(loaded_chunks) + 1)].upper())
            loaded_chunks.append(True)
        else:  # just elongate the list without loading anything
            loaded_sequences.append(None)
            loaded_chunks.append(False)
    # check for too many loaded chunks
    if sum(loaded_chunks) > max_loaded:
        remove_n = sum(loaded_chunks) - max_loaded
        # unload chunk with most neighboring not loaded entries; if tied, unload leftmost of these chunks
        for _i in range(remove_n):
            unloaded_neighbors = [0] * len(loaded_chunks)
            inverted_loaded_chunks = [not entry for entry in loaded_chunks]
            for pos in range(len(loaded_chunks)):
                unloaded_neighbors[pos] = sum(inverted_loaded_chunks[pos + 1:]) + sum(inverted_loaded_chunks[:pos])
            idx_chunk_to_remove = unloaded_neighbors.index(max(unloaded_neighbors))
            # unload leftmost chunk with the most not-loaded neighbors
            loaded_sequences[idx_chunk_to_remove] = None
            loaded_chunks[idx_chunk_to_remove] = False
    if trigger_garbage_collection:
        gc.collect()
    return loaded_sequences, loaded_chunks


def load_specific_chunk(target_index: int, scaffold_length: int, scaffold_name: str,
                        loaded_sequences: List[Optional[str]], loaded_chunks: List[bool],
                        chrom_handle: TwoBitSequence, parent_logger: logging.Logger, max_loaded=10,
                        chunk_loading_size=5000000, trigger_garbage_collection=False) \
        -> Tuple[List[Optional[str]], List[bool]]:
    """
    Use only if target_index lies within len(laoded_seqeunces) - 1!
    :param parent_logger:
    :param scaffold_name:
    :param target_index:
    :param scaffold_length:
    :param loaded_sequences:
    :param loaded_chunks:
    :param chrom_handle:
    :param max_loaded:
    :param chunk_loading_size:
    :param trigger_garbage_collection:
    :return:
    """
    if target_index >= len(loaded_sequences):
        log(message=f"Cannot load chunk {target_index} for scaffold {scaffold_name}. Terminating..",
            log_level=logging.CRITICAL, flush=True, close_handlers=True, i_log_with=parent_logger)
        raise AttributeError
    if chunk_loading_size*(target_index + 1) > scaffold_length:
        loaded_sequences[target_index] = chrom_handle[chunk_loading_size * target_index:scaffold_length].upper()
    else:
        loaded_sequences[target_index] = chrom_handle[chunk_loading_size * target_index:
                                                      chunk_loading_size * (target_index + 1)].upper()
    loaded_chunks[target_index] = True
    # check for too many loaded chunks
    if sum(loaded_chunks) > max_loaded:
        remove_n = sum(loaded_chunks) - max_loaded
        # unload chunk with most neighboring not-loaded-entries; if tied, unload leftmost of these chunks
        for _i in range(remove_n):
            unloaded_neighbors = [0] * len(loaded_chunks)
            inverted_loaded_chunks = [not entr for entr in loaded_chunks]
            for pos in range(len(loaded_chunks)):
                unloaded_neighbors[pos] = sum(inverted_loaded_chunks[pos + 1:]) + sum(inverted_loaded_chunks[:pos])
            idx_chunk_to_remove = unloaded_neighbors.index(max(unloaded_neighbors))
            if target_index == idx_chunk_to_remove:  # must not be target index!
                dummy = unloaded_neighbors[:]
                dummy[idx_chunk_to_remove] = -1  # set to lowest
                idx_chunk_to_remove = dummy.index(max(dummy))  # find next maximum
            loaded_sequences[idx_chunk_to_remove] = None
            loaded_chunks[idx_chunk_to_remove] = False
    if trigger_garbage_collection:
        gc.collect()
    return loaded_sequences, loaded_chunks


def unaligned_bam_worker(bam_path: Union[str, Path], output_path: Union[str, Path], tag_name: str):
    with AlignmentFile(bam_path, mode='rb') as f_in:
        with AlignmentFile(output_path, header=f_in.header, mode='wb') as f_unaligned_tagged:
            for aln in f_in.fetch(until_eof=True, multiple_iterators=True):
                if not aln.is_mapped:
                    aln.set_tag(tag=tag_name, value=0.)  # unaligned reads get a GC correction weight of 0 because the
                    # fragment sequence cannot be safely inferred from the read have
                    f_unaligned_tagged.write(aln)


def bam_tagging_worker_single_chunk(bam_path: str, correction_weights: np.array, temp_dir: str, gc_bases_offset: int,
                                    fragment_length_range: range, two_bit_reference_path: str,
                                    tagging_chunks_list: List[str], reference_lengths: Dict[str, int],
                                    sender_connection: mp_connection.Connection, parent_logger: logging.Logger,
                                    tag_name=DEFAULT_TAG_NAME, ref_chunk_loading_size=500000, annotation_str=None,
                                    use_multithreading=True):
    """

    :param tagging_chunks_list: (chromosome, chunk_start, chunk_end, chunk_size)
    :param gc_bases_offset:
    :param parent_logger:
    :param use_multithreading:
    :param two_bit_reference_path:
    :param annotation_str:
    :param bam_path:
    :param correction_weights:
    :param temp_dir:
    :param fragment_length_range:
    :param reference_lengths: length of reference chromosomes and other scaffolds
    :param sender_connection:
    :param tag_name:
    :param ref_chunk_loading_size:
    :return:
    """
    # get chromosome sequence handle
    reference_handle = TwoBitFile(two_bit_reference_path)
    # (chunk is sliced again subsequently for faster access)
    min_frag_len = fragment_length_range.start
    # compute corrected alignments -> use buffer of 200000 entries
    tagged_bam_files = []
    with silently_open_alignment_file(bam_path, mode='rb', threads=2 if use_multithreading else 1) as input_bam_file:
        for c_idx, (chromosome, start_coord, stop_coord, _ch_len) in enumerate(tagging_chunks_list):
            scaffold_length = reference_lengths[chromosome]
            chromosome_handle = reference_handle[chromosome]  # is large; just slice for chunk sequence retrieval
            tagged_bam_file = str(Path(temp_dir) /
                                  '.'.join(Path(bam_path).name.split('.')[:-2] +
                                           [f"{Path(bam_path).name.split('.')[-2]}"  # anno str will be None
                                            f"+{chromosome}-{start_coord}-{stop_coord}.GCcorr"] +
                                           ([annotation_str, 'bam'] if annotation_str else ['bam'])))
            with silently_open_alignment_file(tagged_bam_file, mode='wb', template=input_bam_file,
                                              threads=2 if use_multithreading else 1) as f_gc_tagged:
                # preload first 5 Mbp of reference sequence; use 2 lists, one stores the sequences, the other stores
                # loading status of each chunk as boolean value (chunks are loaded consecutively without gaps)
                if ref_chunk_loading_size > scaffold_length:
                    ref_chunk_sequences = [chromosome_handle[0:scaffold_length].upper()]
                    loaded_ref_chunks = [True]
                elif scaffold_length < 2 * ref_chunk_loading_size:
                    ref_chunk_sequences = [chromosome_handle[0:ref_chunk_loading_size].upper(),
                                           chromosome_handle[ref_chunk_loading_size:scaffold_length].upper()]
                    loaded_ref_chunks = [True, True]
                else:
                    ref_chunk_sequences = [chromosome_handle[0:ref_chunk_loading_size].upper(),
                                           chromosome_handle[ref_chunk_loading_size:2 * ref_chunk_loading_size].upper()]
                    loaded_ref_chunks = [True, True]
                # fastforward until no N-contigs are in ref_contig_chunks-deque any more
                try:
                    ref_chunk_sequences, loaded_ref_chunks = load_until_leftmost_not_poly_n(
                        loaded_sequences=ref_chunk_sequences, loaded_chunks=loaded_ref_chunks, max_loaded=10,
                        chrom_seq=chromosome_handle, chunk_loading_size=ref_chunk_loading_size)
                except AttributeError:
                    log(message=f"An error occurred when trying to load initial chunks for tagging. Cannot continue "
                                f"processing scaffold '{chromosome}'. Exiting BAM tagging worker..",
                        log_level=logging.ERROR, i_log_with=parent_logger)
                    return -1
                # iterate over alignments from specified reference contigs
                aln_buffer = []
                for aln_seg in input_bam_file.fetch(chromosome, start_coord, stop_coord, multiple_iterators=True):
                    if not aln_seg.is_mapped or aln_seg.pos < start_coord or aln_seg.pos >= stop_coord:
                        continue  # skip alns outside of chunk borders and unaligned reads
                    try:
                        frag_start_scaffold = min(aln_seg.reference_start, aln_seg.next_reference_start)
                        frag_size = abs(aln_seg.template_length)
                        correction_weights_row = correction_weights[frag_size - min_frag_len]
                    except TypeError:  # unaligned: '<' not supported between instances of 'NoneType' and 'int'
                        aln_seg.set_tag(tag_name, value=0.,  # give unaligned reads a GC weight of 0.
                                        value_type="f", replace=True)  # should not occur when fetching chunks
                        aln_buffer.append(aln_seg)  # no need to check aln buffer; will be written in non-default case
                        continue
                    except IndexError:  # fragment length not in reduced weight matrix -> use default value
                        aln_seg.set_tag(tag_name, value=1.,  # give unaligned reads a GC weight of 0.
                                        value_type="f", replace=True)
                        aln_buffer.append(aln_seg)  # no need to check aln buffer; will be written in non-default case
                        continue
                    if frag_size <= abs(aln_seg.query_alignment_length):  # get sequence from aligned read portion
                        frag_seq = aln_seg.query_alignment_sequence
                    else:  # get fragment sequence from reference genome
                        target_chunk_index = frag_start_scaffold // ref_chunk_loading_size
                        if target_chunk_index >= len(loaded_ref_chunks):  # required chunk was not loaded -> load it!
                            try:
                                ref_chunk_sequences, loaded_ref_chunks = extend_chunks_for_index(
                                    target_index=target_chunk_index, loaded_sequences=ref_chunk_sequences,
                                    parent_logger=parent_logger, loaded_chunks=loaded_ref_chunks,
                                    chrom_handle=chromosome_handle, scaffold_name=chromosome, max_loaded=10,
                                    chunk_loading_size=ref_chunk_loading_size, scaffold_length=scaffold_length)
                            except AttributeError:
                                log(message=f"Error occurred when trying to load a downstream chunk. Cannot continue "
                                            f"processing scaffold '{chromosome}'. Exiting BAM tagging worker..",
                                    log_level=logging.ERROR, i_log_with=parent_logger)
                                return -1
                        elif not loaded_ref_chunks[target_chunk_index]:  # chunk in range but was unloaded/not loaded
                            try:
                                ref_chunk_sequences, loaded_ref_chunks = load_specific_chunk(  # shouldn't be necessary
                                    target_index=target_chunk_index, loaded_sequences=ref_chunk_sequences,
                                    loaded_chunks=loaded_ref_chunks, chrom_handle=chromosome_handle, max_loaded=10,
                                    scaffold_name=chromosome, chunk_loading_size=ref_chunk_loading_size,
                                    scaffold_length=scaffold_length, parent_logger=parent_logger)
                            except AttributeError:
                                log(message="Error occurred when trying to load an unloaded/intermediate chunk. Cannot"
                                            f"continue processing scaffold '{chromosome}'. "
                                            "Exiting BAM tagging worker..",
                                    log_level=logging.ERROR, i_log_with=parent_logger)
                                return -1
                        frag_start_chunk = frag_start_scaffold % ref_chunk_loading_size
                        frag_seq = ref_chunk_sequences[target_chunk_index][
                                          frag_start_chunk:frag_start_chunk + frag_size]
                    try:  # to retrieve correction weight
                        corr_factor = correction_weights_row[frag_seq.count('G') + frag_seq.count('C') +
                                                             gc_bases_offset]
                    except IndexError:  # also: for not-in-proper-pair alns..
                        corr_factor = 1.  # ..(across breakpoints) or alns mapping far apart due to deletion
                        # also for fragment lengths + GC bases pointing to corr-factor of 1. (weights matrix is trimmed)
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
        log(message=message, log_level=logging.ERROR, i_log_with=LOGGER)
    if Path(tmp_dir).is_dir():
        shutil.rmtree(tmp_dir, ignore_errors=True)  # remove full directory with content -> use shutil rm_tree
    if LOGGER:
        for hdlr in LOGGER.handlers:
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
    Remove non-informative weights so that numpy matrix is smaller (faster access?)
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


def get_unaligned_reads(bam_path: Union[str, Path], output_dir: Union[str, Path], tag_name: str) \
        -> Tuple[Optional[Path], mp.Process]:
    sample_id = Path(bam_path).stem
    output_bam_unaligned = Path(output_dir) / f'{sample_id}.unaligned.bam'
    output_dir.mkdir(parents=True, exist_ok=True)  # ensure output directory exists
    unaligned_process_handle = mp.Process(target=unaligned_bam_worker,
                                          kwargs={'bam_path': bam_path, 'output_path': output_bam_unaligned,
                                                  'tag_name': tag_name})
    unaligned_process_handle.start()
    return output_bam_unaligned, unaligned_process_handle


def samtools_cat_bams(list_of_bams: List[str], samtools_path: Union[str, Path],
                      tmp_dir: Union[str, Path], output_bam: Path, keep_input=False):
    concatenation_command = [str(samtools_path), 'cat', '-o', output_bam, '--no-PG', '--threads', '4'] + list_of_bams
    sp.call(concatenation_command)
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


def bring_bams_in_order(bam_list: List[str]) -> List[str]:
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


def get_genomic_chunks_for_tagging(bam_for_tagging: Union[str, Path], chunk_size=TAGGING_CHUNK_SIZE, offset=0) \
        -> list:
    ref_lengths = get_reference_tuples(bam=bam_for_tagging)
    whole_genome_regions = [(chrm, offset, r_len) for chrm, r_len in ref_lengths]
    genomic_chunks = []
    for chrom, _strt, stop in whole_genome_regions:
        n_splits = (stop - offset) // chunk_size
        cum_size = offset
        for _split_idx in range(0, n_splits, 1):
            genomic_chunks.append((chrom, cum_size, cum_size + chunk_size, chunk_size))
            cum_size += chunk_size
        resudual_bases = stop - n_splits * chunk_size
        if resudual_bases * 10 <= chunk_size and n_splits:  # remaining part of scaffold is <= 10% of chunk size
            # -> just add to last! BUT: we need to have at least one full chunk. Otherwise, jsut add the entire sequence
            #    because it is smaller than one chunk
            last_chrom, last_start, last_end, _chk_size = genomic_chunks[-1]
            genomic_chunks[-1] = (last_chrom, last_start, stop, chunk_size + resudual_bases)
        else:  # just add as separate chunk otherwise
            genomic_chunks.append((chrom, cum_size, stop, stop-cum_size))
    return genomic_chunks


def tag_bam_with_correction_weights_parallel(sample_output_dir: str, two_bit_genome_file: str, threads: int,
                                             correction_matrix: np.array, frag_len_range: range, bam_path: str,
                                             ref_lengths: Dict[str, int], temporary_directory_sample: str,
                                             gc_base_limits: range, multithread_access=True, output_unaligned=False,
                                             tag_name=DEFAULT_TAG_NAME, samtools_path=DEFAULT_SAMTOOLS_PATH):
    """
    Size increase of BAM file: 6.8 Gb to 6.9 Gb ~= 1.5%
    Test on the 22/11/2022: duration of BAM file tagging was 0:11:30 (h:mm:ss)
    :param gc_base_limits:
    :param frag_len_range:
    :param threads:
    :param multithread_access:
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
    chunks_per_proc = []
    chunk_lengths_per_proc = []
    proc_list_done = []
    for _idx in range(n_processes):
        all_worker_kwargs.append([])
        recv, sender = mp.Pipe(duplex=False)
        receivers.append(recv)
        senders.append(sender)
        chunks_per_proc.append([])
        chunk_lengths_per_proc.append([])
        proc_list_done.append(False)
    # reference_contigs -> ref_lengths! compute total number of bases;
    # add contigs until over threshold, then record base+/-; regard in next iteration (add fewer/more bases)
    # -> as soon as the next added contig would overrun the threshold,
    #    look for the most fitting contig and add that one instead
    genomic_chunks = get_genomic_chunks_for_tagging(bam_for_tagging=bam_path, chunk_size=TAGGING_CHUNK_SIZE, offset=0)
    target_base_sum_per_process = 1 + sum(map(lambda c: c[3], genomic_chunks)) // n_processes
    for scaff_idx, (scaff, strt, stp, scaff_length) in enumerate(genomic_chunks):
        target_proc_idx = scaff_idx % n_processes
        # check if adding to list is possible; if so, just add to current (base_offset does not change)
        if sum(chunk_lengths_per_proc[target_proc_idx]) + scaff_length <= target_base_sum_per_process:
            chunk_lengths_per_proc[target_proc_idx].append(scaff_length)
            chunks_per_proc[target_proc_idx].append((scaff, strt, stp, scaff_length))
        else:  # the current contig should be added somewhere else and choose a better suited list
            bp_overshot = [sum(chunk_lengths_per_proc[list_search_index]) + scaff_length - target_base_sum_per_process
                           for list_search_index in range(n_processes)]
            final_index = bp_overshot.index(min(bp_overshot))
            chunk_lengths_per_proc[final_index].append(scaff_length)
            chunks_per_proc[final_index].append((scaff, strt, stp, scaff_length))
    # feedback to user
    log(message='\n'.join([f"process ID {cid} processes contigs: {', '.join(map(lambda e: e[0], chrms))}\n"
                           f"This corresponds to {sum(chunk_lengths_per_proc[cid]):,} bp (bp compared to target of "
                           f"{target_base_sum_per_process:,}bp = "
                           f"{sum(chunk_lengths_per_proc[cid]) / target_base_sum_per_process:.1%})"
                          for cid, chrms in enumerate(chunks_per_proc)]), log_level=logging.DEBUG, i_log_with=LOGGER)
    # start unaligned reads extraction
    worker_output_path = Path(temporary_directory_sample) / 'scaffold_BAMs_pre-merging'
    worker_output_path.mkdir(parents=True, exist_ok=True)
    unaligned_bam = None
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
                                    'tagging_chunks_list': chunks_per_proc[w_idx],
                                    'sender_connection': senders[w_idx],
                                    'reference_lengths': ref_lengths,
                                    'use_multithreading': multithread_access,
                                    'parent_logger': LOGGER,
                                    'gc_bases_offset': gc_start}
    # create worker processes
    bam_tagging_workers = [mp.Process(target=bam_tagging_worker_single_chunk, kwargs=worker_kwargs)
                           for worker_kwargs in all_worker_kwargs]
    # start all workers, receive results, wait until finished, and close them
    _ = deque(map(lambda w: w.start(), bam_tagging_workers), maxlen=0)  # ~1GB RAM usage using 12 dual-thread processes
    tagged_scaffold_bam_files = []
    try:
        _ = [tagged_scaffold_bam_files.extend(incoming_result.recv())
             for incoming_result in receivers]  # get returns from workers
    except TypeError:  # " 'int' object is not iterable" occurs if one of the bam tagging workers returns with error
        log(message="Error in BAM tagging worker detected. Cannot proceed. Terminating..",
            log_level=logging.ERROR, close_handlers=True, i_log_with=LOGGER)
        sys.exit(3)
    _ = deque(map(lambda w: w.join(), bam_tagging_workers), maxlen=0)
    _ = deque(map(lambda w: w.close(), bam_tagging_workers), maxlen=0)
    # sanity check: do we have any scaffold/contig BAM files for merging?
    if len(tagged_scaffold_bam_files) == 0:
        log(message="Did not get any scaffold/contig-wise BAM files for merging. Cannot proceed. Terminating..",
            log_level=logging.ERROR, close_handlers=True, i_log_with=LOGGER)
        sys.exit(3)
    # define final output path
    tagged_bam_file_path = Path(temporary_directory_sample) / \
        '.'.join(Path(bam_path).name.split('.')[:-2] +
                 [f"{Path(bam_path).name.split('.')[-2]}", "GCtagged", "bam"])
    tagged_bam_file = tagged_bam_file_path
    # concatenate BAM files and index
    tagged_scaffold_bam_files_in_order = bring_bams_in_order(bam_list=tagged_scaffold_bam_files)
    if output_unaligned:
        unaligned_extraction_handle.join(timeout=600)  # wait for max. 10 minutes
        unaligned_extraction_handle.close()
    # move from temporary sample dir into output dir after checking that the latter does not exist
    target_path = Path(sample_output_dir)  # this is the sample output dir
    target_path.mkdir(parents=True, exist_ok=True)  # ensure parent path exists to be able to move sample dir
    log(message=f"Moving GC-bias computation output from '{temporary_directory_sample}' to '{sample_output_dir}' ..",
        log_level=logging.INFO, i_log_with=LOGGER)
    if unaligned_bam is not None and unaligned_bam.is_file():
        create_bam_index(bam_path=unaligned_bam, samtools_path=samtools_path, check_success=True)
        # check number of total reads in unaligned BAM file
        with silently_open_alignment_file(unaligned_bam, mode='rb') as f_ubam:
            index_statistics = f_ubam.get_index_statistics()
        if sum([total for _cont, _mp, _ump, total in index_statistics]) == 0:
            log(message=f"There were no unaligned reads detected in the input file. Nothing will be output.",
                log_level=logging.INFO, i_log_with=LOGGER)
            unaligned_bam.unlink()  # delete empty uBAM
            Path(f'{unaligned_bam}.bai').unlink()  # delete empty index
        else:  # there were unaligned reads - add unaligned reads at the end of the BAM file!
            if output_unaligned:
                _ = shutil.move(unaligned_bam, target_path)  # copy uBAM to target dir
                _ = shutil.move(f'{unaligned_bam}.bai', target_path)  # move index for index-free samtools cat
            else:
                unaligned_bam.unlink()  # delete empty uBAM
                Path(f'{unaligned_bam}.bai').unlink()  # delete empty index
    samtools_cat_bams(list_of_bams=tagged_scaffold_bam_files_in_order, samtools_path=samtools_path,
                      keep_input=False, tmp_dir=temporary_directory_sample, output_bam=tagged_bam_file)  # also indexes
    for f in Path(temporary_directory_sample).iterdir():  # handle all files/dirs separately
        # -> don't just batch delete existing target dir! Will exist probably due to GC-bias computation step
        if f.is_dir():  # subdir -> can be completely moved & target completely deleted
            # (will only concern the bam parts if kept)
            target_dir_path = target_path / f.name
            if target_dir_path.exists():
                log(message=f"Target path for directory '{f.name}' exists and will be completely deleted! "
                            f"Deleting '{target_dir_path}' ..", log_level=logging.WARNING, i_log_with=LOGGER)
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
                               weights_flen_range=None, mask_flen_range=None) -> Tuple[np.array, range, range]:
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
    if mask_matrix is None:  # try to load from mask file if no mask matrix is given
        if mask_path is None or not mask_path.is_file():
            try:  # to find mask based on correction_matrix name
                candidate_matrix = sorted(list(Path(weights_path.parent).glob(
                    f"{sample_id}*_gc_bias_computation_mask.txt*")), key=lambda c: len(c.name), reverse=True)[0]
                mask_matrix, mask_flen_range = load_txt_to_matrix_with_meta(filename=candidate_matrix, to_dtype=bool,
                                                                            loading_logger=LOGGER)
            except (IndexError, AttributeError):   # IndexError out of range; 'NoneType' object no attribute 'parent'
                pass  # there is no weights file -> only matrix; behaviour not optimal!
                # Will potentially use correction weights that are non-default just because of blurring
    if mask_flen_range is not None and \
            mask_flen_range == weights_flen_range:  # sanity check dimensions - discard mask if different
        reduced_mask, (deleted_rows, deleted_columns) = reduce_matrix(matrix_to_trim=mask_matrix, border_elements=0,
                                                                      trim_dimensions_exclusively_containing=[False])
    if reduced_mask is None:  # workaround: reduce based on correction weights themselves
        reduced_weights, (deleted_rows, deleted_columns) = reduce_matrix(
            matrix_to_trim=correction_matrix, trim_dimensions_exclusively_containing=[1., 0.], border_elements=0)
        resulting_flen_range = range(weights_flen_range.start + deleted_rows.start,
                                     weights_flen_range.stop - deleted_rows.stop)
    else:  # default behavior
        resulting_flen_range = range(mask_flen_range.start + deleted_rows.start,
                                     mask_flen_range.stop - deleted_rows.stop)
        reduced_weights = trim_2d_matrix(matrix=correction_matrix, rows=deleted_rows, columns=deleted_columns)
    resulting_gc_bases_range = range(deleted_columns.start, correction_matrix.shape[1] - 1 - deleted_columns.stop)
    return reduced_weights, resulting_flen_range, resulting_gc_bases_range


def get_reference_tuples(bam: Union[str, Path]) -> Tuple[Tuple[Any, Any]]:
    with silently_open_alignment_file(bam, mode='rb', threads=1) as input_bam_file:
        # allow for any error output resulting from BAM access to be printed here
        return tuple(zip(input_bam_file.references, input_bam_file.lengths))


def get_reference_contig_lengths(bam: Union[str, Path]):
    reference_tuples = get_reference_tuples(bam=bam)
    # create chromosome/scaffold length dictionary
    reference_lengths = {}.fromkeys([r_n for r_n, _r_l in reference_tuples])
    for ref_name, ref_len in reference_tuples:
        reference_lengths[ref_name] = ref_len
    return reference_lengths


def sufficient_aligned_reads_available(bam_path: str, target_fragments_processed: int,
                                       chunks_for_processing: ChunksList) -> Tuple[bool, int]:
    """
    Only standard-non-sex-contigs are counted. Number of mappable GRCh38 bases (without blacklist): 2,848,552,403
    :param bam_path:
    :param target_fragments_processed:
    :param chunks_for_processing:
    :return: bool if enough fragments can be expected from te computation
    """
    std_chroms = [f'chr{idx}' for idx in range(1, 23, 1)]
    with silently_open_alignment_file(bam_path, mode='rb') as f_in_bam:
        index_statistics = f_in_bam.get_index_statistics()
    mapped_std_reads = sum(filter(lambda x: x is not None,
                                  [idx_stat.mapped if idx_stat.contig.lower() in std_chroms else None
                                   for idx_stat in index_statistics]))
    # get number of bases in all chunks = max. portion of reference genome we can compute over
    bases_in_all_std_chunks = sum(filter(lambda x: x is not None,
                                         [stop - start if chrom in std_chroms else None
                                          for (chrom, start, stop, _score) in chunks_for_processing]))
    fraction_chunked_genome = bases_in_all_std_chunks / 2848552403  # second number is "usable" bases of hg38 reference
    expected_number_fragments_in_chunks = int(round(mapped_std_reads / 2 * fraction_chunked_genome, ndigits=0))
    # paired-end data -> divide by 2
    # expected_number_fragments_in_chunks = 0  # coding instruction: if exception, catch it here and set value to 0
    return expected_number_fragments_in_chunks >= target_fragments_processed, expected_number_fragments_in_chunks


def create_bam_index(bam_path: Union[Path, str], samtools_path: Union[Path, str], check_success=True):
    bam_path = Path(bam_path)
    indexing_command = [str(samtools_path), 'index', str(bam_path)]
    ran_indexing_subp = sp.run(indexing_command)
    try:
        ran_indexing_subp.check_returncode()
    except sp.CalledProcessError:
        log(message=f"BAM indexing process for file '{Path(bam_path).name}' returned with error. Cannot proceed."
                    f"Terminating ..",
            log_level=logging.ERROR, close_handlers=True, i_log_with=LOGGER)
        sys.exit(2)
    if check_success:
        # check for existence of putative BAM index:
        if not any([Path(potential_index).is_file()
                    for potential_index in (f'{bam_path}.bai',
                                            Path(bam_path.parent) / f'{bam_path.stem}.bai')]):
            log(message=f"No BAM index found for file '{bam_path}'. Cannot proceed. Terminating ..",
                log_level=logging.ERROR, close_handlers=True, i_log_with=LOGGER)
            sys.exit(2)


def fix_bam_index(bam_path: Union[Path, str], samtools_path: str, silent=True):
    with AlignmentFile(bam_path, 'rb') as f_aln:
        try:
            f_aln.get_index_statistics()  # check_index(self) could also be used
        except ValueError:  # ValueError: mapping information not recorded in index or index not available
            # create missing index
            if not silent:
                log(message=f"GCparagon requires index statistics to be present in BAM index file. Creating missing "
                            f"index file for bam '{bam_path}' ..", log_level=logging.INFO, i_log_with=LOGGER)
            create_bam_index(bam_path=bam_path, samtools_path=samtools_path)
        except AttributeError:  # if htsfile is SAM formatted and thus has no index
            log(message=f"input BAM file is actually a SAM file. Code requires a BAM file. Terminating ..",
                log_level=logging.ERROR, close_handlers=True, i_log_with=LOGGER)
            sys.exit(2)


def manage_bad_chunks(bad_chunks_bed: Optional[str]) \
        -> Tuple[Optional[Dict[Tuple[str, int, int],
                               List[int]]],
                 Optional[str]]:
    if bad_chunks_bed is not None:  # figure out most recent one!
        candidate_bad_library_files = Path(bad_chunks_bed).parent.glob(
            Path(bad_chunks_bed).name)  # yields generator object -> no iteration if empty/nothing found
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
            new_bad_chunks_bed_path = Path(
                time.strftime(str(Path(bad_chunks_bed).parent / BAD_CHUNKS_FILE_PATTERN_AS_PATH),
                              most_recent_lib_date))
            if new_bad_chunks_bed_path.is_file():
                log(message=f"Using bad chunks library file '{new_bad_chunks_bed_path}'",
                    log_level=logging.INFO, i_log_with=LOGGER)
                bad_chunks_bed = new_bad_chunks_bed_path
        except ValueError:  # max() arg is an empty sequence
            pass  # leave bad_chunks_bed as was input
    if bad_chunks_bed is None or not Path(bad_chunks_bed).exists():  # if none was found/path is not valid
        return None, None
    bad_chunks = read_bad_chunks_bed_file(bed_path=bad_chunks_bed)  # requires score in col 4 and integer field in col 5
    return bad_chunks, bad_chunks_bed


def main() -> int:
    global LOGGER

    cmd_args = get_cmdline_args()
    # input options
    input_bam = cmd_args.input_bam
    two_bit_reference_file = cmd_args.two_bit_reference_file
    chunks_bed_file = cmd_args.chunks_bed_file
    exclude_chunks_bed_file = cmd_args.exclude_chunks_bed_file
    correction_weights_matrix_path = cmd_args.correction_weights
    mask_path = cmd_args.weights_mask
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
    use_multithreading = cmd_args.multithread_bam_access
    min_unclipped_aln_fracton = cmd_args.min_unclipped_aln_fracton
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
    keep_chunk_data = cmd_args.keep_chunk_data
    output_corrected_bam = cmd_args.output_corrected_bam
    floating_point_precision = cmd_args.floating_point_precision
    gc_tag_name = cmd_args.gc_tag_name
    write_updated_bad_chunks_library = cmd_args.write_updated_bad_chunks_library
    focus_plots = not cmd_args.dont_focus_plots
    output_unaligned_reads = cmd_args.output_unaligned_reads
    # processing settings
    if correction_weights_matrix_path is not None:
        correction_weights_matrix_path = Path(correction_weights_matrix_path)
    if mask_path is not None:
        mask_path = Path(mask_path)
    np.seterr(all='raise')
    # CHECK INPUT:
    # check unfixable parameters:
    if not os.path.isfile(input_bam):
        raise AttributeError(f"input BAM file '{input_bam}' does not exist!")
    if not os.path.isfile(two_bit_reference_file):
        raise AttributeError(f"2bit reference genome file '{two_bit_reference_file}' does not exist!")
    if not samtools_path or not Path(samtools_path).exists() or not Path(samtools_path).is_file():
        raise AttributeError("path to samtools executable either not found or not accessible. Please provide a valid "
                             "and accessible path using '-sp' or '--samtools-path'.")
    # check and fix correctable parameters
    exit_after_warnings = 0
    print_warnings = []
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
    if not (0. <= min_unclipped_aln_fracton <= 1.):
        print_warnings.append(f"Minimum unclipped alignment fraction was set to {min_unclipped_aln_fracton} but must "
                              "be a floating point value between 0 and 1. Setting to default value of "
                              f"{DEFAULT_MIN_UNCLIPPED_ALN_FRACTION} instead ..")
        min_unclipped_aln_fracton = DEFAULT_MIN_UNCLIPPED_ALN_FRACTION
    # set preset parameters if defined
    if preset_number:  # 1, 2, or 3
        match preset_number:
            case 1:
                min_frag_occurs = 2
                process_n_fragments = 5*10**6
                n_simulations = 6
                smooth_weights = True
                smoothing_kernel = 'gauss'
                smoothing_intensity = 5
                detect_outliers = True
                outliers_method = 'IQR'
                outlier_stringency = 2
            case 2:
                min_frag_occurs = 10
                process_n_fragments = 50*10**6
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
        print_warnings.append(f"Number of simulations per chunk was set to {n_simulations} but must be at least 1. "
                              f"Setting to default of {DEFAULT_SIMULATION_REPETITIONS} instead ..")
        n_simulations = DEFAULT_SIMULATION_REPETITIONS
    # find cpu count boundaries (asserts hyper-threading architecture)
    available_logical_cores = len(os.sched_getaffinity(0))
    max_efficiently_usable_physical_cores = available_logical_cores // 2 if available_logical_cores > 1 else 1
    max_efficiently_usable_threads = max_efficiently_usable_physical_cores * (2 if use_multithreading else 1)
    if max_efficiently_usable_physical_cores * 2 < 4:  # general HW check
        print_warnings.append('GCparagon requires at least 4 logical cores for being able to run. Only '
                              f'{max_efficiently_usable_physical_cores * 2} were estimated to be available. Exiting..')
        exit_after_warnings = 1
    if total_number_threads > max_efficiently_usable_threads:
        print_warnings.append(f'CPUs to use in multiprocessing operations was set to {total_number_threads} but number '
                              f'of available efficiently usable logical cores estimated to be only '
                              f'{max_efficiently_usable_threads}. Setting to '
                              f'{max_efficiently_usable_threads}.')
        total_number_threads = max_efficiently_usable_threads
    # manage imago parameters
    sample_id = os.path.basename(input_bam).split('.')[0]
    compute_bias = not only_tag_bam
    if only_tag_bam and not Path(correction_weights_matrix_path).is_file():
        print_warnings.append('input argument --correction-weights missing. Tag-only-mode not possible. Exiting..')
        exit_after_warnings = 1
    if not output_directory or output_directory == str(Path(input_bam).parent):
        print_warnings.append('Output directory is either input BAM parent directory or was None. Setting '
                              "it to subdirectory of input BAM parent directory: 'GC_correction_output'")
        output_directory = str(Path(input_bam).parent / 'GC_correction_output')
    # choose most recent bad chunks library version if multiple are present in parent directory
    bad_chunks, exclude_chunks_bed_file = manage_bad_chunks(bad_chunks_bed=exclude_chunks_bed_file)
    # set up target output directory and logfile
    start_time = time.localtime()
    sample_out_dir_path = Path(output_directory) / sample_id
    if sample_out_dir_path.exists() and compute_bias:
        print_warnings.append(f"Output path for GC bias computation exists. Deleting completely: "
                              f"'{sample_out_dir_path}'")
        shutil.rmtree(sample_out_dir_path)  # will fail if readonly files are present!
    try:
        sample_out_dir_path.mkdir(parents=True, exist_ok=True)  # ensure output directory for sample exists
    except FileExistsError:  # path is a file -> delete it!
        print_warnings.append(f"Output directory for sample {sample_id} path is a file. Deleting file ..")
        os.remove(sample_out_dir_path)
        sample_out_dir_path.mkdir(parents=True, exist_ok=True)
    sample_out_dir = str(sample_out_dir_path)
    # set up logging
    LOGGER = set_up_logging(logfile_path=str(sample_out_dir_path / f"{sample_id}_GCbiasCorrection_"
                                                                   f"{time.strftime('%d-%m-%Y_%H-%M', start_time)}"
                                                                   ".log"), logger_name='GCparagon', verbose=verbose)
    if print_warnings:
        for warning_message in print_warnings:
            log(message=warning_message, log_level=logging.WARNING, i_log_with=LOGGER)
    if exit_after_warnings:
        return exit_after_warnings
    if temporary_directory:
        if temporary_directory == Path(input_bam).parent:
            log(message="Temporary directory is identical to input BAM parent directory. Setting it to "
                        "subdirectory of input BAM parent directory: 'GC_correction_tmp'",
                log_level=logging.WARNING, i_log_with=LOGGER)
            sample_temp_dir_path = Path(input_bam).parent / 'GC_correction_tmp'
        else:
            sample_temp_dir_path = Path(temporary_directory) / sample_id
    else:  # define a temporary directory
        sample_temp_dir_path = Path(sample_out_dir) / 'GC_correction_tmp'  # required for BAM merging (samtools)
    sample_temp_dir = str(sample_temp_dir_path)
    # check if index is there. if not, create it!
    fix_bam_index(bam_path=input_bam, samtools_path=samtools_path, silent=False)
    # get reference contigs and lengths
    reference_contig_lengths = get_reference_contig_lengths(bam=input_bam)  # stderr enabled for AlignmentFile
    # TASK 1: compute the GC bias present in the BAM file of the current sample
    #         WARNING: don't use multi-sample/run BAM files!)
    log(message=f"Number of available physical cores: {max_efficiently_usable_physical_cores:,}. "
                f"Estimate of available efficiently usable logical cores: {max_efficiently_usable_threads}. "
                f"Will use {total_number_threads:,} threads.",
        log_level=logging.INFO, i_log_with=LOGGER)
    correction_weights_matrix = None
    weights_mask = None
    if compute_bias:
        # read chunk list (should be >=150Mbp in total)
        all_processable_chunks = read_scored_regions_bed_file(bed_path=chunks_bed_file)
        # remove all bad chunks based on BED library file that occurred for previous samples based on target fragments
        # processed as lower boundary of how many fragments are in the dataset
        expect_sufficient_fragment_count, n_expected_fragments = sufficient_aligned_reads_available(
            bam_path=input_bam, target_fragments_processed=process_n_fragments,
            chunks_for_processing=all_processable_chunks)
        log(message='GCparagon (GC-bias computation) Started.\n' +
                    f"|---------------------------------------------------------------------------------\n"
                    f"|   Configuration for processing sample {sample_id} was:\n"
                    f"|   ++++++++++++++++++++++++++++++++++++++++{'+' * len(sample_id)}\n" +
                    (f'|   Using parameter preset {preset_number} for analysis setup\n'
                     if preset_number in (1, 2, 3) else '') +
                    f"|   Minimum fragment length: {lower_limit_fragment_length:,}bp\n"
                    f"|   Maximum fragment length: {upper_limit_fragment_length:,}bp\n"
                    f"|   Minimum number of specific fragment attribute combination occurrences: {min_frag_occurs:,}\n"
                    f"|   Target number of fragments to process: " + ('all' if process_n_fragments == 9999999999 else
                                                                      f'{process_n_fragments:,}') +
                    f" ({'' if expect_sufficient_fragment_count else 'not '}expected to be reached)\n" +
                    ('' if expect_sufficient_fragment_count or not n_expected_fragments else
                     f"|   Estimated number of fragments that will be processed: {n_expected_fragments:,}\n") +
                    f"|   Repetitions for simulation of expected GC-content per chunk: {n_simulations:,}\n"
                    f"|   Random seed was: {random_seed}\n"
                    f"|   Temporary data will be written to: {sample_temp_dir}\n"
                    f"|   Final results will be moved from temporary path to directory: {sample_out_dir}\n"
                    f"|---------------------------------------------------------------------------------",
            log_level=logging.INFO, i_log_with=LOGGER)
        if not expect_sufficient_fragment_count:
            if n_expected_fragments == 0:
                log(message="Number of expected aligned fragments from all chunks could not be determined."
                            "The BAM index might not contain information about mapped reads per contig/scaffold or "
                            "you likely have provided an unaligned BAM file (uBAM).",
                    log_level=logging.WARNING, i_log_with=LOGGER)
            else:
                log(message=f"Number of aligned fragments expected from all chunks (={n_expected_fragments:,}) "
                            "over the genome was lower than the target fragment count!",
                    log_level=logging.WARNING, i_log_with=LOGGER)
        # sort chunks based on blacklist overlap
        # check all predefined chunks against bad chunks library (remove exactly matching chunks)
        target_chunks = sort_chunks_by_blacklist_overlap(
            all_chunks=all_processable_chunks, bad_chunks=bad_chunks, expected_dataset_fragments=n_expected_fragments,
            max_overlap_percentage=DEFAULT_MAX_CHUNK_PERCENTAGE_BLACKLIST_OVERLAP)
        log(message=f"Predefined genomic chunks loaded: {len(target_chunks):,} chunks available for processing. Will "
                    "stop either after processing " +
                    ('all ' if process_n_fragments == 9999999999 else f'{process_n_fragments:,} ') +
                    "fragments or having exhausted list of genomic chunks, whichever occurs first.",
            log_level=logging.INFO, i_log_with=LOGGER)
        (correction_weights_matrix_path, correction_weights_matrix), \
            (mask_path, weights_mask) = compute_gc_bias_parallel(
            visualize_matrices=plot_result, output_all=output_simulation_results, in_bam=input_bam,
            out_dir_sample=sample_out_dir, use_multithreading=use_multithreading,
            max_flen=upper_limit_fragment_length, chunks_to_process=target_chunks, min_flen=lower_limit_fragment_length,
            simulation_count=n_simulations, threads=total_number_threads, tmp_dir_sample=sample_temp_dir,
            keep_chunk_data=keep_chunk_data, random_seed=random_seed, sample_name=sample_id,
            chrom_sizes=reference_contig_lengths, float_precision=floating_point_precision,
            two_bit_reference_file=two_bit_reference_file, min_frag_occurs=min_frag_occurs,
            target_fragments_processed=process_n_fragments, expected_yield=n_expected_fragments,
            write_updated_bad_chunks_library=write_updated_bad_chunks_library,
            bad_chunks_library_file=exclude_chunks_bed_file, strict_n_base_exclusion=strict_n_base_exclusion,
            detect_outliers=detect_outliers, outlier_detection_method=outliers_method,
            outlier_detection_stringency=outlier_stringency, smooth_weights=smooth_weights,
            smoothing_kernel=smoothing_kernel, smoothing_intensity=smoothing_intensity,
            plot_focus_border=10 if focus_plots else None, min_unclipped_aln_fracton=min_unclipped_aln_fracton)
        # compute end time and give feedback
        log(message=f"Correction weights matrix averaged from {n_simulations:,} simulations written to file: "
                    f"'{correction_weights_matrix_path}'", log_level=logging.INFO, i_log_with=LOGGER)
        computation_end_time = time.localtime()
        elapsed_time = datetime.timedelta(seconds=time.mktime(computation_end_time) - time.mktime(start_time))
        log(message='GCparagon (GC-bias computation) Finished Successfully. '
                    f"Elapsed time: {elapsed_time} (h:mm:ss)", log_level=logging.INFO, i_log_with=LOGGER)
    # TASK 2: create a GC-weights-tagged version of the input BAM file
    if output_corrected_bam or only_tag_bam:
        start_time = time.localtime()
        if compute_bias:  # add separation line in this case for prettier output
            log(message="---------------------------------------------------------------------------------",
                log_level=logging.INFO, i_log_with=LOGGER)
        log(message='GCparagon (adding weights to BAM file) Tagging Started.',
            log_level=logging.INFO, i_log_with=LOGGER)
        default_flen_range = range(lower_limit_fragment_length, upper_limit_fragment_length)
        weights_matrix_for_tagging, flen_range, gc_bas_range = reduce_weights_for_tagging(
            weights_path=correction_weights_matrix_path, mask_path=mask_path, mask_matrix=weights_mask,
            sample_id=sample_id, correction_matrix=correction_weights_matrix, weights_flen_range=default_flen_range,
            mask_flen_range=default_flen_range)
        tag_bam_with_correction_weights_parallel(
            bam_path=input_bam, tag_name=gc_tag_name, multithread_access=use_multithreading,
            correction_matrix=weights_matrix_for_tagging, gc_base_limits=gc_bas_range, samtools_path=samtools_path,
            output_unaligned=output_unaligned_reads, two_bit_genome_file=two_bit_reference_file,
            temporary_directory_sample=sample_temp_dir, threads=total_number_threads, sample_output_dir=sample_out_dir,
            ref_lengths=reference_contig_lengths, frag_len_range=flen_range)
        # compute duration of execution and message user
        computation_end_time = time.localtime()
        elapsed_time = datetime.timedelta(seconds=time.mktime(computation_end_time) - time.mktime(start_time))
        log(message='GCparagon (adding weights to BAM file) Finished Successfully. '
                    f"Elapsed time: {elapsed_time} (h:mm:ss)", log_level=logging.INFO, i_log_with=LOGGER)
    # close logging handlers
    for hdlr in LOGGER.handlers:
        hdlr.flush()
        hdlr.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
