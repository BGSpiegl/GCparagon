#!/usr/bin/env python3

import sys
import numpy as np
from pathlib import Path
from GCparagon.correct_GC_bias import save_matrix_to_txt, DEFAULT_FLOAT_PRECISION

SOURCE_ROOT_DIR = Path(__file__).parent.parent  # ../src/GCparagon
if SOURCE_ROOT_DIR not in sys.path:
    sys.path.append(str(SOURCE_ROOT_DIR))

INVERT = True  # the values given by Griffin are actually the bias, not the correction weight!

MIN_FRAG_LENGTH = 20
MAX_FRAG_LENGTH = 550

samples = ('B01', 'C01', 'H01', 'P01')

########################################################################################################################
# TODO: SET YOUR GCparagon SOURCE DIRECTORY HERE IF NEEDED !!!!!!!!
MY_ABSOLUTE_SOURCE_DIRECTORY_PATH = SOURCE_ROOT_DIR  # write_your_GCparagon_directory_absolute_path_here
# TODO: SET YOUR GCparagon SOURCE DIRECTORY HERE IF NEEDED !!!!!!!!
########################################################################################################################
# check it:
if not MY_ABSOLUTE_SOURCE_DIRECTORY_PATH.is_dir():
    raise FileNotFoundError(f"Your GCparagon source ('GCparagon') directory does not exist!"
                            f"You specified: '{MY_ABSOLUTE_SOURCE_DIRECTORY_PATH}'")

griffin_correction_weights = {}.fromkeys(samples)
gc_bias_files_parent_path = MY_ABSOLUTE_SOURCE_DIRECTORY_PATH / ('validation/01_transformed_Griffin_bias_matrices/'
                                                                 'Griffin_bias_matrices')
bias_files = tuple(gc_bias_files_parent_path.glob('*.GC_bias.txt'))
if not bias_files:
    raise FileNotFoundError("could not find any GC bias files (Griffin output) in directory "
                            f"'{gc_bias_files_parent_path}'")
for gc_bias_file_path in bias_files:
    sample_name = gc_bias_file_path.name.split('_')[0]
    assert griffin_correction_weights.get(sample_name) is None
    griffin_correction_weights[sample_name] = gc_bias_file_path

assert None not in griffin_correction_weights.values()  # all bias file found and assigned

transformed_matrices_output_path = MY_ABSOLUTE_SOURCE_DIRECTORY_PATH / 'validation/01_transformed_Griffin_bias_matrices'
transformed_matrices_output_path.mkdir(parents=True, exist_ok=True)

# IMPORTANT INFORMATION FOR PROCESSING:
# THE OUTPUT REPRESENTS THE ALREADY FILTERED WEIGHTS, OBSERVATIONS AND TARGET OBSERVATION COUNTS!!!
#
# 1) Griffin uses the smoothed_GC_bias values for correction (smoothed via a median_filter)
#    -> since GCparagon also smoothes the matrix, the smoothing is fine for the outcome comparison
# 2) It should be noted that in the griffin GC correction, only bias values over 0.05 are considered and that fragments
#    without matching bias value are not considered. "these fragments are extremely rare, so it is difficult to get a
#    good estimate of GC bias"
#    !! Griffin does not compute bias estimate for rare fragments and excludes them !!
# code line 213ff in griffin_coverage.py: """
#     # get rid of extremely low GC bias values
#     # these fragments will now be excluded
#     # these fragments are extremely rare so it is difficult to get a good estimate of GC bias
#     GC_bias['smoothed_GC_bias'] = np.where(GC_bias['smoothed_GC_bias']<0.05,np.nan,GC_bias['smoothed_GC_bias'])"""
#
# 3) Empty cells for GC_bias and smoothed_GC_bias: these are NaN values, which Griffin uses if no fragments were
# observed
# code line 402ff in griffin_coverage.py:"""
# (https://github.com/adoebley/Griffin/blob/6dfe39d7aa4172f372311c1ce3e6efbb831f2e2e/scripts/ [..]
#  [..] griffin_coverage.py#L397C1-L403C64)
#         #count the fragment weighted by GC bias
#         if not np.isnan(read_GC_bias):
#             GC_cov_dict[midpoint]+=(1/read_GC_bias)
# test for nan values: GC_bias = pd.read_csv(GC_bias_path, sep='\t') -> these lines exist
# (extreme combinations of fragment length and GC bases)
# NOTE concernign filtering:
# ########################
# #count coverage
# ########################
# for read in fetched:
#     #filter out reads
#     if abs(read.template_length)>=sz_range[0] and abs(read.template_length)<=sz_range[1] \
#        and read.is_paired==True and read.mapping_quality>=map_q and read.is_duplicate==False and read.is_qcfail==False:
#         #only use fw reads with positive fragment lengths (negative indicates an abnormal pair)
#         #all paired end reads have a fw and rv read so we don't need the rv read to find the midpoint.
#         if read.is_reverse==False and read.template_length>0:


if __name__ == '__main__':
    # fill with qualified values according to the Griffin procedure
    # (This does not include any coverage processing, just the GC bias matrix processing!
    #  This also does not include any mappability correction!)
    for sample, griffin_bias_matrix in griffin_correction_weights.items():
        # create empty matrix: rows = fragment length (starting at MIN_FRAG_LENGTH); columns is GC base count
        sample_weights_matrix = np.zeros((MAX_FRAG_LENGTH - MIN_FRAG_LENGTH + 1, MAX_FRAG_LENGTH + 1), dtype=np.float64)
        sample_observed_matrix = np.zeros((MAX_FRAG_LENGTH - MIN_FRAG_LENGTH + 1, MAX_FRAG_LENGTH + 1), dtype=np.uint64)
        # sample_target_matrix = np.zeros((MAX_FRAG_LENGTH - MIN_FRAG_LENGTH, MAX_FRAG_LENGTH + 1), dtype=np.float64)
        with open(griffin_bias_matrix, 'rt') as f_gb:
            header_columns = f_gb.readline().strip().split('\t')
            for line in f_gb.readlines():
                try:
                    length, num_GC, number_of_fragments, GC_content, number_of_positions, GC_bias, smoothed_GC_bias = \
                        line.strip().split('\t')
                except ValueError:  # ValueError: not enough values to unpack (expected 7, got 5) -> nan values
                    weight = .0  # leave at default
                    continue
                # cast to correct types
                length = int(length)
                num_GC = int(num_GC)
                number_of_fragments = int(number_of_fragments)
                GC_content = float(GC_content)
                number_of_positions = int(number_of_positions)
                # skip insane records
                if num_GC > length:
                    continue
                    # Griffin: # get rid of values where the num_GC is greater than the length
                    # (included due to the way I made the dict)
                try:
                    GC_bias = float(GC_bias)
                    smoothed_GC_bias = float(smoothed_GC_bias)
                except ValueError:  # ValueError: could not convert string to float: ''
                    # weight = .0  # leave at default
                    continue
                if smoothed_GC_bias <= 0.05:
                    # low correction weight entry exist; number of fragments for it may be zero
                    # and GC bias is also zero; SMOOTHED_GC_BIAS MAY NOT BE ZERO!
                    if str(smoothed_GC_bias) != '0.0':
                        print(f"extremely low bias encountered: {smoothed_GC_bias:.3f} which is below expected lower "
                              f"limit of 0.05")
                    continue
                else:
                    if INVERT:
                        weight = 1 / smoothed_GC_bias
                    else:
                        weight = smoothed_GC_bias
                # update matrix
                try:
                    assert sample_weights_matrix[length-MIN_FRAG_LENGTH, num_GC] == 0.0
                    sample_weights_matrix[length - MIN_FRAG_LENGTH, num_GC] = weight
                    sample_observed_matrix[length - MIN_FRAG_LENGTH, num_GC] += number_of_fragments
                except AssertionError:
                    print(f"Logic Error - Existing weight warning - Griffin weight for fragment length {length:,} bp "
                          f"and GC base count {num_GC:,} already existed: "
                          f"{sample_weights_matrix[length-MIN_FRAG_LENGTH, num_GC]}.\nWanted to set it to: {weight}")
        # save transformed Griffin weights matrix in GCparagon style:
        weights_matrix_output_path = transformed_matrices_output_path / \
            f"{sample}_gc_weights.Griffin.txt.gz"
        save_matrix_to_txt(matrix=sample_weights_matrix, output_dir=str(weights_matrix_output_path.parent),
                           max_frag_length=MAX_FRAG_LENGTH, min_frag_length=MIN_FRAG_LENGTH, gzipped=True,
                           filename=weights_matrix_output_path.name, float_data_precision=DEFAULT_FLOAT_PRECISION)
        # save observed counts matrix
        observations_matrix_output_path = transformed_matrices_output_path / \
            f"{sample}_observed_attributes_matrix.Griffin.txt.gz"
        save_matrix_to_txt(matrix=sample_observed_matrix, output_dir=str(observations_matrix_output_path.parent),
                           max_frag_length=MAX_FRAG_LENGTH, min_frag_length=MIN_FRAG_LENGTH, gzipped=True,
                           filename=observations_matrix_output_path.name, float_data_precision=DEFAULT_FLOAT_PRECISION)
        # save artificial mask file
        smoothed_gc_bias_based_mask = sample_observed_matrix.astype(bool)
        mask_matrix_output_path = transformed_matrices_output_path / f"{sample}_gc_bias_computation_mask.Griffin.txt.gz"
        save_matrix_to_txt(matrix=smoothed_gc_bias_based_mask, output_dir=str(mask_matrix_output_path.parent),
                           max_frag_length=MAX_FRAG_LENGTH, min_frag_length=MIN_FRAG_LENGTH, gzipped=True,
                           filename=mask_matrix_output_path.name, float_data_precision=DEFAULT_FLOAT_PRECISION)
        # save artificial target counts matrix
        # is identical to one which is continuously aggregated (as should be!)
        sample_target_matrix = np.multiply(sample_weights_matrix, sample_observed_matrix)
        # former "target_fragment_count_matrix"
        target_matrix_output_path = transformed_matrices_output_path / \
            f"{sample}_target_attributes_matrix.Griffin.txt.gz"
        save_matrix_to_txt(matrix=sample_target_matrix, output_dir=str(target_matrix_output_path.parent),
                           max_frag_length=MAX_FRAG_LENGTH, min_frag_length=MIN_FRAG_LENGTH, gzipped=True,
                           filename=target_matrix_output_path.name, float_data_precision=DEFAULT_FLOAT_PRECISION)
        # give feedback:
        observed_fragment_count = sample_observed_matrix.sum()
        target_fragment_count = sample_target_matrix.sum()
        griffin_count_difference = target_fragment_count-observed_fragment_count
        print(f"Finished processing bias matrix '{griffin_bias_matrix.name}' from sample '{sample}'.\n"
              f"Total number of observed fragments: {observed_fragment_count:,} (standard fragment counting)\n"
              f"Total number of target fragments: {target_fragment_count:,} "
              f"(fragments counted by summing their weights)\n"
              f"Total difference observed - target fragments: {griffin_count_difference:,} "
              f"(= {griffin_count_difference / observed_fragment_count:%})")
        print(f"Resulting weights matrix for sample '{sample}'\n{'-'*100}\n"
              f"All weights mean: {np.mean(sample_weights_matrix):.4f}\n"
              f"All weights median: {np.median(sample_weights_matrix[sample_weights_matrix != 0.0]):.4f}\n"
              f"Only non-zero weights mean: {np.mean(sample_weights_matrix, where=sample_weights_matrix != 0.0):.4f}\n"
              f"Only non-zero weights median: "
              f"{np.median(sample_weights_matrix[sample_weights_matrix != 0.0]):.4f}\n{'-'*100}\n")
        #
        # RESULTING OUTPUT:
        # +++++++++++++++++
        # Finished processing bias matrix 'B01_griffin.GC_bias.txt' from sample 'B01'.
        # Total number of observed fragments: 403,562,610 (standard fragment counting)
        # Total number of target fragments: 222,429,331.91008586 (fragments counted by summing their weights)
        # Total difference observed - target fragments: -181,133,278.08991414 (= -44.883563%)
        # Resulting weights matrix for sample 'B01'
        # ----------------------------------------------------------------------------------------------------
        # All weights mean: 0.2560
        # All weights median: 0.5820
        # Only non-zero weights mean: 0.8819
        # Only non-zero weights median: 0.5820
        # ----------------------------------------------------------------------------------------------------
        #
        # Finished processing bias matrix 'C01_griffin.GC_bias.txt' from sample 'C01'.
        # Total number of observed fragments: 471,840,152 (standard fragment counting)
        # Total number of target fragments: 355,282,845.5548582 (fragments counted by summing their weights)
        # Total difference observed - target fragments: -116,557,306.44514179 (= -24.702710%)
        # Resulting weights matrix for sample 'C01'
        # ----------------------------------------------------------------------------------------------------
        # All weights mean: 0.0894
        # All weights median: 0.6992
        # Only non-zero weights mean: 1.4444
        # Only non-zero weights median: 0.6992
        # ----------------------------------------------------------------------------------------------------
        #
        # Finished processing bias matrix 'H01_griffin.GC_bias.txt' from sample 'H01'.
        # Total number of observed fragments: 521,544,609 (standard fragment counting)
        # Total number of target fragments: 336,025,093.4060935 (fragments counted by summing their weights)
        # Total difference observed - target fragments: -185,519,515.59390652 (= -35.571169%)
        # Resulting weights matrix for sample 'H01'
        # ----------------------------------------------------------------------------------------------------
        # All weights mean: 0.3826
        # All weights median: 0.6568
        # Only non-zero weights mean: 1.1163
        # Only non-zero weights median: 0.6568
        # ----------------------------------------------------------------------------------------------------
        #
        # Finished processing bias matrix 'P01_griffin.GC_bias.txt' from sample 'P01'.
        # Total number of observed fragments: 806,122,625 (standard fragment counting)
        # Total number of target fragments: 703,940,645.3128626 (fragments counted by summing their weights)
        # Total difference observed - target fragments: -102,181,979.68713737 (= -12.675736%)
        # Resulting weights matrix for sample 'P01'
        # ----------------------------------------------------------------------------------------------------
        # All weights mean: 0.3356
        # All weights median: 0.8728
        # Only non-zero weights mean: 2.2731
        # Only non-zero weights median: 0.8728
        # ----------------------------------------------------------------------------------------------------
