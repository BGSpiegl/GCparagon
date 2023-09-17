#!/usr/bin/env python3

import numpy as np
from typing import Union, Optional, Tuple
from pathlib import Path
import gzip

########################################################################################################################
# TODO: SET YOUR GCparagon SOURCE DIRECTORY HERE !!!!!!!!
MY_ABSOLUTE_SOURCE_DIRECTORY_PATH = Path('write_your_GCparagon_directory_absolute_path_here')
# TODO: SET YOUR GCparagon SOURCE DIRECTORY HERE !!!!!!!!
########################################################################################################################
# check it:
if not MY_ABSOLUTE_SOURCE_DIRECTORY_PATH.is_dir():
    raise FileNotFoundError(f"Your GCparagon source ('GCparagon') directory does not exist!"
                            f"You specified: '{MY_ABSOLUTE_SOURCE_DIRECTORY_PATH}'")

griffin_matrices_path = MY_ABSOLUTE_SOURCE_DIRECTORY_PATH / 'validation/transformed_Griffin_bias_matrices'
gcparagon_matrices_parent_path = MY_ABSOLUTE_SOURCE_DIRECTORY_PATH / 'GCparagon_matrices_ALL_PRESETS'

correction_algos = ('Griffin', 'GCparagon')
all_samples = ('B01', 'C01', 'H01', 'P01')

matrices = {'O_gc': {'Griffin': {'B01': griffin_matrices_path / 'B01_observed_attributes_matrix.Griffin.txt.gz',
                                 'C01': griffin_matrices_path / 'C01_observed_attributes_matrix.Griffin.txt.gz',
                                 'H01': griffin_matrices_path / 'H01_observed_attributes_matrix.Griffin.txt.gz',
                                 'P01': griffin_matrices_path / 'P01_observed_attributes_matrix.Griffin.txt.gz'},
                     'GCparagon': {'B01': 'preset{}/B01_observed_attributes_matrix.txt.gz',
                                   'C01': 'preset{}/C01_observed_attributes_matrix.txt.gz',
                                   'H01': 'preset{}/H01_observed_attributes_matrix.txt.gz',
                                   'P01': 'preset{}/P01_observed_attributes_matrix.txt.gz'}
                     # format these with preset taken from (1, 2, 3)
                     },
            'W_gc': {'Griffin': {'B01': griffin_matrices_path / 'B01_gc_weights.Griffin.txt.gz',
                                 'C01': griffin_matrices_path / 'C01_gc_weights.Griffin.txt.gz',
                                 'H01': griffin_matrices_path / 'H01_gc_weights.Griffin.txt.gz',
                                 'P01': griffin_matrices_path / 'P01_gc_weights.Griffin.txt.gz'},
                     'GCparagon': {'B01': 'preset{}/'
                                          'B01_gc_weights_{}simsMean.2IQRoutliersRemoved.{}IgaussSmoothed.txt.gz',
                                   'C01': 'preset{}/'
                                          'C01_gc_weights_{}simsMean.2IQRoutliersRemoved.{}IgaussSmoothed.txt.gz',
                                   'H01': 'preset{}/'
                                          'H01_gc_weights_{}simsMean.2IQRoutliersRemoved.{}IgaussSmoothed.txt.gz',
                                   'P01': 'preset{}/'
                                          'P01_gc_weights_{}simsMean.2IQRoutliersRemoved.{}IgaussSmoothed.txt.gz'}
                     # format these with preset taken from (1, 2, 3) AND with postprocessing parameters
                     }
            }

preset_postprocessing_parameters = {1: [6, 5],
                                    2: [4, 2],
                                    3: [4, 2]}


def load_txt_to_matrix_with_meta(filename: Union[str, Path], loading_logger: Optional[str] = None,
                                 to_dtype=np.float64) -> Tuple[np.array, range]:
    """
    :param loading_logger:
    :param filename:
    :param to_dtype:
    :return:
    """
    statistic_matrix = np.loadtxt(filename, delimiter='|', skiprows=0, dtype=to_dtype)
    with gzip.open(str(filename), 'rt') as f_mat:
        hdr = f_mat.readline()
    elements = hdr.split('# rows representing fragment lengths (')[1].split()  # split on whitespace + trim empty fields
    fragment_lengths = range(int(elements[0]), int(elements[3]))  # 20 bp to 550 bp (non-Pythonic in header)
    return statistic_matrix, fragment_lengths


if __name__ == '__main__':
    for correction_algo in correction_algos:
        print(f"Processing now results for '{correction_algo}' correction ..")
        if correction_algo == 'GCparagon':  # process GCparagon matrices
            for preset in range(1, 4, 1):
                print(f"Processing now results for preset{preset} ..")
                across_samples = []
                for sample in all_samples:
                    observed_attributes, flength_range_o = load_txt_to_matrix_with_meta(
                        gcparagon_matrices_parent_path / matrices['O_gc'][correction_algo][sample].format(preset))
                    size_original_fragment_pool = observed_attributes.sum()
                    attribute_weights, flength_range_w = load_txt_to_matrix_with_meta(
                        gcparagon_matrices_parent_path / matrices['W_gc'][correction_algo][sample].format(
                            preset, *preset_postprocessing_parameters[preset]))
                    size_weighted_fragment_pool = (observed_attributes * attribute_weights).sum()
                    discrepancy_fraction = size_weighted_fragment_pool / size_original_fragment_pool - 1.
                    across_samples.append(discrepancy_fraction)
                    print(f"sample '{sample}': discrepancy pool sizes post {correction_algo}-preset{preset} "
                          f"weighting vs. original size: "
                          f"{discrepancy_fraction:.3%}")
                # summarize
                print(f"The average discrepancy for '{correction_algo}-preset{preset}' correction was: "
                      f"{np.mean(across_samples):.3%}")
        else:  # process Griffin matrices
            across_samples = []
            for sample in all_samples:
                observed_attributes, flength_range_o = load_txt_to_matrix_with_meta(
                    matrices['O_gc'][correction_algo][sample])
                size_original_fragment_pool = observed_attributes.sum()
                attribute_weights, flength_range_w = load_txt_to_matrix_with_meta(
                    matrices['W_gc'][correction_algo][sample])
                size_weighted_fragment_pool = (observed_attributes * attribute_weights).sum()
                discrepancy_fraction = size_weighted_fragment_pool / size_original_fragment_pool - 1.
                across_samples.append(discrepancy_fraction)
                print(f"sample '{sample}': discrepancy pool sizes post {correction_algo} weighting vs. original "
                      f"size: {discrepancy_fraction:.3%}")
            # summarize
            print(f"The average discrepancy for '{correction_algo}' correction across all samples was: "
                  f"{np.mean(across_samples):.3%}")


# cmd output:
#
# Processing now results for 'Griffin' correction ..
# sample 'B01': discrepancy pool sizes post Griffin weighting vs. original size: -44.884%
# sample 'C01': discrepancy pool sizes post Griffin weighting vs. original size: -24.703%
# sample 'H01': discrepancy pool sizes post Griffin weighting vs. original size: -35.571%
# sample 'P01': discrepancy pool sizes post Griffin weighting vs. original size: -6.474%
#
# The average discrepancy for 'Griffin' correction across all samples was: -27.908%
# Processing now results for 'GCparagon' correction ..

# Processing now results for preset1 ..
# sample 'B01': discrepancy pool sizes post GCparagon-preset1 weighting vs. original size: 0.447%
# sample 'C01': discrepancy pool sizes post GCparagon-preset1 weighting vs. original size: 0.210%
# sample 'H01': discrepancy pool sizes post GCparagon-preset1 weighting vs. original size: 0.851%
# sample 'P01': discrepancy pool sizes post GCparagon-preset1 weighting vs. original size: 0.822%
#
# The average discrepancy for 'GCparagon-preset1' correction was: 0.583%

# Processing now results for preset2 ..
# sample 'B01': discrepancy pool sizes post GCparagon-preset2 weighting vs. original size: -0.114%
# sample 'C01': discrepancy pool sizes post GCparagon-preset2 weighting vs. original size: -0.045%
# sample 'H01': discrepancy pool sizes post GCparagon-preset2 weighting vs. original size: -0.044%
# sample 'P01': discrepancy pool sizes post GCparagon-preset2 weighting vs. original size: 0.043%
# The average discrepancy for 'GCparagon-preset2' correction was: -0.040%

# Processing now results for preset3 ..
# sample 'B01': discrepancy pool sizes post GCparagon-preset3 weighting vs. original size: -0.128%
# sample 'C01': discrepancy pool sizes post GCparagon-preset3 weighting vs. original size: -0.064%
# sample 'H01': discrepancy pool sizes post GCparagon-preset3 weighting vs. original size: -0.122%
# sample 'P01': discrepancy pool sizes post GCparagon-preset3 weighting vs. original size: -0.027%
#
# The average discrepancy for 'GCparagon-preset3' correction was: -0.085%
