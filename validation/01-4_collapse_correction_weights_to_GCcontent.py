#!/usr/bin/env python3

import re
import sys
import logging
import gzip as gz
import numpy as np
import pandas as pd
from math import ceil
from pathlib import Path
from plotly import graph_objs as go
import plotly.express as px
from typing import Union, Tuple, Optional, Dict

SOURCE_ROOT_DIR = Path(__file__).parent.parent  # ../src/GCparagon
if SOURCE_ROOT_DIR not in sys.path:
    sys.path.append(str(SOURCE_ROOT_DIR))

########################################################################################################################
# TODO: SET YOUR GCparagon SOURCE DIRECTORY HERE IF NEEDED !!!!!!!!
MY_ABSOLUTE_SOURCE_DIRECTORY_PATH = SOURCE_ROOT_DIR  # write_your_GCparagon_directory_absolute_path_here
# TODO: SET YOUR GCparagon SOURCE DIRECTORY HERE IF NEEDED !!!!!!!!
########################################################################################################################

# TODO use new results!
parent_paths = {'Griffin': MY_ABSOLUTE_SOURCE_DIRECTORY_PATH / 'validation/01_transformed_Griffin_bias_matrices',
                'GCparagon': MY_ABSOLUTE_SOURCE_DIRECTORY_PATH / 'validation/01_transformed_Griffin_bias_matrices/'
                                                                 'GCparagon_matrices_per_preset/preset{}'}
# Mind that for expected fragment length distributions (which represents a theoretical distribution) we will use results
# from preset 3 as a first approximation of the observed fragment length distribution.


IMAGE_FORMATS = ('png', 'svg', 'pdf')


correction_weights_matrices = {'Griffin': {'B01': parent_paths['Griffin'] / 'B01_gc_weights.Griffin.txt.gz',
                                           'C01': parent_paths['Griffin'] / 'C01_gc_weights.Griffin.txt.gz',
                                           'H01': parent_paths['Griffin'] / 'H01_gc_weights.Griffin.txt.gz',
                                           'P01': parent_paths['Griffin'] / 'P01_gc_weights.Griffin.txt.gz'}}

GRIFFIN_OBSERVATION_MATRICES_PATH = parent_paths['Griffin']
GCPARAGON_OBSERVATION_MATRICES_PATH_STR = str(parent_paths['GCparagon'])

# Griffin matrix paths:
observed_fragments = {'Griffin': {'H01': GRIFFIN_OBSERVATION_MATRICES_PATH /
                                  'H01_observed_attributes_matrix.Griffin.txt.gz',
                                  'P01': GRIFFIN_OBSERVATION_MATRICES_PATH /
                                  'P01_observed_attributes_matrix.Griffin.txt.gz',
                                  'C01': GRIFFIN_OBSERVATION_MATRICES_PATH /
                                  'C01_observed_attributes_matrix.Griffin.txt.gz',
                                  'B01': GRIFFIN_OBSERVATION_MATRICES_PATH /
                                  'B01_observed_attributes_matrix.Griffin.txt.gz'}}

# GCparagon matrix paths:
sims = {1: 6,
        2: 4,
        3: 4}
smoothing = {1: 5,
             2: 2,
             3: 2}

for preset in range(1, 3, 1):  # no preset 3 validation yet
    observed_fragments[f'GCparagon-preset{preset}'] = {'B01': Path(
        GCPARAGON_OBSERVATION_MATRICES_PATH_STR.format(preset)) / 'B01_observed_attributes_matrix.txt.gz',
                                                       'C01': Path(
        GCPARAGON_OBSERVATION_MATRICES_PATH_STR.format(preset)) / 'C01_observed_attributes_matrix.txt.gz',
                                                       'H01': Path(
        GCPARAGON_OBSERVATION_MATRICES_PATH_STR.format(preset)) / 'H01_observed_attributes_matrix.txt.gz',
                                                       'P01': Path(
        GCPARAGON_OBSERVATION_MATRICES_PATH_STR.format(preset)) / 'P01_observed_attributes_matrix.txt.gz'}
    correction_weights_matrices[f'GCparagon-preset{preset}'] = {'B01': Path(
        str(parent_paths['GCparagon']).format(preset)) / f'B01_gc_weights_{sims[preset]}simsMean.2IQRoutliersRemoved.'
                                                         f'{smoothing[preset]}IgaussSmoothed.txt.gz',
                                                                'C01': Path(
        str(parent_paths['GCparagon']).format(preset)) / f'C01_gc_weights_{sims[preset]}simsMean.2IQRoutliersRemoved.'
                                                         f'{smoothing[preset]}IgaussSmoothed.txt.gz',
                                                                'H01': Path(
        str(parent_paths['GCparagon']).format(preset)) / f'H01_gc_weights_{sims[preset]}simsMean.2IQRoutliersRemoved.'
                                                         f'{smoothing[preset]}IgaussSmoothed.txt.gz',
                                                                'P01': Path(
        str(parent_paths['GCparagon']).format(preset)) / f'P01_gc_weights_{sims[preset]}simsMean.2IQRoutliersRemoved.'
                                                         f'{smoothing[preset]}IgaussSmoothed.txt.gz'}

output_path = MY_ABSOLUTE_SOURCE_DIRECTORY_PATH / ('validation/01_transformed_Griffin_bias_matrices/'
                                                   'weights_per_fragment_GC_content_only')
output_path.mkdir(parents=True, exist_ok=True)

combine_percentages = 2  # like in genome-wide FGCD plots
if combine_percentages < 1:
    raise ValueError(f"combine_percentages must be a positive integer but was {combine_percentages}!")


def log(message: str, log_level: int, logger_name: str, flush=True, close_handlers=False):
    """
    :param message:
    :param log_level:
    :param logger_name:
    :param flush:
    :param close_handlers:
    :return:
    """
    current_logger = logging.getLogger(logger_name)
    match log_level:
        case logging.NOTSET:
            current_logger.info(message)
        case logging.DEBUG:
            current_logger.debug(message)
        case logging.INFO:
            current_logger.info(message)
        case logging.WARNING:
            current_logger.warning(message)
        case logging.ERROR:
            current_logger.error(message)
        case logging.CRITICAL:
            current_logger.critical(message)
    if flush and isinstance(logger_name, logging.Logger):  # do nothing if logger_name is undefined
        for hdlr in current_logger.handlers:
            hdlr.flush()
    if close_handlers:
        for hdlr in current_logger.handlers:
            hdlr.flush()
            hdlr.close()


def load_txt_to_matrix_with_meta(filename: Union[str, Path], loading_logger: Optional[str] = None,
                                 to_dtype=np.float64) -> Tuple[np.array, range]:
    """
    :param loading_logger:
    :param filename:
    :param to_dtype:
    :return:
    """
    if loading_logger is not None:
        log(message=f"Loading statistic matrix from {filename}", log_level=logging.INFO, logger_name=loading_logger)
    else:
        print(f"Loading statistic matrix from {filename}")
    statistic_matrix = np.loadtxt(filename, delimiter='|', skiprows=0, dtype=float).astype(to_dtype)  # is float matrix
    with gz.open(str(filename), 'rt') as f_mat:
        hdr = f_mat.readline()
    elements = hdr.split('# rows representing fragment lengths (')[1].split()  # split on whitespace + trim empty fields
    fragment_lengths = range(int(elements[0]), int(elements[3]))  # 20 bp to 550 bp (non-Pythonic in header)
    return statistic_matrix, fragment_lengths


def plot_gc_content_curve(correction_weights_gccontent_dict: Dict[str, np.array], output_dir: Union[str, Path],
                          annotation: Optional[str] = None, fig_width=1600, fig_height=1000, fig_fontsize=34,
                          show_figure=False, normalize_total_weights=False, reduced_bins=False,
                          spline_interpolation=False, image_formats=('png',)):
    """
    WARNING - function expects already collapsed data if collapsed is wanted for printing
    :param correction_weights_gccontent_dict:
    :param annotation:
    :param normalize_total_weights:
    :param reduced_bins:
    :return:
    """
    # create total count sum for all samples
    trace_ids = sorted(correction_weights_gccontent_dict.keys())
    output_dir = Path(output_dir)
    # normalize datasets if requested
    if normalize_total_weights:  # divide each dataset by number of total fragments and multiply with 100
        raise NotImplementedError
    plot_data_list = []
    columns = ['sample', 'GC content', 'GC content weight (weighted average)']
    # ASSERTS all bins present!
    for trace_id in trace_ids:
        plot_data_list.extend([trace_id,
                               gc_idx * (2 if reduced_bins else 1),
                               correction_weights_gccontent_dict[trace_id][gc_idx]]
                              for gc_idx in range(len(correction_weights_gccontent_dict[trace_id])))
    plot_data = pd.DataFrame(plot_data_list, columns=columns)
    gc_content_fig = px.line(plot_data, x='GC content', y='GC content weight (weighted average)',
                             color='sample',  # line_dash_sequence=['longdash', 'solid'],
                             line_shape='spline' if spline_interpolation else 'linear',
                             template="simple_white", width=fig_width, height=fig_height,
                             title=f"Fragment-length-weighted Correction Weights per GC Content" +
                                   (' (' + annotation + ')') if annotation is not None else '',
                             color_discrete_map={'B01': 'rgba(60, 120, 255, 1.)',
                                                 'C01': 'rgba(255, 125, 0, 1.)',
                                                 'H01': 'rgba(40, 188, 40, 1.)',
                                                 'P01': 'rgba(255, 0, 0, 1.)',
                                                 'GCparagon': 'rgba(40, 188, 40, 0.75)',
                                                 'Griffin': 'rgba(255, 125, 0, 0.75)'})
    # add expected max at 167 bp an 316 bp:
    gc_content_fig.add_trace(go.Scatter(x=(0., 100.), y=(0., 0.),
                                        mode='lines', line={'color': 'rgba(33, 33, 33, 0.5)', 'width': 1.5}))
    gc_content_fig.add_trace(go.Scatter(x=(0., 100.), y=(1., 1.),
                                        mode='lines', line={'color': 'rgba(133, 133, 133, 0.5)', 'width': 1.}))
    # ADJUST ATTRIBUTES OF INDIVIDUAL TRACES
    for dat_idx in range(len(gc_content_fig.data)):
        gc_content_fig.data[dat_idx].line.width = 3
    gc_content_fig.update_layout(font_family="Ubuntu", font_size=fig_fontsize,
                                 xaxis_title='fragment GC content / %',
                                 yaxis_title='GC content weight (weighted average) / 1',
                                 legend={'orientation': 'h', 'xanchor': 'center', 'yanchor': 'top',
                                         'x': 0.5, 'y': -0.2, 'title': ''})
    if show_figure:
        gc_content_fig.show()
    output_dir.mkdir(parents=True, exist_ok=True)
    for image_format in image_formats:
        out_file = output_dir / f"GCparagon_GC-content_WeightedAverageAccordingToRelFlengthFreq" \
                                f"{'_SPLINE' if spline_interpolation else '_LINEAR'}" \
                                f"{('_' + re.sub(', ', '_', annotation)) if annotation else ''}_cfDNAref.{image_format}"
        gc_content_fig.write_image(out_file)


if __name__ == '__main__':
    number_bins = ceil(100. / combine_percentages)
    all_collapsed_weights = {}
    algos = set()
    samples = set()
    for correction_algo, weight_matrix_dict in correction_weights_matrices.items():
        correction_algo_complete = correction_algo
        all_collapsed_weights.update({correction_algo: {}})
        algos.update((correction_algo,))
        for sample_id, weight_matrix_path in correction_weights_matrices[correction_algo].items():
            samples.update((sample_id,))
            all_collapsed_weights[correction_algo_complete].update({sample_id: None})
            weights_matrix, weights_flength_range = load_txt_to_matrix_with_meta(filename=weight_matrix_path)
            # load observations matrix
            observed_matrix, observations_flength_range = load_txt_to_matrix_with_meta(
                filename=observed_fragments[correction_algo_complete][sample_id])
            flength_frequencies = observed_matrix.sum(axis=1)
            relative_flength_frequencies = np.divide(flength_frequencies, flength_frequencies.sum(),
                                                     where=lambda x: x != 0)
            print(f"the sum of all relative fragment length occurrences was: "
                  f"{relative_flength_frequencies.sum():.4f} (1.0 expected)")
            cumulative_gc_content_sample = np.zeros(number_bins, dtype=float)
            cumulated_attribute_weights = np.zeros(number_bins, dtype=float)
            for row_idx in range(weights_matrix.shape[0]):  # iter over rows = fragment lengths - min. length
                if flength_frequencies[row_idx] == 0:  # fragment length was not observed
                    continue  # -> ignore though unlikely
                relative_flength_frequency = relative_flength_frequencies[row_idx]  # for weighted mean
                current_fragment_length = weights_flength_range.start + row_idx
                current_weights_per_gc_bases = weights_matrix[row_idx, :]
                # create border start values for (binned) GC content counting (always use either 50 or 100 bins)
                bin_borders = np.arange(0.,
                                        current_fragment_length*(1+1/number_bins),
                                        current_fragment_length/number_bins)
                try:  # slice start which is last expected non-zero value (can be zero if no observations)
                    bin_borders = bin_borders[:[b >= current_fragment_length-0.001 for b in bin_borders].index(True)]
                except ValueError:  # True is not in list -> nothing to slice
                    pass
                for gc_base_count in range(current_fragment_length+1):  # only consider possible attribute combinations
                    if current_weights_per_gc_bases[gc_base_count] != 0.:  # ignored fragments are not counted
                        try:
                            gc_content_bin_index = [b_start > gc_base_count for b_start in bin_borders].index(True) - 1
                        except ValueError:  # no True inside -> raises this error; only occurs if last bin is correct
                            gc_content_bin_index = len(bin_borders) - 1  # assign to last bin
                        cumulated_attribute_weights[gc_content_bin_index] += relative_flength_frequency
                        cumulative_gc_content_sample[gc_content_bin_index] += \
                            current_weights_per_gc_bases[gc_base_count] * relative_flength_frequency
            # compute average - weighted by relative occurrence of fragment length
            average_weight_gc_content = np.divide(cumulative_gc_content_sample, cumulated_attribute_weights,
                                                  where=(cumulated_attribute_weights > 0.000000001) |
                                                        (cumulative_gc_content_sample > 0.000000001))
            # set zero values to zero because these tend to be faulty/extreme after division
            average_weight_gc_content[(cumulated_attribute_weights < 0.000000001) |
                                      (cumulative_gc_content_sample < 0.000000001)] = 0.
            print(f"{correction_algo} - {sample_id}: average weight across all GC content percentages: "
                  f"{np.mean(average_weight_gc_content):.2f}")  # TODO: this shows a singularity for B01
            all_collapsed_weights[correction_algo][sample_id] = average_weight_gc_content
    # visualize per sample:
    for sample_id in samples:
        plotting_data = {}
        for algo in algos:
            plotting_data.update({algo: all_collapsed_weights[algo][sample_id]})
        plot_gc_content_curve(correction_weights_gccontent_dict=plotting_data, output_dir=output_path,
                              annotation=f"{sample_id}", reduced_bins=True, show_figure=True,
                              image_formats=IMAGE_FORMATS)

# CMD output was:
#
# Loading statistic matrix from /mnt/NVMeScratch/PycharmProjects/GCparagon_public/validation/01_transformed_Griffin_bias_matrices/B01_gc_weights.Griffin.txt.gz
# Loading statistic matrix from /mnt/NVMeScratch/PycharmProjects/GCparagon_public/validation/01_transformed_Griffin_bias_matrices/B01_observed_attributes_matrix.Griffin.txt.gz
# the sum of all relative fragment length occurrences was: 1.0000 (1.0 expected)
# Griffin - B01: average weight across all GC content percentages: 2.14
# Loading statistic matrix from /mnt/NVMeScratch/PycharmProjects/GCparagon_public/validation/01_transformed_Griffin_bias_matrices/C01_gc_weights.Griffin.txt.gz
# Loading statistic matrix from /mnt/NVMeScratch/PycharmProjects/GCparagon_public/validation/01_transformed_Griffin_bias_matrices/C01_observed_attributes_matrix.Griffin.txt.gz
# the sum of all relative fragment length occurrences was: 1.0000 (1.0 expected)
# Griffin - C01: average weight across all GC content percentages: 1.98
# Loading statistic matrix from /mnt/NVMeScratch/PycharmProjects/GCparagon_public/validation/01_transformed_Griffin_bias_matrices/H01_gc_weights.Griffin.txt.gz
# Loading statistic matrix from /mnt/NVMeScratch/PycharmProjects/GCparagon_public/validation/01_transformed_Griffin_bias_matrices/H01_observed_attributes_matrix.Griffin.txt.gz
# the sum of all relative fragment length occurrences was: 1.0000 (1.0 expected)
# Griffin - H01: average weight across all GC content percentages: 1.70
# Loading statistic matrix from /mnt/NVMeScratch/PycharmProjects/GCparagon_public/validation/01_transformed_Griffin_bias_matrices/P01_gc_weights.Griffin.txt.gz
# Loading statistic matrix from /mnt/NVMeScratch/PycharmProjects/GCparagon_public/validation/01_transformed_Griffin_bias_matrices/P01_observed_attributes_matrix.Griffin.txt.gz
# the sum of all relative fragment length occurrences was: 1.0000 (1.0 expected)
# Griffin - P01: average weight across all GC content percentages: 2.98
# Loading statistic matrix from /mnt/NVMeScratch/PycharmProjects/GCparagon_public/validation/01_transformed_Griffin_bias_matrices/GCparagon_matrices_per_preset/preset1/B01_gc_weights_6simsMean.2IQRoutliersRemoved.5IgaussSmoothed.txt.gz
# Loading statistic matrix from /mnt/NVMeScratch/PycharmProjects/GCparagon_public/validation/01_transformed_Griffin_bias_matrices/GCparagon_matrices_per_preset/preset1/B01_observed_attributes_matrix.txt.gz
# the sum of all relative fragment length occurrences was: 1.0000 (1.0 expected)
# GCparagon-preset1 - B01: average weight across all GC content percentages: 1.30
# Loading statistic matrix from /mnt/NVMeScratch/PycharmProjects/GCparagon_public/validation/01_transformed_Griffin_bias_matrices/GCparagon_matrices_per_preset/preset1/C01_gc_weights_6simsMean.2IQRoutliersRemoved.5IgaussSmoothed.txt.gz
# Loading statistic matrix from /mnt/NVMeScratch/PycharmProjects/GCparagon_public/validation/01_transformed_Griffin_bias_matrices/GCparagon_matrices_per_preset/preset1/C01_observed_attributes_matrix.txt.gz
# the sum of all relative fragment length occurrences was: 1.0000 (1.0 expected)
# GCparagon-preset1 - C01: average weight across all GC content percentages: 1.14
# Loading statistic matrix from /mnt/NVMeScratch/PycharmProjects/GCparagon_public/validation/01_transformed_Griffin_bias_matrices/GCparagon_matrices_per_preset/preset1/H01_gc_weights_6simsMean.2IQRoutliersRemoved.5IgaussSmoothed.txt.gz
# Loading statistic matrix from /mnt/NVMeScratch/PycharmProjects/GCparagon_public/validation/01_transformed_Griffin_bias_matrices/GCparagon_matrices_per_preset/preset1/H01_observed_attributes_matrix.txt.gz
# the sum of all relative fragment length occurrences was: 1.0000 (1.0 expected)
# GCparagon-preset1 - H01: average weight across all GC content percentages: 1.08
# Loading statistic matrix from /mnt/NVMeScratch/PycharmProjects/GCparagon_public/validation/01_transformed_Griffin_bias_matrices/GCparagon_matrices_per_preset/preset1/P01_gc_weights_6simsMean.2IQRoutliersRemoved.5IgaussSmoothed.txt.gz
# Loading statistic matrix from /mnt/NVMeScratch/PycharmProjects/GCparagon_public/validation/01_transformed_Griffin_bias_matrices/GCparagon_matrices_per_preset/preset1/P01_observed_attributes_matrix.txt.gz
# the sum of all relative fragment length occurrences was: 1.0000 (1.0 expected)
# GCparagon-preset1 - P01: average weight across all GC content percentages: 1.22
# Loading statistic matrix from /mnt/NVMeScratch/PycharmProjects/GCparagon_public/validation/01_transformed_Griffin_bias_matrices/GCparagon_matrices_per_preset/preset2/B01_gc_weights_4simsMean.2IQRoutliersRemoved.2IgaussSmoothed.txt.gz
# Loading statistic matrix from /mnt/NVMeScratch/PycharmProjects/GCparagon_public/validation/01_transformed_Griffin_bias_matrices/GCparagon_matrices_per_preset/preset2/B01_observed_attributes_matrix.txt.gz
# the sum of all relative fragment length occurrences was: 1.0000 (1.0 expected)
# GCparagon-preset2 - B01: average weight across all GC content percentages: 1.44
# Loading statistic matrix from /mnt/NVMeScratch/PycharmProjects/GCparagon_public/validation/01_transformed_Griffin_bias_matrices/GCparagon_matrices_per_preset/preset2/C01_gc_weights_4simsMean.2IQRoutliersRemoved.2IgaussSmoothed.txt.gz
# Loading statistic matrix from /mnt/NVMeScratch/PycharmProjects/GCparagon_public/validation/01_transformed_Griffin_bias_matrices/GCparagon_matrices_per_preset/preset2/C01_observed_attributes_matrix.txt.gz
# the sum of all relative fragment length occurrences was: 1.0000 (1.0 expected)
# GCparagon-preset2 - C01: average weight across all GC content percentages: 1.18
# Loading statistic matrix from /mnt/NVMeScratch/PycharmProjects/GCparagon_public/validation/01_transformed_Griffin_bias_matrices/GCparagon_matrices_per_preset/preset2/H01_gc_weights_4simsMean.2IQRoutliersRemoved.2IgaussSmoothed.txt.gz
# Loading statistic matrix from /mnt/NVMeScratch/PycharmProjects/GCparagon_public/validation/01_transformed_Griffin_bias_matrices/GCparagon_matrices_per_preset/preset2/H01_observed_attributes_matrix.txt.gz
# the sum of all relative fragment length occurrences was: 1.0000 (1.0 expected)
# GCparagon-preset2 - H01: average weight across all GC content percentages: 1.10
# Loading statistic matrix from /mnt/NVMeScratch/PycharmProjects/GCparagon_public/validation/01_transformed_Griffin_bias_matrices/GCparagon_matrices_per_preset/preset2/P01_gc_weights_4simsMean.2IQRoutliersRemoved.2IgaussSmoothed.txt.gz
# Loading statistic matrix from /mnt/NVMeScratch/PycharmProjects/GCparagon_public/validation/01_transformed_Griffin_bias_matrices/GCparagon_matrices_per_preset/preset2/P01_observed_attributes_matrix.txt.gz
# the sum of all relative fragment length occurrences was: 1.0000 (1.0 expected)
# GCparagon-preset2 - P01: average weight across all GC content percentages: 1.25
#
# Process finished with exit code 0