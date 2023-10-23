#!/usr/bin/env python3

import sys
import numpy as np
import pandas as pd
from pathlib import Path
from GCparagon.utilities.plot_distributions import load_txt_to_matrix_with_meta, plot_fragment_length_dists
from GCparagon.utilities.plot_GC_matrices import plot_statistic_matrices, limit_extreme_outliers
from GCparagon.utilities.gc_logging import gib_cmd_logger
from GCparagon.correct_GC_bias import reduce_matrix, trim_2d_matrix, \
    DEFAULT_OUTLIER_DETECTION_STRINGENCY, DEFAULT_OUTLIER_DETECTION_METHOD

SOURCE_ROOT_DIR = Path(__file__).parent.parent  # ../src/GCparagon
if SOURCE_ROOT_DIR not in sys.path:
    sys.path.append(str(SOURCE_ROOT_DIR))
# CONTENT_ROOT_DIR = SOURCE_ROOT_DIR.parent.parent

IMAGE_FORMATS = ('png', 'svg', 'pdf')

# ANALYSIS DEFINITIONS
MIN_FRAG_LENGTH = 20
MAX_FRAG_LENGTH = 550

########################################################################################################################
# TODO: SET YOUR GCparagon SOURCE DIRECTORY HERE !!!!!!!!
MY_ABSOLUTE_SOURCE_DIRECTORY_PATH = SOURCE_ROOT_DIR  # write_your_GCparagon_directory_absolute_path_here
# TODO: SET YOUR GCparagon SOURCE DIRECTORY HERE !!!!!!!!
########################################################################################################################
# check it:
if not MY_ABSOLUTE_SOURCE_DIRECTORY_PATH.is_dir():
    raise FileNotFoundError(f"Your GCparagon source ('GCparagon') directory does not exist!"
                            f"You specified: '{MY_ABSOLUTE_SOURCE_DIRECTORY_PATH}'")

matrices_parent_dir = MY_ABSOLUTE_SOURCE_DIRECTORY_PATH / \
                                          'validation/01_transformed_Griffin_bias_matrices'
plot_output_path = matrices_parent_dir / 'Griffin_matrix_plots'
plot_output_path_limited = matrices_parent_dir / 'Griffin_matrix_plots-OutliersLimited'

# process paths - glob matrices
plot_output_path.mkdir(parents=True, exist_ok=True)
plot_output_path_limited.mkdir(parents=True, exist_ok=True)
sample_bias_matrices = tuple(matrices_parent_dir.glob("*gc_weights.Griffin.txt.gz"))
sample_names_original_order = [bias_path.name.split('_')[0] for bias_path in sample_bias_matrices]
sample_bias_matrices_dict = {}.fromkeys(sample_names_original_order)
sample_observation_matrices_dict = {}.fromkeys(sample_names_original_order)
sample_target_matrices_dict = {}.fromkeys(sample_names_original_order)
sample_observation_masks_dict = {}.fromkeys(sample_names_original_order)
for bias_path in sample_bias_matrices:
    sample_bias_matrices_dict[bias_path.name.split('_')[0]] = bias_path
sample_observation_matrices = tuple(matrices_parent_dir.glob("*observed_attributes_matrix.Griffin.txt.gz"))
for observation_path in sample_observation_matrices:
    sample_observation_matrices_dict[observation_path.name.split('_')[0]] = observation_path
sample_observation_masks = tuple(matrices_parent_dir.glob("*gc_bias_computation_mask.Griffin.txt.gz"))
for mask_path in sample_observation_masks:
    sample_observation_masks_dict[mask_path.name.split('_')[0]] = mask_path
sample_target_matrices = tuple(matrices_parent_dir.glob("*target_attributes_matrix.Griffin.txt.gz"))
for target_matrix_path in sample_target_matrices:
    sample_target_matrices_dict[target_matrix_path.name.split('_')[0]] = target_matrix_path


if __name__ == '__main__':
    cmd_logger = gib_cmd_logger()
    observation_matrices = {}
    # 1) plot weight matrices
    for sample_bias_matrix in sample_bias_matrices:
        sample_name = sample_bias_matrix.name.split('_')[0]
        sample_np_bias_matrix, flen_range_weights = load_txt_to_matrix_with_meta(filename=sample_bias_matrix,
                                                                                 loading_logger=cmd_logger,
                                                                                 to_dtype=np.float64)
        sample_mask, flen_range_mask = load_txt_to_matrix_with_meta(filename=sample_observation_masks_dict[sample_name],
                                                                    loading_logger=cmd_logger, to_dtype=np.float64)
        sample_observations_matrix, flen_range_obs = load_txt_to_matrix_with_meta(
            filename=sample_observation_matrices_dict[sample_name], loading_logger=cmd_logger, to_dtype=np.uint64)
        sample_target_attributes_matrix, flen_range_trg = load_txt_to_matrix_with_meta(
            filename=sample_target_matrices_dict[sample_name], loading_logger=cmd_logger, to_dtype=np.uint64)
        deleted_rows, deleted_columns = range(0), range(0)
        # create focused Griffin_matrix_plots
        sample_mask_focused, (deleted_rows, deleted_columns) = reduce_matrix(
            matrix_to_trim=sample_mask, trim_dimensions_exclusively_containing=[False], border_elements=10)
        assert flen_range_weights == flen_range_mask and flen_range_mask == flen_range_obs and \
               flen_range_obs == flen_range_trg
        # concerning border_elements: 10 default "pixels" distance drawn
        sample_np_bias_matrix_focused = trim_2d_matrix(matrix=sample_np_bias_matrix,
                                                       rows=(deleted_rows.start, deleted_rows.stop),
                                                       columns=(deleted_columns.start, deleted_columns.stop))
        sample_observations_matrix_focused = trim_2d_matrix(matrix=sample_observations_matrix,
                                                            rows=(deleted_rows.start, deleted_rows.stop),
                                                            columns=(deleted_columns.start, deleted_columns.stop))
        observation_matrices.update({sample_name: (pd.DataFrame(sample_observations_matrix_focused),
                                                   flen_range_obs.start + deleted_rows.start)})
        sample_target_attributes_matrix_focused = trim_2d_matrix(matrix=sample_target_attributes_matrix,
                                                                 rows=(deleted_rows.start, deleted_rows.stop),
                                                                 columns=(deleted_columns.start, deleted_columns.stop))
        frq_data = {'W_gc': pd.DataFrame(sample_np_bias_matrix_focused),
                    'O_gc': pd.DataFrame(sample_observations_matrix_focused),
                    'S_gc': pd.DataFrame(sample_target_attributes_matrix_focused),
                    'Mask': pd.DataFrame(sample_mask_focused)}
        data_for_limited_outliers = {'W_gc': sample_np_bias_matrix_focused,
                                     'O_gc': sample_observations_matrix_focused,
                                     'S_gc': sample_target_attributes_matrix_focused,
                                     'Mask': sample_mask_focused}
        for data_id in frq_data.keys():
            plot_statistic_matrices(frq_data=frq_data, data_id_to_show=data_id, image_formats=IMAGE_FORMATS,
                                    y_tick_label_offset=deleted_rows.start + MIN_FRAG_LENGTH,
                                    x_tick_label_offset=deleted_columns.start, show_figure=False,
                                    in_file=str(matrices_parent_dir / f'{sample_name}.bam'),
                                    output_dir=plot_output_path, sample_id=sample_name, fig_width=1800,
                                    fig_height=2000, fig_fontsize=50, parent_logger=cmd_logger)
            # create special weights matrix plots
            if data_id == 'W_gc':
                ZERO_REPLACEMENT_VALUE = 25
                # create version with weights differences close to zero = very positive
                weights_focused_zero_marked = np.where(
                    (-0.001 < sample_np_bias_matrix_focused) & (sample_np_bias_matrix_focused < 0.001),
                    ZERO_REPLACEMENT_VALUE, sample_np_bias_matrix_focused)
                plot_statistic_matrices(frq_data={'W_gc': pd.DataFrame(weights_focused_zero_marked)},
                                        data_id_to_show=data_id, image_formats=IMAGE_FORMATS,
                                        y_tick_label_offset=deleted_rows.start + MIN_FRAG_LENGTH,
                                        x_tick_label_offset=deleted_columns.start, show_figure=False,
                                        in_file=str(matrices_parent_dir / f'{sample_name}.bam'),
                                        output_dir=plot_output_path, sample_id=f'{sample_name}-close2zeroMarked',
                                        fig_width=1800, fig_height=2000, fig_fontsize=50, parent_logger=cmd_logger)
                # create version with weights differences limited to a max. of 4
                weights_focused_4_limited = np.where(sample_np_bias_matrix_focused > 4.,
                                                     4., sample_np_bias_matrix_focused)
                plot_statistic_matrices(frq_data={'W_gc': pd.DataFrame(weights_focused_4_limited)},
                                        data_id_to_show=data_id, image_formats=IMAGE_FORMATS,
                                        y_tick_label_offset=deleted_rows.start + MIN_FRAG_LENGTH,
                                        x_tick_label_offset=deleted_columns.start, show_figure=False,
                                        in_file=str(matrices_parent_dir / f'{sample_name}.bam'),
                                        output_dir=plot_output_path, sample_id=f'{sample_name}-Diff4limited',
                                        fig_width=1800, fig_height=2000, fig_fontsize=50, parent_logger=cmd_logger)
                # create version with weights differences limited to a max. of 2
                weights_focused_2_limited = np.where(sample_np_bias_matrix_focused > 2.,
                                                     2., sample_np_bias_matrix_focused)
                plot_statistic_matrices(frq_data={'W_gc': pd.DataFrame(weights_focused_2_limited)},
                                        data_id_to_show=data_id, image_formats=IMAGE_FORMATS,
                                        y_tick_label_offset=deleted_rows.start + MIN_FRAG_LENGTH,
                                        x_tick_label_offset=deleted_columns.start, show_figure=False,
                                        in_file=str(matrices_parent_dir / f'{sample_name}.bam'),
                                        output_dir=plot_output_path, sample_id=f'{sample_name}-Diff2limited',
                                        fig_width=1800, fig_height=2000, fig_fontsize=50, parent_logger=cmd_logger)
            # create limited outliers plot:
            # plot limited distributions (Griffin limits everything that has a bias lower than 0.05
            # (extremely underrepresented attribute tuples) but to make it comparable in terms of scale, let's create
            # the GCparagon-limited samples!
            if data_id != 'Mask':  # don't create limited mask -> nonsense
                sample_np_bias_matrix_focused_limited = limit_extreme_outliers(
                    outliers_matrix=data_for_limited_outliers[data_id],
                    outliers_factor=10 - DEFAULT_OUTLIER_DETECTION_STRINGENCY,
                    detection_method=DEFAULT_OUTLIER_DETECTION_METHOD, parent_logger=None)
                plot_statistic_matrices(frq_data={data_id: pd.DataFrame(sample_np_bias_matrix_focused_limited)},
                                        data_id_to_show=data_id, image_formats=IMAGE_FORMATS,
                                        y_tick_label_offset=deleted_rows.start + MIN_FRAG_LENGTH,
                                        x_tick_label_offset=deleted_columns.start, show_figure=False,
                                        in_file=str(matrices_parent_dir / f'{sample_name}-OutliersLimited.bam'),
                                        output_dir=plot_output_path_limited, sample_id=f'{sample_name}-OutliersLimited',
                                        fig_width=1800, fig_height=2000, fig_fontsize=50, parent_logger=cmd_logger)
    # 2) plot fragment length distributions
    for sample_name, (pd_matrix, flength_offset) in observation_matrices.items():
        fragment_occurrences_per_sample = pd.DataFrame(np.zeros(MAX_FRAG_LENGTH), columns=[sample_name])
        fragment_occurrences_per_sample[sample_name] = pd.Series(
            list(pd_matrix.sum(axis=1)), index=range(flength_offset, flength_offset+pd_matrix.shape[0]))
        fragment_occurrences_per_sample = fragment_occurrences_per_sample.fillna(0)
        plot_fragment_length_dists(matrix_data_frame=fragment_occurrences_per_sample,
                                   matrix_file_list=None, image_formats=IMAGE_FORMATS,
                                   out_dir_path=plot_output_path, normalize_to_dataset_size=True, show_figure=False,
                                   strip_xaxis_end_zeros=True, parent_logger=cmd_logger, sample_id=sample_name)
