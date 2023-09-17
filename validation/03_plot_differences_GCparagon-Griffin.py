#!/usr/bin/env python3

import numpy as np
import pandas as pd
from pathlib import Path
from GCparagon.utilities.plot_distributions import load_txt_to_matrix_with_meta, plot_fragment_length_dists
from GCparagon.utilities.plot_GC_matrices import plot_statistic_matrices
from GCparagon.utilities.gc_logging import gib_cmd_logger
from GCparagon.correct_GC_bias import reduce_matrix

########################################################################################################################
# TODO: SET YOUR GCparagon SOURCE DIRECTORY HERE !!!!!!!!
MY_ABSOLUTE_SOURCE_DIRECTORY_PATH = Path('write_your_GCparagon_directory_absolute_path_here')
# TODO: SET YOUR GCparagon SOURCE DIRECTORY HERE !!!!!!!!
########################################################################################################################
# check it:
if not MY_ABSOLUTE_SOURCE_DIRECTORY_PATH.is_dir():
    raise FileNotFoundError(f"Your GCparagon source ('GCparagon') directory does not exist!"
                            f"You specified: '{MY_ABSOLUTE_SOURCE_DIRECTORY_PATH}'")

# ANALYSIS DEFINITIONS
MIN_FRAG_LENGTH = 20
MAX_FRAG_LENGTH = 550

transformed_griffin_matrices_parent_dir = MY_ABSOLUTE_SOURCE_DIRECTORY_PATH / \
                                          'validation/transformed_Griffin_bias_matrices'

gcparagon_matrices_parent_dir = MY_ABSOLUTE_SOURCE_DIRECTORY_PATH / 'GCparagon_matrices_ALL_PRESETS'

# process paths - glob matrices
# Griffin matrices:

sample_bias_matrices = tuple(transformed_griffin_matrices_parent_dir.glob("*gc_weights.Griffin.txt.gz"))
sample_names_original_order = [bias_path.name.split('_')[0] for bias_path in sample_bias_matrices]
sample_griffin_bias_matrices_dict = {}.fromkeys(sample_names_original_order)
sample_griffin_observation_matrices_dict = {}.fromkeys(sample_names_original_order)
sample_griffin_observation_masks_dict = {}.fromkeys(sample_names_original_order)
sample_griffin_target_attributes_matrix_dict = {}.fromkeys(sample_names_original_order)
# target matrix is in silico construct from GC weights and fragment frequencies!
# is the equivalent of GCParagon's simulated_attribute_matrix


for bias_path in sample_bias_matrices:
    sample_griffin_bias_matrices_dict[bias_path.name.split('_')[0]] = bias_path
sample_observation_matrices = tuple(transformed_griffin_matrices_parent_dir.glob(
    "*observed_attributes_matrix.Griffin.txt.gz"))
for observation_path in sample_observation_matrices:
    sample_griffin_observation_matrices_dict[observation_path.name.split('_')[0]] = observation_path
sample_observation_masks = tuple(transformed_griffin_matrices_parent_dir.glob(
    "*gc_bias_computation_mask.Griffin.txt.gz"))
for mask_path in sample_observation_masks:
    sample_griffin_observation_masks_dict[mask_path.name.split('_')[0]] = mask_path
sample_target_matrices = tuple(transformed_griffin_matrices_parent_dir.glob(
    "*target_attributes_matrix.Griffin.txt.gz"))
for target_matrix_path in sample_target_matrices:
    sample_griffin_target_attributes_matrix_dict[target_matrix_path.name.split('_')[0]] = target_matrix_path


# GCparagon matrices:
observed_attribute_matrices_gcparagon = {}
weight_matrices_gcparagon = {}
mask_matrices_gcparagon = {}
simulated_attribute_matrices_gcparagon = {}
gcparagon_matrices_collection = {}

for preset in range(1, 4, 1):
    preset_label = f'preset{preset}'  # preset_dir.name
    assert preset_label not in gcparagon_matrices_collection
    gcparagon_matrices_collection.update(
        {preset_label:
             {'observations': {},
              'simulations': {},
              'mask': {},
              'weights': {}}})
    weights_matrices = tuple((gcparagon_matrices_parent_dir / f'preset{preset}').glob(
        '*_gc_weights_*simsMean.*IQRoutliersRemoved.*IgaussSmoothed.txt.gz'))
    for w_mat in weights_matrices:
        sample_name = w_mat.stem.split('_')[0]
        gcparagon_matrices_collection[preset_label]['weights'].update(
            {sample_name: w_mat})
    print(gcparagon_matrices_collection)
    observations_matrices = tuple((gcparagon_matrices_parent_dir / f'preset{preset}').glob(
        '*_observed_attributes_matrix.txt.gz'))
    for o_mat in observations_matrices:
        sample_name = o_mat.stem.split('_')[0]
        gcparagon_matrices_collection[preset_label]['observations'].update(
            {sample_name: o_mat})
    print(gcparagon_matrices_collection)
    simulations_matrices = tuple((gcparagon_matrices_parent_dir / f'preset{preset}').glob(
        '*_simulated_attributes_matrix.txt.gz'))
    for s_mat in simulations_matrices:
        sample_name = s_mat.stem.split('_')[0]
        gcparagon_matrices_collection[preset_label]['simulations'].update(
            {sample_name: s_mat})
    masks = tuple((gcparagon_matrices_parent_dir / f'preset{preset}').glob('*_gc_bias_computation_mask.txt.gz'))
    print(gcparagon_matrices_collection)
    for m_mat in masks:
        sample_name = m_mat.stem.split('_')[0]
        gcparagon_matrices_collection[preset_label]['mask'].update(
            {sample_name: m_mat})
    print(f"\nINFO - results for preset {preset}:")
    print(gcparagon_matrices_collection)


skip_samples = []  # all results received from Sebastian


if __name__ == '__main__':
    # ASSUMPTION: all matrices have the same fragment length range!!!
    cmd_logger = gib_cmd_logger()
    # 1) plot differences of weight, observations, simulations, and mask matrices (Griffin - GCparagon)
    for sample_bias_matrix in sample_bias_matrices:
        sample_name = str(sample_bias_matrix.name.split('_')[0])
        if sample_name in skip_samples:
            continue
        sample_np_bias_matrix, flen_range_weights = load_txt_to_matrix_with_meta(filename=sample_bias_matrix,
                                                                                 loading_logger=cmd_logger,
                                                                                 to_dtype=np.float64)
        sample_mask, flen_range_mask = load_txt_to_matrix_with_meta(
            loading_logger=cmd_logger, filename=sample_griffin_observation_masks_dict[sample_name], to_dtype=np.float64)
        sample_observations_matrix, flen_range_obs = load_txt_to_matrix_with_meta(
            loading_logger=cmd_logger, filename=sample_griffin_observation_matrices_dict[sample_name],
            to_dtype=np.uint64)
        sample_target_attributes_matrix, flen_range_trg = load_txt_to_matrix_with_meta(
            loading_logger=cmd_logger, filename=sample_griffin_target_attributes_matrix_dict[sample_name],
            to_dtype=np.uint64)
        # load GCparagon matrices
        for preset_label in gcparagon_matrices_collection.keys():
            plot_output_path = (transformed_griffin_matrices_parent_dir /
                                f'difference_plots_GCparagon-Griffin/{preset_label}')
            plot_output_path.mkdir(parents=True, exist_ok=True)
            gcparagon_mask, flen_range_mask_gcp = load_txt_to_matrix_with_meta(
                loading_logger=cmd_logger, filename=gcparagon_matrices_collection[preset_label]['mask'][sample_name],
                to_dtype=np.float64)
            gcparagon_observations_matrix, flen_range_obs_gcp = load_txt_to_matrix_with_meta(
                loading_logger=cmd_logger, filename=gcparagon_matrices_collection[preset_label]['observations'][sample_name],
                to_dtype=np.uint64)
            gcparagon_simulated_attributes_matrix, flen_range_trg_gcp = load_txt_to_matrix_with_meta(
                loading_logger=cmd_logger, filename=gcparagon_matrices_collection[preset_label]['simulations'][sample_name],
                to_dtype=np.uint64)
            gcparagon_weights_matrix, flen_range_weights_gcp = load_txt_to_matrix_with_meta(
                loading_logger=cmd_logger, filename=gcparagon_matrices_collection[preset_label]['weights'][sample_name],
                to_dtype=np.uint64)
            assert flen_range_mask_gcp == flen_range_obs_gcp == flen_range_trg_gcp == flen_range_weights_gcp
            # check ranges against Griffin:
            assert flen_range_mask == flen_range_mask_gcp
            # compute differences; but: always extend shorter matrix to larger one!
            # MASK:
            assert sample_mask.shape == gcparagon_mask.shape
            mask_differences = sample_mask.astype(bool).astype(int) - gcparagon_mask.astype(bool).astype(int)
            # OBSERVATIONS: (scale to 100k fragments dataset size -> prevent glitching)
            theoretical_possibilities_gcp = ((flen_range_obs_gcp[1] - flen_range_obs_gcp[0]) *  # fragment lengths
                                             (flen_range_obs_gcp[1] + 1) -  # max. GC base count
                                             ((flen_range_obs_gcp[1] - flen_range_obs_gcp[0]) ** 2) / 2)
            print(f"INFO - scaling up to {theoretical_possibilities_gcp:,} fragments for sample '{sample_name}'")
            observations_differences = sample_observations_matrix * (theoretical_possibilities_gcp /
                                                                     sample_observations_matrix.sum()) - \
                                       gcparagon_observations_matrix * (theoretical_possibilities_gcp /
                                                                        gcparagon_observations_matrix.sum())
            target_differences = sample_target_attributes_matrix * (theoretical_possibilities_gcp /
                                                                    sample_target_attributes_matrix.sum()) - \
                                 gcparagon_simulated_attributes_matrix * (theoretical_possibilities_gcp /
                                                                          gcparagon_simulated_attributes_matrix.sum())
            # WEIGHTS difference:
            weights_differences = np.subtract(sample_np_bias_matrix, gcparagon_weights_matrix)
            # A) process difference matrices (CODE BLOCK ADAPTED FROM GCparagon) - focus on non-default values!
            deleted_rows_dict = {'observations': None, 'simulations': None, 'mask': None, 'weights': None}
            # create focused plots -> different focus for different plots possible
            mask_differences_focused, (deleted_mask_rows, deleted_mask_columns) = reduce_matrix(
                matrix_to_trim=mask_differences, trim_dimensions_exclusively_containing=[False], border_elements=10)
            deleted_rows_dict['Mask'] = (deleted_mask_rows, deleted_mask_columns)
            observations_differences_focused, (deleted_obs_rows, deleted_obs_columns) = reduce_matrix(
                matrix_to_trim=observations_differences, trim_dimensions_exclusively_containing=[0.],
                border_elements=10)
            deleted_rows_dict['O_gc'] = (deleted_obs_rows, deleted_obs_columns)
            target_differences_focused, (deleted_target_rows, deleted_target_columns) = reduce_matrix(
                matrix_to_trim=target_differences, trim_dimensions_exclusively_containing=[0.],
                border_elements=10)
            deleted_rows_dict['S_gc'] = (deleted_target_rows, deleted_target_columns)
            weights_differences_focused, (deleted_weight_rows, deleted_weight_columns) = reduce_matrix(
                matrix_to_trim=weights_differences, trim_dimensions_exclusively_containing=[0.], border_elements=10)
            deleted_rows_dict['W_gc'] = (deleted_weight_rows, deleted_weight_columns)
            frq_data_diff = {'W_gc': pd.DataFrame(weights_differences_focused),
                             'O_gc': pd.DataFrame(observations_differences_focused),
                             'S_gc': pd.DataFrame(target_differences_focused),
                             'Mask': pd.DataFrame(mask_differences_focused)}

            # CODE BLOCK END; Plot matrices!
            for data_id in frq_data_diff.keys():
                plot_statistic_matrices(frq_data=frq_data_diff, data_id_to_show=data_id,
                                        y_tick_label_offset=deleted_rows_dict[data_id][0].start + MIN_FRAG_LENGTH,
                                        x_tick_label_offset=deleted_rows_dict[data_id][1].start, show_figure=False,
                                        in_file=str(transformed_griffin_matrices_parent_dir / f'{sample_name}.bam'),
                                        fig_width=1800, output_dir=plot_output_path, sample_id=sample_name,
                                        fig_height=2000, fig_fontsize=50, parent_logger=cmd_logger)
                if data_id == 'W_gc':
                    ZERO_REPLACEMENT_VALUE = 25
                    # create version with weights differences close to zero = very positive
                    weights_differences_focused_zero_marked = np.where(
                        (-0.001 < weights_differences_focused) & (weights_differences_focused < 0.001),
                        ZERO_REPLACEMENT_VALUE, weights_differences_focused)
                    plot_statistic_matrices(frq_data={'W_gc': pd.DataFrame(weights_differences_focused_zero_marked)},
                                            data_id_to_show=data_id,
                                            y_tick_label_offset=deleted_rows_dict[data_id][0].start + MIN_FRAG_LENGTH,
                                            x_tick_label_offset=deleted_rows_dict[data_id][1].start, show_figure=False,
                                            in_file=str(transformed_griffin_matrices_parent_dir / f'{sample_name}.bam'),
                                            output_dir=plot_output_path,
                                            sample_id=f'{sample_name}-close2zeroMarked',
                                            fig_width=1800, fig_height=2000, fig_fontsize=50, parent_logger=cmd_logger)
                    # create version with weights differences limited to a max. of 4
                    weights_differences_focused_4_limited = np.where(weights_differences_focused > 4.,
                                                                     4., weights_differences_focused)
                    plot_statistic_matrices(frq_data={'W_gc': pd.DataFrame(weights_differences_focused_4_limited)},
                                            data_id_to_show=data_id,
                                            y_tick_label_offset=deleted_rows_dict[data_id][0].start+MIN_FRAG_LENGTH,
                                            x_tick_label_offset=deleted_rows_dict[data_id][1].start, show_figure=False,
                                            in_file=str(transformed_griffin_matrices_parent_dir / f'{sample_name}.bam'),
                                            output_dir=plot_output_path,
                                            sample_id=f'{sample_name}-Diff4limited',
                                            fig_width=1800, fig_height=2000, fig_fontsize=50, parent_logger=cmd_logger)
                    # create version with weights differences limited to a max. of 2
                    weights_differences_focused_2_limited = np.where(weights_differences_focused > 2.,
                                                                     2., weights_differences_focused)
                    plot_statistic_matrices(frq_data={'W_gc': pd.DataFrame(weights_differences_focused_2_limited)},
                                            data_id_to_show=data_id,
                                            y_tick_label_offset=deleted_rows_dict[data_id][0].start + MIN_FRAG_LENGTH,
                                            x_tick_label_offset=deleted_rows_dict[data_id][1].start, show_figure=False,
                                            in_file=str(transformed_griffin_matrices_parent_dir / f'{sample_name}.bam'),
                                            output_dir=plot_output_path,
                                            sample_id=f'{sample_name}-Diff2limited',
                                            fig_width=1800, fig_height=2000, fig_fontsize=50, parent_logger=cmd_logger)
            fragment_occurrences_per_sample = pd.DataFrame(np.zeros(MAX_FRAG_LENGTH), columns=[sample_name])
            fragment_occurrences_per_sample[sample_name] = pd.Series(
                list(frq_data_diff['O_gc'].sum(axis=1)),
                index=range(MIN_FRAG_LENGTH + deleted_rows_dict['O_gc'][0].start,
                            MIN_FRAG_LENGTH + deleted_rows_dict['O_gc'][0].start + frq_data_diff['O_gc'].shape[0]))
            fragment_occurrences_per_sample = fragment_occurrences_per_sample.fillna(0)
            plot_fragment_length_dists(matrix_data_frame=fragment_occurrences_per_sample, fig_fontsize=32,
                                       expected_fragment_lengths=None,
                                       matrix_file_list=None,
                                       out_dir_path=plot_output_path, normalize_to_dataset_size=False,
                                       show_figure=False, strip_xaxis_end_zeros=True, parent_logger=cmd_logger,
                                       sample_id=f'{sample_name}_Differences_GriffinVsGCparagon')
