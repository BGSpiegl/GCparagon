#!/usr/bin/env python3
import math
import re
import sys
import gzip
import logging
import argparse
import numpy as np
import pandas as pd
from pathlib import Path
import plotly.express as px
from plotly import graph_objs as go
from plotly.subplots import make_subplots
from collections import defaultdict
from scipy.stats import kurtosis, skew
from typing import List, Tuple, Optional, Union as OneOf, Dict

# fix for /tmp environment not being available to user (e.g., on very restrictive HPC environments) ----------
import plotly.io as pio
pio.kaleido.scope.chromium_args = tuple(
    [arg for arg in pio.kaleido.scope.chromium_args if arg != "--disable-dev-shm-usage"])  # remove flag
# ------------------------------------------ END FIX ---------------------------------------------------------


code_root = Path(__file__).parent.parent
if code_root not in sys.path:
    sys.path.append(str(code_root))


from utilities.gc_logging import log, gib_cmd_logger


def get_cmdline_args():
    """
    :return:
    """
    commandline_parser = argparse.ArgumentParser()
    # define argument groups
    input_args = commandline_parser.add_argument_group('Input (required)')
    output_args = commandline_parser.add_argument_group('Output options')
    processing_args = commandline_parser.add_argument_group('Processing options')
    input_args.add_argument('-m', '--observation-matrices', dest='observation_matrices', required=True, nargs='+',
                            help="Path(s) to one or more '*_observations_matrix.txt.gz' file(s).")
    output_args.add_argument('-o', '--output-directory', dest='output_dir',
                             help="Path to an output directory. If not specified, the parent directory of the "
                                  "correction matrix will be used. Will be created if does not exist. Output file will "
                                  "be overwritten if it already exists.")
    processing_args.add_argument('-av', '--absolute-values', dest='normalize_datasets', action='store_false',
                                 help='Flag which deactivates normalizing datasets based on their size if specified. '
                                      'Plotted values would then be absolute counts instead of the default percentage '
                                      'of overall dataset.')
    processing_args.add_argument('-dsz', '--dont-strip-zero-counts', dest='strip_zero_counts', action='store_false',
                                 help='Flag which lets the fragment length x-axis as received. By default, x-axis '
                                      'is reduced to lengths for which any received sample exhibits a non-zero count.')
    return commandline_parser.parse_args()


def load_txt_to_dataframe(file_list: List[str]) -> pd.DataFrame:
    fragment_occurrences_per_sample = None
    for matrix_file in file_list:
        sample_id = '_'.join(Path(matrix_file).stem.split('_')[:2])
        try:
            if sample_id in fragment_occurrences_per_sample.columns:
                print(f"WARNING: sample '{sample_id}' already present in data frame! Existing sample will be skipped..")
                continue
        except AttributeError:  # first sample -> fragment_occurrences_per_sample = None
            pass
        pd_matrix, frag_len_range = load_txt_to_matrix_with_meta(filename=matrix_file, loading_logger=None)
        len_max = frag_len_range.stop  # 550 bp default
        try:
            total_max = max(len_max, fragment_occurrences_per_sample.index.stop)
        except AttributeError:  # first sample
            total_max = len_max
        # check if first sample
        if fragment_occurrences_per_sample is None:
            fragment_occurrences_per_sample = pd.DataFrame(np.zeros(len_max), columns=[sample_id])
            fragment_occurrences_per_sample[sample_id] = pd.Series(
                list(pd_matrix.sum(axis=1)), index=range(frag_len_range.start, len_max + 1))  # include max. f-length!
            fragment_occurrences_per_sample = fragment_occurrences_per_sample.fillna(0)
            continue
        # check if new sample has more fragment length entries than current DataFrame and append values if necessary
        if fragment_occurrences_per_sample.index.stop + 1 < len_max:  # update
            additional_flengths = {}.fromkeys(fragment_occurrences_per_sample.columns)
            for sample_name in additional_flengths.keys():
                additional_flengths[sample_name] = np.zeros(len_max-fragment_occurrences_per_sample.index.stop)
            fragment_occurrences_per_sample.append(pd.DataFrame(additional_flengths), ignore_index=True)  # + zero rows
        # add additional sample values (new column)
        new_sample_values = pd.DataFrame(np.zeros(total_max), columns=[sample_id])
        new_sample_values[sample_id] = pd.Series(
            list(pd_matrix.sum(axis=1)), index=range(frag_len_range.start, len_max + 1))  # include max. f-length!
        new_sample_values = new_sample_values.fillna(0)
        fragment_occurrences_per_sample[sample_id] = new_sample_values
    return fragment_occurrences_per_sample


def load_txt_to_matrix(filename: OneOf[str, Path], loading_logger: Optional[str],
                       to_dtype=np.float64) -> np.array:
    """
    :param filename:
    :param to_dtype:
    :return:
    """
    if loading_logger is not None:
        log(message=f"Loading statistic matrix from {filename}", log_level=logging.INFO, logger_name=loading_logger)
    else:
        print(f"Loading statistic matrix from {filename}")
    statistic_matrix = np.loadtxt(filename, delimiter='|', skiprows=0, dtype=to_dtype)
    return statistic_matrix


def load_txt_to_matrix_with_meta(filename: OneOf[str, Path], loading_logger: Optional[str] = None,
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
    statistic_matrix = np.loadtxt(filename, delimiter='|', skiprows=0, dtype=to_dtype)
    with gzip.open(str(filename), 'rt') as f_mat:
        hdr = f_mat.readline()
    elements = hdr.split('# rows representing fragment lengths (')[1].split()  # split on whitespace + trim empty fields
    fragment_lengths = range(int(elements[0]), int(elements[3]))  # 20 bp to 550 bp (non-Pythonic in header)
    assert statistic_matrix.shape[0] == fragment_lengths.stop - fragment_lengths.start + 1  # check rows
    assert statistic_matrix.shape[1] == fragment_lengths.stop + 1  # check columns
    return statistic_matrix, fragment_lengths


def plot_fragment_length_dists(matrix_data_frame: Optional[pd.DataFrame], sample_id: Optional[str],
                               matrix_file_list: Optional[List[OneOf[str, Path]]], out_dir_path: Path, image_formats=('png',),
                               parent_logger: Optional[str] = None, fig_width=1500, fig_height=1000, fig_fontsize=24,
                               expected_fragment_lengths=(167, 167+149), normalize_to_dataset_size=True,
                               strip_xaxis_end_zeros=True, show_figure=False):
    if matrix_data_frame is None:  # load from file list!
        if matrix_file_list is None:
            log(message="either a matrix in pd.DataFrame format or a list of matrix file paths must be "
                        "supplied.", log_level=logging.ERROR, logger_name=parent_logger, flush=True, close_handlers=True)
            raise AttributeError
        matrix_data_frame = load_txt_to_dataframe(file_list=matrix_file_list)
    samples = matrix_data_frame.columns
    single_sample = len(samples) == 1
    if strip_xaxis_end_zeros:
        # remove starting all-zero lines
        while matrix_data_frame.iloc[0, :].sum() == 0:
            matrix_data_frame = matrix_data_frame.iloc[1:, :]
        # remove terminal all-zero lines
        while matrix_data_frame.iloc[-1, :].sum() == 0:
            matrix_data_frame = matrix_data_frame.iloc[:-1, :]
    # create total count sum for all samples
    sample_total_counts = {}.fromkeys(samples)
    for sample in samples:
        sample_total_counts[sample] = [matrix_data_frame[sample].sum()]
    # normalize datasets if requested
    if normalize_to_dataset_size:  # divide each dataset by number of total fragments and multiply with 100
        for sample in samples:
            matrix_data_frame[sample] /= matrix_data_frame[sample].sum() / 100.
    y_label = f"{'relative frequency / ' if normalize_to_dataset_size else 'count / '}" + \
              ('%' if normalize_to_dataset_size else '1')
    length_distribution_fig = px.line(matrix_data_frame, template="simple_white", width=fig_width, height=fig_height,
                                      labels={'index': 'fragment length / bp',
                                              'value': y_label})
    length_distribution_fig.update_layout(showlegend=False, font_family="Ubuntu", font_size=fig_fontsize,
                                          title_font_size=fig_fontsize+2,
                                          title={'text': 'Observed Fragment Lengths for Sample' +
                                                 (': ' if single_sample else 's: ') + ', '.join(samples) +
                                                 (f' ({int(tuple(sample_total_counts.values())[0][0]):,} fragments)'
                                                  if single_sample else ''),
                                                 'font': {'family': 'Ubuntu', 'size': 28, 'color': 'rgb(20, 20, 20)'},
                                                 'xanchor': 'center', 'yanchor': 'middle', 'x': 0.5})
    length_distribution_fig.data[0].update(line={'width': 3})
    length_distribution_fig.update_xaxes(showgrid=True, dtick=50, gridwidth=1.5, gridcolor='rgb(220, 220, 220)',
                                         minor_showgrid=True, minor_dtick=25, minor_gridwidth=1, minor_griddash='dash',
                                         minor_gridcolor='rgb(220, 220, 220)')
    # add zero line
    x_min = length_distribution_fig.data[0].x.min()
    x_max = length_distribution_fig.data[0].x.max()
    y_min = length_distribution_fig.data[0].y.min()
    y_max = length_distribution_fig.data[0].y.max()
    length_distribution_fig.add_trace(go.Scatter(x=(x_min, x_max), y=(y_min, y_min),
                                                 line={'color': 'rgb(45, 45, 45)', 'width': 1.5}, mode='lines'))
    added_annotations = 0
    if expected_fragment_lengths is not None:
        for expected_fragment_length in expected_fragment_lengths:
            try:
                upper_line_end = length_distribution_fig.data[0].y[
                    length_distribution_fig.data[0].x == expected_fragment_length][0]
                if upper_line_end < (y_max-y_min)*0.06:
                    raise ValueError  # line would be too small
            except (ValueError, IndexError):
                continue  # skip this annotation - not in figure range
            # add expected max at 167 bp an 316 bp:
            length_distribution_fig.add_trace(go.Scatter(x=(expected_fragment_length, expected_fragment_length),
                                                         y=(y_min, upper_line_end), mode='lines',
                                                         line={'color': 'rgb(255, 45, 45)', 'width': 1.5}))
            # add annotation at 167 bp and 316 bp:
            length_distribution_fig.add_annotation(x=expected_fragment_length+15, y=(y_max-y_min)*0.1,
                                                   text=f'{expected_fragment_length} bp', showarrow=False,
                                                   font={'color': 'rgb(255, 45, 45)', 'family': 'Ubuntu', 'size': 22})
            added_annotations += 1
    # swap traces to have zero line on bottom
    length_distribution_fig.data = tuple([length_distribution_fig.data[1]] +
                                         [length_distribution_fig.data[2+i]
                                          for i in range(added_annotations)] +
                                         [length_distribution_fig.data[0]])
    if show_figure:
        length_distribution_fig.show()
    for image_format in image_formats:
        out_file = (out_dir_path /
                    f"{((samples[0] + '_') if single_sample else '') if sample_id is None else sample_id}"
                    f".fragment_length_distribution.{image_format}")
        length_distribution_fig.write_image(out_file)


def plot_gc_dists(original_gc_data: Dict[str, np.array], corrected_gc_data: Dict[str, np.array],
                  out_dir_path: Path, fig_width=1500, fig_height=1000, fig_fontsize=24,
                  spline_interpolation=True, normalize_to_dataset_size=True, annotation=None, reduced_bins=True,
                  reads='both', reference_dist=OneOf[Dict[str, float], defaultdict[float], None],
                  reference_normalized=True, show_figure=False, image_formats=('png',)):
    reads = reads.lower()
    if reads not in ('r1', 'r2', 'both'):
        raise AttributeError(f"invalid use of reads: '{reads}'. Attribute reads must be one of 'r1', 'r2', or 'both'.")
    orig_samples = sorted(list(original_gc_data.keys()))
    corrected_samples = sorted(list(corrected_gc_data.keys()))
    assert orig_samples == corrected_samples
    if reads != 'both':
        orig_samples = list(filter(lambda x: x is not None,
                                   [smp if reads in smp.lower() else None for smp in orig_samples]))
    # create total count sum for all samples
    sample_orig_total_counts = {}.fromkeys(orig_samples)
    sample_corr_total_counts = {}.fromkeys(orig_samples)
    # normalize datasets if requested
    if normalize_to_dataset_size:  # divide each dataset by number of total fragments and multiply with 100
        for sample in orig_samples:
            sample_orig_total_counts[sample] = original_gc_data[sample].sum()
            sample_corr_total_counts[sample] = corrected_gc_data[sample].sum()
            original_gc_data[sample] /= sample_orig_total_counts[sample]
            original_gc_data[sample] *= 100.
            corrected_gc_data[sample] /= sample_corr_total_counts[sample]
            corrected_gc_data[sample] *= 100.
        if reference_dist is not None and not reference_normalized:  # normalize if is not normalized!
            sample_corr_total_counts['expected GC'] = sum(list(reference_dist.values))
            for gc in reference_dist.keys():
                reference_dist[gc] /= sample_corr_total_counts['expected GC']
    plot_data_list = []
    columns = ['sample', 'status', 'GC', 'relative frequency percent']
    gc_values = range(0, len(original_gc_data[orig_samples[0]]), 2) \
        if reduced_bins else range(0, len(original_gc_data[orig_samples[0]]), 1)  # ASSERTS all bins present!
    if reference_dist is not None:
        if reduced_bins:  # left included; right excluded
            plot_data_list.extend([['hg38 expected GC', 'original', gc,
                                    reference_dist[gc] + (reference_dist[gc + 1] if gc != 100 else 0)]
                                   for gc in gc_values])
        else:
            plot_data_list.extend([['hg38 expected GC', 'original', gc, reference_dist[gc]] for gc in gc_values])
    for status in ('original', 'corrected'):  # plot corrected on top of original
        for sample in orig_samples:
            if reads != 'both':
                if reads not in sample.lower():  # skip 'sample' of wrong read
                    continue
            if reduced_bins:  # left included; right excluded
                plot_data_list.extend([[sample, status, gc,
                                        (original_gc_data[sample][gc] +
                                         (original_gc_data[sample][gc+1] if gc != 100 else 0))
                                        if 'original' == status else
                                        (corrected_gc_data[sample][gc] +
                                         (corrected_gc_data[sample][gc+1] if gc != 100 else 0))]
                                       for gc in gc_values])
            else:
                plot_data_list.extend([[sample, status, gc,
                                        original_gc_data[sample][gc]
                                        if 'original' == status else corrected_gc_data[sample][gc]]
                                       for gc in gc_values])
    plot_data = pd.DataFrame(plot_data_list, columns=columns)
    length_distribution_fig = px.line(plot_data, x='GC', y='relative frequency percent', line_dash='status',
                                      color='sample', line_dash_sequence=['dash', 'solid'],
                                      line_shape='spline' if spline_interpolation else 'linear',
                                      template="simple_white", width=fig_width, height=fig_height,
                                      labels={'index': 'GC content / %', 'value': 'relative frequency / 1'},
                                      title=f"Original vs. Corrected GC Content" +
                                            (' (' + annotation + ')') if annotation is not None else '',
                                      color_discrete_map={
                                          'hg38 expected GC': '#252525',
                                          'B01_R1': '#3362ff', 'B01_R2': '#3362ff',
                                          'C01_R1': '#ff730c', 'C01_R2': '#ff730c',
                                          'H01_R1': '#00b212', 'H01_R2': '#00b212',
                                          'P01_R1': '#ff0000', 'P01_R2': '#ff0000'})
    for dat_idx in range(len(length_distribution_fig.data)):
        trace_name = length_distribution_fig.data[dat_idx].name  # .split(', ')[0]  # only if single preset per plot
        cur_status = length_distribution_fig.data[dat_idx].name.split(', ')[1]
        avg_gc = (plot_data[(plot_data['sample'] == trace_name) * (plot_data['status'] == cur_status)]['GC'] *
                  plot_data[(plot_data['sample'] == trace_name) *
                            (plot_data['status'] == cur_status)]['relative frequency percent']).sum() / 100.
        if reference_dist is not None and \
                length_distribution_fig.data[dat_idx].name == 'hg38 expected GC, original':  # expected distribution
            length_distribution_fig.data[dat_idx].line.color = 'black'
            length_distribution_fig.data[dat_idx].line.width = 6
            length_distribution_fig.data[dat_idx].line.dash = 'dot'
            length_distribution_fig.data[dat_idx].name = f"GRCh38, 150 bp reads ({avg_gc:.1f}% avg. GC)"
        elif length_distribution_fig.data[dat_idx].line.dash == 'solid':
            length_distribution_fig.data[dat_idx].line.width = 3
            length_distribution_fig.data[dat_idx].name = f"{re.sub('_R2', '', re.sub('_R1', '', trace_name))} " \
                                                         f"({avg_gc:.1f}% GC)"
        else:  # original, uncorrected values; dashed
            length_distribution_fig.data[dat_idx].line.width = 2
            length_distribution_fig.data[dat_idx].name = f"{re.sub('_R2', '', re.sub('_R1', '', trace_name))} " \
                                                         f"({avg_gc:.1f}% GC)"
    length_distribution_fig.update_layout(font_family="Ubuntu", font_size=fig_fontsize,
                                          legend={'orientation': 'h', 'xanchor': 'center', 'yanchor': 'top',
                                                  'x': 0.5, 'y': -0.2, 'title': ''})
    if show_figure:
        length_distribution_fig.show()
    out_dir_path.mkdir(parents=True, exist_ok=True)
    for image_format in image_formats:
        out_file = out_dir_path / (f"GCparagon_GC_content_comparison_pre-post_correction_"
                                   f"{re.sub(', ', '_', annotation)}.{image_format}")
        try:
            length_distribution_fig.write_image(out_file)
        except:
            print(f"WARNING - could not save figure '{out_file.stem}' in '{image_format}' image format! Continuing ..")


def plot_ref_gc_content(data_to_plot: Dict[str, Dict[str, np.array]], transparencies: Dict[str, float],
                        figure_title: str, signal_colors: Dict[str, Tuple[int, int, int]],
                        output_file_path: Path, fig_width=1500, fig_height=1000, fig_fontsize=24,
                        y_is_percentage=True, show_figure=False):
    color_map = {}
    data_frame_lines = []
    for data_id, plot_data_signals in data_to_plot.items():
        for signal_id, signal in plot_data_signals.items():
            midpoint_idx = int(len(signal) // 2)
            data_frame_lines.extend([[f'{signal_id}, {data_id}', abs_pos-midpoint_idx, val]
                                     for abs_pos, val in enumerate(signal)])
            sig_col = signal_colors[signal_id]
            sig_r, sig_g, sig_b = sig_col
            try:
                sig_alph = transparencies[data_id]
            except KeyError:
                sig_alph = 1.
            color_map.update({f'{signal_id}, {data_id}': f'rgba({sig_r}, {sig_g}, {sig_b}, {sig_alph})'})
    if data_frame_lines is None:
        return 1
    figure_data = pd.DataFrame(data_frame_lines,
                               columns=['gene group, processing', 'relative position / bp', 'GC percentage / %'])
    loci_group_gc_fig = px.line(figure_data,
                                x='relative position / bp', color='gene group, processing',
                                y='GC percentage / %' if y_is_percentage else 'GC content / 1',
                                template="simple_white", width=fig_width, height=fig_height,
                                title=figure_title, color_discrete_map=color_map)
    # change details
    for dat_idx in range(len(loci_group_gc_fig.data)):
        trace_name = loci_group_gc_fig.data[dat_idx].name  # .split(', ')[0]  # only if single preset per plot
        if 'original' in trace_name:
            loci_group_gc_fig.data[dat_idx].line.width = 2
        loci_group_gc_fig.data[dat_idx].name = re.sub('hamming', 'Hamming',
                                                      re.sub(', original', '',
                                                             re.sub('_', ' ', trace_name)))
    loci_group_gc_fig.update_layout(showlegend=True, font_family="Ubuntu", font_size=fig_fontsize,
                                    legend={'orientation': 'h', 'xanchor': 'center', 'yanchor': 'top',
                                            'x': 0.5, 'y': -0.2, 'title': ''})
    loci_group_gc_fig.update_xaxes(showgrid=True, dtick=250, gridwidth=2)
    loci_group_gc_fig.update_yaxes(showgrid=True, gridwidth=2)
    if show_figure:
        loci_group_gc_fig.show()
    output_file_path.parent.mkdir(parents=True, exist_ok=True)
    loci_group_gc_fig.write_image(output_file_path)
    return 0


def plot_fragment_gc_dists(original_gc_data: Dict[str, np.array], corrected_gc_data: Dict[str, np.array],
                           out_dir_path: Path,
                           reference_dists: defaultdict[str, OneOf[Dict[int, float], defaultdict[float], None]],
                           markers=True, fig_width=1200, fig_height=800, fig_fontsize=30, show_figure=False,
                           normalize_to_dataset_size=True, annotation=None, reduced_bins=True,
                           spline_interpolation=True, reference_normalized=True):
    # sample_R1\tDATA
    orig_samples = sorted(list(original_gc_data.keys()))
    corrected_samples = sorted(list(corrected_gc_data.keys()))
    # assert orig_samples == corrected_samples  # not valid for all presets in 1 plot
    # create total count sum for all samples
    sample_orig_total_counts = {}.fromkeys(orig_samples)
    sample_corr_total_counts = {}.fromkeys(corrected_samples)
    # normalize datasets if requested
    if normalize_to_dataset_size:  # divide each dataset by number of total fragments and multiply with 100
        for sample in orig_samples:
            sample_orig_total_counts[sample] = original_gc_data[sample].sum()
            original_gc_data[sample] /= sample_orig_total_counts[sample]
            original_gc_data[sample] *= 100.
        for sample in corrected_samples:
            sample_corr_total_counts[sample] = corrected_gc_data[sample].sum()
            corrected_gc_data[sample] /= sample_corr_total_counts[sample]
            corrected_gc_data[sample] *= 100.
        if reference_dists is not None and not reference_normalized:  # normalize if is not normalized!
            for sample in reference_dists.keys():
                sample_corr_total_counts[f'{sample} simulated'] = sum(list(reference_dists[sample].values()))
                for gc in reference_dists[sample].keys():
                    reference_dists[sample][gc] /= sample_corr_total_counts[f'{sample} simulated']
    plot_data_list = []
    columns = ['sample', 'status', 'GC', 'relative frequency percent']
    gc_values = range(0, len(original_gc_data[orig_samples[0]]), 2) \
        if reduced_bins else range(0, len(original_gc_data[orig_samples[0]]), 1)  # ASSERTS all bins present!
    # ADD SIMULATED REFERENCE LINES TO PLOT DATAFRAME PRECURSOR
    if reference_dists is not None:
        if reduced_bins:  # left included; right excluded
            for sample in reference_dists.keys():
                plot_data_list.extend([[f'{sample} simulated', 'original', gc, reference_dists[sample][gc] +
                                        (reference_dists[sample][gc + 1] if gc != 100 else 0)]
                                       for gc in gc_values])
        else:
            for sample in reference_dists.keys():
                plot_data_list.extend([[f'{sample} simulated', 'original', gc,
                                        reference_dists[sample][gc]] for gc in gc_values])
    for status in ('original', 'corrected'):  # plot corrected on top of original
        for sample in orig_samples:
            if reduced_bins:  # left included; right excluded
                plot_data_list.extend([[(', '.join(sample.split(', ')[:-1]) if 'original' == status else sample),
                                        status, gc,
                                        ((original_gc_data[sample][gc] +
                                         (original_gc_data[sample][gc+1] if gc != 100 else 0))
                                         if sample == orig_samples[0] else None)
                                        if 'original' == status else
                                        (corrected_gc_data[sample][gc] +
                                         (corrected_gc_data[sample][gc+1] if gc != 100 else 0))]
                                       for gc in gc_values])
                plot_data_list = list(filter(lambda e: e is not None, plot_data_list))
            else:
                plot_data_list.extend([[sample, status, gc,
                                        original_gc_data[sample][gc]
                                        if 'original' == status else corrected_gc_data[sample][gc]]
                                       for gc in gc_values])
    plot_data = pd.DataFrame(plot_data_list, columns=columns)
    length_distribution_fig = px.line(plot_data, x='GC', y='relative frequency percent', line_dash='status',
                                      color='sample', line_dash_sequence=['longdash', 'solid'],
                                      line_shape='spline' if spline_interpolation else 'linear',
                                      template="simple_white", width=fig_width, height=fig_height,
                                      title=f"Original vs. Corrected GC Content" +
                                            (' (' + annotation + ')') if annotation is not None else '',
                                      color_discrete_map={'B01 simulated': 'rgba(37, 37, 37, 1)',
                                                          'H01 simulated': 'rgba(37, 37, 37, 1)',
                                                          'C01 simulated': 'rgba(37, 37, 37, 1)',
                                                          'P01 simulated': 'rgba(37, 37, 37, 1)',
                                                          'B01, Preset 1': 'rgba(51, 98, 255, 0.7)',
                                                          'B01, Preset 2': 'rgba(0, 178, 18, 0.75)',
                                                          'B01, Preset 3': 'rgba(255, 0, 0, 0.7)',
                                                          'C01, Preset 1': 'rgba(51, 98, 255, 0.7)',
                                                          'C01, Preset 2': 'rgba(0, 178, 18, 0.75)',
                                                          'C01, Preset 3': 'rgba(255, 0, 0, 0.7)',
                                                          'H01, Preset 1': 'rgba(51, 98, 255, 0.7)',
                                                          'H01, Preset 2': 'rgba(0, 178, 18, 0.75)',
                                                          'H01, Preset 3': 'rgba(255, 0, 0, 0.7)',
                                                          'P01, Preset 1': 'rgba(51, 98, 255, 0.7)',
                                                          'P01, Preset 2': 'rgba(0, 178, 18, 0.75)',
                                                          'P01, Preset 3': 'rgba(255, 0, 0, 0.7)',
                                                          'B01': 'rgba(51, 98, 255, 0.7)',
                                                          'C01': 'rgba(255, 115, 12, 0.75)',
                                                          'H01': 'rgba(0, 178, 18, 0.75)',
                                                          'P01': 'rgba(255, 0, 0, 0.7)'})

    ref_iter = 0
    ref_markers = (134, 133, 101, 102)  # ('x-thin-open','cross-thin-open',  'square-open', 'diamond-open')
    ref_dashes = ('solid', 'dot', 'dash', 'dashdot')
    for dat_idx in range(len(length_distribution_fig.data)):
        trace_name = ', '.join(length_distribution_fig.data[dat_idx].name.split(', ')[:-1])
        cur_status = length_distribution_fig.data[dat_idx].name.split(', ')[-1]
        avg_gc = (plot_data[(plot_data['sample'] == trace_name) * (plot_data['status'] == cur_status)]['GC'] *
                  plot_data[(plot_data['sample'] == trace_name) *
                            (plot_data['status'] == cur_status)]['relative frequency percent']).sum() / 100.
        if reference_dists is not None and \
                'simulated' in length_distribution_fig.data[dat_idx].name:  # ADJUST ATTRIBUTES OF SIMULATED LINES
            if markers:
                length_distribution_fig.data[dat_idx].line.dash = 'solid'
                length_distribution_fig.data[dat_idx].mode = 'lines+markers'
                length_distribution_fig.data[dat_idx].marker.symbol = ref_markers[ref_iter]
                length_distribution_fig.data[dat_idx].marker.size = 10
                length_distribution_fig.data[dat_idx].marker.line.width = 1.5
            else:
                length_distribution_fig.data[dat_idx].line.dash = ref_dashes[ref_iter]
            length_distribution_fig.data[dat_idx].line.width = 2
            length_distribution_fig.data[dat_idx].name = f"{trace_name} ({avg_gc:.1f}% GC)"
            ref_iter += 1
        elif length_distribution_fig.data[dat_idx].line.dash == 'solid':  # ADJUST ATTRIBUTES OF GC BIAS-CORRECTED LINES
            length_distribution_fig.data[dat_idx].line.width = 2.5
            length_distribution_fig.data[dat_idx].name = f"{trace_name} corrected ({avg_gc:.1f}% GC)"
        else:  # ADJUST ATTRIBUTES OF ORIGINAL, UNCORRECTED LINES
            length_distribution_fig.data[dat_idx].line.width = 1.5
            length_distribution_fig.data[dat_idx].name = f"{trace_name} ({avg_gc:.1f}% GC)"
    length_distribution_fig.update_layout(font_family="Ubuntu", font_size=fig_fontsize,
                                          xaxis_title='fragment GC content / %', yaxis_title='relative frequency / %',
                                          legend={'orientation': 'h', 'xanchor': 'center', 'yanchor': 'top',
                                                  'x': 0.5, 'y': -0.2, 'title': ''})
    length_distribution_fig.update_xaxes(range=[10, 85])
    if show_figure:
        length_distribution_fig.show()
    out_dir_path.mkdir(parents=True, exist_ok=True)
    out_file = out_dir_path / \
        (('individual_samples_and_presets_separate/' if len(orig_samples) == 1 else '') +
        (f'{orig_samples[0]}_' if len(orig_samples) == 1 else '') +
         f"GCparagon_GC-content-comparison_GC-bias-correction{'_SPLINE' if spline_interpolation else '_LINEAR'}"
         f"{'_MARKERS' if markers else ''}{('_' + re.sub(', ', '_', annotation)) if annotation else ''}"
         f"_cfDNAref.png")
    length_distribution_fig.write_image(out_file)


def visualize_weights(region_weights: np.array, out_dir: OneOf[str, Path], sample_label: str,
                      show_figure=False, image_formats=('png',), compute_skew=True, compute_curtosis=True,
                      fig_width=1800, fig_height=1000, fig_fontsize=30, max_weight: float = 10.):
    plot_data = pd.DataFrame({'interval weights': region_weights / max(region_weights) * max_weight})  # transform to plot-able object
    weights_fig = px.histogram(plot_data, x="interval weights",  # create histogram
                               title='Weights of Genomic Interval to approximate Genomic GC Content',
                               histnorm='percent',  #  If “percent” / “probability”, the span of each bar corresponds to
                               # the percentage / fraction of occurrences with respect to the total number of sample
                               # points (here, the sum of all bin HEIGHTS equals 100% / 1)
                               width=fig_width, height=fig_height, nbins=20)
    x_max_value = weights_fig.data[0].x.max()
    n_weights = len(region_weights)
    region_weights *= n_weights
    if compute_skew:
        dat_skew = skew(region_weights)
        weights_fig.add_annotation(dict(font=dict(size=fig_fontsize - fig_fontsize / 8),
                                        x=x_max_value * 0.5,
                                        y=n_weights / 10 * 0.1, showarrow=False,
                                        text=f"skew = {dat_skew:.2e}",
                                        xanchor='center'))
    if compute_curtosis:
        dat_kurt = kurtosis(region_weights)
        weights_fig.add_annotation(dict(font=dict(size=fig_fontsize - fig_fontsize/8),
                                        x=x_max_value * 0.9,
                                        y=n_weights / 10 * 0.1, showarrow=False,
                                        text=f"kurtosis = {dat_kurt:.2e}",
                                        xanchor='center'))
    weights_fig.update_layout(font_family="Ubuntu", font_size=fig_fontsize, showlegend=False,
                              xaxis_title='genomic interval weights for weighted mean of weight matrices / 1',
                              yaxis_title='interval count / 1',
                              title={'text': 'Weights of Genomic Interval to approximate Genomic GC Content',
                                     'font': {'family': 'Ubuntu', 'size': fig_fontsize + 6,
                                              'color': 'rgb(20, 20, 20)'},
                                     'xanchor': 'center', 'yanchor': 'middle', 'x': 0.5}
                              )
    if show_figure:
        weights_fig.show()
    out_dir_path = Path(out_dir)
    out_dir_path.mkdir(parents=True, exist_ok=True)  # ensure existence of output directory
    for image_format in image_formats:
        out_file = out_dir_path / f"{sample_label}.preselectedRegionWeights_refGCreconstruction.{image_format}"
        try:
            weights_fig.write_image(out_file)
        except:
            print(f"WARNING - could not save figure '{out_file.stem}' in '{image_format}' image format! Continuing ..")


def visualize_reconstruction_result(out_dir: OneOf[str, Path], target_dist: np.array, reconstructed_dist: np.array,
                                    original_dist: np.array, abs_err: float, sample_id: str, fig_width=1800,
                                    fig_height=1100, fig_fontsize=30, image_formats=('png',), show: bool = False,
                                    reduced_bins: bool = True, spline_interpolation: bool = True,
                                    residual_dists: Optional[OneOf[List[np.array], np.array]] = None):
    if not isinstance(residual_dists, list):
        residual_dists = [residual_dists]
    plot_data_list = []
    columns = ['name', 'fragment GC content / %',
               'relative frequency / %']
    # create GC content range integers (either in steps of 1 or 2) according to reduced_bins setting!
    gc_values = range(0, len(target_dist), 2) \
        if reduced_bins else range(0, len(target_dist), 1)  # ASSERTS all bins present!
    # add reference distribution
    if target_dist is not None:
        # insanity check: create normalized target distribution if is insane
        if not (0.9999 < target_dist.sum() < 1.0001):
            target_dist /= target_dist.sum()
        if reduced_bins:  # left included; right excluded
            plot_data_list.extend(
                [['expected reference GC', gc,
                  (target_dist[gc] + (target_dist[gc + 1] if gc != 100 else target_dist[gc]))]
                 for gc in gc_values])
        else:
            plot_data_list.extend([['expected reference GC', gc, target_dist[gc]] for gc in gc_values])
    # add intervals data
    if reduced_bins:  # left included; right excluded
        plot_data_list.extend([['original', gc,
                                original_dist[gc] + (original_dist[gc + 1] if gc != 100 else 0)]
                               for gc in gc_values])
        plot_data_list.extend([['reconstructed', gc,
                                reconstructed_dist[gc] + (reconstructed_dist[gc + 1] if gc != 100 else 0)]
                               for gc in gc_values])
    else:
        plot_data_list.extend([['original', gc, original_dist[gc]]
                               for gc in gc_values])
        plot_data_list.extend([['reconstructed', gc, reconstructed_dist[gc]]
                               for gc in gc_values])
    # create plot
    gc_dist_fig = make_subplots(rows=1, cols=1, vertical_spacing=0.03,
                                shared_xaxes=False, shared_yaxes=False,
                                subplot_titles=('',),
                                start_cell='top-left',
                                specs=[[{"secondary_y": True, 'colspan': 1, 'rowspan': 1}]])
    # SECONDARY Y-AXIS TRACE FIRST:
    # add zero line to secondary axis
    gc_dist_fig.add_trace(go.Scatter(x=[0, 100], y=[0, 0], mode='lines', name=''), row=1, col=1, secondary_y=True)
    # create residual distributions for secondary axis (actually, only 1 expected)
    residual_fig = px.line(pd.DataFrame(([[f"residual error{' ' + str(res_idx) if len(residual_dists) > 1 else ''}",
                                           gc, residual_dist[gc] + (residual_dist[gc + 1] if gc != 100 else 0)]
                                          for res_idx, residual_dist in enumerate(residual_dists) for gc in gc_values]
                                         if reduced_bins else
                                         [[f"residual error{' ' + str(res_idx) if len(residual_dists) > 1 else ''}",
                                           gc, residual_dist[gc]]
                                          for res_idx, residual_dist in enumerate(residual_dists) for gc in gc_values]
                                         ),
                                        columns=['name', 'fragment GC content / %',
                                                 'relative frequency difference / %']),
                           x='fragment GC content / %', color='name', width=fig_width, height=fig_height,
                           line_shape='spline' if spline_interpolation else 'linear',
                           y='relative frequency difference / %', template="simple_white")
    for residual_trace in residual_fig.data:  # add all residual traces
        gc_dist_fig.add_trace(residual_trace, row=1, col=1, secondary_y=True)
    # PRIMARY Y-AXIS TRACES:
    for line_trace in px.line(pd.DataFrame(plot_data_list, columns=columns), x='fragment GC content / %',
                              y='relative frequency / %', template="simple_white", width=fig_width, height=fig_height,
                              line_shape='spline' if spline_interpolation else 'linear', color='name',
                              title=f"Fragment GC Content - Reconstructed from Genomic Intervals vs. Genome").data:
        gc_dist_fig.add_trace(line_trace, row=1, col=1, secondary_y=False)
    # customize FGCD lines:
    for dat_idx in range(len(gc_dist_fig.data)):
        if target_dist is not None and \
                gc_dist_fig.data[dat_idx].name == 'expected reference GC':  # expected distribution
            gc_dist_fig.data[dat_idx].line.color = 'red'
            gc_dist_fig.data[dat_idx].line.width = 5
        elif gc_dist_fig.data[dat_idx].name == 'original':
            gc_dist_fig.data[dat_idx].line.color = 'rgba(80, 80, 80, 0.8)'
            gc_dist_fig.data[dat_idx].line.width = 3
        elif (gc_dist_fig.data[dat_idx].name != '' and   # ignore the secondary axis trace(s)
              'residual error' not in gc_dist_fig.data[dat_idx].name):  # vanilla residuals on the secondary axis
            gc_dist_fig.data[dat_idx].line.color = 'rgb(0, 0, 0)'
            gc_dist_fig.data[dat_idx].line.width = 4
    # add residual error annotation:
    gc_dist_fig.add_annotation(dict(font=dict(
        size=fig_fontsize - fig_fontsize / 8,
        color=gc_dist_fig.data[[trc.name == 'residual error' for trc in gc_dist_fig.data].index(True)].line.color),
        x=gc_values.stop * 0.8, y=0.08 if reduced_bins else 0.04, showarrow=False, text=f"AES = {abs_err:.4f}",
        xanchor='center'))
    # add mean GC content annotations:
    mean_target_gc = sum([gc_pc * rel_freq for gc_pc, rel_freq in enumerate(target_dist)]) / len(target_dist)
    gc_dist_fig.add_annotation(dict(font=dict(
        size=fig_fontsize - fig_fontsize / 8,
        color=gc_dist_fig.data[[trc.name == 'expected reference GC' for trc in gc_dist_fig.data].index(True)
                               ].line.color),
        x=gc_values.stop * 0.8, y=0.07 if reduced_bins else 0.035, showarrow=False,
        text=f"avg. reference GC% = {mean_target_gc:.3%}", xanchor='center'))
    mean_original_gc = sum([gc_pc * rel_freq for gc_pc, rel_freq in enumerate(original_dist)]) / len(original_dist)
    gc_dist_fig.add_annotation(dict(font=dict(
        size=fig_fontsize - fig_fontsize / 8,
        color=gc_dist_fig.data[[trc.name == 'original' for trc in gc_dist_fig.data].index(True)
                               ].line.color), x=gc_values.stop * 0.8, y=0.065 if reduced_bins else 0.0325,
        showarrow=False, text=f"avg. original GC% = {mean_original_gc:.3%}", xanchor='center'))
    mean_reconstructed_gc = sum([gc_pc * rel_freq
                                 for gc_pc, rel_freq in enumerate(reconstructed_dist)]) / len(reconstructed_dist)
    gc_dist_fig.add_annotation(dict(font=dict(
        size=fig_fontsize - fig_fontsize / 8,
        color=gc_dist_fig.data[[trc.name == 'reconstructed' for trc in gc_dist_fig.data].index(True)
                               ].line.color), x=gc_values.stop * 0.8, y=0.06 if reduced_bins else 0.03,
        showarrow=False, text=f"avg. reconstructed GC% = {mean_reconstructed_gc:.3%}", xanchor='center'))

    # Update figure
    gc_dist_fig.update_layout(font_family="Ubuntu", font_size=fig_fontsize, showlegend=True,
                              width=fig_width, height=fig_height,
                              title_font_size=fig_fontsize + 2, template='simple_white',
                              title={'text': f"Fragment GC Content - Reconstructed from Genomic Intervals vs. Genome",
                                     'font': {'family': 'Ubuntu', 'size': fig_fontsize + 6, 'color': 'rgb(20, 20, 20)'},
                                     'xanchor': 'center', 'yanchor': 'middle', 'x': 0.5},
                              legend=dict(orientation='h', xanchor='center', yanchor='top', x=0.5, y=-0.1,
                                          traceorder='reversed', font={'size': fig_fontsize - 2}))
    # add axis labels
    gc_dist_fig.update_yaxes(title_text=f"relative frequency{' (2% GC bins)' if reduced_bins else ''} / fraction",
                             row=1, col=1)
    gc_dist_fig.update_yaxes(title_text='residual error / 1', row=1, col=1, secondary_y=True)
    gc_dist_fig.update_xaxes(title_text='fragment GC content / %', row=1, col=1)
    # adjust axis limits
    gc_limits = [None, None]
    res_limits = [None, None]
    for trc_idx, trace_dat in enumerate(gc_dist_fig.data):
        if trace_dat.name == '':  # customize the zero line
            gc_dist_fig.data[trc_idx].line.color = 'rgb(0,0,0)'
            gc_dist_fig.data[trc_idx].line.width = 1
            gc_dist_fig.data[trc_idx].showlegend = False
            continue
        elif 'residual error' in trace_dat.name:  # signals plotted on primary y-axis
            if res_limits[0] is None or res_limits[0] > trace_dat.y.min():
                res_limits = [trace_dat.y.min(), res_limits[1]]
            if res_limits[1] is None or res_limits[1] < trace_dat.y.max():
                res_limits = [res_limits[0], trace_dat.y.max()]
        else:  # find limits for secondary y-axis
            if gc_limits[0] is None or gc_limits[0] > trace_dat.y.min():
                gc_limits = [trace_dat.y.min(), gc_limits[1]]
            if gc_limits[1] is None or gc_limits[1] < trace_dat.y.max():
                gc_limits = [gc_limits[0], trace_dat.y.max()]
    # set y-ranges for primary and secondary axes
    if None not in res_limits:
        res_limits_range = res_limits[1] - res_limits[0]
        res_limits = (res_limits[0] - 0.01 * res_limits_range, res_limits[1] + 0.01 * res_limits_range)
        gc_dist_fig.update_yaxes(range=res_limits, row=1, col=1, secondary_y=True)
    if None not in gc_limits:
        if gc_limits[0] < 0:  # make non-negative minimum
            gc_limits[0] = 0
        gc_limits_range = gc_limits[1] - gc_limits[0]
        gc_limits = (gc_limits[0] - 0.01 * gc_limits_range, gc_limits[1] + 0.01 * gc_limits_range)
        gc_dist_fig.update_yaxes(range=gc_limits, row=1, col=1, secondary_y=False)
    if show:
        gc_dist_fig.show()
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    for image_format in image_formats:
        out_file = out_dir / f"{sample_id}.GCcontent_RefVsReconstructed.{image_format}"
        try:
            gc_dist_fig.write_image(out_file)
        except:
            print(f"WARNING - could not save figure '{out_file.stem}' in '{image_format}' image format! Continuing ..")


def main() -> int:
    cmd_args = get_cmdline_args()
    observations_matrix_paths = cmd_args.observation_matrices
    output_dir = cmd_args.output_dir
    normalize_datasets = cmd_args.normalize_datasets
    strip_zero_counts = cmd_args.strip_zero_counts
    test_logger = gib_cmd_logger()
    # read input and care for existence of output directory
    observation_matrix = load_txt_to_dataframe(file_list=observations_matrix_paths)
    if output_dir:
        output_dir = Path(output_dir)
    else:
        output_dir = Path(observations_matrix_paths[0]).parent
        log(f"output will be written to path '{output_dir}'", log_level=logging.INFO, logger_name=test_logger)
    output_dir.mkdir(parents=True, exist_ok=True)
    # create line plot
    plot_fragment_length_dists(matrix_data_frame=observation_matrix, matrix_file_list=None, out_dir_path=output_dir,
                               normalize_to_dataset_size=normalize_datasets, sample_id=None, show_figure=False,
                               strip_xaxis_end_zeros=strip_zero_counts, parent_logger=test_logger, fig_fontsize=32)
    return 0


if __name__ == '__main__':
    sys.exit(main())
