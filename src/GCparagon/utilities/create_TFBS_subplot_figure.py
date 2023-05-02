#!/usr/bin/env python3

import re
import sys
import math
import numpy as np
import polars as pl
import pandas as pd
from pathlib import Path
import plotly_express as px
from collections import deque
import plotly.graph_objects as go
from twobitreader import TwoBitFile
from plotly.subplots import make_subplots
from typing import Union, Dict, Tuple

SOURCE_ROOT_DIR = Path(__file__).parent.parent  # ../src/GCparagon
if SOURCE_ROOT_DIR not in sys.path:
    sys.path.append(str(SOURCE_ROOT_DIR))
CONTENT_ROOT_DIR = SOURCE_ROOT_DIR.parent.parent

# ref GC content definitions:
PERCENTAGE_PLOT = True
two_bit_reference_file = SOURCE_ROOT_DIR / '2bit_reference/hg38.analysisSet.2bit'
target_regions_per_class = {'LYL1': CONTENT_ROOT_DIR / 'accessory_files/TFBSs/'
                            'LYL1.sorted.gtrd_version_21_12.10000_sites.chrY-chrM-chrMT_excluded.hg38.bed',
                            'GRHL2': CONTENT_ROOT_DIR / 'accessory_files/TFBSs/'
                            'GRHL2.sorted.gtrd_version_21_12.10000_sites.chrY-chrM-chrMT_excluded.hg38.bed'}
TFBS_COLORS = {'LYL1': (110, 0, 110), 'GRHL2': (0, 110, 110)}
SAMPLE_COLORS = {'B01': 'rgba(51, 98, 255, 1)',
                 'C01': 'rgba(255, 115, 12, 1)',
                 'H01': 'rgba(0, 178, 18, 1)',
                 'P01': 'rgba(255, 0, 0, 1)'}

# TFBSs DoC definitions:
tfbs_coverage_csv_paths = {'stem': CONTENT_ROOT_DIR / 'accessory_files/TFBSs/coverages',
                           'relative': ('GRHL2_original', 'GRHL2_corrected', 'LYL1_corrected', 'LYL1_original')}
output_path = Path(CONTENT_ROOT_DIR / 'accessory_files/TFBSs')
output_file_path = output_path / 'DoC_bias_correction_effect_TFBSs.png'
REF_GC_DUMP = output_path / 'TFBS_ref_gc_content.npydct'
# created in the first iteration; MUST be deleted if the reference regions in tfbs_coverage_csv_paths change!


def load_bed_regions(region_list_dict: Dict[str, Union[str, Path]]) -> Dict[str, Tuple[str, int, int]]:
    gene_data = {}.fromkeys(region_list_dict)
    for list_name, list_path in region_list_dict.items():
        with open(list_path, mode='rt') as f_gene_list:
            gene_data.update({list_name: tuple(map(lambda t: (t[0], int(t[1]), int(t[2])),
                                                   filter(lambda e: e != [''],
                                                          map(lambda l: l.strip().split('\t')[:3],
                                                              f_gene_list.readlines()))))})
    return gene_data


def get_region_gc_content(reference_regions_dict: Dict[str, Tuple[str, int, int]],
                          window_size=10001, as_percentage=True) -> Dict[str, np.array]:
    gc_arrays = {}.fromkeys(reference_regions_dict)
    actual_window_size = int(window_size // 2) * 2 + 1
    if actual_window_size != window_size:
        print(f"-> WARNING :  window size specified was {window_size:,} bp but will use "
              f"{actual_window_size:,} bp instead.")
    for target_id, region_tuple in reference_regions_dict.items():
        print(f" i :  processing region group '{target_id}'\n------------------------------------")
        chromosome_wise_loci = {}
        for chrom, start, stop in region_tuple:
            try:
                chromosome_wise_loci[chrom].append((start, stop))
            except KeyError:
                chromosome_wise_loci.update({chrom: [(start, stop)]})
        # sort chromosome-wise lists
        for chrom, loci_list in chromosome_wise_loci.items():
            chromosome_wise_loci[chrom] = sorted(loci_list, key=lambda t: t[0], reverse=False)
        # extract GC per position
        reference_handle = TwoBitFile(two_bit_reference_file)
        chrom_handles = {}
        held_handles = deque()
        max_handle_num = 2
        total_loci = 0
        gc_window = np.zeros(actual_window_size)
        for chrom, loci_list in chromosome_wise_loci.items():
            try:
                chrom_handle = chrom_handles[chrom]
            except KeyError:
                chrom_handles.update({chrom: reference_handle[chrom]})
                chrom_handle = chrom_handles[chrom]
                held_handles.append(chrom)
                if len(held_handles) > max_handle_num:
                    del chrom_handles[held_handles.popleft()]
            for region_start, region_end in loci_list:
                center_pos = region_start + int((region_end - region_start) / 2.)
                # get uppercase reference sequence
                tss_seq = chrom_handle[center_pos - int(window_size // 2):
                                       center_pos + int(window_size // 2) + 1].upper()
                # update numpy array
                gc_window += tuple(map(lambda c: int(c in 'GC'), tss_seq))  # adds either 0 or 1 for every position
                total_loci += 1
        # compute GC content per relative position:
        gc_window /= total_loci
        if as_percentage:
            gc_window *= 100.
        gc_arrays[target_id] = gc_window
    # write to file
    np_default_print_thresh = np.get_printoptions()['threshold']
    np.set_printoptions(threshold=100000)  # otherwise, the repr() will use summarisation
    with open(REF_GC_DUMP, mode='w') as f_ref_dump:
        f_ref_dump.write(repr(gc_arrays))
    np.set_printoptions(threshold=np_default_print_thresh)  # set to default again
    return gc_arrays


def smooth_1d_data(x=None, window_len=5, window='hamming'):
    # adapted from Scipy cookbook: https://scipy-cookbook.readthedocs.io/items/SignalSmooth.html; accessed 05-01-2021
    # Smoothing of a 1D signal
    # Date:	2017-07-13 (last modified), 2006-10-31 (created)
    # Section author: Unknown[1], GaelVaroquaux, Unknown[142], Unknown[143], Unknown[144], Unknown[145],
    # Unknown[146], Unknown[147], WesTurner, Christian Gagnon, clecocel
    """smooth the data using a window with requested size.
        This method is based on the convolution of a scaled window with the signal.
        The signal is prepared by introducing reflected copies of the signal
        (with the window size) in both ends so that transient parts are minimized
        in the begining and end part of the output signal.

        input:
            x: the input signal
            window_len: the dimension of the smoothing window; should be an odd integer
            window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
                flat window will produce a moving average smoothing.

        output:
            the smoothed signal

        example:
        t=linspace(-2,2,0.1)
        x=sin(t)+randn(len(t))*0.1
        y=smooth(x)

        see also:
        numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
        scipy.signal.lfilter

        NOTE: length(output) != length(input), to correct this:
              return y[floor(window_len/2-1):-ceil(window_len/2)] instead of just y.
    """
    if not isinstance(x, np.ndarray):
        x = np.array(x)
    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")

    if window_len < 3:
        return x
    if not window_len % 2:  # force odd window length
        window_len += 1

    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")

    s = np.r_[x[window_len - 1:0:-1], x, x[-2:-window_len - 1:-1]]
    if window == 'flat':  # moving average
        w = np.ones(window_len, 'd')
    else:
        w = eval('np.' + window + '(window_len)')

    y = np.convolve(w / w.sum(), s, mode='valid')
    return y[math.floor(window_len / 2 - 1):-math.ceil(window_len / 2)]


def get_multi_signals_plot(plot_data_signals: Dict[str, np.array], data_category: str, x_axis_data_label: str,
                           y_axis_data_label: str, signal_colors=None, annotation=None) -> Union[go.Figure, int]:
    """
    Assumes centered data
    """
    if signal_colors is None:
        signal_colors = {}
    color_map = {}
    data_frame_lines = []
    signal_id_ordered = sorted(list(plot_data_signals.keys()), reverse=True)
    for signal_id in signal_id_ordered:  # original, corrected
        signal = plot_data_signals[signal_id]
        midpoint_idx = int(len(signal) // 2)
        data_frame_lines.extend([[f'{data_category} {signal_id}', abs_pos - midpoint_idx, val]
                                 for abs_pos, val in enumerate(signal)])
        if color_map.get(annotation) is None:
            try:
                sig_col = signal_colors[annotation]
            except KeyError:
                sig_col = (185, 185, 185)  # light grey
            color_map.update({annotation: sig_col})
    if data_frame_lines is None:
        return 1
    figure_data = pd.DataFrame(data_frame_lines,
                               columns=['TF correction status', x_axis_data_label, y_axis_data_label])
    loci_group_gc_fig = px.line(figure_data, x=x_axis_data_label, y=y_axis_data_label, color='TF correction status',
                                template="simple_white", line_shape='linear')
    # change details
    for dat_idx in range(len(loci_group_gc_fig.data)):
        tf, correction_status = loci_group_gc_fig.data[dat_idx].name.split(' ')  # only if single preset per plot
        match correction_status:
            case 'original':
                loci_group_gc_fig.data[dat_idx].line.width = 2
            case 'corrected':
                loci_group_gc_fig.data[dat_idx].line.width = 3.5
        loci_group_gc_fig.data[dat_idx].line.color = color_map[annotation]  # use sample for color
        loci_group_gc_fig.data[dat_idx].name = f'{annotation} {loci_group_gc_fig.data[dat_idx].name}'
    return loci_group_gc_fig


def ref_gc_trace_name_replacemenet_func(trc_nm: str) -> str:
    return re.sub('hamming', 'Hamming', re.sub(', original', '', re.sub('_', ' ', trc_nm)))


def plot_ref_gc_content(data_to_plot: Dict[str, Dict[str, np.array]], transparencies: Dict[str, float],
                        signal_colors: Dict[str, Tuple[int, int, int]], fig_fontsize=24,
                        fig_width=1500, fig_height=1000, y_is_percentage=True) -> Union[int, go.Figure]:
    color_map = {}
    data_frame_lines = []
    for data_id, plot_data_signals in data_to_plot.items():
        for signal_id, signal in plot_data_signals.items():
            midpoint_idx = int(len(signal) // 2)
            data_frame_lines.extend([[f'{signal_id}, {data_id}', abs_pos - midpoint_idx, val]
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
    loci_group_gc_fig = px.line(figure_data, x='relative position / bp', color='gene group, processing',
                                y='GC percentage / %' if y_is_percentage else 'GC content / 1',
                                template="simple_white", width=fig_width, height=fig_height,
                                color_discrete_map=color_map)
    # change details (needs to be done outside)
    for dat_idx in range(len(loci_group_gc_fig.data)):
        trace_name = loci_group_gc_fig.data[dat_idx].name
        if 'original' in trace_name:
            loci_group_gc_fig.data[dat_idx].line.width = 2
        loci_group_gc_fig.data[dat_idx].name = re.sub('hamming', 'Hamming', re.sub(', original', '',
                                                                                   re.sub('_', ' ', trace_name)))
    loci_group_gc_fig.update_layout(showlegend=True, font_family="Ubuntu", font_size=fig_fontsize,
                                    legend={'orientation': 'h', 'xanchor': 'center', 'yanchor': 'top',
                                            'x': 0.5, 'y': -0.2, 'title': ''})
    loci_group_gc_fig.update_xaxes(showgrid=True, dtick=250, gridwidth=2)
    loci_group_gc_fig.update_yaxes(showgrid=True, gridwidth=2)
    return loci_group_gc_fig


def main() -> int:
    # TFBS ref seq GC content subplot:
    # ------------------------------------------------------------------------------------------------------------------
    # read target genes
    hg38_ref_regions = load_bed_regions(region_list_dict=target_regions_per_class)
    # load all_reference_genes; accumulate GC content
    view_window = 2001
    if REF_GC_DUMP.is_file():  # fast
        with open(REF_GC_DUMP, 'r') as f_gc_dump:
            gc_per_group_raw = f_gc_dump.read()
            gc_per_group_raw = gc_per_group_raw.replace('array', 'np.array')
            gc_per_group = eval(gc_per_group_raw)
    else:  # slow
        gc_per_group = get_region_gc_content(reference_regions_dict=hg38_ref_regions,
                                             window_size=view_window, as_percentage=PERCENTAGE_PLOT)
    # create hamming window filtered versions
    hamming_window_length = 15
    transparent_signals = {'original': 0.15}
    data_to_plot = {'original': gc_per_group}
    data_to_plot.update({f'{hamming_window_length}bp_hamming_filtered': None})
    transparent_signals.update({f'{hamming_window_length}bp_hamming_filtered': 0.8})
    filtered_data = {}
    for group, np_arr in gc_per_group.items():
        filtered_data.update({group: smooth_1d_data(x=np_arr, window='hamming',
                                                    window_len=hamming_window_length)})
    data_to_plot[f'{hamming_window_length}bp_hamming_filtered'] = filtered_data

    # TFBS DoC subplot data prep:
    # ------------------------------------------------------------------------------------------------------------------
    stem_path = Path(tfbs_coverage_csv_paths['stem'])
    cov_data = {}
    for csv in tfbs_coverage_csv_paths['relative']:
        gene, signal_type = csv.split('_')
        csv_path = stem_path / csv
        csv_candidates = tuple(csv_path.glob('*.csv'))
        if len(csv_candidates) == 0:
            print(f"WARNING: no csv file found under path '{csv_path}'. Ignoring this directory ..")
            continue
        for csv_candidate in csv_candidates:
            sample_name = csv_candidate.stem
            if 'X' in sample_name:
                continue
            covs = pl.read_csv(csv_candidate)
            avg_covs = np.array(covs[:, 3:-1].mean(axis=0))[0]  # -1001 to 1000 and 1220 -> -1000, +1000
            if cov_data.get(sample_name) is None:
                cov_data.update({sample_name: {gene: {signal_type: None}}})
            elif cov_data[sample_name].get(gene) is None:
                cov_data[sample_name].update({gene: {signal_type: None}})
            elif cov_data[sample_name][gene].get(signal_type) is None:
                cov_data[sample_name][gene].update({signal_type: None})
            cov_data[sample_name][gene][signal_type] = avg_covs  # add avg. signal

    # ACTUALLY CREATE PLOTS:
    # ------------------------------------------------------------------------------------------------------------------
    # create Multiplot:
    fig = make_subplots(rows=7, cols=1, shared_xaxes=False, shared_yaxes=False, vertical_spacing=0.03,
                        subplot_titles=('average reference sequence GC content',
                                        'B01 TFBS coverage LYL1 (-2.2% GC)',
                                        'H01 TFBS coverage LYL1 (+1.1% GC)',
                                        'P01 TFBS coverage LYL1 (+5.0% GC)',
                                        'B01 TFBS coverage GRHL2 (-2.2% GC)',
                                        'H01 TFBS coverage GRHL2 (+1.1% GC)',
                                        'P01 TFBS coverage GRHL2 (+5.0% GC)'),
                        start_cell='top-left',
                        specs=[[{"secondary_y": False, 'colspan': 1, 'rowspan': 1}],
                               [{"secondary_y": True, 'colspan': 1, 'rowspan': 1}],
                               [{"secondary_y": True, 'colspan': 1, 'rowspan': 1}],
                               [{"secondary_y": True, 'colspan': 1, 'rowspan': 1}],
                               [{"secondary_y": True, 'colspan': 1, 'rowspan': 1}],
                               [{"secondary_y": True, 'colspan': 1, 'rowspan': 1}],
                               [{"secondary_y": True, 'colspan': 1, 'rowspan': 1}]])
    # 0.) LYL1/GRHL2 reference GC content
    ref_gc_content_traces = plot_ref_gc_content(data_to_plot=data_to_plot, transparencies=transparent_signals,
                                                signal_colors=TFBS_COLORS, y_is_percentage=PERCENTAGE_PLOT,
                                                fig_fontsize=24, fig_width=1000, fig_height=800)
    if ref_gc_content_traces == 1:
        raise ValueError("something went wrong during plotting of P01 LYL1")
    for ref_trace in ref_gc_content_traces.data:
        fig.add_trace(ref_trace, row=1, col=1)  # add all traces to first subplot
        # add smooth signals to DoC subplots
        if 'LYL1' in ref_trace.name and 'filtered' in ref_trace.name:
            for subp_idx in range(2, 5):
                fig.add_trace(ref_trace, row=subp_idx, col=1, secondary_y=True)
        elif 'GRHL2' in ref_trace.name and 'filtered' in ref_trace.name:
            for subp_idx in range(5, 8):
                fig.add_trace(ref_trace, row=subp_idx, col=1, secondary_y=True)
    # 1.) B01 LYL1
    b01_lyl1_traces = get_multi_signals_plot(plot_data_signals=cov_data['B01']['LYL1'], data_category='LYL1',
                                             signal_colors=SAMPLE_COLORS, x_axis_data_label='relative position / bp',
                                             annotation='B01', y_axis_data_label='average relative coverage / 1')
    if b01_lyl1_traces == 1:
        raise ValueError("something went wrong during plotting of P01 LYL1")
    for doc_trace in b01_lyl1_traces.data:
        fig.add_trace(doc_trace, row=2, col=1)
    # 2.) H01 LYL1
    h01_lyl1_traces = get_multi_signals_plot(plot_data_signals=cov_data['H01']['LYL1'], data_category='LYL1',
                                             signal_colors=SAMPLE_COLORS, x_axis_data_label='relative position / bp',
                                             annotation='H01', y_axis_data_label='average relative coverage / 1')
    if h01_lyl1_traces == 1:
        raise ValueError("something went wrong during plotting of P01 LYL1")
    for doc_trace in h01_lyl1_traces.data:
        fig.add_trace(doc_trace, row=3, col=1)
    # 3.) P01 LYL1
    p01_lyl1_traces = get_multi_signals_plot(plot_data_signals=cov_data['P01']['LYL1'], data_category='LYL1',
                                             signal_colors=SAMPLE_COLORS, x_axis_data_label='relative position / bp',
                                             annotation='P01', y_axis_data_label='average relative coverage / 1')
    if p01_lyl1_traces == 1:
        raise ValueError("something went wrong during plotting of P01 LYL1")
    for doc_trace in p01_lyl1_traces.data:
        fig.add_trace(doc_trace, row=4, col=1)
    # 4.) B01 GRHL2
    b01_grhl2_traces = get_multi_signals_plot(plot_data_signals=cov_data['B01']['GRHL2'], data_category='GRHL2',
                                              signal_colors=SAMPLE_COLORS, x_axis_data_label='relative position / bp',
                                              annotation='B01', y_axis_data_label='average relative coverage / 1')
    if b01_grhl2_traces == 1:
        raise ValueError("something went wrong during plotting of B01 GRHL2")
    for doc_trace in b01_grhl2_traces.data:
        fig.add_trace(doc_trace, row=5, col=1)
    # 5.) H01 GRHL2
    h01_grhl2_traces = get_multi_signals_plot(plot_data_signals=cov_data['H01']['GRHL2'], data_category='GRHL2',
                                              signal_colors=SAMPLE_COLORS, x_axis_data_label='relative position / bp',
                                              annotation='H01', y_axis_data_label='average relative coverage / 1')
    if h01_grhl2_traces == 1:
        raise ValueError("something went wrong during plotting of H01 GRHL2")
    for doc_trace in h01_grhl2_traces.data:
        fig.add_trace(doc_trace, row=6, col=1)
    # 6.) P01 GRHL2
    p01_grhl2_traces = get_multi_signals_plot(plot_data_signals=cov_data['P01']['GRHL2'], data_category='GRHL2',
                                              signal_colors=SAMPLE_COLORS, x_axis_data_label='relative position / bp',
                                              annotation='P01', y_axis_data_label='average relative coverage / 1')
    if p01_grhl2_traces == 1:
        raise ValueError("something went wrong during plotting of P01 GRHL2")
    for doc_trace in p01_grhl2_traces.data:
        fig.add_trace(doc_trace, row=7, col=1)
    # formatting of ref GC plot: change details
    lyl1_limits = [None, None]
    grhl2_limits = [None, None]
    ref_limits = [None, None]
    for dat_idx in range(len(fig.data)):
        trace_name = fig.data[dat_idx].name
        fig.data[dat_idx].line.dash = 'solid'
        if 'filtered' in trace_name:
            fig.data[dat_idx].showlegend = False
            if dat_idx in (0, 1):
                fig.data[dat_idx].line.width = 2
            else:  # filtered reference GC content signals
                fig.data[dat_idx].line.width = 1.5
                if dat_idx in (2, 6):
                    fig.data[dat_idx].showlegend = True
                else:
                    tf_name = trace_name.split(',')[0]
                    r, g, b = TFBS_COLORS[tf_name]
                    fig.data[dat_idx].line.color = f"rgba({r}, {g}, {b}, {0.5})"
                fig.data[dat_idx].name = f"reference GC content {fig.data[dat_idx].name.split(',')[0]} genes"
                if ref_limits[0] is None or ref_limits[0] > min(fig.data[dat_idx].y):
                    ref_limits[0] = min(fig.data[dat_idx].y)
                if ref_limits[1] is None or ref_limits[1] < max(fig.data[dat_idx].y):
                    ref_limits[1] = max(fig.data[dat_idx].y)
        elif 'original' in trace_name:
            fig.data[dat_idx].line.width = 2
            r, g, b, _o = fig.data[dat_idx].line.color.strip('rgba(').strip(')').split(', ')
            fig.data[dat_idx].line.color = f'rgba({r}, {g}, {b}, 0.5)'  # make transparent
            # update limits
            if 'LYL1' in trace_name:
                if lyl1_limits[0] is None or lyl1_limits[0] > fig.data[dat_idx].y.min():
                    lyl1_limits[0] = fig.data[dat_idx].y.min()
                if lyl1_limits[1] is None or lyl1_limits[1] < fig.data[dat_idx].y.max():
                    lyl1_limits[1] = fig.data[dat_idx].y.max()
            elif 'GRHL2' in trace_name:
                if grhl2_limits[0] is None or grhl2_limits[0] > fig.data[dat_idx].y.min():
                    grhl2_limits[0] = fig.data[dat_idx].y.min()
                if grhl2_limits[1] is None or grhl2_limits[1] < fig.data[dat_idx].y.max():
                    grhl2_limits[1] = fig.data[dat_idx].y.max()
            if 'B01 LYL1' in trace_name:
                fig.data[dat_idx].name = 'B01 original'
            elif 'H01 LYL1' in trace_name:
                fig.data[dat_idx].name = 'H01 original'
            elif 'P01 LYL1' in trace_name:
                fig.data[dat_idx].name = 'P01 original'
            else:
                fig.data[dat_idx].showlegend = False
        elif 'corrected' in trace_name:
            fig.data[dat_idx].line.width = 3.5
            # update limits
            if 'LYL1' in trace_name:
                if lyl1_limits[0] is None or lyl1_limits[0] > fig.data[dat_idx].y.min():
                    lyl1_limits[0] = fig.data[dat_idx].y.min()
                if lyl1_limits[1] is None or lyl1_limits[1] < fig.data[dat_idx].y.max():
                    lyl1_limits[1] = fig.data[dat_idx].y.max()
            elif 'GRHL2' in trace_name:
                if grhl2_limits[0] is None or grhl2_limits[0] > fig.data[dat_idx].y.min():
                    grhl2_limits[0] = fig.data[dat_idx].y.min()
                if grhl2_limits[1] is None or grhl2_limits[1] < fig.data[dat_idx].y.max():
                    grhl2_limits[1] = fig.data[dat_idx].y.max()
            if 'B01 LYL1' in trace_name:
                fig.data[dat_idx].name = 'B01 corrected'
            elif 'H01 LYL1' in trace_name:
                fig.data[dat_idx].name = 'H01 corrected'
            elif 'P01 LYL1' in trace_name:
                fig.data[dat_idx].name = 'P01 corrected'
            else:
                fig.data[dat_idx].showlegend = False
        elif dat_idx in (0, 1):
            fig.data[dat_idx].showlegend = False
    # update overall figure formatting
    fig.update_layout(height=3400, width=1200, template='simple_white', showlegend=True,
                      font={'family': 'Ubuntu', 'size': 22, 'color': 'rgb(20, 20, 20)'},
                      title={'text': 'TFBS GC Bias',
                             'font': {'family': 'Ubuntu', 'size': 30, 'color': 'rgb(20, 20, 20)'},
                             'xanchor': 'center', 'yanchor': 'middle', 'x': 0.48},
                      legend=dict(orientation='h', xanchor='center', yanchor='top', x=0.5, y=-0.025,
                                  traceorder='reversed'))
    fig.update_annotations(font={'family': 'Ubuntu', 'size': 26, 'color': 'rgb(20, 20, 20)'})  # format subplot titles
    fig.update_yaxes(showgrid=True, gridwidth=2, row=1, col=1)
    for row in range(1, 8):  # turn on y-grid
        fig.update_xaxes(showgrid=True, dtick=250, gridwidth=2, row=row, col=1)
        if row == 1:
            fig.update_yaxes(title_text='average reference GC content / %', row=row, col=1)
        else:
            fig.update_yaxes(title_text='normalized average DoC / ratio', row=row, col=1)
            fig.update_yaxes(title_text='average reference GC content / %', row=row, col=1, secondary_y=True)
            if row == 7:
                fig.update_xaxes(title_text='relative position to TFBS / bp', row=row, col=1)
    # update y ranges
    ref_limits = (ref_limits[0]-1, ref_limits[1]+1)
    lyl1_limits_range = lyl1_limits[1] - lyl1_limits[0]
    lyl1_limits = (lyl1_limits[0]-0.01*lyl1_limits_range, lyl1_limits[1]+0.01*lyl1_limits_range)
    grhl2_limits_range = grhl2_limits[1] - grhl2_limits[0]
    grhl2_limits = (grhl2_limits[0]-0.01*grhl2_limits_range, grhl2_limits[1]+0.01*grhl2_limits_range)
    for sec_idx in range(1, 8):
        if sec_idx in range(2, 5):
            fig.update_yaxes(range=lyl1_limits, row=sec_idx, col=1, secondary_y=False)
        elif sec_idx in range(5, 8):
            fig.update_yaxes(range=grhl2_limits, row=sec_idx, col=1, secondary_y=False)
        if sec_idx != 1:
            fig.update_yaxes(range=ref_limits, row=sec_idx, col=1, secondary_y=True)
    fig.show()
    output_file_path.parent.mkdir(parents=True, exist_ok=True)
    fig.write_image(output_file_path)
    return 0


if __name__ == '__main__':
    sys.exit(main())
