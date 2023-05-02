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
from typing import Union, Dict, Tuple, List

StrandedLocus = Tuple[str, int, Union[int, bool]]

SOURCE_ROOT_DIR = Path(__file__).parent.parent  # ../src/GCparagon
if SOURCE_ROOT_DIR not in sys.path:
    sys.path.append(str(SOURCE_ROOT_DIR))
CONTENT_ROOT_DIR = SOURCE_ROOT_DIR.parent.parent

# ref GC content definitions:
PERCENTAGE_PLOT = True
two_bit_reference_file = SOURCE_ROOT_DIR / '2bit_reference/hg38.analysisSet.2bit'
target_region_library_bed_path = Path(CONTENT_ROOT_DIR / 'accessory_files/TSSs/hg38_MANE.bed')
target_regions_per_class = {'HK': CONTENT_ROOT_DIR / 'accessory_files/TSSs/HK.txt',
                            'PAU': CONTENT_ROOT_DIR / 'accessory_files/TSSs/PAU.txt'}

TSS_COLORS = {'HK': 'rgba(190, 50, 0, 1)',
              'PAU': 'rgba(150, 150, 150, 1)'}
SAMPLE_COLORS = {'B01': 'rgba(51, 98, 255, 1)',
                 'C01': 'rgba(255, 115, 12, 1)',
                 'H01': 'rgba(0, 178, 18, 1)',
                 'P01': 'rgba(255, 0, 0, 1)'}

# TSSs DoC definitions:
tss_coverage_csv_paths = {'stem': CONTENT_ROOT_DIR / 'accessory_files/TSSs/coverages',
                          'relative': ('TSS_corrected', 'TSS_original')}
output_path = Path(CONTENT_ROOT_DIR / 'accessory_files/TSSs')
output_file_path = output_path / 'DoC_bias_correction_effect_TSSs.png'
REF_GC_DUMP = output_path / 'TSS_ref_gc_content.npydct'
# created in the first iteration;  MUST be deleted if the reference regions in tss_coverage_csv_paths change!


def load_bed_regions(region_list_dict: Dict[str, Union[str, Path]], library_file=None, stranded=False) \
        -> Tuple[Dict[str, Tuple[Tuple[str, int, Union[int, bool]]]], bool]:
    gene_data = {}
    if library_file is None:  # load regions directly from paths as defined in region_list_dict
        for list_name, list_path in region_list_dict.items():
            with open(list_path, mode='rt') as f_gene_list:
                gene_data.update({list_name: tuple(map(lambda t: (t[0], int(t[1]), int(t[2])),
                                                       filter(lambda e: e != [''],
                                                              map(lambda l: l.strip().split('\t')[:3],
                                                                  f_gene_list.readlines()))))})
    else:  # need to load library and then filter for the regions in region_list_dict
        gene_library = {}
        with open(library_file, mode='rt') as f_gene_list:
            if stranded:
                _ = deque(map(lambda t: gene_library.update({t[3].upper():
                                                            ((t[0], int(t[1]), False) if t[5] == '+' else
                                                             (t[0], int(t[2]), True))}),  # needs to be flipped later?
                              filter(lambda e: e != [''],
                                     map(lambda l: l.strip().split('\t')[:6], f_gene_list.readlines()))), maxlen=0)
            else:
                _ = deque(map(lambda t: gene_library.update({t[3].upper(): (t[0], int(t[1]), int(t[2]))}),
                              filter(lambda e: e != [''],
                                     map(lambda l: l.strip().split('\t')[:4], f_gene_list.readlines()))), maxlen=0)
        if len(gene_library) == 0:
            raise ValueError("-> ERROR :  gene library is empty! Exiting..")
        # read the gene lists and retrieve genomic coordinates from the library
        for list_name, list_path in region_list_dict.items():
            with open(list_path, mode='rt') as f_gene_list:
                gene_data[list_name] = tuple(map(lambda g: gene_library[g],
                                                               filter(lambda g: g != '' and gene_library.get(g.upper()),
                                                                      map(lambda l: l.strip(),
                                                                          f_gene_list.readlines()))))
            with open(list_path, mode='rt') as f_gene_list:
                genes_not_found = sum(tuple(map(lambda g: 0 if gene_library.get(g) else 1,
                                                filter(lambda s: s != '',
                                                       map(lambda l: l.strip(), f_gene_list.readlines())))))
            print(f"info: {len(gene_data[list_name]):,} genes found for list '{list_name}'")
            total_genes = len(gene_data[list_name]) + genes_not_found
            print(f"info: {genes_not_found:,} out of {total_genes:,} genes (= {genes_not_found/total_genes:%}) NOT "
                  f"found for list '{list_name}'")
    return gene_data, stranded


def get_region_gc_content_around_center(reference_regions_dict: Dict[str, Tuple[StrandedLocus]],
                                        stranded: bool, window_size=10001, as_percentage=True) -> Dict[str, np.array]:
    gc_arrays = {}.fromkeys(reference_regions_dict)
    actual_window_size = int(window_size // 2) * 2 + 1
    if actual_window_size != window_size:  # if even window size was specified, odd one is used
        print(f"-> WARNING :  window size specified was {window_size:,} bp but will use "
              f"{actual_window_size:,} bp instead.")
    half_window_size_without_center_base = actual_window_size // 2
    for target_id, region_tuple in reference_regions_dict.items():
        print(f" i :  processing region group '{target_id}'\n------------------------------------")
        chromosome_wise_loci = {}
        if stranded:
            for chrom, start, need_to_flip in region_tuple:
                region_to_add = (start-half_window_size_without_center_base,
                                 start+half_window_size_without_center_base+1, need_to_flip)
                try:
                    chromosome_wise_loci[chrom].append(region_to_add)
                except KeyError:
                    chromosome_wise_loci.update({chrom: [region_to_add]})
        else:
            for chrom, start, stop in region_tuple:
                try:
                    chromosome_wise_loci[chrom].append((start, stop))
                except KeyError:
                    chromosome_wise_loci.update({chrom: [(start, stop)]})
        # sort chromosome-wise lists
        for chrom, loci_list in chromosome_wise_loci.items():
            chromosome_wise_loci[chrom] = sorted(loci_list, key=lambda t: t[0], reverse=False)  # sort ascending
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
            if stranded:
                for region_start, region_end, flip in loci_list:
                    center_pos = region_start + int((region_end - region_start) / 2.)
                    # get uppercase reference sequence
                    tss_seq = chrom_handle[center_pos - int(window_size // 2):
                                           center_pos + int(window_size // 2) + 1].upper()
                    if flip:
                        tss_seq = tss_seq[::-1]
                    # update numpy array
                    gc_window += tuple(map(lambda c: int(c in 'GC'), tss_seq))  # adds either 0 or 1 for every position
                    total_loci += 1
            else:
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
    np.set_printoptions(threshold=100000)  # otherwise, the repr() will use summarization
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


def get_multi_signals_plot(plot_data_signals: Dict[str, Dict[str, np.array]],
                           x_axis_data_label: str, y_axis_data_label: str, annotation=None, signal_colors=None,
                           signal_dash=None, transparencies=None, figure_width=1500, figure_height=1000) \
        -> Union[List[go.Scattergl], int]:
    """
    Assumes centered data
    """
    if transparencies is None:
        transparencies = {}
    if signal_colors is None:
        signal_colors = {}
    if signal_dash is None:
        signal_dash = {}
    color_map = {}
    data_frame_lines = []
    signal_id_ordered = sorted(list(plot_data_signals.keys()), reverse=True)
    for signal_id in signal_id_ordered:  # HK, PAU
        for status_id in plot_data_signals[signal_id]:  # original, corrected
            signal = plot_data_signals[signal_id][status_id]  # np.array
            # change name of trace using gene group name (signal_id)
            if signal_dash.get(status_id) is None:
                signal_dash[status_id] = 'solid'
            midpoint_idx = int(len(signal) // 2)
            data_frame_lines.extend([[annotation, signal_id, status_id, abs_pos - midpoint_idx, val]
                                     for abs_pos, val in enumerate(signal)])
            if color_map.get(signal_id) is None:
                try:
                    sig_col = signal_colors[signal_id]
                except KeyError:
                    sig_col = (31, 119, 180)
                sig_r, sig_g, sig_b = sig_col
                try:
                    sig_alph = transparencies[signal_id]
                except KeyError:
                    sig_alph = 1.
                color_map.update({signal_id: f'rgba({sig_r}, {sig_g}, {sig_b}, {sig_alph})'})
    if data_frame_lines is None:
        return 1
    figure_data = pd.DataFrame(data_frame_lines,
                               columns=['sample', 'signal', 'status', x_axis_data_label, y_axis_data_label])
    current_figure = px.line(figure_data, x=x_axis_data_label, y=y_axis_data_label, width=figure_width,
                             color='signal', color_discrete_map=color_map, template="simple_white",
                             height=figure_height, line_shape='linear', title=annotation,
                             line_dash='status', line_dash_map=signal_dash)
    return [trc for trc in current_figure.data]


def ref_gc_trace_name_replacemenet_func(trc_nm: str) -> str:
    return re.sub('hamming', 'Hamming', re.sub(', original', '', re.sub('_', ' ', trc_nm)))


def plot_ref_gc_content(data_to_plot: Dict[str, Dict[str, np.array]], transparencies: Dict[str, float],
                        signal_colors: Dict[str, str], fig_fontsize=24,
                        fig_width=1500, fig_height=1000, y_is_percentage=True) -> Union[int, go.Figure]:
    color_map = {}
    data_frame_lines = []
    for data_id, plot_data_signals in data_to_plot.items():
        for signal_id, signal in plot_data_signals.items():
            midpoint_idx = int(len(signal) // 2)
            data_frame_lines.extend([[f'{signal_id}, {data_id}', abs_pos - midpoint_idx, val]
                                     for abs_pos, val in enumerate(signal)])
            sig_col = signal_colors[signal_id]
            try:
                sig_alph = transparencies[data_id]
                sig_col = f"{', '.join(sig_col.split(', ')[:3])}, {sig_alph})"
            except KeyError:
                pass
            color_map.update({f'{signal_id}, {data_id}': sig_col})
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


def sort_regions(regions_dict: Dict[str, Tuple[StrandedLocus]]) -> Dict[str, Tuple[StrandedLocus]]:
    for region_class, regions in regions_dict.items():
        regions_dict[region_class] = sorted(regions, key=lambda x: (x[0], x[1]))
    return regions_dict


def main() -> int:
    # TSS ref seq GC content subplot:
    # ------------------------------------------------------------------------------------------------------------------
    if REF_GC_DUMP.is_file():  # fast
        with open(REF_GC_DUMP, 'r') as f_gc_dump:
            gc_per_group_raw = f_gc_dump.read()
            gc_per_group_raw = gc_per_group_raw.replace('array', 'np.array')
            gc_per_group = eval(gc_per_group_raw)
    else:  # slow
        # read target genes
        hg38_ref_regions, stranded_data = load_bed_regions(region_list_dict=target_regions_per_class,
                                                           library_file=target_region_library_bed_path, stranded=True)
        hg38_ref_regions = sort_regions(regions_dict=hg38_ref_regions)
        view_window = 2001
        # load all_reference_genes; accumulate GC content
        gc_per_group = get_region_gc_content_around_center(reference_regions_dict=hg38_ref_regions,
                                                           stranded=stranded_data, window_size=view_window,
                                                           as_percentage=PERCENTAGE_PLOT)
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

    # TSS DoC subplot data prep:
    # ------------------------------------------------------------------------------------------------------------------
    stem_path = Path(tss_coverage_csv_paths['stem'])
    # get genes per class for signal averaging
    genes_per_class = {}
    for region_class, gene_list_path in target_regions_per_class.items():
        with open(gene_list_path, 'rt') as f_genes:
            genes_per_class.update({region_class: tuple(filter(lambda s: s != '',
                                                               map(lambda l: l.strip(), f_genes.readlines()
                                                                   )))})
    # target_regions_per_class
    cov_data = {}
    type_mapping = {'TSS_corrected': 'corrected', 'TSS_original': 'original'}
    for csv in tss_coverage_csv_paths['relative']:   # iter over correction state
        # gene, signal_type = csv.split('_')
        csv_path = stem_path / csv
        csv_candidates = tuple(csv_path.glob('*.csv'))
        if len(csv_candidates) == 0:
            print(f"WARNING: no csv file found under path '{csv_path}'. Ignoring this directory ..")
            continue
        for csv_candidate in csv_candidates:   # iter over samples
            sample_name = csv_candidate.stem
            if 'X' in sample_name:
                continue
            signal_type = csv_candidate.parent.stem
            signal_type = type_mapping[signal_type]
            covs = pl.read_csv(csv_candidate)
            # get list of genes per class and filter data for it:
            for region_class, gene_list in genes_per_class.items():   # iter over region groups
                if cov_data.get(sample_name) is None:
                    cov_data.update({sample_name: {region_class: {signal_type: None}}})
                elif cov_data[sample_name].get(region_class) is None:
                    cov_data[sample_name].update({region_class: {signal_type: None}})
                elif cov_data[sample_name][region_class].get(signal_type) is None:
                    cov_data[sample_name][region_class].update({signal_type: None})
                avg_cov = np.array(covs.filter(covs['Gene'].is_in(gene_list))[:, 3:-1].mean(axis=0))[0]
                cov_data[sample_name][region_class][signal_type] = avg_cov

    # ACTUALLY CREATE PLOTS:
    # ------------------------------------------------------------------------------------------------------------------
    # create Multiplot:
    fig = make_subplots(rows=7, cols=1, shared_xaxes=False, shared_yaxes=False, vertical_spacing=0.03,
                        subplot_titles=('average reference sequence GC content',
                                        'B01 TSS coverage HK genes (-2.2% GC)',
                                        'H01 TSS coverage HK genes (+1.1% GC)',
                                        'P01 TSS coverage HK genes (+5.0% GC)',
                                        'B01 TSS coverage PAU genes (-2.2% GC)',
                                        'H01 TSS coverage PAU genes (+1.1% GC)',
                                        'P01 TSS coverage PAU genes (+5.0% GC)'),
                        start_cell='top-left',
                        specs=[[{"secondary_y": False, 'colspan': 1, 'rowspan': 1}],
                               [{"secondary_y": True, 'colspan': 1, 'rowspan': 1}],
                               [{"secondary_y": True, 'colspan': 1, 'rowspan': 1}],
                               [{"secondary_y": True, 'colspan': 1, 'rowspan': 1}],
                               [{"secondary_y": True, 'colspan': 1, 'rowspan': 1}],
                               [{"secondary_y": True, 'colspan': 1, 'rowspan': 1}],
                               [{"secondary_y": True, 'colspan': 1, 'rowspan': 1}]])
    # 0.) HK/PAU reference GC content
    ref_gc_content_traces = plot_ref_gc_content(data_to_plot=data_to_plot, transparencies=transparent_signals,
                                                signal_colors=TSS_COLORS, y_is_percentage=PERCENTAGE_PLOT,
                                                fig_fontsize=24, fig_width=1000, fig_height=1000)
    if ref_gc_content_traces == 1:
        raise ValueError("something went wrong during plotting of P01 LYL1")
    for ref_trace in ref_gc_content_traces.data:
        fig.add_trace(ref_trace, row=1, col=1)  # add all traces to first subplot
        # add smooth signals to DoC subplots
        if 'filtered' in ref_trace.name:
            if 'HK' in ref_trace.name:
                fig.add_trace(ref_trace, row=2, col=1, secondary_y=True)
                fig.add_trace(ref_trace, row=3, col=1, secondary_y=True)
                fig.add_trace(ref_trace, row=4, col=1, secondary_y=True)
            elif 'PAU' in ref_trace.name:
                fig.add_trace(ref_trace, row=5, col=1, secondary_y=True)
                fig.add_trace(ref_trace, row=6, col=1, secondary_y=True)
                fig.add_trace(ref_trace, row=7, col=1, secondary_y=True)
    # 1.) B01 HK+PAU genes
    b01_hkpau_traces = get_multi_signals_plot(plot_data_signals=cov_data['B01'], annotation='B01',
                                              signal_colors=SAMPLE_COLORS, x_axis_data_label='relative position / bp',
                                              y_axis_data_label='average relative coverage / 1',
                                              signal_dash={'corrected': 'solid', 'original': 'dot'},
                                              figure_width=1000, figure_height=380)
    if b01_hkpau_traces == 1:
        raise ValueError("something went wrong during plotting of B01 HK+PAU genes")
    for doc_trace in b01_hkpau_traces:
        if 'HK' in doc_trace.name:
            if 'corrected' in doc_trace.name:
                doc_trace.name = 'B01 HK corrected'
            elif 'original' in doc_trace.name:
                doc_trace.name = 'B01 HK original'
            fig.add_trace(doc_trace, row=2, col=1)
        elif 'PAU' in doc_trace.name:
            if 'corrected' in doc_trace.name:
                doc_trace.name = 'B01 PAU corrected'
            elif 'original' in doc_trace.name:
                doc_trace.name = 'B01 PAU original'
            fig.add_trace(doc_trace, row=5, col=1)
    # 2.) H01  HK+PAU genes
    h01_hkpau_traces = get_multi_signals_plot(plot_data_signals=cov_data['H01'], annotation='H01',
                                              signal_colors=SAMPLE_COLORS, x_axis_data_label='relative position / bp',
                                              y_axis_data_label='average relative coverage / 1',
                                              signal_dash={'corrected': 'solid', 'original': 'dot'},
                                              figure_width=1000, figure_height=380)
    if h01_hkpau_traces == 1:
        raise ValueError("something went wrong during plotting of H01 HK+PAU genes")
    for doc_trace in h01_hkpau_traces:
        if 'HK' in doc_trace.name:
            if 'corrected' in doc_trace.name:
                doc_trace.name = 'H01 HK corrected'
            elif 'original' in doc_trace.name:
                doc_trace.name = 'H01 HK original'
            fig.add_trace(doc_trace, row=3, col=1)
        elif 'PAU' in doc_trace.name:
            if 'corrected' in doc_trace.name:
                doc_trace.name = 'H01 PAU corrected'
            elif 'original' in doc_trace.name:
                doc_trace.name = 'H01 PAU original'
            fig.add_trace(doc_trace, row=6, col=1)
    # 3.) P01 HK+PAU genes
    p01_hkpau_traces = get_multi_signals_plot(plot_data_signals=cov_data['P01'], annotation='P01',
                                              signal_colors=SAMPLE_COLORS, x_axis_data_label='relative position / bp',
                                              y_axis_data_label='average relative coverage / 1',
                                              signal_dash={'corrected': 'solid', 'original': 'dot'},
                                              figure_width=1000, figure_height=380)
    if p01_hkpau_traces == 1:
        raise ValueError("something went wrong during plotting of P01 HK+PAU genes")
    for doc_trace in p01_hkpau_traces:
        if 'HK' in doc_trace.name:
            if 'corrected' in doc_trace.name:
                doc_trace.name = 'P01 HK corrected'
            elif 'original' in doc_trace.name:
                doc_trace.name = 'P01 HK original'
            fig.add_trace(doc_trace, row=4, col=1)
        elif 'PAU' in doc_trace.name:
            if 'corrected' in doc_trace.name:
                doc_trace.name = 'P01 PAU corrected'
            elif 'original' in doc_trace.name:
                doc_trace.name = 'P01 PAU original'
            fig.add_trace(doc_trace, row=7, col=1)
    # customize formatting and compute shared axis limits
    hk_limits = [None, None]
    pau_limits = [None, None]
    ref_limits = [None, None]
    plotting_order = range(len(fig.data))
    for dat_idx in plotting_order:
        trace_name = fig.data[dat_idx].name
        fig.data[dat_idx].line.dash = 'solid'
        if 'P01' in trace_name:
            fig.data[dat_idx].line.color = SAMPLE_COLORS['P01']
        elif 'H01' in trace_name:
            fig.data[dat_idx].line.color = SAMPLE_COLORS['H01']
        elif 'B01' in trace_name:
            fig.data[dat_idx].line.color = SAMPLE_COLORS['B01']
        if 'filtered' in trace_name:
            fig.data[dat_idx].line.width = 1.5
            fig.data[dat_idx].name = f"reference GC content {fig.data[dat_idx].name.split(',')[0]} genes"
            if dat_idx not in (2, 6):
                fig.data[dat_idx].line.color = \
                    f"{', '.join(TSS_COLORS['HK' if 'HK' in trace_name else 'PAU'].split(', ')[:3])}, " \
                    f"{0.5 if 'HK' in trace_name else 0.65})"
            if dat_idx not in (2, 6):
                fig.data[dat_idx].line.width = 2
                fig.data[dat_idx].showlegend = False
            else:  # show filtered reference genome GC content signals in legend & rename
                if ref_limits[0] is None or ref_limits[0] > fig.data[dat_idx].y.min():
                    ref_limits = [fig.data[dat_idx].y.min(), ref_limits[1]]
                if ref_limits[1] is None or ref_limits[1] < fig.data[dat_idx].y.max():
                    ref_limits = [ref_limits[0], fig.data[dat_idx].y.max()]
                fig.data[dat_idx].showlegend = True
        elif 'original' in trace_name:
            fig.data[dat_idx].line.width = 2
            r, g, b, _o = fig.data[dat_idx].line.color.strip('rgba(').strip(')').split(', ')
            fig.data[dat_idx].line.color = f'rgba({r}, {g}, {b}, 0.5)'  # make transparent
            # update limits
            if 'HK' in trace_name:
                fig.data[dat_idx].showlegend = True
                while 'HK' in fig.data[dat_idx].name:
                    fig.data[dat_idx].name = fig.data[dat_idx].name.replace('HK ', '')
                if hk_limits[0] is None or hk_limits[0] > fig.data[dat_idx].y.min():
                    hk_limits = [fig.data[dat_idx].y.min(), hk_limits[1]]
                if hk_limits[1] is None or hk_limits[1] < fig.data[dat_idx].y.max():
                    hk_limits = [hk_limits[0], fig.data[dat_idx].y.max()]
            elif 'PAU' in trace_name:
                fig.data[dat_idx].showlegend = False
                if pau_limits[0] is None or pau_limits[0] > fig.data[dat_idx].y.min():
                    pau_limits = [fig.data[dat_idx].y.min(), pau_limits[1]]
                if pau_limits[1] is None or pau_limits[1] < fig.data[dat_idx].y.max():
                    pau_limits = [pau_limits[0], fig.data[dat_idx].y.max()]
        elif 'corrected' in trace_name:
            fig.data[dat_idx].line.width = 3.5
            # update limits
            if 'HK' in trace_name:
                fig.data[dat_idx].showlegend = True
                while 'HK' in fig.data[dat_idx].name:
                    fig.data[dat_idx].name = fig.data[dat_idx].name.replace('HK ', '')
                if hk_limits[0] is None or hk_limits[0] > fig.data[dat_idx].y.min():
                    hk_limits = [fig.data[dat_idx].y.min(), hk_limits[1]]
                if hk_limits[1] is None or hk_limits[1] < fig.data[dat_idx].y.max():
                    hk_limits = [hk_limits[0], fig.data[dat_idx].y.max()]
            elif 'PAU' in trace_name:
                fig.data[dat_idx].showlegend = False
                if pau_limits[0] is None or pau_limits[0] > fig.data[dat_idx].y.min():
                    pau_limits = [fig.data[dat_idx].y.min(), pau_limits[1]]
                if pau_limits[1] is None or pau_limits[1] < fig.data[dat_idx].y.max():
                    pau_limits = [pau_limits[0], fig.data[dat_idx].y.max()]
        elif dat_idx in (0, 1):
            fig.data[dat_idx].showlegend = False
    # update overall figure formatting
    fig.update_layout(height=3400, width=1200, template='simple_white', showlegend=True,
                      font={'family': 'Ubuntu', 'size': 22, 'color': 'rgb(20, 20, 20)'},
                      title={'text': 'TSS GC Bias',
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
            if row == 7:
                fig.update_xaxes(title_text='relative position to TSS / bp', row=row, col=1)
            fig.update_yaxes(title_text='normalized average DoC / ratio', row=row, col=1)
            fig.update_yaxes(title_text='average reference GC content / %', row=row, col=1, secondary_y=True)
    # set same y-range for PAU and for HK gene DoC plots (separately)
    ref_limits = (ref_limits[0] - 1, ref_limits[1] + 1)
    hk_limits_range = hk_limits[1] - hk_limits[0]
    hk_limits = (hk_limits[0] - 0.01 * hk_limits_range, hk_limits[1] + 0.01 * hk_limits_range)
    pau_limits_range = pau_limits[1] - pau_limits[0]
    pau_limits = (pau_limits[0] - 0.01 * pau_limits_range, pau_limits[1] + 0.01 * pau_limits_range)
    for subplot_idx in range(2, 5):  # y-lims for HK
        fig.update_yaxes(range=hk_limits, row=subplot_idx, col=1, secondary_y=False)
    for subplot_idx in range(5, 8):  # y-lims for PAU
        fig.update_yaxes(range=pau_limits, row=subplot_idx, col=1, secondary_y=False)
    for subplot_idx in range(2, 8):  # y-lims for reference GC content
        fig.update_yaxes(range=ref_limits, row=subplot_idx, col=1, secondary_y=True)
    fig.show()
    output_file_path.parent.mkdir(parents=True, exist_ok=True)
    fig.write_image(output_file_path)
    return 0


if __name__ == '__main__':
    sys.exit(main())
