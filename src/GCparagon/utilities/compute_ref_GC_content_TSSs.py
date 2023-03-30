#!/usr/bin/env python3
import sys
import math
import numpy as np
from pathlib import Path
from collections import deque
from twobitreader import TwoBitFile
from typing import Union, Dict, Tuple
from plot_distributions import plot_ref_gc_content

CODE_ROOT_DIR = Path(__file__).parent.parent
CONTENT_ROOT_DIR = CODE_ROOT_DIR.parent.parent
PERCENTAGE_PLOT = True
two_bit_reference_file = CODE_ROOT_DIR / '2bit_reference/hg38.analysisSet.2bit'
hg38_genes_file = CONTENT_ROOT_DIR / 'accessory_files/TSSs/hg38_MANE_hg38.bed'
target_genes_per_class = {'HK': CONTENT_ROOT_DIR / 'accessory_files/TSSs/housekeeping_genes_hg38.txt',
                          'PAU': CONTENT_ROOT_DIR / 'accessory_files/TSSs/PAU_genes_hg38.txt'}


def load_gene_lists(gene_list_dict: Dict[str, Union[str, Path]]) -> Dict[str, Tuple[str]]:
    gene_data = {}.fromkeys(gene_list_dict)
    for list_name, list_path in gene_list_dict.items():
        with open(list_path, mode='rt') as f_gene_list:
            gene_data.update({list_name: tuple(filter(lambda e: e != '',
                                                      map(lambda l: l.strip(),
                                                          f_gene_list.readlines())))})
    return gene_data


def load_reference_genes(reference_genes_tsv: Union[str, Path]) -> Dict[str, Tuple[str, str, int, int]]:
    ref_gene_tuples = {}
    # chr1	33306765	33321098	ENSG00000184389	.	-
    with open(reference_genes_tsv, mode='rt') as f_ref_genes:
        _ = [ref_gene_tuples.update({ens_id: data})
             for ens_id, data in
             map(lambda c: (c[3], (c[0], c[5], int(c[1]), int(c[2]))),
                 filter(lambda c: c != [],
                        [gene_line.strip().split('\t') for gene_line in f_ref_genes.readlines()]))]
    return ref_gene_tuples


def get_gc_content(target_gene_ids: Dict[str, Tuple[str]], reference_gene_dict: Dict[str, Tuple[str, str, int, int]],
                   window_size=10001, as_percentage=True) -> Dict[str, np.array]:
    gc_arrays = {}.fromkeys(target_gene_ids)
    actual_window_size = int(window_size//2) * 2 + 1
    if actual_window_size != window_size:
        print(f"-> WARNING :  window size specified was {window_size:,} bp but will use {actual_window_size:,} bp instead.")
    for target_id, gene_ids in target_gene_ids.items():
        print(f" i :  processing gene group '{target_id}'\n------------------------------------")
        chromosome_wise_loci = {}
        for target_gene in gene_ids:
            try:
                chrom, strand, cds_leftmost, cds_rightmost = reference_gene_dict[target_gene]
            except KeyError:  # gene not in ensembl list
                print(f"gene '{target_gene}' not found in MANE transcripts. Continuing ..")
                continue
            try:
                chromosome_wise_loci[chrom].append((cds_leftmost, cds_rightmost, strand))
            except KeyError:
                chromosome_wise_loci.update({chrom: [(cds_leftmost, cds_rightmost, strand)]})
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
            for cds_leftmost, cds_rightmost, strand in loci_list:
                # get uppercase reference sequence
                if strand == '+':
                    tss_seq = chrom_handle[cds_leftmost-int(window_size//2):
                                           cds_leftmost+int(window_size//2)+1].upper()
                elif strand == '-':
                    # reverse slice:
                    # if slice_me = 'abcdefghijklmnop'; then
                    # slice_me[4:10] == 'efghij' AND
                    # slice_me[4:10][::-1] == 'jihgfe'
                    tss_seq = chrom_handle[cds_rightmost - int(window_size//2):
                                           cds_rightmost + int(window_size//2)+1].upper()[::-1]
                else:
                    raise ValueError(f'unknown strand {strand} for entry '
                                     f'({chrom}, {cds_leftmost, }, {cds_rightmost}, {strand})')
                # update numpy array
                gc_window += tuple(map(lambda c: int(c in 'GC'), tss_seq))  # adds either 0 or 1 for every position
                total_loci += 1
        # compute GC content per relative position:
        gc_window /= total_loci
        if as_percentage:
            gc_window *= 100.
        gc_arrays[target_id] = gc_window
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
    return y[math.floor(window_len/2-1):-math.ceil(window_len/2)]


def main() -> int:
    # read target genes
    target_gene_ensembl_ids = load_gene_lists(gene_list_dict=target_genes_per_class)
    # load all_reference_genes
    hg38_ref_genes = load_reference_genes(reference_genes_tsv=hg38_genes_file)
    # accumulate GC content
    for view_window in (2001, ):
        gc_per_group = get_gc_content(target_gene_ids=target_gene_ensembl_ids, reference_gene_dict=hg38_ref_genes,
                                      window_size=view_window, as_percentage=PERCENTAGE_PLOT)
        # create hamming window filtered versions
        for hamming_window_length in (15, ):
            transparent_signals = {'original': 0.25}
            data_to_plot = {'original': gc_per_group}
            data_to_plot.update({f'{hamming_window_length}bp_hamming_filtered': None})
            transparent_signals.update({f'{hamming_window_length}bp_hamming_filtered': 1})
            filtered_data = {}
            for group, np_arr in gc_per_group.items():
                filtered_data.update({group: smooth_1d_data(x=np_arr, window='hamming',
                                                            window_len=hamming_window_length)})
            data_to_plot[f'{hamming_window_length}bp_hamming_filtered'] = filtered_data
            # create GC plot plotly express scatter with lines
            figure_output_path = CONTENT_ROOT_DIR / 'accessory_files/TSSs/' \
                                                    f'gene_groups_ref_gc_content_{view_window}bp_' \
                                                    f'{hamming_window_length}bpHammingSmoothed.png'
            plot_ref_gc_content(data_to_plot=data_to_plot, transparencies=transparent_signals, fig_fontsize=24,
                                signal_colors={'HK': (31, 119, 180), 'PAU': (255, 127, 14)}, fig_height=1000,
                                output_file_path=figure_output_path, fig_width=1000, y_is_percentage=PERCENTAGE_PLOT,
                                figure_title='TSSs Group Average Reference GC Content')
    return 0


if __name__ == "__main__":
    sys.exit(main())
