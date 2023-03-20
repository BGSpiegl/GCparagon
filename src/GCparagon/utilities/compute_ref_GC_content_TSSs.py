#!/usr/bin/env python3
import sys
import math
import numpy as np
from pathlib import Path
from collections import deque
from twobitreader import TwoBitFile
from typing import Union, Dict, Tuple
from plot_distributions import plot_ref_gc_content

CONTENT_ROOT_DIR = Path(__file__).parent.parent
PERCENTAGE_PLOT = True
two_bit_reference_file = CONTENT_ROOT_DIR / '2bit_reference/hg38.analysisSet.2bit'
hg38_genes_file = CONTENT_ROOT_DIR / 'accessory_files/gene_lists/hg38_MANE_hg38.bed'
target_genes_per_class = {'HK': CONTENT_ROOT_DIR / 'accessory_files/gene_lists/housekeeping_genes_hg38.txt',
                          'PAU': CONTENT_ROOT_DIR / 'accessory_files/gene_lists/PAU_genes_hg38.txt'}


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
    for view_window in (4001, 2001):  # 10001, 6001, 3001,
        gc_per_group = get_gc_content(target_gene_ids=target_gene_ensembl_ids, reference_gene_dict=hg38_ref_genes,
                                      window_size=view_window)
        # create hamming window filtered versions
        for hamming_window_length in (15, 17):  # 11, 15,
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
            figure_output_path = CONTENT_ROOT_DIR / 'accessory_files/gene_lists/' \
                                                    f'gene_groups_ref_gc_content_{view_window}bp_' \
                                                    f'{hamming_window_length}bpHammingSmoothed.png'
            plot_ref_gc_content(data_to_plot=data_to_plot, transparencies=transparent_signals, fig_fontsize=24,
                                signal_colors={'HK': (31, 119, 180), 'PAU': (255, 127, 14)}, fig_height=1000,
                                output_file_path=figure_output_path, fig_width=1000, y_is_percentage=PERCENTAGE_PLOT)
    return 0


if __name__ == "__main__":
    sys.exit(main())


# 74 GENES NOT FOUND for housekeeping genes (1267 remaining):
#  i :  processing gene group 'HK'
# ------------------------------------
# gene 'ENSG00000275837' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000276701' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000278494' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000131100' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000142168' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000206208' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000236490' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000093010' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000099901' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000128185' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000206281' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000112493' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000099940' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000099942' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000185651' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000100030' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000099956' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000240972' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000100028' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000133422' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000285302' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000110651' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000145494' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000126934' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000100348' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000100353' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000154767' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000187866' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000100142' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000100201' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000100216' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000128272' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000187051' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000100359' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000100387' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000100401' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000100410' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000198911' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000100243' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000242247' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000100266' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000100347' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000127445' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000254505' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000100422' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000100258' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000085721' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000120948' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000231925' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000179218' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000158290' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000146223' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000119392' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000156603' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000156467' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000119655' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000169032' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000135956' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000189241' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000131495' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000173933' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000042753' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000182768' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000108433' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000027697' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000141378' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000265241' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000179295' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000108946' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000115524' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000143420' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000159352' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000116030' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000153827' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000058673' not found in MANE transcripts. Continuing ..
#
# 370 GENES NOT FOUND for group protein atlas unexpressed (1421 remaining):
# i :  processing gene group 'PAU'
# -------------------------------------------------------------------
# gene 'ENSG00000277339' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000277630' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000276011' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000277725' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000267728' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000262594' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000276256' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000277400' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000273748' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000278384' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000274412' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000276816' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000276420' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000272712' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000275499' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000281296' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000275960' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000273496' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000276033' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000276848' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000277390' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000277316' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000278760' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000275918' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000224336' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000183704' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000226187' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000279081' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000235769' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000236210' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000279783' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000282175' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000282424' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000224922' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000237482' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000224708' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000206401' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000130538' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000283150' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000283767' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000283949' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000229728' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000203618' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000284588' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000267793' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000283834' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000288542' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000288389' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000226872' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000261456' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000205740' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000169548' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000220891' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000285084' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000230212' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000285029' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000230479' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000288159' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000288431' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000125879' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000176515' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000258555' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000196431' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000132671' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000273838' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000173335' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000265818' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000237378' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000225344' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000236637' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000253896' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000177596' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000162006' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000266990' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000261720' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000261071' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000230801' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000272968' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000230944' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000255310' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000186190' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000177710' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000224876' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000268798' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000279083' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000279029' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000182545' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000235837' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000263970' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000173976' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000268818' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000128313' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000183562' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000269125' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000225655' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000273514' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000260476' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000259780' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000268670' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000259214' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000232406' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000260739' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000230952' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000228327' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000230699' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000224969' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000171987' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000128310' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000225337' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000180913' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000180909' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000122548' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000229116' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000267571' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000212302' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000158077' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000170846' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000284633' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000101074' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000228919' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000249695' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000250913' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000260642' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000271806' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000223509' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000234396' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000262456' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000103426' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000240131' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000232135' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000286401' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000234965' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000214248' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000255309' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000255552' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000100341' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000270617' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000243207' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000175746' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000244255' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000271010' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000261216' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000272836' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000262089' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000260814' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000136275' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000262211' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000280401' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000225705' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000269210' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000124091' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000205086' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000233727' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000276945' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000248871' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000116726' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000179172' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000229990' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000270601' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000204480' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000257529' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000223889' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000230386' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000180269' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000251105' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000179342' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000288611' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000258377' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000262560' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000103310' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000176659' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000282885' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000278002' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000179002' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000255521' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000267633' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000269970' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000231102' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000121314' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000206549' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000125522' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000280393' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000258592' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000227782' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000283667' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000258682' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000229717' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000206592' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000275186' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000255921' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000205176' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000275928' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000276540' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000278743' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000214643' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000008197' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000204049' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000228510' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000197665' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000236022' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000276380' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000279847' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000177359' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000253175' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000265041' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000258517' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000228634' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000273791' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000273769' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000198283' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000284873' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000225850' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000231170' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000222604' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000174417' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000287625' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000227416' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000258234' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000227958' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000228768' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000201384' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000277511' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000228335' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000197734' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000259827' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000262075' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000197322' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000213590' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000213605' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000220130' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000280055' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000185933' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000278126' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000278342' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000260122' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000112238' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000257829' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000161850' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000257500' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000260146' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000250733' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000226490' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000262691' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000272861' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000235319' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000173612' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000187812' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000231505' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000198083' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000171360' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000204961' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000278901' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000268975' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000171815' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000259474' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000223720' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000274884' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000258663' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000179038' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000259332' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000268643' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000184224' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000234997' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000214285' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000270394' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000269177' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000283653' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000274606' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000226699' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000286789' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000280409' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000109471' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000258001' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000198711' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000248903' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000267496' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000267211' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000204850' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000259992' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000279512' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000259009' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000279980' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000238244' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000248458' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000253630' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000227197' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000183643' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000226180' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000270655' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000258971' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000279551' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000280291' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000242295' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000251209' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000222031' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000205018' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000253445' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000259112' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000213121' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000268170' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000288547' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000279742' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000258044' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000257127' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000233184' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000213065' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000243415' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000281652' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000281058' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000280799' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000258834' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000281627' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000280457' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000260596' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000222033' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000104970' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000278998' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000204960' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000229365' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000254702' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000240652' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000223872' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000280414' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000214919' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000255292' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000268472' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000269794' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000242107' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000232098' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000230381' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000174325' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000279862' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000273240' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000231384' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000252119' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000236334' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000264491' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000255537' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000222017' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000203804' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000255535' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000197915' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000197364' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000279046' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000234238' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000287064' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000267624' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000200063' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000267432' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000261924' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000189030' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000238934' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000163081' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000233538' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000213057' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000227141' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000222007' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000222022' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000258984' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000225057' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000286307' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000196758' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000220256' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000219159' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000162888' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000234863' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000229930' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000177800' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000270106' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000270188' not found in MANE transcripts. Continuing ..
# gene 'ENSG00000229960' not found in MANE transcripts. Continuing ..
