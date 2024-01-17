#!/usr/bin/env python3

import re
import numpy as np
import pandas as pd
import plotly.express as px
from typing import Union as OneOf, Dict, Tuple
from pathlib import Path
# from correct_GC_bias import create_region_label

# DEFINITIONS
# intervals (opaque -> 0.01; we have ~1,700 of these for hg38; use already available functions!)

SOURCE_CODE_PATH = Path(__file__).parent.parent.parent

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# ADAPT THE FOLLOWING PARAMETERS ACCORDING TO YOUR NEEDS! MUST BE IDENTICAL TO THE ONES FROM THE PREVIOUS SCRIPTS !!!!!!
GENOME_BUILD = 'hg19'
INTERVAL_SIZE = 10**6  # 1Mbp -> change according to input BED files!
SHIFT_N_TIMES = 16  # you might want to select a higher number
REGION_OVERLAP_PERCENTAGE_THRESHOLD = 33  # percent max. bad region overlap!
genomic_interval_fgcd = (SOURCE_CODE_PATH /
                         f'accessory_files/'
                         f'{GENOME_BUILD}_minimalExclusionListOverlap_{INTERVAL_SIZE//10**6}Mbp_intervals_'
                         f'{REGION_OVERLAP_PERCENTAGE_THRESHOLD}pcOverlapLimited.FGCD.bed')
#                           ^----- compute with 03-GI-preselection_compute_genomic_interval_fragment_GC_content_hg19.py
reference_fgcd = SOURCE_CODE_PATH / f'accessory_files/{GENOME_BUILD}_reference_GC_content_distribution.tsv'
#                           ^----- compute with 04-GI-preselection_simulate_genomewide_reference_FGCD_hg19.py

output_directory = (SOURCE_CODE_PATH /
                    f'accessory_files/'
                    f'genomic_interval_preselection-shifted{SHIFT_N_TIMES}x_{GENOME_BUILD}/'
                    f'{GENOME_BUILD}_FGCD_intervals_vs_reference')
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


# from correct_GC_bias.py
def create_region_label(chrm: str, start: int, end: int):
    return f"{chrm}_{start:,}-{end:,}"


def read_ref_gc_dist(dist_table_path: OneOf[str, Path]) -> np.array:
    with open(dist_table_path, 'rt') as f_ref:
        _hdr = f_ref.readline()
        counts = np.array([float(line.strip().split('\t')[1]) for line in f_ref.readlines()], dtype=float)
    return counts / counts.sum()


def read_interval_fgcd_table(table_path: OneOf[str, Path]) -> Dict[str, Tuple[int, np.array]]:
    gc_content_dists = {}
    with open(table_path, 'rt') as f_fgcd:
        _hdr = f_fgcd.readline()  # "chromosome	start	end	exclusion_marked_bases	0	1	2	..." etc.
        for line in f_fgcd:
            chrom, start, end, bad_score, *gc_fragment_counts = line.strip().split('\t')
            interval_start = int(start)
            interval_stop = int(end)
            fragment_counts = np.array(list(map(lambda s: int(s), gc_fragment_counts)), dtype=float)
            gc_content_dists[create_region_label(chrm=chrom, start=interval_start, end=interval_stop)] = (
                int(bad_score), fragment_counts / fragment_counts.sum())
    return gc_content_dists


def plot_interval_gc_dists_vs_ref(interval_gc_data: Dict[str, Tuple[int, np.array]], reference_gc_dist: np.array,
                                  out_dir_path: Path, fig_width=1500, fig_height=1000, fig_fontsize=24,
                                  spline_interpolation=True, normalize_to_dataset_size=True, annotation=None,
                                  reference_normalized=True, show_figure=False, image_formats=('png',),
                                  reduced_bins=True):
    intervals = sorted(list(interval_gc_data.keys()))
    # create total count sum for all samples
    interval_gc_content_total_fragments = {}.fromkeys(intervals)
    # normalize datasets if requested
    if normalize_to_dataset_size:  # divide each dataset by number of total fragments and multiply with 100
        for ntvl in intervals:
            fgcd = interval_gc_data[ntvl][1]
            interval_gc_content_total_fragments[ntvl] = fgcd.sum()
            if interval_gc_content_total_fragments[ntvl] != 1.0:
                fgcd /= interval_gc_content_total_fragments[ntvl]
            interval_gc_data[ntvl] = (interval_gc_data[ntvl][0], fgcd * 100.)
        if reference_gc_dist is not None and not reference_normalized:  # normalize if is not normalized!
            reference_gc_dist = reference_gc_dist / reference_gc_dist.sum() * 100.
    plot_data_list = []
    columns = ['interval', 'non-exclusion marked bases fraction', 'fragment GC content / %',
               'relative frequency / %']
    # create GC content range integers (either in steps of 1 or 2)
    gc_values = range(0, len(interval_gc_data[intervals[0]][1]), 2) \
        if reduced_bins else range(0, len(interval_gc_data[intervals[0]][1]), 1)  # ASSERTS all bins present!
    # add intervals data
    for ntvl in intervals:
        if reduced_bins:  # left included; right excluded
            plot_data_list.extend([[ntvl, interval_gc_data[ntvl][0], gc,
                                    (interval_gc_data[ntvl][1][gc] +
                                     (interval_gc_data[ntvl][1][gc+1] if gc != 100 else 0))]
                                   for gc in gc_values])
        else:
            plot_data_list.extend([[ntvl, interval_gc_data[ntvl][0], gc, interval_gc_data[ntvl][1][gc]]
                                   for gc in gc_values])
    # add reference distribution
    if reference_gc_dist is not None:
        if reduced_bins:  # left included; right excluded
            plot_data_list.extend([['hg38 expected GC', '1.', gc,
                                    reference_gc_dist[gc] + (reference_gc_dist[gc + 1] if gc != 100 else 0)]
                                   for gc in gc_values])
        else:
            plot_data_list.extend([['hg38 expected GC', 'original', gc, reference_gc_dist[gc]] for gc in gc_values])
    # create DataFrame for plotting
    plot_data = pd.DataFrame(plot_data_list, columns=columns)
    gc_dist_fig = px.line(plot_data, x='fragment GC content / %', y='relative frequency / %',
                          line_shape='spline' if spline_interpolation else 'linear',
                          color='interval',
                          template="simple_white", width=fig_width, height=fig_height,
                          title=f"cfDNA GC Content - 1Mbp intervals vs. genome" +
                                (' (' + annotation + ')' if annotation is not None else ''))
    for dat_idx in range(len(gc_dist_fig.data)):
        if reference_gc_dist is not None and \
                gc_dist_fig.data[dat_idx].name == 'hg38 expected GC':  # expected distribution
            gc_dist_fig.data[dat_idx].line.color = 'red'
            gc_dist_fig.data[dat_idx].line.width = 4
        else:
            gc_dist_fig.data[dat_idx].line.width = 2
            interval_name = gc_dist_fig.data[dat_idx].name
            start_str, end_str = interval_name.split('_')[1].split('-')
            int_start, int_end = int(re.sub(',', '', start_str)), int(re.sub(',', '', end_str))
            non_exclusion_fraction = 1. - (interval_gc_data[ntvl][0] / (int_end - int_start))
            # compute transparency -> will have hundreds of lines!
            opacity = min(1.0, 0.05 * (non_exclusion_fraction - 0.6666666) * 3)
            # non_exclusion_fraction in [1, 0.6666...]s -> scale to a range between 1.0 and 0.0.
            color = f'rgba(10, 10, 10, {round(opacity, ndigits=4)})'
            gc_dist_fig.data[dat_idx].line.color = color
            gc_dist_fig.data[dat_idx].name = re.sub('_', ':', gc_dist_fig.data[dat_idx].name)
    # set legend options
    gc_dist_fig.update_layout(font_family="Ubuntu", font_size=fig_fontsize, showlegend=False,
                              title_font_size=fig_fontsize + 2,
                              title={'text': f"cfDNA GC Content - {len(interval_gc_data):,} 1Mbp intervals vs. genome" +
                                             (' (' + annotation + ')' if annotation is not None else ''),
                                     'font': {'family': 'Ubuntu', 'size': 28, 'color': 'rgb(20, 20, 20)'},
                                     'xanchor': 'center', 'yanchor': 'middle', 'x': 0.5})
    if show_figure:
        gc_dist_fig.show()
    out_dir_path.mkdir(parents=True, exist_ok=True)
    for image_format in image_formats:
        out_file = out_dir_path / (f"GC_content_reference_genome_vs_preselected_intervals"
                                   f"{re.sub(', ', '_', annotation) if annotation else ''}."
                                   f"{len(gc_values)}bins.{image_format}")
        try:
            gc_dist_fig.write_image(out_file)
        except:
            print(f"WARNING - could not save figure '{out_file.stem}' in '{image_format}' image format! Continuing ..")


if __name__ == '__main__':
    interval_gc_distributions = read_interval_fgcd_table(table_path=genomic_interval_fgcd)
    ref_fgcd = read_ref_gc_dist(dist_table_path=reference_fgcd)
    # create expensive plot ... (hundreds of lines expected!)
    plot_interval_gc_dists_vs_ref(interval_gc_data=interval_gc_distributions, reference_gc_dist=ref_fgcd,
                                  out_dir_path=output_directory, fig_width=1500, fig_height=1000, fig_fontsize=24,
                                  spline_interpolation=False,  # does not work with so many lines apparently:
                                  # ValueError:
                                  #     Invalid value of type 'builtins.str' received for the 'shape' property of
                                  #     scattergl.line
                                  #         Received value: 'spline'
                                  #
                                  #     The 'shape' property is an enumeration that may be specified as:
                                  #       - One of the following enumeration values:
                                  #             ['linear', 'hv', 'vh', 'hvh', 'vhv']
                                  # (NOTE: plotly.graph_objs.scattergl.Line
                                  normalize_to_dataset_size=True, annotation=None,
                                  reference_normalized=False, show_figure=False, image_formats=('png', 'svg', 'pdf'),
                                  reduced_bins=True)  # smoothes a bit
