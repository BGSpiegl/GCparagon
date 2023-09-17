#!/usr/bin/env python3

import re
import sys
import numpy as np
import pandas as pd
from pathlib import Path
import plotly.express as px
from collections import defaultdict
from typing import Optional, Union, Dict, Tuple

SCRIPT_ROOT_PATH = Path(__file__).parent
SCRIPT_ROOT_DIR = str(SCRIPT_ROOT_PATH)

output_dir = SCRIPT_ROOT_PATH / 'genome_wide_correction_fidelity_plots'

IMAGE_FORMATS = ('png', 'svg', 'pdf')

all_samples_dists_file = SCRIPT_ROOT_PATH / \
                         'COMBINED-CORRECTION-STATISTICS_allSamples_GC-percentage_perFragmentSequence_28-08-2023.tsv'

reference_gc_path = SCRIPT_ROOT_PATH / 'simulated_hg38_ref_gc_fromEntireBAM'
reference_gc_dists_samples = {'B01': reference_gc_path /
                                     'hg38.analysisSet.2bit_B01-entireBAM-fragLenDist_ref_GC_dist_500Mreads.tsv',
                              'H01': reference_gc_path /
                                     'hg38.analysisSet.2bit_H01-entireBAM-fragLenDist_ref_GC_dist_500Mreads.tsv',
                              'C01': reference_gc_path /
                                     'hg38.analysisSet.2bit_C01-entireBAM-fragLenDist_ref_GC_dist_500Mreads.tsv',
                              'P01': reference_gc_path /
                                     'hg38.analysisSet.2bit_P01-entireBAM-fragLenDist_ref_GC_dist_500Mreads.tsv'}

SPLINE_INTERPOLATION = True


def read_stats(stats_path: Union[str, Path], return_as_percentage=False) \
        -> Optional[Tuple[Dict[str, Dict[str, np.array]]]]:
    with open(stats_path, 'rt') as f_stats:
        stats_header = f_stats.readline()
        stats_data = f_stats.readlines()
    header_content = stats_header.strip().split('\t')
    leading_class_columns = header_content[:[ct.isdigit() for ct in header_content].index(True)]
    try:
        status_column = leading_class_columns.index('status')
    except ValueError:
        try:
            status_column = leading_class_columns.index('Status')
        except ValueError:
            print(f"'preset' column not found in header. Exiting ...")
            return None
    try:
        algorithm_column = leading_class_columns.index('algorithm')
    except ValueError:
        try:
            algorithm_column = leading_class_columns.index('Algorithm')
        except ValueError:
            print(f"'preset' column not found in header. Exiting ...")
            return None
    try:
        sample_column = leading_class_columns.index('sample')
    except ValueError:
        try:
            sample_column = leading_class_columns.index('Sample')
        except ValueError:
            print(f"'sample' column not found in header. Exiting ...")
            return None
    # gather corrected and uncorrected of specific presets for all 4 samples
    stats_data_content = list(map(lambda l: l.strip().split('\t'), stats_data))
    original_data = {}
    corrected_data = {}
    return_data = {'original': original_data, 'corrected': corrected_data}
    for line_content in stats_data_content:
        if len(line_content) <= 1:
            continue
        cur_stt = line_content[status_column]
        cur_smpl = line_content[sample_column]
        # change names
        if cur_smpl not in ('B01', 'C01', 'H01', 'P01'):
            match cur_smpl:
                case 'B58_3':
                    cur_smpl = 'B01'
                case 'P244_7':
                    cur_smpl = 'P01'
                case 'C219_5':
                    cur_smpl = 'C01'
                case 'NPH_011':
                    cur_smpl = 'H01'
        cur_algo = line_content[algorithm_column]  # for original counts (= uncorrected), the algorithm is irrelevant
        cur_data = np.array(line_content[len(leading_class_columns):], dtype=float)
        # must be float - will be divided later! avoid error:
        # "numpy.core._exceptions._UFuncOutputCastingError: Cannot cast ufunc 'divide' output from dtype('float64') to
        #  dtype('int64') with casting rule 'same_kind'"
        if return_as_percentage:
            data_sum = cur_data.sum()
            scaled_array = cur_data * 100. / data_sum
            try:
                return_data[cur_stt][cur_algo].update({f'{cur_smpl}': scaled_array})
            except KeyError:
                return_data[cur_stt].update({cur_algo: {f'{cur_smpl}': scaled_array}})
        else:  # return native as read in
            try:
                return_data[cur_stt][cur_algo].update({f'{cur_smpl}': cur_data})
            except KeyError:
                return_data[cur_stt].update({cur_algo: {f'{cur_smpl}': cur_data}})
    if return_data == {'original': {}, 'corrected': {}}:
        return None
    return original_data, corrected_data


def read_tsv_data_to_dict(table_path: Path, normalize_to_1=False, normalize_to_100=False) -> defaultdict[int, Union[float, int]]:
    """

    :param table_path:
    :param normalize_to_1:
    :param normalize_to_100: takes precedence over 'normalize_to_1'
    :return:
    """
    with open(table_path, 'rt') as f_tab:
        _header = f_tab.readline()
        file_content = [line.strip().split('\t') for line in f_tab.readlines()]
    content_dict = defaultdict(float)
    for gc_perc, frequency_value in file_content:
        content_dict[int(gc_perc)] = float(frequency_value)
    if normalize_to_1 or normalize_to_100:
        all_values_sum = sum(content_dict.values())
        normalized_content_dict = defaultdict(float)
        for gc_perc, abs_freq in content_dict.items():
            normalized_content_dict[gc_perc] = abs_freq / all_values_sum
        content_dict = normalized_content_dict
    if normalize_to_100:
        normalized_content_dict = defaultdict(float)
        for gc_perc, rel_freq in content_dict.items():
            normalized_content_dict[gc_perc] = rel_freq * 100.
        content_dict = normalized_content_dict
    return content_dict


def plot_fragment_gc_dists_griffin(original_gc_data: Dict[str, np.array], corrected_gc_data: Dict[str, np.array],
                                   out_dir_path: Path,
                                   reference_dists: Dict[str, Union[defaultdict[int, float], defaultdict[float], None]],
                                   sample_name: Optional = None,
                                   markers=True, fig_width=1200, fig_height=800, fig_fontsize=30, show_figure=False,
                                   normalize_to_dataset_size=True, annotation=None, reduced_bins=True,
                                   spline_interpolation=True, reference_normalized=True, image_formats=('png',)):
    orig_samples = sorted(list(original_gc_data.keys()))  # no algorithms in name
    corrected_samples = sorted(list(corrected_gc_data.keys()))  # algorithms in name
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
                    reference_dists[sample][gc] *= 100.
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
        for algorithm in ('Griffin', 'GCparagon'):
            for sample in corrected_gc_data:
                if algorithm not in sample:
                    continue  # print Griffin first and then GCparagon on top
                sample_no_alg = sample.split()[0]
                if reduced_bins:  # left included; right excluded
                    plot_data_list.extend([
                        [(f'{sample_no_alg} original' if 'original' == status else f'{sample} corrected'),
                         status, gc, ((original_gc_data[sample_no_alg][gc] +
                                       (original_gc_data[sample_no_alg][gc+1] if gc != 100 else 0))
                                      if sample_no_alg == orig_samples[0] and algorithm == 'Griffin' else None)
                         # write original series only once!
                         if 'original' == status else
                         (corrected_gc_data[sample][gc] + (corrected_gc_data[sample][gc+1] if gc != 100 else 0))]
                        for gc in gc_values])
                    plot_data_list = list(filter(lambda e: e is not None, plot_data_list))
                else:
                    plot_data_list.extend([[(f'{sample} original' if 'original' == status else f'{sample} corrected'),
                                            status, gc, original_gc_data[sample][gc]
                                            if 'original' == status else corrected_gc_data[sample][gc]]
                                           for gc in gc_values])
    plot_data = pd.DataFrame(plot_data_list, columns=columns)
    opaque_blank = 1.
    opaque_griffin = 1.
    opaque_gcparagon = .93
    opaque_original = 1.
    opaque_simulated = 1.
    color_str_b01 = '95, 145, 255'
    color_str_c01 = '255, 125, 0,'
    color_str_h01 = '47, 207, 33'
    color_str_p01 = '255, 0, 0,'
    griffin_str = '180, 180, 180'
    gc_content_figure = px.line(plot_data, x='GC', y='relative frequency percent', line_dash='status',
                                color='sample', line_dash_sequence=['longdash', 'solid'],
                                line_shape='spline' if spline_interpolation else 'linear',
                                template="simple_white", width=fig_width, height=fig_height,
                                color_discrete_map={'B01 (Griffin) corrected': f'rgba({griffin_str} {opaque_griffin})',
                                                    'B01 (GCparagon) corrected': f'rgba({color_str_b01} {opaque_gcparagon})',
                                                    'C01 (Griffin) corrected': f'rgba({griffin_str} {opaque_griffin})',
                                                    'C01 (GCparagon) corrected': f'rgba({color_str_c01} {opaque_gcparagon})',
                                                    'H01 (Griffin) corrected': f'rgba({griffin_str} {opaque_griffin})',
                                                    'H01 (GCparagon) corrected': f'rgba({color_str_h01} {opaque_gcparagon})',
                                                    'P01 (Griffin) corrected': f'rgba({griffin_str} {opaque_griffin})',
                                                    'P01 (GCparagon) corrected': f'rgba({color_str_p01} {opaque_gcparagon})',
                                                    'B01 original': f'rgba({color_str_b01} {opaque_original})',
                                                    'C01 original': f'rgba({color_str_c01} {opaque_original})',
                                                    'H01 original': f'rgba({color_str_h01} {opaque_original})',
                                                    'P01 original': f'rgba({color_str_p01} {opaque_original})',
                                                    'B01 simulated': f'rgba(0, 0, 0, {opaque_simulated})',
                                                    'H01 simulated': f'rgba(0, 0, 0, {opaque_simulated})',
                                                    'C01 simulated': f'rgba(0, 0, 0, {opaque_simulated})',
                                                    'P01 simulated': f'rgba(0, 0, 0, {opaque_simulated})',
                                                    'B01': f'rgba({color_str_b01} {opaque_blank})',
                                                    'C01': f'rgba({color_str_c01} {opaque_blank})',
                                                    'H01': f'rgba({color_str_h01} {opaque_blank})',
                                                    'P01': f'rgba({color_str_p01} {opaque_blank})'})

    ref_iter = 0
    ref_markers = (134, 133, 101, 102)  # ('x-thin-open','cross-thin-open',  'square-open', 'diamond-open')
    ref_dashes = ('solid', 'dot', 'dash', 'dashdot')
    for dat_idx in range(len(gc_content_figure.data)):
        trace_name = ', '.join(gc_content_figure.data[dat_idx].name.split(', ')[:-1])
        cur_status = gc_content_figure.data[dat_idx].name.split(', ')[-1]
        avg_gc = (plot_data[(plot_data['sample'] == trace_name) * (plot_data['status'] == cur_status)]['GC'] *
                  plot_data[(plot_data['sample'] == trace_name) *
                            (plot_data['status'] == cur_status)]['relative frequency percent']).sum() / 100.
        if reference_dists is not None and \
                'simulated' in gc_content_figure.data[dat_idx].name:  # ADJUST ATTRIBUTES OF SIMULATED LINES
            if markers:  # not used
                gc_content_figure.data[dat_idx].line.dash = 'solid'
                gc_content_figure.data[dat_idx].mode = 'lines+markers'
                gc_content_figure.data[dat_idx].marker.symbol = ref_markers[ref_iter]
                gc_content_figure.data[dat_idx].marker.size = 10
                gc_content_figure.data[dat_idx].marker.line.width = 1.5
            else:
                gc_content_figure.data[dat_idx].line.dash = ref_dashes[ref_iter]
            # The 'dash' property is an enumeration that may be specified as:
            #       - One of the following dash styles:
            #             ['solid', 'dot', 'dash', 'longdash', 'dashdot', 'longdashdot']
            #       - A string containing a dash length list in pixels or percentages
            #             (e.g. '5px 10px 2px 2px', '5, 10, 2, 2', '10% 20% 40%', etc.)
            # default definitions: (reference line)
            gc_content_figure.data[dat_idx].line.width = 20  # 20
            gc_content_figure.data[dat_idx].name = f"{trace_name} ({avg_gc:.1f}% GC)"
            ref_iter += 1
            # ADJUST ATTRIBUTES OF GCPARAGON-CORRECTED LINES
        elif gc_content_figure.data[dat_idx].line.dash == 'solid':
            gc_content_figure.data[dat_idx].line.width = 8  # 7
            gc_content_figure.data[dat_idx].name = f"{trace_name} ({avg_gc:.1f}% GC)"
            # gc_content_figure.data[dat_idx].line.color = ','.join(gc_content_figure.data[dat_idx].
            #                                                             line.color.split(',')[:-1]) + ', 0.6)'
        else:  # ADJUST ATTRIBUTES OF ORIGINAL, UNCORRECTED LINES
            gc_content_figure.data[dat_idx].line.width = 6  # 6
            gc_content_figure.data[dat_idx].line.dash = 'dot'
            gc_content_figure.data[dat_idx].name = f"{trace_name} ({avg_gc:.1f}% GC)"
        if 'Griffin' in gc_content_figure.data[dat_idx].name:  # ADJUST ATTRIBUTES OF GRIFFIN-CORRECTED LINES
            gc_content_figure.data[dat_idx].line.dash = '1.1%'
            gc_content_figure.data[dat_idx].line.width = 8  # 6
    # could shift ordering of elements in legened by using "legendrank=5" in trace definition
    gc_content_figure.update_layout(font_family="Ubuntu", font_size=fig_fontsize-2,
                                    xaxis_title='fragment GC content / %', yaxis_title='relative frequency / %',
                                    legend={'orientation': 'h', 'xanchor': 'center', 'yanchor': 'top',
                                            'x': 0.5, 'y': -0.2, 'title': '', 'itemsizing': 'trace'})
    gc_content_figure.update_layout(font_family="Ubuntu",
                                    title={'text': f"{sample_name if sample_name is not None else ''} "
                                                   f"Original vs. Corrected GC Content" +
                                                   ((' (' + annotation + ')') if annotation is not None else ''),
                                           'font': {'family': 'Ubuntu', 'size': fig_fontsize + 8,
                                                    'color': 'rgb(20, 20, 20)'},
                                           'xanchor': 'center', 'yanchor': 'middle', 'x': 0.5})
    gc_content_figure.update_xaxes(range=[10, 85])
    gc_content_figure.update_layout(margin=dict(l=5, r=5, t=40, b=10))
    if show_figure:
        gc_content_figure.show()
    out_dir_path.mkdir(parents=True, exist_ok=True)
    if annotation is not None:
        transformed_anno = re.sub('--', '-',
                                  re.sub(':', '-',
                                         re.sub(',', '_',
                                                re.sub(' ', '-', annotation))))
    else:
        transformed_anno = None
    for image_format in image_formats:
        out_file = out_dir_path / \
                   ((f'{orig_samples[0]}_' if len(orig_samples) == 1 else '') +
                    f"GCparagon_GC-content-comparison_GC-bias-correction"
                    f"{'_SPLINE' if spline_interpolation else '_LINEAR'}"
                    f"{'_MARKERS' if markers else ''}{('_' + transformed_anno) if annotation else ''}"
                    f"_cfDNAref.{image_format}")
        gc_content_figure.write_image(out_file)


def main():
    # read all-samples distributions file
    # NOTE: original distributions should be EQUAL between algorithms since the entire BAM file is parsed.
    original_gc_distributions_all_samples, corrected_gc_distributions_all_samples_algorithms = read_stats(
        stats_path=all_samples_dists_file, return_as_percentage=True)
    if not original_gc_distributions_all_samples or not corrected_gc_distributions_all_samples_algorithms:
        return 1
    # iterate over target distribution creation algorithms (different reference fragment length distributions)
    for sample, ref_dist_path_sample in reference_gc_dists_samples.items():
        ref_algo = 'entireBAM'
        reference_gc_dists = {}
        reference_gc_dists[f'{sample}'] = read_tsv_data_to_dict(
            table_path=ref_dist_path_sample, normalize_to_100=False)
        # algorithm observations (GCparagon: 500 M reads simulated from reference DNA
        # -> universally applicable reference distribution!)
        for plot_markers in (False, ):  # (True, False):
            residual_cumulative_deviation = {}
            corrected_single_sample_all_algorithms = {}
            for algo in corrected_gc_distributions_all_samples_algorithms.keys():  # add corrected lines for algorithms
                corrected_single_sample_all_algorithms.update(
                    {f'{sample} ({algo})': corrected_gc_distributions_all_samples_algorithms[algo][sample]})
                # compute cumulative deviation (in percentage)
                residual_cumulative_deviation[algo] = sum(  # in "sum of relative percentage" units (%)
                    [abs(reference_gc_dists[f'{sample}'][gc_idx] - cv)
                     for gc_idx, cv in enumerate(corrected_single_sample_all_algorithms[f'{sample} ({algo})'])])
            print(f"{sample} - cumulative residual bias (CRB):\n"
                  f"   Griffin:   {residual_cumulative_deviation['Griffin']:.2f}%\n"
                  f"   GCparagon: {residual_cumulative_deviation['GCparagon']:.2f}%")
            plot_fragment_gc_dists_griffin(
                original_gc_data={f'{sample}': original_gc_distributions_all_samples['GCparagon'][sample]},
                fig_width=1200, fig_height=800, fig_fontsize=32,
                corrected_gc_data=corrected_single_sample_all_algorithms,
                reference_normalized=True, reduced_bins=True, sample_name=sample,
                # for single preset Griffin_matrix_plots only!
                out_dir_path=output_dir / f"reference_created_from_{ref_algo}",
                reference_dists=reference_gc_dists, show_figure=False, normalize_to_dataset_size=True,
                spline_interpolation=SPLINE_INTERPOLATION, markers=plot_markers, image_formats=IMAGE_FORMATS)
    return 0


if __name__ == '__main__':
    sys.exit(main())


# CMLINE OUTPUT:
# -------------------------------------------
# B01 - cumulative residual bias (CRB):
#    Griffin:   2.17%
#    GCparagon: 4.22%
# H01 - cumulative residual bias (CRB):
#    Griffin:   1.82%
#    GCparagon: 2.32%
# C01 - cumulative residual bias (CRB):
#    Griffin:   1.62%
#    GCparagon: 5.83%
# P01 - cumulative residual bias (CRB):
#    Griffin:   1.44%
#    GCparagon: 9.44%
