#!/usr/bin/env python3
import sys
import numpy as np
from pathlib import Path
from collections import defaultdict
from typing import Optional, Union, Dict, Tuple
from src.GCparagon.utilities import plot_fragment_gc_dists

SCRIPT_ROOT_PATH = Path(__file__).parent
SCRIPT_ROOT_DIR = str(SCRIPT_ROOT_PATH)

output_dir = SCRIPT_ROOT_PATH / 'fragment_GCcontent_plots'  # output for final plots
all_samples_dists_file = SCRIPT_ROOT_PATH / \
                         'STATISTICS_allSamples_GC-percentage_perFragmentSequence_FULL_24-02-2023.tsv'
dist_files = {'B01': SCRIPT_ROOT_PATH / 'FRAGMENT_STATISTICS_B01_GC-percentage_perFragmentSequence_FULL_07-03-2023.tsv',
              'C01': SCRIPT_ROOT_PATH / 'FRAGMENT_STATISTICS_C01_GC-percentage_perFragmentSequence_FULL_07-03-2023.tsv',
              'H01': SCRIPT_ROOT_PATH / 'FRAGMENT_STATISTICS_H01_GC-percentage_perFragmentSequence_FULL_07-03-2023.tsv',
              'P01': SCRIPT_ROOT_PATH / 'FRAGMENT_STATISTICS_P01_GC-percentage_perFragmentSequence_FULL_07-03-2023.tsv'}
reference_gc_hg38_cfDNA = {
    'B01': SCRIPT_ROOT_PATH / 'hg38_ref_dists/B01-preset3-fragLenDist_ref_GC_dist_500Mreads.tsv',
    'H01': SCRIPT_ROOT_PATH / 'hg38_ref_dists/H01-preset3-fragLenDist_ref_GC_dist_500Mreads.tsv',
    'C01': SCRIPT_ROOT_PATH / 'hg38_ref_dists/C01-preset3-fragLenDist_ref_GC_dist_500Mreads.tsv',
    'P01': SCRIPT_ROOT_PATH / 'hg38_ref_dists/P01-preset3-fragLenDist_ref_GC_dist_500Mreads.tsv'}

SPLINE_INTERPOLATION = True


def read_stats(stats_path: Union[str, Path]) -> Optional[Tuple[Dict[str, Dict[str, np.array]]]]:
    # output: {PRESET: {sample: 2D_numpy_array(matrix)}
    with open(stats_path, 'rt') as f_stats:
        stats_header = f_stats.readline()
        stats_data = f_stats.readlines()
    header_content = stats_header.strip().split('\t')
    leading_class_columns = header_content[:[ct.isdigit() for ct in header_content].index(True)]  # positive integer?
    # sample{any_string}, preset{0,1,2,3}, status(original, corrected)
    try:
        status_column = leading_class_columns.index('status')
    except ValueError:
        try:
            status_column = leading_class_columns.index('Status')
        except ValueError:
            print(f"'preset' column not found in header. Exiting ...")
            return None
    try:
        preset_column = leading_class_columns.index('preset')
    except ValueError:
        try:
            preset_column = leading_class_columns.index('Preset')
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
        match cur_smpl:
            case 'B58_3':
                cur_smpl = 'B01'
            case 'P244_7':
                cur_smpl = 'P01'
            case 'C219_5':
                cur_smpl = 'C01'
            case 'NPH_011':
                cur_smpl = 'H01'
        cur_prst = line_content[preset_column]
        cur_data = np.array(line_content[3:], dtype=float)
        assert len(cur_data) == 101
        try:
            return_data[cur_stt][cur_prst].update({f'{cur_smpl}': cur_data})
        except KeyError:
            return_data[cur_stt].update({cur_prst: {f'{cur_smpl}': cur_data}})
    if return_data == {'original': {}, 'corrected': {}}:
        return None
    return original_data, corrected_data


def read_normalized_tsv_data_to_dict(table_path: Path) -> defaultdict[int, Union[float, int]]:
    with open(table_path, 'rt') as f_tab:
        _header = f_tab.readline()
        file_content = [line.strip().split('\t') for line in f_tab.readlines()]
    content_dict = defaultdict(float)
    for gc_perc, frequency_value in file_content:
        content_dict[int(gc_perc)] = float(frequency_value)
    return content_dict


def main():
    for sample, dists_file in dist_files.items():
        original_gc_distributions_per_preset, corrected_gc_distributions_per_preset = read_stats(stats_path=dists_file)
        if not original_gc_distributions_per_preset or not corrected_gc_distributions_per_preset:
            return 1
        # read in reference dist
        # reference_gc_dists = defaultdict(None)
        # # required for all samples in a preset plot:
        # for sample_id, hg38_cfDNA_ref_dist in reference_gc_hg38_cfDNA.items():
        #     reference_gc_dists[sample_id] = read_normalized_tsv_data_to_dict(table_path=Path(hg38_cfDNA_ref_dist))
        # # required for single sample, per preset plots and all presets per sample plots:
        reference_gc_dists = defaultdict()
        reference_gc_dists[sample] = read_normalized_tsv_data_to_dict(table_path=Path(reference_gc_hg38_cfDNA[sample]))
        # reduce number of bins by factor 2:
        # # required for single sample per preset plots:
        # for preset in sorted(list(original_gc_distributions_per_preset.keys())):  # 'corrected' and 'original'
        for plot_markers in (False, ):  # (True, False):
            original_all_presets = {}
            corrected_all_presets = {}
            for preset in sorted(list(original_gc_distributions_per_preset.keys())):  # 'corrected' and 'original'
                original_all_presets.update({f'{sample}, Preset {preset}':
                                             original_gc_distributions_per_preset[preset][sample]})
                corrected_all_presets.update({f'{sample}, Preset {preset}':
                                              corrected_gc_distributions_per_preset[preset][sample]})
            plot_fragment_gc_dists(original_gc_data=original_all_presets,
                                   # single/multiple sample(s) per preset: original_gc_distributions_per_preset[preset],
                                   corrected_gc_data=corrected_all_presets,
                                   # corrected_gc_distributions_per_preset[preset],  # for single preset plots only
                                   out_dir_path=Path(output_dir), fig_width=1200, fig_height=800, fig_fontsize=24,
                                   normalize_to_dataset_size=True, annotation=f'{sample}-allPresets',
                                   # annotation=f'Preset{preset}',
                                   reduced_bins=True, reference_dists=reference_gc_dists,
                                   spline_interpolation=SPLINE_INTERPOLATION, markers=plot_markers)
    return 0


if __name__ == '__main__':
    sys.exit(main())
