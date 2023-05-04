#!/usr/bin/env python3

import sys
from pathlib import Path

REPO_ROOT_DIR = Path(__file__).parent.parent.parent.parent
PRESET_DIR = REPO_ROOT_DIR / 'preset_computation'
benchmark_results_path = PRESET_DIR / 'benchmark_results'

OUTPUT_PATH = benchmark_results_path / 'benchmarking_summary.tsv'


def main() -> int:
    candidate_dirs = tuple(benchmark_results_path.glob('mprof_*'))
    execution_data = {}
    if len(candidate_dirs) == 0:
        raise ValueError(f"nothing found starting with 'mprof_' in directory {benchmark_results_path}!")
    for profiling_result in candidate_dirs:
        try:
            results_file = tuple(profiling_result.glob('benchmark_results_*'))[0]
        except IndexError:
            print(f"no profiling results found in {profiling_result}! Ignorng directory ..")
            continue
        with open(results_file, 'rt') as f_br:
            cur_line = f_br.readline()
            command_line_content = cur_line.strip().split(' ')
            try:
                sample_name = Path(command_line_content[command_line_content.index('--bam') + 1]).stem
            except ValueError:
                sample_name = Path(command_line_content[command_line_content.index('-b') + 1]).stem
            try:
                preset = int(command_line_content[command_line_content.index('--use-parameter-preset') + 1])
            except ValueError:
                preset = int(command_line_content[command_line_content.index('-up') + 1])
            if execution_data.get(preset) is None:  # would be dict otherwise
                execution_data.update({preset: {}})
            while cur_line[:len('iteration')] != 'iteration':
                cur_line = f_br.readline()
            _hdr = cur_line.strip()
            iter_data = f_br.readline().strip()
            while iter_data != '':  # data is followed by two empty lines
                iteration, max_mem_usage, seconds_duration, _timestamp = iter_data.split('\t')
                try:
                    execution_data[preset][sample_name].append([float(max_mem_usage), float(seconds_duration), None])
                except KeyError:
                    execution_data[preset].update({sample_name: [[float(max_mem_usage),
                                                                  float(seconds_duration), None]]})
                iter_data = f_br.readline().strip()
        # iteration	max. memory usage (MiB)	duration (s)	timestamp (d-m-Y,H:M:S)
        # 1	2490.64	110.597	02-05-2023,17:22:34
        # 2	2430.17	110.504	02-05-2023,17:24:28
    # gather processed fragments statistics
    preset_candidate_dirs = tuple(PRESET_DIR.glob('*'))
    if len(preset_candidate_dirs) == 0:
        raise ValueError(f"nothing found in profiling results directory {benchmark_results_path}!")
    for preset_dir in preset_candidate_dirs:
        try:
            preset = int(preset_dir.stem)
        except ValueError:
            # not a preset directory; will be ignored
            continue
        sample_preset_dirs = tuple(preset_dir.glob('*'))
        if len(sample_preset_dirs) == 0:
            raise ValueError(f"nothing found in preset computation output directory {preset_dir}!")
        for sample_dir in sample_preset_dirs:
            sample_name = sample_dir.stem
            try:
                log_file = tuple(sample_dir.glob('*GCbiasCorrection*.log'))[0]
            except IndexError:
                print(f"no log file found in directory {sample_dir}")
            log_line = ''
            with open(log_file, 'rt') as f_log:
                itr = 0
                while 'Target number of fragments to process reached:' not in log_line and \
                        'Could not reach target count of processed fragments: ' not in log_line:
                    log_line = f_log.readline()  # skip header
                    itr += 1
                    if itr == 1000:
                        raise ValueError(f"no processed fragments information found in log file {log_file}!")
                if 'Could not reach target count of processed fragments: ' in log_line:
                    processed_frags_str = log_line.split('Could not reach target count of processed fragments: ') \
                        [1].split('/')[0]
                else:
                    processed_frags_str = log_line.split('Target number of fragments to process reached: ') \
                        [1].split('/')[0]
                while ',' in processed_frags_str:
                    processed_frags_str = processed_frags_str.replace(',', '')
                processed_fragments = int(processed_frags_str)
            for i, (mem, dur, _) in enumerate(execution_data[preset][sample_name]):
                execution_data[preset][sample_name][i][2] = processed_fragments
    # write tabular output
    out_lines = ['sample\tpreset\titeration\tmemory_consumption\tduration\tprocessed_fragments\n']
    for sample in sorted(list(execution_data[1].keys())):
        for preset in sorted(list(execution_data.keys())):
            for i, (mem, dur, frags) in enumerate(execution_data[preset][sample]):
                out_lines.append(f'{sample}\t{preset}\t{i+1}\t{mem}\t{dur}\t{frags}\n')
    with open(OUTPUT_PATH, 'wt') as f_out:
        f_out.writelines(out_lines)
    return 0


if __name__ == '__main__':
    sys.exit(main())
