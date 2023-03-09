#!/usr/bin/python3
import os
import sys
import time
import datetime
import numpy as np
from pathlib import Path
from argparse import ArgumentParser, Namespace


# ----------- DEFAULT SETUP ------------
DEFAULT_SAMPLING_INTERVAL = 5  # seconds
DEFAULT_REPETITIONS = 3
# --------------------------------------


def get_cmdline_args() -> tuple[Namespace, list[str]]:
    prsr = ArgumentParser()
    prsr.add_argument('-s', '--script', dest='script_path', required=True,
                      help='Path to a Python script that should be benchmarked. [ REQUIRED ]')
    prsr.add_argument('-i', '--iter', dest='n_iteration', default=DEFAULT_REPETITIONS,
                      help='Repeat benchmarking n times and build the average, min, max and median. '
                           f'[ DEFAULT: {DEFAULT_REPETITIONS} repetitions ]')
    prsr.add_argument('-t', '--sampling-frequency', dest='sample_freq', default=DEFAULT_SAMPLING_INTERVAL,
                      help=f'Repeat sampling every t seconds. [ DEFAULT: {DEFAULT_SAMPLING_INTERVAL} seconds ]')
    cwd = Path().absolute()
    prsr.add_argument('-o', '--output-path', dest='abs_path', default=cwd,
                      help='Optional: Absolute path to the output directory. If none is specified, the current '
                           f'working directory will be used. [ DEFAULT: {cwd} ]')
    return prsr.parse_known_args()


def main() -> int:
    args = get_cmdline_args()
    execute_time = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    execute_time_print = datetime.datetime.now().strftime("%d-%m-%Y, %H:%M:%S")
    script_name = Path(args[0].script_path).name
    output_dir = Path(args[0].abs_path) / f'mprof_{script_name}_{execute_time}'
    output_dir.mkdir(parents=True, exist_ok=True)
    arguments = ' '.join(sys.argv[sys.argv.index(args[0].script_path):])
    info_message = f'Command called: {arguments} \n' \
                   f'Date and time of execution: {execute_time_print}\n' \
                   f'Repeat computation {args[0].n_iteration} times, sampling every {args[0].sample_freq} seconds\n'
    stats_file_out = f'{output_dir}/benchmark_results_{script_name}_{execute_time}.txt'
    with open(stats_file_out, 'wt') as f_bench:
        f_bench.write(info_message + f'\niteration\tmax. memory usage (MiB)\tduration (s)\ttimestamp (d-m-Y,H:M:S)\n')
    statistic_matrix = np.zeros((int(args[0].n_iteration), 2), )  # 1st column is time, 2nd column is memory usage
    print(info_message, flush=True)
    for ind in range(int(args[0].n_iteration)):
        print(f"------------------------------------------------------------------------------\n\n"
              f"Starting iteration: {ind + 1}", flush=True)
        dat_file = f'{output_dir}/mprofile_{script_name}_{ind + 1}.dat'
        os.system('mprof clean')
        iter_time = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        iter_time_print = datetime.datetime.now().strftime("%d-%m-%Y,%H:%M:%S")
        start_iter = time.perf_counter()
        os.system(f'mprof run -o "{dat_file}" -T {args[0].sample_freq} {arguments}')
        end_iter = time.perf_counter()
        statistic_matrix[ind, 0] = round(end_iter - start_iter, 3)
        statistic_matrix[ind, 1] = round(max((np.loadtxt(dat_file, dtype=str, usecols=1)[1:]).astype(float)), ndigits=2)
        print(f'Finished in {statistic_matrix[ind, 0]} seconds', flush=True)
        title = f'{script_name}, {iter_time}\n' \
                f' Maximum memory usage in MiB: {statistic_matrix[ind, 1]}'
        if os.system(f'mprof plot -o {output_dir}/plot_{ind + 1}.png --title "{title}" {dat_file}') != 0:
            return 1  # os.system should return int
        with open(stats_file_out, 'at') as f_bench:
            f_bench.write(f'{ind + 1}\t{statistic_matrix[ind, 1]}\t{statistic_matrix[ind, 0]}\t{iter_time_print}\n')
    with open(stats_file_out, 'at') as f_bench:
        f_bench.write(f"\n\nOverall time statistic (in seconds):\n Average execution time: "
                      f"{np.mean(statistic_matrix, 0)[0]}\n"
                      f" Min execution time: {np.min(statistic_matrix, 0)[0]}\n"
                      f" Max execution time: {np.max(statistic_matrix, 0)[0]}"
                      f"\n\nAverage memory consumption (in MiB): {np.mean(statistic_matrix, 0)[1]}\n"
                      f" Min memory usage: {np.min(statistic_matrix, 0)[1]}\n"
                      f" Max memory usage: {np.max(statistic_matrix, 0)[1]}")
    print(f"\n------------------------------------------------------------------------------\n"
          f"Profiling finished. Output is in: {output_dir}", flush=True)
    return 0


if __name__ == '__main__':
    sys.exit(main())
