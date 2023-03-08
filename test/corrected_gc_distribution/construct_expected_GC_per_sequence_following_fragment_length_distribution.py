#!/usr/bin/python3
import sys
import gzip
import numpy as np
from pathlib import Path
from random import randint
import multiprocessing as mp
from collections import defaultdict, deque
from twobitreader import TwoBitFile
from typing import Tuple, Dict, Union, List

SOURCE_CODE_ROOT_PATH = Path(__file__).parent.parent.parent
SOURCE_CODE_ROOT_DIR = str(SOURCE_CODE_ROOT_PATH)

mp.set_start_method('spawn', force=True)  # for safety reasons, don't fork (don't screw up file handles)

hg38_2bit_ref_genome = SOURCE_CODE_ROOT_PATH / '2bit_reference/hg38.analysisSet.2bit'  # needs to be downloaded using
# the EXECUTE_reference_download.sh script
hg38_ref_genome_chrom_sizes = SOURCE_CODE_ROOT_PATH / '2bit_reference/hg38.analysisSet.chrom.sizes'
draw_n_fragments = 500*10**6  # five hundred million
n_processes = 24  # actually 1 per stdchrom -> chr1 will take the longest

SCRIPT_ROOT_PATH = Path(__file__).parent
out_path = SCRIPT_ROOT_PATH / 'simulated_hg38_ref_dists'

observed_fragments = {'preset3': {'H01': SOURCE_CODE_ROOT_PATH /
                                  'preset_computation/3/H01/H01_observed_attributes_matrix.txt.gz',
                                  'P01': SOURCE_CODE_ROOT_PATH /
                                  'preset_computation/3/P01/P01_observed_attributes_matrix.txt.gz',
                                  'C01': SOURCE_CODE_ROOT_PATH /
                                  'preset_computation/3/C01/C01_observed_attributes_matrix.txt.gz',
                                  'B01': SOURCE_CODE_ROOT_PATH /
                                  'preset_computation/3/B01/B01_observed_attributes_matrix.txt.gz'}}


def load_txt_to_matrix_with_meta(filename: Union[str, Path], to_dtype=np.float64) -> Tuple[np.array, range]:
    """

    :param loading_logger:
    :param filename:
    :param to_dtype:
    :return:
    """
    print(f"Loading statistic matrix from {filename}")
    statistic_matrix = np.loadtxt(filename, delimiter='|', skiprows=0, dtype=to_dtype)
    with gzip.open(str(filename), 'rt') as f_mat:
        hdr = f_mat.readline()
    elements = hdr.split('# rows representing fragment lengths (')[1].split()  # split on whitespace + trim empty fields
    fragment_lengths = range(int(elements[0]), int(elements[3])+1)  # end exclusive
    return statistic_matrix, fragment_lengths


def load_reference_fragment_length_dist(total_simulate_fragments: int,
                                        observed_fragments_matrix_file: Union[str, Path]) -> defaultdict[int, int]:
    observed_flen_gc_np_mat, flen_range = load_txt_to_matrix_with_meta(filename=observed_fragments_matrix_file)
    fragment_lengths_absolute = np.sum(observed_flen_gc_np_mat, axis=1)
    fragment_lengths_relative = fragment_lengths_absolute / fragment_lengths_absolute.sum()
    fragment_length_shares = defaultdict(float)
    for f_idx, frag_len in enumerate(range(flen_range.start, flen_range.stop)):
        fragment_length_shares[frag_len] = fragment_lengths_relative[f_idx]
    sample_target_fragment_lengths = defaultdict(int)
    for f_idx, (frag_len, share) in enumerate(fragment_length_shares.items()):
        sample_target_fragment_lengths[frag_len] = int(round(share * total_simulate_fragments, ndigits=0))
    return sample_target_fragment_lengths


def compute_sampling_numbers_per_chrom(chrom_sizes_path: Path, target_frags_per_flen: defaultdict[int, int]) \
        -> Tuple[defaultdict[str, Dict[int, int]], Dict[str, int]]:
    # read chrom sizes file
    relevant_chrom_sizes = {}
    relevant_chrom_samples = defaultdict(dict)
    total_size = 0
    with open(chrom_sizes_path, 'rt') as f_chrms:
        for line in f_chrms:
            chrom, size = line.strip().split('\t')
            # if chrom != 'chr21':
            if '_' in chrom or 'chrY' in chrom or 'chrEBV' in chrom or 'chrM' in chrom or 'random' in chrom.lower() or \
                    'unk' in chrom.lower():  # ignore line -> chrY is too short
                # -> cannot sample 3.2k fragments from that thing using 5000 random positions!
                continue
            total_size += int(size)
            relevant_chrom_sizes[chrom] = int(size)
    # compute number of fragments per chromosomes
    for chrm, size in relevant_chrom_sizes.items():
        for frag_len, n_samples in target_frags_per_flen.items():
            relevant_chrom_samples[chrm].update({frag_len: int(round(n_samples * size / total_size, ndigits=0))})
    print("Checking target fragment samples balance now (max +-1% of total observed fragment countn allowed)..")
    assert sum(target_frags_per_flen.values())*0.99 < \
           sum([sum(chrom_samples.values()) for chrom, chrom_samples in relevant_chrom_samples.items()]) < \
           sum(target_frags_per_flen.values())*1.01  # within +-1% deviation of target value
    return relevant_chrom_samples, relevant_chrom_sizes


def sample_worker(samples_per_flen_per_scaffold: List[Tuple[str, Dict[int, int]]], ref_genome_handle: TwoBitFile,
                  chromosome_sizes: Dict[str, int], sender_connection: mp.Pipe):
    gathered_gc_samples = defaultdict(int)  # counts of GC% samples
    random_number_generator = np.random.default_rng(seed=randint(0, 999))  # use random seed for reproducibility here!
    for chrom, sampling_data in samples_per_flen_per_scaffold:
        if chromosome_sizes[chrom] < 20000 + max(sampling_data.keys()):
            print(f"can't sample scaffold '{chrom}' - too small! ({chromosome_sizes[chrom]:,} bp). Skipping ..")
            continue
        chrom_handle = ref_genome_handle[chrom]
        for frag_len, n_samples in sampling_data.items():
            rand_samp_pos = random_number_generator.integers(low=10000, high=chromosome_sizes[chrom] - frag_len - 9999,
                                                             size=max(int(n_samples * 2), 5000))
            # skip 10k telomeres and create at least 5000 random integers
            if int(n_samples * 1.5) > 5000:
                rand_samp_pos.sort()  # for sweeping-like file pointer placement (if that even matters for twobitreader)
            chrom_sampling_offset = 0
            try:
                for sample_idx in range(n_samples):
                    sample_pos = rand_samp_pos[sample_idx + chrom_sampling_offset]
                    cur_fseq = chrom_handle[sample_pos:sample_pos + frag_len].upper()
                    while 'N' in cur_fseq:
                        chrom_sampling_offset += 1
                        sample_pos = rand_samp_pos[sample_idx + chrom_sampling_offset]
                        cur_fseq = chrom_handle[sample_pos:sample_pos + frag_len].upper()
                    gathered_gc_samples[int(round(
                        (cur_fseq.count('G') + cur_fseq.count('C')) * 100. / frag_len, ndigits=0))] += 1
            except IndexError as e:
                print(f"could not sample {n_samples:,} fragments from two-bit reference genome file for chromosome "
                      f"'{chrom}'. Adapt!\nError was:\n{e}\nExiting now..", flush=True)
                sender_connection.send(None)
                sender_connection.close()
                return
    sender_connection.send(gathered_gc_samples)
    sender_connection.close()


def sample_ref_genome(ref_genome_file: Union[str, Path], chom_sizes: Dict[str, int],
                      sample_n_from_chroms_per_flen: defaultdict[str, Dict[int, int]]) \
        -> Tuple[defaultdict[int, int], bool]:
    samplers = []
    receivers = []
    chroms_per_proc = []
    cum_lengths = []
    for _i in range(n_processes):
        chroms_per_proc.append([])
        cum_lengths.append(0)
    target_sum_bases = sum(chom_sizes.values()) // n_processes
    process_offset = 0
    for chrom_idx, (chrom, bases) in enumerate(chom_sizes.items()):
        added = False
        process_index = chrom_idx % n_processes
        if cum_lengths[process_index + process_offset] >= target_sum_bases:  # increase the offset
            while cum_lengths[process_index + process_offset] >= target_sum_bases:
                if process_offset == n_processes:  # all have reached the target base sum -> assign to smalles list
                    smalles_scaffold_list_process_index = cum_lengths.index(min(cum_lengths))
                    chroms_per_proc[smalles_scaffold_list_process_index].append(chrom)
                    cum_lengths[smalles_scaffold_list_process_index] += bases
                    added = True
                    break
                process_offset += 1
            if added:
                continue
        chroms_per_proc[process_index + process_offset].append(chrom)
        cum_lengths[process_index + process_offset] += bases
    # create workers
    for p_idx, chrom_list in enumerate(chroms_per_proc):
        if not chrom_list:  # don't create workers for empty chromosome lists! (number parallel processes > chromosomes)
            continue
        receiver_end, sender_end = mp.Pipe(duplex=False)
        receivers.append(receiver_end)
        worker_ref_handle = TwoBitFile(str(ref_genome_file))
        samplers.append(mp.Process(target=sample_worker,
                                   kwargs={'samples_per_flen_per_scaffold': [(chrm, sample_n_from_chroms_per_flen[chrm])
                                                                             for chrm in chrom_list],
                                           'ref_genome_handle': worker_ref_handle,
                                           'chromosome_sizes': chom_sizes,
                                           'sender_connection': sender_end}))
    # start processes
    _ = deque([smplr.start() for smplr in samplers], maxlen=0)
    # receive all results
    results = []
    for rcvr in receivers:
        results.append(rcvr.recv())
    # check if any sampling procedure failed:
    combined_gc_frequencies = defaultdict(int)
    broken = False
    if None in results:
        broken = True
        print(f"ERROR - received None from some subprocess. Could not sample refrence genome successfully.")
    for res in results:
        if res is None:  # skip buggy chromosome list
            continue
        for gc_pc, sample_count in res.items():
            combined_gc_frequencies[gc_pc] += sample_count
    return combined_gc_frequencies, broken


def write_gc_percentages(gc_percentages: defaultdict, sample_name: str, preset: int, was_failure: bool):
    ref_name = Path(hg38_2bit_ref_genome).name
    out_file = out_path / f'{ref_name}_{sample_name}-preset{preset}-fragLenDist_ref_GC_dist_{draw_n_fragments//10**6}' \
                          f"Mreads{'-FAILED' if was_failure else ''}.tsv"
    out_path.mkdir(parents=True, exist_ok=True)
    out_lines = ['gc_percentage\trelative_frequency\n']
    total_sum = sum(gc_percentages.values())  # must not be 0 but can't in a non-empty dataset reality
    out_lines.extend([f'{gc_perc}\t{gc_percentages[gc_perc] * 100. /total_sum}\n'
                      for gc_perc in sorted(list(gc_percentages.keys()))])
    with open(out_file, 'wt') as f_out:
        f_out.writelines(out_lines)


def main():
    for sample_id, observed_frags in observed_fragments['preset3'].items():  # USE PRESET 3 distributions as these are
        # the most complete approximations of the entire dataset! (the fragment length distribution is not to be
        # corrected) -> Nonsense question: "what would a universally normal fragment length distribution be that we
        # could expect?" The GC bias does not change the length of the fragments but the relative observed frequency!
        sample_fragments_per_fragment_lengths = load_reference_fragment_length_dist(
            total_simulate_fragments=draw_n_fragments, observed_fragments_matrix_file=observed_frags)
        # read sizes of standard ref contigs and compute sampling proportions
        sample_n_from_chroms, chom_sizes = compute_sampling_numbers_per_chrom(
            chrom_sizes_path=Path(hg38_ref_genome_chrom_sizes),
            target_frags_per_flen=sample_fragments_per_fragment_lengths)
        # sample using multiprocessing!
        gc_percs, failure = sample_ref_genome(ref_genome_file=hg38_2bit_ref_genome, chom_sizes=chom_sizes,
                                              sample_n_from_chroms_per_flen=sample_n_from_chroms)
        if gc_percs is None:
            return 1
        # write!
        write_gc_percentages(gc_percentages=gc_percs, was_failure=failure, sample_name=sample_id, preset=3)


if __name__ == '__main__':
    sys.exit(main())
