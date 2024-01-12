#!/usr/bin/env python3

import sys
import gzip
import numpy as np
from pathlib import Path

from random import randint
import multiprocessing as mp
from collections import defaultdict, deque
from twobitreader import TwoBitFile
from typing import Tuple, Dict, Union as OneOf, List

SOURCE_CODE_ROOT_PATH = Path(__file__).parent.parent.parent  # ../src/GCparagon
SOURCE_CODE_ROOT_DIR = str(SOURCE_CODE_ROOT_PATH)

mp.set_start_method('spawn', force=True)  # for safety reasons, don't fork (don't screw up file handles)

# GENOME BUILD FILE DEFINITIONS - set according to your target genome build!
# ----------------------------------------------------------------------------------------------------------------------
GENOME_BUILD = 'hg19'
two_bit_ref_genome = (SOURCE_CODE_ROOT_PATH /
                      f'src/GCparagon/2bit_reference/{GENOME_BUILD}.2bit')
if not two_bit_ref_genome.is_file():
    raise FileNotFoundError(f"required 2bit reference file not found at: '{two_bit_ref_genome}'")
# ^--- needs to be downloaded using the EXECUTE_reference_download.sh script or a commented-out command inside
ref_genome_chrom_sizes = (SOURCE_CODE_ROOT_PATH /
                          f'src/GCparagon/2bit_reference/{GENOME_BUILD}.chrom.sizes')
draw_n_fragments = 500*10**6  # five hundred million
n_processes = 4  #24  # if 24, actually 1 used per std-chrom -> chr1 will take the longest

SCRIPT_ROOT_PATH = Path(__file__).parent

out_path = SCRIPT_ROOT_PATH.parent.parent  # ../GCparagon/accessory_files
# !!!!!
# READ THE FOLLOWING
# !!!!!
# CURRENT ASSAY: PlasmaSeq
# CURRENT SAMPLE TYPE: blood plasma
# SEQUENCING TECHNOLOGY: Illumina (paired-end short reads)
# WARNING: if your sample type and/or wet lab procedure (and/or the sequencing technology) differes (significantly),
#          it is STRONGLY recommended to recompute the fragment length distribution percentages! (20-800bp fragments)
# WARNING: GCparagon was designed for untargeted, PAIRED-END sequencing read data of ccfDNA from plasma samples!
# ----------------------------------------------------------------------------------------------------------------------
observed_fragments_tsv_file = SCRIPT_ROOT_PATH.parent / 'plasmaSeq_ccfDNA_reference_fragment_length_distribution.tsv'
# ^--- USE THIS REFERENCE DISTRIBUTION ONLY FOR cfDNA from PLASMA and PlasmaSeq library prep!
# IF YOU USE ANOTHER WET LAB PROCEDURE AND/OR ANALYTE, CREATE THIS REFERENCE DISTRIBUTION YOURSELF FROM AT LEAST 4
# SAMPLES WITH MATCHING WET LAB PROCEDURE AND IDENTICAL ANALYTE TYPE (urin, CSF, sputum, etc.) !!!!!!!!!!!!!!!!!!!!!!!!!
# ----------------------------------------------------------------------------------------------------------------------


def load_txt_to_matrix_with_meta(filename: OneOf[str, Path], to_dtype=np.float64) -> Tuple[np.array, range]:
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


def load_fragment_length_tsv(tsv_path: Path) -> Tuple[Dict[int, float], range]:
    all_flength_counts = defaultdict(float)
    with open(tsv_path, 'rt') as f_flength_stats:
        hdr_content = f_flength_stats.readline().strip().split('\t')
        if any([not s.isnumeric() for s in hdr_content]):
            raise ValueError(f"unexpected content of reference fragment length distribution! "
                             f"An element in the header line was not numeric")
        hdr_content = [int(s) for s in hdr_content]
        sorted_hdr_content = sorted(hdr_content)
        min_flength = sorted_hdr_content[0]
        max_flength = sorted_hdr_content[-1]
        flength_range = range(min_flength, max_flength, 1)
        if len(sorted_hdr_content) != max_flength - min_flength + 1:
            raise ValueError(f"unexpected fragment length step - assumed each fragment length between {min_flength:,} "
                             f"bp and {max_flength:,} bp was present but "
                             f"{max_flength - min_flength + 1 - len(sorted_hdr_content):,} fragment lengths are "
                             f"missing!")
        length_frequencies = f_flength_stats.readline().strip().split('\t')
        for flength_idx in range(len(sorted_hdr_content)):
            all_flength_counts[min_flength + flength_idx] = float(hdr_content[flength_idx])
    return all_flength_counts, flength_range


def load_reference_fragment_length_dist(fragment_number_to_draw: int, observed_fragments_tsv_file: OneOf[str, Path]) \
        -> defaultdict[int, int]:
    flength_counts, flength_range = load_fragment_length_tsv(tsv_path=observed_fragments_tsv_file)
    sample_these_flengths = defaultdict(int)
    fragment_count_sum = sum(flength_counts.values())
    for f_idx, (frag_len, count) in enumerate(flength_counts.items()):
        sample_these_flengths[frag_len] = round(
            count * fragment_number_to_draw / fragment_count_sum)  # round without 'ndigits=' returns int
    return sample_these_flengths


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
                    'unk' in chrom.lower():  # ignore line
                # -> chrY is too short: cannot sample 3.2k fragments from that thing using 5000 random positions!
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


def sample_ref_genome(ref_genome_file: OneOf[str, Path], chom_sizes: Dict[str, int],
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
    _ = deque(map(lambda w: w.start(), samplers), maxlen=0)  # deque to consume map iterator
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


def write_gc_percentages(gc_percentages: defaultdict, was_failure: bool):
    ref_name = Path(two_bit_ref_genome).stem
    out_file = out_path / f"{ref_name}_reference_GC_content_distribution{'-FAILED' if was_failure else ''}.tsv"
    out_path.mkdir(parents=True, exist_ok=True)  # should exist
    out_lines = ['gc_percentage\trelative_frequency\n']
    total_sum = sum(gc_percentages.values())  # must not be 0 but can't in a non-empty dataset reality
    out_lines.extend([f'{gc_perc}\t{gc_percentages[gc_perc] * 100. /total_sum}\n'
                      for gc_perc in sorted(list(gc_percentages.keys()))])
    with open(out_file, 'wt') as f_out:
        f_out.writelines(out_lines)


def main():
    fragment_lengths_to_draw = load_reference_fragment_length_dist(
        fragment_number_to_draw=draw_n_fragments, observed_fragments_tsv_file=observed_fragments_tsv_file)
    # read sizes of standard ref contigs and compute sampling proportions
    sample_n_from_chroms, chom_sizes = compute_sampling_numbers_per_chrom(
        chrom_sizes_path=Path(ref_genome_chrom_sizes),
        target_frags_per_flen=fragment_lengths_to_draw)
    # sample using multiprocessing!
    print(f"Will sample the reference genome now evenly for {draw_n_fragments:,} fragments. "
          "High IOPS expected! This will take some time (30 minutes+) depending on your hardware setup. "
          "Please be patient..")
    gc_percs, failure = sample_ref_genome(ref_genome_file=two_bit_ref_genome, chom_sizes=chom_sizes,
                                          sample_n_from_chroms_per_flen=sample_n_from_chroms)
    if gc_percs is None:
        return 1
    # write!
    print(f"writing reference fragment GC content distribution percentages ...")
    write_gc_percentages(gc_percentages=gc_percs, was_failure=failure)


if __name__ == '__main__':
    sys.exit(main())

# approx. computation time per sample: 30 min
