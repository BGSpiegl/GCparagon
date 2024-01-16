#!/usr/bin/env python3
import argparse
import sys
import math
import random
import numpy as np
from pathlib import Path
from math import log10
import multiprocessing as mp
from twobitreader import TwoBitFile
from typing import Union as OneOf, Optional, Tuple, Dict, List

# get root path
REPO_ROOT_DIR = Path(__file__).parent.parent.parent

# from previous script
SHIFT_N_TIMES = 4  # you might want to select a higher number

# table path definitions
# the following two paths are expected to be right -> change otherwise!
GENOME_BUILD = 'hg19'
default_2bit_reference_genome_path = (REPO_ROOT_DIR /
                                      f'src/GCparagon/2bit_reference/{GENOME_BUILD}.2bit')  # <---- required input!
default_predefined_genomic_regions = \
    REPO_ROOT_DIR / (f'accessory_files/genomic_interval_preselection-shifted{SHIFT_N_TIMES}x_{GENOME_BUILD}/'
                     f'{GENOME_BUILD}_minimalExclusionListOverlap_1Mbp_intervals_33pcOverlapLimited.bed')
# ^------- OUTPUT TABLE PATH!
default_output_directory = REPO_ROOT_DIR / 'accessory_files'
# below not required if
default_putative_ref_flength_dist_table = REPO_ROOT_DIR / 'accessory_files/plasmaSeq_ccfDNA_reference_fragment_length_distribution.tsv'
SAVE_COUNTS = True  # saves relative frequency otherwise

# sample definitions
# use_average_across_these_samples = ('B01', 'H01', 'P01')

# analysis setup
sample_n_fragments_per_mbp_default = 1 * 10**6
# OUTPUT PRECISION
DEFAULT_FLOAT_PRECISION_DIGITS = round(log10(sample_n_fragments_per_mbp_default)) + 3  # 9 digits after the comma
# ^--- the parameter above is only used if fractions are stored instead of counts (i.e., SAVE_COUNTS = False)
# Only has an impact if SAVE_COUNTS == False
default_parallel_processes = 3  # 24
exclude_n_bases_default = True
# random numbers:
RANDOM_SEED = random.randint(1, 999)
DEFAULT_FRAGMENT_N_CONTENT_THRESHOLD = 0.3  # from correct_GC_bias.py

# definitions of specific exceptions
class OutOfGenomicBoundsError(Exception):
    """
    To be raised in conditions where a genomic locus (= chromosome/scaffold/contig + coordinate) does not exist in a
    given reference genome build, whichever that may be.
    """
    pass


# from correct_GC_bias.py
def rhui(num: OneOf[int, float]):
    """

    :param num:
    :return:
    """
    return int(math.floor(num + 0.5))


# from correct_GC_bias.py
def safe_gc_base_count_inference_thresh(f_seq, f_len, threshold=0.3):
    """
    Execution takes around 30-60 micro seconds
    :param threshold:
    :param f_seq:
    :param f_len:
    :return:
    :raises: nothing but may lead to IndexError
    """
    n_count = f_seq.count('N')
    if n_count / f_len >= threshold:  # too many N-bases compared to fragment length
        return 99999999  # will lead to an IndexError
    try:
        return rhui(f_len / (f_len - n_count) * (f_seq.count('G') + f_seq.count('C')))
    except ZeroDivisionError:  # if all bases returned are 'N's
        return 99999999  # will lead to an IndexError


# from correct_GC_bias.py
def gc_count_rejecting_n_containing(f_seq):
    n_count = f_seq.count('N')
    if n_count:  # too many N-bases compared to fragment length
        return 99999999  # will lead to an IndexError
    return f_seq.count('G') + f_seq.count('C')


def load_table_with_flength_hdr(table_path: OneOf[str, Path]) -> Tuple[np.array, range]:
    table_path = Path(table_path)
    if not table_path.is_file():
        raise FileNotFoundError(f"the provided path does not exist: '{table_path}'")
    with (open(table_path, 'rt') as f_tbl):
        hdr_cont = f_tbl.readline().strip().split('\t')
        if hdr_cont == ['']:
            raise ValueError(f"Table '{table_path.name}' was empty!")
        return np.array(list(map(lambda s: float(s), f_tbl.readline().strip().split('\t')))), \
            range(int(hdr_cont[0]), int(hdr_cont[-1]))


def load_expected_fgcds_and_average(tables_dict: Dict[str, OneOf[str, Path]], output_dir: OneOf[str, Path],
                                    output_table: bool = True, genome_build_label: Optional[str] = None) \
        -> np.array:
    all_target_fgcds = []
    pctgs = None
    for sample_id, table_path in tables_dict.items():
        percentages = []
        with open(table_path, 'rt') as f_fgcd:
            hdr = f_fgcd.readline()  # header for optional output
            current_fgcd = []
            for line in f_fgcd.readlines():
                percentage, f_rel = line.strip().split('\t')
                percentages.append(percentage)
                current_fgcd.append(float(f_rel))
        # create numpy array
        current_fgcd = np.array(current_fgcd)
        try:
            assert 0.99 < current_fgcd.sum() < 1.01  # should sum up to 1
        except AssertionError:
            current_fgcd = current_fgcd / current_fgcd.sum()
        assert 0.99 < current_fgcd.sum() < 1.01  # should sum up to 1
        all_target_fgcds.append(current_fgcd)
        if pctgs is None:
            pctgs = percentages
    # compute average
    avg_fgcd = np.mean(all_target_fgcds, axis=0)
    if output_table:
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        reference_gc_content_distribution_path = Path(output_dir) / \
            f"{genome_build_label if genome_build_label else 'UNKNOWN'}_reference_GC_content_distribution.tsv"
        with open(reference_gc_content_distribution_path, 'wt') as f_fgcd_ideal:
            f_fgcd_ideal.writelines([hdr] + [f"{pc}\t{fr}\n" for pc, fr in zip(pctgs, avg_fgcd)])
    # return average otherwise
    return avg_fgcd


def get_ref_flength_rel_frequencies(corrected_fraglengths_table_path: Path, normalize_sample_counts: bool = True,
                                    output_table: bool = False, output_dir: Optional[OneOf[Path, str]] = None) \
        -> Tuple[np.array, range]:
    # read multi-sample FLD table
    with open(corrected_fraglengths_table_path, 'rt') as f_flengths:
        hdr_content = f_flengths.readline().strip().split('\t')
        sample_column = hdr_content.index('sample')
        status_column = hdr_content.index('status')
        sample_flengths = {}.fromkeys(use_average_across_these_samples)
        # find first column with number inside ASSERTS only number afterward
        is_int_callable = False
        first_integer_column_index = None
        for first_num_idx, cell in enumerate(hdr_content):
            try:
                _ = int(cell)
                if not is_int_callable:
                    first_integer_column_index = first_num_idx
                is_int_callable = True
            except ValueError:  # invalid literal for int() with base 10:
                if is_int_callable:
                    print(f"LOGIC ERROR - EXITING because the element '{cell}' was expected to be an integer.")
                    sys.exit(3)
                continue
        # insanity check
        if first_integer_column_index is None:
            raise ValueError(f"no integer callable cell found in: '" + ", ".join(hdr_content) + "'")
        # read distributions
        for line in f_flengths.readlines():
            line_content = line.strip().split('\t')
            if line_content[status_column].lower() != 'corrected':  # skip original, biased fragment length counts
                continue
            if line_content[sample_column] not in use_average_across_these_samples:  # use specified samples
                continue
            # update
            sample_flengths[line_content[sample_column]] = np.array(list(map(lambda s: float(s),
                                                                         line_content[first_integer_column_index:])))
    # output fragment length boundaries:
    min_flength = int(hdr_content[first_integer_column_index])
    max_flength = int(hdr_content[-1])
    print(f"fragment length counts found for fragments between {min_flength:,} bp and {max_flength:,} bp.")
    # combine counts
    if normalize_sample_counts:  # transform to relative frequencies
        for sample in sample_flengths.keys():
            sample_flengths[sample] /= sample_flengths[sample].sum()
    # compute the average across all relative frequencies
    all_dists = np.array(list(sample_flengths.values()))
    reference_flength_dist = np.sum(all_dists, axis=0) / \
        np.sum(all_dists, axis=(0, 1))  # divide by total sum across all (should be ~ num samples)
    if output_table:
        # output reference table
        output_dir.mkdir(exist_ok=True, parents=True)
        reference_flength_dist_path = output_dir / 'plasmaSeq_ccfDNA_reference_fragment_length_distribution.tsv'
        reference_flength_dist_path.parent.mkdir(parents=True, exist_ok=True)
        with open(reference_flength_dist_path, 'wt') as f_ref_fgcd:
            hdr = '\t'.join([str(length) for length in range(min_flength, max_flength + 1)]) + '\n'
            f_ref_fgcd.writelines([hdr, '\t'.join(map(lambda x: str(x), reference_flength_dist))])
    return reference_flength_dist, range(min_flength, max_flength)


def load_predefined_genomic_regions(genomic_regions_table: OneOf[str, Path]) -> List[Tuple[str, int, int, int]]:
    genomic_intervals = []
    with open(genomic_regions_table, 'rt') as f_reg:
        for line in f_reg.readlines():
            chrom, start, end, overlap, *_ = line.strip().split('\t')
            genomic_intervals.append((chrom, int(start), int(end), int(overlap)))
    return genomic_intervals


def simulate_fgcd_worker(genomic_intervals: List[Tuple[str, int, int, int]], two_bit_reference_path: OneOf[str, Path],
                         fragment_length_range: range, sample_n_per_mb: int, parent_connection: mp.Pipe,
                         fld: List[float], random_seed=RANDOM_SEED,
                         strict_n_ref_bases_handling=True,
                         frag_n_cont_thresh=DEFAULT_FRAGMENT_N_CONTENT_THRESHOLD):
    """
    Simulate the fragment GC content distribution (FGCD) based on a reference fragment length distribution (FLD) using a
    target number of simulated fragments.
    :param genomic_intervals: list of genomic intervals + their exclusion list overlap
    :param two_bit_reference_path: the path to a reference genome sequence file in 2bit format (NOT FASTA!)
    :param fragment_length_range: shortest and longest fragment length in bp; must match with received fld data!
    :param sample_n_per_mb: target number of simulated fragments per interval Mega base
    :param parent_connection: multiprocessing.Pipe connection; the parent process waits until this process has finished
                              and then receives all data put on the pipe by all child processes.
    :param fld: reference fragment length distribution which is used to draw fragments from the 2bit reference sequence
    :param random_seed: (might be useless)
    :param strict_n_ref_bases_handling: whether to reject all N-base(s) containing drawn fragment sequence;
                                        uses a fraction-of-bases threshold otherwise
    :param frag_n_cont_thresh: the threshold fraction of N bases for repeating fragment sequence drawing if
                               strict_n_ref_bases_handling is True. Not used otherwise.
    :return: nothing (sends tuple of: List[Tuple[str, int, int, int, List[float]]] (= list of good intervals) and
             List[Tuple[str, int, int, int]] (= bad intervals) back to parent process which contains the expected FGCD
             for all processed intervals.
    """
    ref_genome_handle = TwoBitFile(two_bit_reference_path)
    held_chromosome_hanldes = {}
    fgcd_results = []
    bad_intervals = []
    min_frag_len, max_frag_len = fragment_length_range.start, fragment_length_range.stop
    total_bad_intervals = 0
    for chrom, start, end, ovrlp in genomic_intervals:
        n_fragments = (end - start) / 10**6 * sample_n_per_mb
        try:
            chromosome_handle = held_chromosome_hanldes[chrom]
        except KeyError:
            chromosome_handle = ref_genome_handle[chrom]
            held_chromosome_hanldes[chrom] = chromosome_handle
        # check if interval coordinates are valid
        try:
            if start < 0 or end > ref_genome_handle.sequence_sizes()[chrom]:
                raise OutOfGenomicBoundsError(f"Check your coordinates! Your provided {chrom}:{start}-{end}")
        except KeyError:
            print(f"The following contig/scaffold was not found in the 2bit reference genome file: {chrom}")
            sys.exit(2)
        ref_seq_chunk_slice = chromosome_handle[start:end].upper()
        fragment_gc_content_observations = [0] * 101
        # CODE BELOW ADAPTED FROM MAIN CODE SIMULATION:
        total_fragments_simulated = 0
        exclude_bad_interval = False
        random_number_generator = np.random.default_rng(seed=random_seed)  # use random seed for reproducibility here!
        # for sim_iter_idx in range(simulation_repetitions): -> simulate only once!
        for fragment_length in range(min_frag_len, max_frag_len + 1):  # simulate each fragment length separately
            rel_freq = fld[fragment_length - min_frag_len]
            amount_fragments = round(n_fragments * rel_freq)
            if not amount_fragments:  # do not simulate fragments of this length -> rounded number was zero
                continue
            total_fragments_simulated += amount_fragments
            # create 2 random values for each fragment
            rand_ints = random_number_generator.integers(low=0,
                                                         high=end - start - fragment_length + 1,
                                                         size=amount_fragments)
            rand_ints.sort()  # for sweeping-like file pointer placement (if that even matters for twobitreader)
            sampling_failure_threshold = max(int(amount_fragments / 3), 55 if strict_n_ref_bases_handling else 33)
            # drawing fragment must fail at least 55 times (default is strict handling; 33 times otherwise) or
            # one third of all required draws, whichever is higher
            unsorted_randoms_index = 0  # index of new fall-back random positions for each fragment length
            unsorted_randoms = random_number_generator.integers(low=0, high=end - start - fragment_length + 1,
                                                                size=sampling_failure_threshold + 1)
            # could be rare fragments -> threshold +1 as minimum
            # 1/3 of target fragment number - whichever is higher - for chunk to be marked for exclusion.
            # Approach if not strict_n_ref_bases_handling: linear extrapolation
            # subtract the N bases from the fragment length, compute GC content of the analyzable portion and
            # extrapolate the number of GC-bases according to the fragment's actual length.
            # Still not perfect but the best we can do for a fast algorithm
            # (implemented in 'safe_gc_base_count_inference_thresh()')
            if strict_n_ref_bases_handling:  # faster
                gc_content_iterator = map(lambda q: rhui(100 * gc_count_rejecting_n_containing(f_seq=q) /
                                                         fragment_length),
                                          map(lambda s: ref_seq_chunk_slice[s:s + fragment_length].upper(),
                                              rand_ints))
            else:
                gc_content_iterator = map(lambda q: rhui(100 * safe_gc_base_count_inference_thresh(
                    f_seq=q[0], f_len=q[1], threshold=frag_n_cont_thresh) / fragment_length),
                                          map(lambda s: (ref_seq_chunk_slice[s:s + fragment_length].upper(),
                                                         fragment_length), rand_ints))
            for gc_content in gc_content_iterator:  # process all simulated fragments, consuming iterator from above
                try:
                    fragment_gc_content_observations[gc_content] += 1
                except IndexError:  # all-Ns-fragment: slow routine here -> get another random position
                    backup_seq = 'N' * fragment_length  # initialize for while loop
                    while backup_seq.count('N') if strict_n_ref_bases_handling else \
                            (backup_seq.count('N') / fragment_length >= frag_n_cont_thresh):  # keep drawing
                        try:
                            cur_start = unsorted_randoms[unsorted_randoms_index]  # get new random fragment start
                        except IndexError:  # all random integers consumed -> bad chunk! (should not occur)
                            print(f"Too many attempts of drawing random fragments ({sampling_failure_threshold} "
                                  f"attempts) for chunk '{chrom}:{start}-{end}' were in vain. Triggered for fragments "
                                  f"of {fragment_length}bp length. Marking interval as failed ..")
                            exclude_bad_interval = True  # stop redrawing once all fall-back random integers are used up
                            break
                        unsorted_randoms_index += 1  # can max. be 33, then clause below should trigger chunk marking
                        if unsorted_randoms_index >= sampling_failure_threshold:  # should always be the reason why a
                            # chunk is marked as bad
                            print(f"Too many attempts of drawing random fragments for chunk '{chrom}:{start}-{end}' "
                                  f"were in vain (threshold was {sampling_failure_threshold} tries). Marking interval "
                                  f"as failed ..")
                            exclude_bad_interval = True  # stop redrawing once the threshold has been reached
                            break
                        backup_seq = ref_seq_chunk_slice[cur_start:cur_start + fragment_length].upper()
                    if exclude_bad_interval:
                        break
                    if strict_n_ref_bases_handling:
                        gc_content = rhui(100 * gc_count_rejecting_n_containing(f_seq=backup_seq) / fragment_length)
                    else:
                        gc_content = rhui(100 * safe_gc_base_count_inference_thresh(f_seq=backup_seq,
                                                                                    f_len=fragment_length,
                                                                                    threshold=frag_n_cont_thresh) /
                                          fragment_length)
                    fragment_gc_content_observations[gc_content] += 1
            if exclude_bad_interval:
                break
        if exclude_bad_interval:
            bad_intervals.append((chrom, start, end, ovrlp))
            total_bad_intervals += 1
        else:
            fgcd_results.append((chrom, start, end, ovrlp, fragment_gc_content_observations))
    if total_bad_intervals:
        print(f"INFO - {total_bad_intervals:,} genomic interval(s) was/were marked for exclusion (i.e., not enough "
              f"fragments passing the N-bases threshold could've been drawn from the reference sequence).")
    else:
        print(f"INFO - no genomic intervals were marked for exclusion.")
    parent_connection.send((fgcd_results, bad_intervals))
    parent_connection.close()


def simulate_expected_fgcd_for_intervals(reference_fld: np.array, genomic_intervals: List[Tuple[str, int, int, int]],
                                         reference_genome: OneOf[str, Path], output_path: OneOf[str, Path],
                                         fragment_length_range: range, sample_n_fragments_per_mbp: int,
                                         processes: int = 8, strict_n_exclusion: bool = True, n_content_threshold=0.3,
                                         precision: int = 9, save_counts: bool = True):
    if reference_genome is None:
        raise ValueError("reference genome was not defined!")
    regions_per_process = len(genomic_intervals) // processes
    intervals_for_processes = []
    for list_idx in range(processes):
        if list_idx == processes - 1:  # last list consumes the rest
            intervals_for_processes.append(genomic_intervals[regions_per_process * list_idx:])  # until end
        else:
            intervals_for_processes.append(genomic_intervals[regions_per_process * list_idx:
                                                             regions_per_process * (list_idx + 1)])
    # create processes
    receivers = []
    worker_processes = []
    for proc_id in range(processes):
        receiver, sender = mp.Pipe(duplex=False)
        receivers.append(receiver)
        worker_processes.append(mp.Process(target=simulate_fgcd_worker,
                                           kwargs={'genomic_intervals': intervals_for_processes[proc_id],
                                                   'two_bit_reference_path': reference_genome,
                                                   'fragment_length_range': fragment_length_range,
                                                   'sample_n_per_mb': sample_n_fragments_per_mbp,
                                                   'parent_connection': sender,
                                                   'fld': reference_fld,
                                                   'random_seed': RANDOM_SEED,
                                                   'strict_n_ref_bases_handling': strict_n_exclusion,
                                                   'frag_n_cont_thresh': n_content_threshold}))
    # consolidate results form child processes
    for sim_proc in worker_processes:
        sim_proc.start()
    # wait until finished
    for sim_proc in worker_processes:
        sim_proc.join()
    # collect results
    interval_results = []
    there_were_bad_intervals = []
    for rcvr in receivers:
        good_intervals_results, bad_intervals_list = rcvr.recv()
        interval_results.extend(good_intervals_results)
        # ^----- a list of: (chrom, start, end, ovrlp, fragment_gc_content_observations)
        there_were_bad_intervals.extend(bad_intervals_list)  # a list of: (chrom, start, end, ovrlp)
    # close processes
    for sim_proc in worker_processes:
        sim_proc.close()
    # output FGCD results table
    output_table_path = (Path(output_path) /
                         f'{GENOME_BUILD}_minimalExclusionListOverlap_1Mbp_intervals_33pcOverlapLimited.FGCD.bed')
    # i:                            ^--- ADAPT NAME IF PARAMETERS WERE CHANGED! (interval size and/or overlap threshold)
    output_table_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_table_path, 'wt') as f_out:
        hdr = '\t'.join(['chromosome', 'start', 'end', 'exclusion_marked_bases'] +
                        [str(gc_pc) for gc_pc in range(0, 101, 1)]) + '\n'
        f_out.write(hdr)
        f_out.writelines(['\t'.join(list(map(
            lambda t: (str(t) if t != interval_result[-1] else  # unpack GC content counts list
                       ('\t'.join(list(map(lambda p: (str(p if save_counts else round(p/sum(t), ndigits=precision))),
                                           t)))
                        )), interval_result))) + '\n' for interval_result in interval_results])
    # output bad intervals exclusion list if -> use a sample yield of 0 such that the interval will always be rejected
    # independent of the sample yield
    if there_were_bad_intervals:
        bad_intervals_bed_path = output_table_path.parent / (f"{GENOME_BUILD}_genomic_intervals_failed_fragment_drawing"
                                                             f"{'-counts' if save_counts else '-f_rel'}.bed")
        bad_intervals_bed_path.parent.mkdir(parents=True, exist_ok=True)
        with open(bad_intervals_bed_path, 'wt') as f_bad:  # always re-write
            # hdr = '\t'.join(['chromosome', 'start', 'end', 'exclusion_marked_bases']) + '\n'
            f_bad.writelines(['\t'.join(list(map(lambda t: (str(t)
                                                            if t != bad_interval[-1] else
                                                            f'0,{t}'),  # add a zero for yield -> always reject interval
                                                 bad_interval))) + '\n'
                              for bad_interval in there_were_bad_intervals])


def get_cmd_args():
    commandline_parser = argparse.ArgumentParser()
    # TYPICAL ARGS:
    commandline_parser.add_argument('-gi', '--genomic-intervals-table', dest='predefined_genomic_intervals_table',
                                    help='Path to a table containing (equal sized) genomic intervals '
                                         'which are intended to be used for GC bias correction later. This table '
                                         "should follow the output format of GCparagon utility scripts "
                                         "'01-GI-preselection_test_Mbp_intervals_against_ExclusionList_hg19.py' and "
                                         "'02-GI-preselection_select_low_score_regions_from_overlapping_hg19.py'.",
                                    default=default_predefined_genomic_regions, metavar='File')
    # EITHER: SAMPLE-FragLenDist (avg. from multiple samples):
    commandline_parser.add_argument('-fldt', '--fragment-length-distributions-table',
                                    dest='samples_flengths_table_path',
                                    help='If the -rfldt is not available yet, this input is required. Contains '
                                         'fragment length distributions for multiple samples. A subset of samples '
                                         'included in the averaging process can be defined via --samples-subset. '
                                         'Initially, the fragment length distributions of the GCparagon validation '
                                         'samples were used. Use a --reference-fragment-length-distribution-table '
                                         'instead if you have one created already for your the specific combination of '
                                         'your experiment conditions.',
                                    metavar='File')
    # OR: REFERENCE FragLenDist
    # (reference distribution for plasmaSeq ccfDNA data from blood plasma samples already present):
    commandline_parser.add_argument('-rfldt', '--reference-fragment-length-distribution-table',
                                    dest='putative_ref_flength_dist_table',
                                    help='If the --fragment-length-distributions-table (above) has already been used '
                                         'to compute a --reference-fragment-length-distribution-table, the path to '
                                         'this --reference-fragment-length-distribution-table can be used to skip the '
                                         'averaging step. The TSV contains a theoretical fragment length distributions '
                                         'matching specific experiment conditions like: cfDNA isolated from a blood '
                                         'plasma sample; cfDNA isolated and sequencing library prepared following '
                                         'specific protocols; DNA sequenced with a specific technology, etc.',
                                    default=default_putative_ref_flength_dist_table, metavar='File')
    commandline_parser.add_argument('-o', '--output-dir', dest='output_directory', default=default_output_directory,
                                    help='Path to which the resulting table(s) will be written. If does not exist, '
                                         f'it will be created. [ DEFAULT: {default_output_directory}]',
                                    metavar='Directory')
    commandline_parser.add_argument('-ss', '--samples-subset', dest='sample_subset',
                                    help='This can be used to select a subset of samples from the '
                                         '--fragment-length-distributions-table. Any input will be reduced to '
                                         'data associated with these sample identifiers. Has no effect in the default '
                                         'case where the --reference-fragment-length-distribution-table for plasmaSeq '
                                         '(paired-end) sequenced ccfDNA from blood plasma samples is used.', nargs='+',
                                    metavar='List[String]')
    # OTHER (OPTIONAL) PATH DEFINITIONS
    commandline_parser.add_argument('-rtb', '--two-bit-reference-genome', dest='ref_genome_path',
                                    help='Path to 2bit version of the reference genome FastA file which was used for '
                                         'read alignment of the input BAM file. If only a FastA version is available, '
                                         "one can create the file using the following command: 'faToTwoBit "
                                         "<PATH_TO_REF_FASTA> -long <PATH_TO_OUT_2BIT>' "
                                         "(see genome.ucsc.edu/goldenPath/help/twoBit.html for more details)"
                                         f"[ DEFAULT: {default_2bit_reference_genome_path} ]",
                                    default=default_2bit_reference_genome_path, metavar='File')
    # FURTHER ANALYSIS PARAMETERS:
    commandline_parser.add_argument('-p', '--processes', dest='n_parallel', type=int,
                                    help='Number of parallel processes used for simultaneous simulation of '
                                         'fragment GC content distributions (FGCD) across multiple genomic intervals.',
                                    default=default_parallel_processes, metavar='Integer')
    commandline_parser.add_argument('-in', '--include-ns', dest='exclude_n_bases', action='store_false',
                                    help='Flag: if set, the presence of N bases in simulated fragments will be handled '
                                         'more lenient - all fragments can contain a fraction of Ns before these must '
                                         'be redrawn. This can also cause less genomic intervals to be rejected '
                                         f'(excluded from the output table. [ DEFAULT: {exclude_n_bases_default} ]',
                                    default=exclude_n_bases_default)
    commandline_parser.add_argument('-mnc', '--max-frag-n-content', dest='n_cont_max',
                                    help='Fragments with N bases fraction above this threshold will be redrawn. This '
                                         'might lead to genomic intervals to be marked for exclusion and not be '
                                         'contained in the output table. Not used if --include-ns flag is set. '
                                         f'[ DEFAULT: {DEFAULT_FRAGMENT_N_CONTENT_THRESHOLD}]',
                                    default=DEFAULT_FRAGMENT_N_CONTENT_THRESHOLD, metavar='Float')
    commandline_parser.add_argument('-fp', '--float-precision', dest='float_precision', type=int,
                                    help='Number of digits after the comma that will be stored for floating point '
                                         'numbers in the final output table. These floats are fragment fractions of a '
                                         'specific GC content. Only used if --save-fractions (instead of fragment '
                                         'counts) was set. Not used for any other output. The connection between '
                                         '--sample-n-fragments and this value should be as follows: '
                                         'round(log10(sample_n_fragments_per_mbp_default)) + 3  (default: 9 digits '
                                         'after the comma', default=DEFAULT_FLOAT_PRECISION_DIGITS, metavar='Integer')
    commandline_parser.add_argument('-snf', '--sample-n-fragments', dest='fragments_per_mbp_sampled', type=int,
                                    help='Number of fragments that will be sampled per megabase genomic interval '
                                         'following the fragment length distribution defined in -rfldt. '
                                         f'[ DEFAULT: {sample_n_fragments_per_mbp_default} ]',
                                    default=sample_n_fragments_per_mbp_default, metavar='Integer')
    commandline_parser.add_argument('-sf', '--save-fractions', dest='save_counts', action='store_false',
                                    help='Flag: if set, the simulated fragment GC contents will be stored as fractions '
                                         'for each fragment length rather than pure counts.')
    return commandline_parser.parse_args()


def gather_tsvs(tsvs_path: OneOf[str, Path], reduce_to_samples: Optional[OneOf[List, Tuple]] = None) \
        -> Dict[str, OneOf[str, Path]]:
    candidate_tsvs = set()
    candidate_tsvs.update(Path(tsvs_path).glob('*.tsv'))
    candidate_tsvs.update(Path(tsvs_path).glob('*.TSV'))
    # filter
    if reduce_to_samples:
        filtered_transformed_tsvs = list(filter(lambda x: x is not None,
                                                [(prt, ct) if prt in reduce_to_samples else None
                                                 for ct in candidate_tsvs
                                                 for prt in ct.stem.split('-')[0].split('_')[:3]]))
    else:
        try:
            filtered_transformed_tsvs = [(ct.stem.split('-')[0].split('_')[1], ct) for ct in candidate_tsvs]
        except IndexError:
            filtered_transformed_tsvs = [(ct.stem.split('-')[0].split('_')[0], ct) for ct in candidate_tsvs]
    tsvs_dict = {}
    for sample_id, tsv_path in filtered_transformed_tsvs:
        tsvs_dict[sample_id] = tsv_path
    return tsvs_dict


if __name__ == '__main__':
    # INPUT:
    cmd_args = get_cmd_args()
    predefined_genomic_regions_path = cmd_args.predefined_genomic_intervals_table
    parallel_processes = cmd_args.n_parallel
    exclude_n_bases = cmd_args.exclude_n_bases
    output_directory = cmd_args.output_directory
    output_path = Path(output_directory)
    reference_genome_path = cmd_args.ref_genome_path
    fragment_n_content_threshold = cmd_args.n_cont_max
    fragment_lengths_table_path = cmd_args.samples_flengths_table_path
    putative_ref_flength_dist_table = cmd_args.putative_ref_flength_dist_table
    use_average_across_these_samples = cmd_args.sample_subset
    fragments_per_mbp_sampled = cmd_args.fragments_per_mbp_sampled
    float_precision = cmd_args.float_precision
    save_counts = cmd_args.save_counts
    # check sanity of cmd line arguments
    if fragment_lengths_table_path is None and putative_ref_flength_dist_table is None:
        if default_putative_ref_flength_dist_table.is_file():
            putative_ref_flength_dist_table = default_putative_ref_flength_dist_table  # rescue
            print(f"WARNING - no reference fragment length distribution table was provided. ")
        else:
            raise ValueError("UNABLE TO PROCEED - path to either a table with fragment length distributions (FLDs) "
                             "nor a table to a consensus FLD was provided via the commandline.")
    # 1) read reference fragment length distribution (= FLD)
    if putative_ref_flength_dist_table is not None and Path(putative_ref_flength_dist_table).is_file():
        use_reference_fld, flength_range = load_table_with_flength_hdr(table_path=putative_ref_flength_dist_table)
    else:  # TODO: test this clause!
        use_reference_fld, flength_range = get_ref_flength_rel_frequencies(
            corrected_fraglengths_table_path=fragment_lengths_table_path,
            normalize_sample_counts=True, output_table=True, output_dir=output_directory)
    # 2) read predefined genomic regions
    predefined_genomic_regions = load_predefined_genomic_regions(
        genomic_regions_table=predefined_genomic_regions_path)
    print(f"INSANITY CHECK - read {len(predefined_genomic_regions):,} lines from BED file "
          f"'{predefined_genomic_regions_path}'")
    # 3) simulate fragment GC content distribution (= FGCD) in genomic intervals using multiprocessing
    simulate_expected_fgcd_for_intervals(genomic_intervals=predefined_genomic_regions,
                                         processes=parallel_processes,
                                         reference_genome=reference_genome_path,
                                         output_path=output_path,
                                         reference_fld=use_reference_fld,
                                         strict_n_exclusion=exclude_n_bases,
                                         fragment_length_range=flength_range,
                                         n_content_threshold=fragment_n_content_threshold,
                                         sample_n_fragments_per_mbp=fragments_per_mbp_sampled,
                                         precision=float_precision)
