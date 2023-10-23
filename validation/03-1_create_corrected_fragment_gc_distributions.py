#!/usr/bin/env python3
# tested with GCparagon Python 3.10 conda environment
import sys
import pysam
import logging
import numpy as np

from json import dumps
from math import floor
import multiprocessing as mp
from collections import deque
from natsort import humansorted
from twobitreader import TwoBitFile
from collections import defaultdict
from time import localtime, strftime
from pathlib import Path
from pysam import AlignmentFile
from typing import Union, Optional, Tuple

from GCparagon.correct_GC_bias import DEFAULT_MIN_UNCLIPPED_ALN_FRACTION  # should be 0.75
from GCparagon.utilities.gc_logging import gib_cmd_logger

SCRIPT_PARENT_PATH = Path(__file__).parent
SOURCE_CODE_ROOT_PATH = SCRIPT_PARENT_PATH.parent  # ../src/GCparagon
SOURCE_CODE_ROOT_DIR = str(SOURCE_CODE_ROOT_PATH)

if SOURCE_CODE_ROOT_PATH not in sys.path:
    sys.path.append(str(SOURCE_CODE_ROOT_PATH))

########################################################################################################################
# TODO: set your Griffin tagged BAM files parent directory here !!!!!!!!
GRIFFIN_TAGGED_BAM_FILES_DIR = Path('write_your_Griffin-tagged_BAM_files_parent_directory_path_here')
# TODO: set your Griffin tagged BAM files parent directory here !!!!!!!!
########################################################################################################################


min_unclipped_aln_fraction = DEFAULT_MIN_UNCLIPPED_ALN_FRACTION  # directly import from main module

hg38_2bit_ref_genome = SOURCE_CODE_ROOT_PATH / 'src/GCparagon/2bit_reference/hg38.analysisSet.2bit'
# above is valid if the EXECUTE_reference_download.sh is used; REPLACE THE PATH ABOVE OTHERWISE!

assert hg38_2bit_ref_genome.is_file()

MIN_FLENGTH = 20
MAX_FLENGTH = 800

# --------------- DEFINE AVAILABLE PRESETS, if GCparagon correction fidelity should be computed -----------
PRESETS = (1, 2)
# --------------------------------------------------------------------------------------------------------
output_path_gcparagon = SCRIPT_PARENT_PATH / ('03_genome-wide_correction_fidelity/'
                                              'original_and_corrected_GC_distribution_GCparagon'
                                              '-preset{}')  # FORMAT WITH PRESET USED!
output_path_griffin = SCRIPT_PARENT_PATH / ('03_genome-wide_correction_fidelity/'
                                            'original_and_corrected_GC_distribution_Griffin')
output_paths = {'GCparagon': output_path_gcparagon,
                'Griffin': output_path_griffin}  # FORMAT GCparagon WITH PRESET USED!
# Griffin correction weights-tagged BAM files
tagged_bams_griffin = (GRIFFIN_TAGGED_BAM_FILES_DIR / 'C01/C01.GCtagged.bam',
                       GRIFFIN_TAGGED_BAM_FILES_DIR / 'H01/H01.GCtagged.bam',
                       GRIFFIN_TAGGED_BAM_FILES_DIR / 'B01/B01.GCtagged.bam',
                       GRIFFIN_TAGGED_BAM_FILES_DIR / 'P01/P01.GCtagged.bam')

# GCparagon
# -------------------------------------------------------------------------
input_path_gcparagon = SOURCE_CODE_ROOT_PATH / 'preset_computation/preset{}'  # FORMAT WITH PRESET USED!
tagged_bams_gcparagon = (input_path_gcparagon / 'C01/C01.GCtagged.bam',
                         input_path_gcparagon / 'H01/H01.GCtagged.bam',
                         input_path_gcparagon / 'B01/B01.GCtagged.bam',
                         input_path_gcparagon / 'P01/P01.GCtagged.bam',)
tagged_bam_lists = {'Griffin': tagged_bams_griffin,
                    'GCparagon': tagged_bams_gcparagon}
correction_weight_tags = {'Griffin': 'GG',
                          'GCparagon': 'GC'}  # HIGHLY IMPORTANT !!! use 'GC' for GCparagon and 'GG' for Griffin


def rhui(num: Union[int, float]):
    return int(floor(num + 0.5))


def get_standard_chromosomes_from_bam(bam_path: str, sort_by_length=False,
                                      remove_scaffolds=None) -> Tuple[list, dict]:
    standard_chromosome_contigs = []
    standard_chromosome_lengths = []
    remove_scaffolds_lowercase = None
    if remove_scaffolds is not None:
        remove_scaffolds_lowercase = [scaf.lower() for scaf in  remove_scaffolds]
    with AlignmentFile(bam_path, "rb") as aln_file:
        reference_names = aln_file.references
        reference_lengths = aln_file.lengths
        for length_index, contig_name in enumerate(reference_names):
            if '_' in contig_name or 'EBV' in contig_name:
                continue
            if remove_scaffolds_lowercase is not None and contig_name.lower() in remove_scaffolds_lowercase:
                continue
            standard_chromosome_contigs.append(contig_name)
            standard_chromosome_lengths.append(reference_lengths[length_index])
    # sort chromosome names based on length in descending order!
    standard_chromosomes = humansorted(list(zip(standard_chromosome_contigs, standard_chromosome_lengths)),
                                       key=lambda d: d[int(sort_by_length)], reverse=False)
    # create genome file for blacklist processing (blacklisted region enlargement for fetch region exclusion)
    # create name to length mapping dictionary
    standard_chromosomes_lengths = {}.fromkeys(standard_chromosome_contigs)
    for chr_n, chr_l in standard_chromosomes:
        standard_chromosomes_lengths[chr_n] = chr_l
    standard_chromosomes = [chrm for chrm, ln in standard_chromosomes]
    return standard_chromosomes, standard_chromosomes_lengths


LOGGER = gib_cmd_logger()  # formatted cmdline output


def log(message: str, log_level: int, logger_name: str, flush=True, close_handlers=False):
    """

    :param message:
    :param log_level:
    :param logger_name:
    :param flush:
    :param close_handlers:
    :return:
    """
    current_logger = logging.getLogger(logger_name)
    match log_level:
        case logging.NOTSET:
            current_logger.info(message)
        case logging.DEBUG:
            current_logger.debug(message)
        case logging.INFO:
            current_logger.info(message)
        case logging.WARNING:
            current_logger.warning(message)
        case logging.ERROR:
            current_logger.error(message)
        case logging.CRITICAL:
            current_logger.critical(message)
    if flush and isinstance(current_logger, logging.Logger):  # do nothing if logger_name is undefined
        for hdlr in current_logger.handlers:
            hdlr.flush()
    if close_handlers:
        for hdlr in current_logger.handlers:
            hdlr.flush()
            hdlr.close()


def scaffold_gc_stats_counter(bam_path: Union[Path, str], tag: str,
                              sender_connection: mp.Pipe,
                              scaffold_to_process: str, sample_data_id: str, max_no_tags=100000):  # ASSERTS pe-reads!
    ref_genome_handle = TwoBitFile(hg38_2bit_ref_genome)
    gc_counter_pre = defaultdict(int)  # max. 101 entries
    fragment_lengths_pre = defaultdict(int)
    fragment_lengths_post = defaultdict(int)
    gc_counter_corr = defaultdict(float)  # sum of GC-tags
    none_type_bases = 0
    tag_not_present = 0
    min_unclipped_aln_fraction = DEFAULT_MIN_UNCLIPPED_ALN_FRACTION
    with pysam.AlignmentFile(bam_path, mode='rb') as f_aln:  # process all alignments passing filters
        # raises ValueError if scaffold_to_process is not in BAM header sequences
        this_ref_length = f_aln.lengths[f_aln.references.index(scaffold_to_process)]
        # binary filter:
        exclude_flags = np.uint32(3852)  # = 256 + 2048 + 512 + 1024 + 4 + 8
        exclude_flags_binary = bin(exclude_flags)
        paired_flag = np.uint32(1)  # only paired!
        paired_flag_binary = bin(paired_flag)
        # -> not necessary to check for "mapped" attribute (filtered out if "read unmapped")

        # complete alignment filter:
        # --------------------------
        # EXCLUDE if:
        # read unmapped = 4                        '0b100'
        # mate unmapped = 8                       '0b1000'
        # not primary = 256                  '0b100000000'
        # vendor/QC fail = 512              '0b1000000000'
        # PCR or optical duplicate = 1024  '0b10000000000'
        # supplementary = 2048            '0b100000000000'
        # = 3852                          '0b111100001100'
        # REQUIRE THAT:
        # (alignment is mapped = 1 '0b1'; requirement removed because covered by "not read unmapped")
        # mates map to different strands
        #    a.is_forward != a.mate_is_forward
        # TLEN column is (positive and) between defined fragment length limits (inclusive)
        #    min_frag_len <= a.template_length <= max_frag_len
        filtered_alignments = filter(lambda a:
                                     bin(np.uint32(a.flag) & paired_flag) == paired_flag_binary and
                                     bin(~np.uint32(a.flag) & exclude_flags) == exclude_flags_binary and
                                     a.is_forward != a.mate_is_forward and
                                     (MIN_FLENGTH <= a.template_length <= MAX_FLENGTH) and
                                     #  ^--- choose theoretically feasible fragment lengths
                                     a.reference_length >= a.query_length * min_unclipped_aln_fraction,
                                     f_aln.fetch(scaffold_to_process, 0, this_ref_length, multiple_iterators=True))
        chrom_handles = {}
        held_chrom_handles = deque()
        for aln in filtered_alignments:
            # if aln.reference_length < int(aln.query_length * min_unclipped_aln_fraction):
            #     continue  # extensively (soft-)clipped read alignment -> ignore
            try:  # to get chromosome handle
                cur_chrom_handle = chrom_handles[aln.reference_name]
            except KeyError:
                chrom_handles[aln.reference_name] = ref_genome_handle[aln.reference_name]
                held_chrom_handles.append(aln.reference_name)
                cur_chrom_handle = chrom_handles[aln.reference_name]
                if len(held_chrom_handles) == 3:
                    del chrom_handles[held_chrom_handles.popleft()]
            # get fragment sequence (estimate)
            try:
                if aln.template_length <= aln.query_alignment_length:  # template_length filtered to be >0; max. 25% clipped
                    bases = aln.query_alignment_sequence.upper()
                    #         ^-- This is a substring of query_sequence, excluding flanking bases that were soft clipped
                else:
                    start_pos = min(aln.reference_start, aln.next_reference_start)
                    end_pos = start_pos + aln.template_length  # is filtered to be positive
                    bases = cur_chrom_handle[start_pos:end_pos].upper()  # coords 0-based, end-open
                gc_pc = rhui(100. * (bases.count('G') + bases.count('C')) / len(bases))  # returns mathematically
                # correct int bin
                gc_counter_pre[gc_pc] += 1
                gc_counter_corr[gc_pc] += aln.get_tag(tag)  # returned value is cast into an appropriate python type
                fragment_lengths_pre[aln.template_length] += 1
                fragment_lengths_post[aln.template_length] += aln.get_tag(tag)
            except AttributeError:  # 'NoneType' object has no attribute 'upper'
                none_type_bases += 1  # no sequence present -> faulty entry?
            except KeyError:
                tag_not_present += 1
                if tag_not_present >= max_no_tags:
                    log(message=f"too many (>={max_no_tags:,}) alignments encountered without identifiable "
                                "GC-tags. Terminating..", log_level=logging.CRITICAL, logger_name=LOGGER)
                    sender_connection.send((False, None, None, None, None, None))
                    sender_connection.close()
                    return
    sender_connection.send((True,  # success
                            fragment_lengths_pre, fragment_lengths_post, gc_counter_pre, gc_counter_corr,  # stats
                            none_type_bases, tag_not_present,  # fragment skipping reasons/debug
                            sample_data_id))  # identification
    sender_connection.close()


def infer_gc_tag(preferred: str, bam_path: Union[Path, str]):
    options = ('YC', 'GG', )
    tag_counts = defaultdict(int)
    infer_from = 1000  # first
    aln_counter = 0
    with pysam.AlignmentFile(bam_path, mode='rb') as f_aln:
        for aln in f_aln.fetch(multiple_iterators=True, until_eof=True):
            aln_counter += 1
            if aln_counter >= infer_from:
                break
            if aln.has_tag(preferred):
                tag_counts[preferred] += 1
                continue
            for opt in options:  # check options otherwise
                if aln.has_tag(opt):
                    tag_counts[opt] += 1
                    continue
            tag_counts['unknown'] += 1
    return sorted(list(tag_counts.items()), key=lambda x: x[1], reverse=True)[0][0]  # returns tag with highest count


def main() -> Optional[int]:
    preset_tables = {}
    for CORRECTION_ALGORITHM in ('Griffin', 'GCparagon'):
        if CORRECTION_ALGORITHM == 'GCparagon':
            preset_tables[CORRECTION_ALGORITHM] = {}.fromkeys(PRESETS)
        else:
            preset_tables[CORRECTION_ALGORITHM] = {}
        print(f"INFO - using analysis definitions for GC bias correction algorithm '{CORRECTION_ALGORITHM}'")
        # select correct definition set
        USE_TAG = correction_weight_tags[CORRECTION_ALGORITHM]
        output_path = output_paths[CORRECTION_ALGORITHM]
        tagged_bams = tagged_bam_lists[CORRECTION_ALGORITHM]
        itr_over = PRESETS if CORRECTION_ALGORITHM == 'GCparagon' else (None,)
        # -------------------------------------------------------------------------
        timestamp_str = strftime('%d-%m-%Y', localtime())
        for PRESET in itr_over:
            # iterate over tuples
            sample_lines_dict = {}
            sample_flength_lines_dict = {}
            percent_part = '\t'.join([str(gc_pc) for gc_pc in range(0, 101, 1)])
            hdr_line = f'sample\talgorithm\tstatus\t{percent_part}\n'
            length_part = '\t'.join([str(length) for length in range(MIN_FLENGTH, MAX_FLENGTH + 1, 1)])
            hdr_line_lengths = f'sample\talgorithm\tstatus\t{length_part}\n'
            formatted_output_path = Path(str(output_path).format(PRESET))  # has no effect if no '{}' present in str
            formatted_output_path.mkdir(parents=True, exist_ok=True)
            output_table_path = \
                Path(formatted_output_path /
                     (f"{CORRECTION_ALGORITHM.upper()}-STATISTICS_allSamples_GC-percentage_perFragmentSequence"
                      f"_{timestamp_str}.tsv"))
            if PRESET is None and CORRECTION_ALGORITHM == 'Griffin':  # 'Griffin'
                preset_tables[CORRECTION_ALGORITHM] = output_table_path
            else:
                if preset_tables[CORRECTION_ALGORITHM] is None:
                    preset_tables[CORRECTION_ALGORITHM] = {PRESET: output_table_path}
                else:
                    preset_tables[CORRECTION_ALGORITHM][PRESET] = output_table_path
            if output_table_path.is_file():
                output_table_path.unlink()  # delete if exists; will be opened as append otherwise!
            output_length_table_path = \
                Path(formatted_output_path /
                     f"{CORRECTION_ALGORITHM.upper()}-STATISTICS_allSamples_fragment_length_{timestamp_str}.tsv")
            if output_length_table_path.is_file():
                output_length_table_path.unlink()  # delete if exists; will be opened as append otherwise
            for bam_idx, bam_path in enumerate(tagged_bams):
                formatted_bam_path = Path(str(bam_path).format(PRESET))
                print(f"INFO - processing BAM file '{formatted_bam_path.name}'")
                std_chroms, std_chrom_lengths = get_standard_chromosomes_from_bam(
                    bam_path=str(formatted_bam_path), remove_scaffolds=['chrY', 'chrM', 'chrMT', 'chrEBV'])
                fragment_workers = []
                receivers = []
                # GC tag for GCparagon tagged files; additional GG for Griffin results-tagged file
                inferred_gc_tag = infer_gc_tag(preferred=USE_TAG, bam_path=formatted_bam_path)
                if inferred_gc_tag == 'unknown':
                    log(message=f"could not infer GC-tag. Code will likely fail with error.", logger_name=LOGGER,
                        log_level=logging.WARNING)
                else:
                    log(message=f"inferred GC-tag: '{inferred_gc_tag}'", logger_name=LOGGER, log_level=logging.INFO)
                sample_id_path = formatted_bam_path.parent
                sample_id = str(sample_id_path.stem.split('.')[0])  # just 'P01' etc. not including '.GCtagged'
                sample_data_id = f'{sample_id}_{CORRECTION_ALGORITHM}'
                sample_lines_dict[sample_data_id] = {'original': f'{sample_id}\t{CORRECTION_ALGORITHM}\toriginal\t',
                                                     'corrected': f'{sample_id}\t{CORRECTION_ALGORITHM}\tcorrected\t'}
                sample_flength_lines_dict[sample_data_id] = {
                    'original': f'{sample_id}\t{CORRECTION_ALGORITHM}\toriginal\t',
                    'corrected': f'{sample_id}\t{CORRECTION_ALGORITHM}\tcorrected\t'}
                # create GC stats using multiprocessing (for first in pair)
                fragment_stats = {'original': [], 'corrected': []}
                fragment_lengths = {'original': [], 'corrected': []}
                for scaffold_to_process in std_chroms:
                    receiver, sender = mp.Pipe(duplex=False)
                    receivers.append(receiver)
                    fragment_workers.append(mp.Process(target=scaffold_gc_stats_counter,
                                                       kwargs={'bam_path': formatted_bam_path,
                                                               'tag': inferred_gc_tag,
                                                               'sender_connection': sender,
                                                               'scaffold_to_process': scaffold_to_process,
                                                               'sample_data_id': sample_data_id}))
                for read_worker in fragment_workers:
                    read_worker.start()
                for rcv in receivers:
                    # vars/vals used in worker: True, gc_counter_pre, gc_counter_corr, none_type_bases, tag_not_present
                    success, flength_dist_pre, flength_dist_post, stats_pre, stats_corr, none_bases, no_tag, \
                        sample_data_id = rcv.recv()
                    # report progress
                    log(message=f"received stats for sample data '{sample_data_id}'",
                        log_level=logging.INFO, logger_name=LOGGER, flush=True)
                    if not success:
                        log(message=f'NOT ENOUGH TAGGED ALNS. TERMINATING..', log_level=logging.CRITICAL,
                            logger_name=LOGGER, close_handlers=True, flush=True)
                        for read_worker in fragment_workers:
                            read_worker.kill()
                        return 1
                    if none_bases:
                        log(message=f"WARNING: unexpected behavior - there were {none_bases:,} None-base alignments in "
                                    f"sample data '{sample_data_id}! NOT EXPECTED",
                            log_level=logging.WARNING, logger_name=LOGGER, flush=True)
                    if no_tag:
                        log(message=f"WARNING: unexpected behavior - there were {no_tag:,} alignments with no "
                                    f"identifiable GC-Tags (GC or YC) in file of sample '{sample_data_id}!",
                            log_level=logging.WARNING, logger_name=LOGGER, flush=True)
                    # record output statistics
                    fragment_stats['original'].append(stats_pre)
                    fragment_stats['corrected'].append(stats_corr)
                    fragment_lengths['original'].append(flength_dist_pre)
                    fragment_lengths['corrected'].append(flength_dist_post)
                # accumulate stats
                print(f"INFO - accumulating statistics for sample '{sample_id}'")
                for cur_status, cur_frag_stats_list in fragment_stats.items():
                    accumulated_stats = defaultdict(float)
                    for frag_stats in cur_frag_stats_list:
                        for gc_key, frag_count in frag_stats.items():
                            accumulated_stats[gc_key] += frag_count
                    sample_lines_dict[sample_data_id][cur_status] += '\t'.join([str(int(accumulated_stats[gc_pc]))
                                                                               for gc_pc in range(0, 101, 1)])
                    sample_lines_dict[sample_data_id][cur_status] += '\n'
                # output accumulated statistics (append if existing at this point))
                with open(output_table_path, 'at') as f_tab:
                    if bam_idx == 0:
                        f_tab.write(hdr_line)
                        f_tab.flush()
                    f_tab.writelines([sample_lines_dict[sample_data_id][cur_status]
                                      for cur_status in ('original', 'corrected')])
                    f_tab.flush()  # continuously write stats per sample to file
                # print result also to cmd
                print(f"GC content data for sample '{sample_id}' written to stats TSV file "
                      f"'{output_table_path.name}':\n" + '----' * 10 + '\n' +
                      dumps(sample_lines_dict, sort_keys=True, indent=3))
                # accumulate fragment length counts
                sample_flength_lines_dict[sample_data_id] = {
                    'original':  f'{sample_id}\t{CORRECTION_ALGORITHM}\toriginal\t',
                    'corrected': f'{sample_id}\t{CORRECTION_ALGORITHM}\tcorrected\t'}
                for cur_status, cur_frag_length_list in fragment_lengths.items():
                    accumulated_lengths = defaultdict(float)
                    for frag_lengths in cur_frag_length_list:
                        for len_key, frag_count in frag_lengths.items():
                            accumulated_lengths[len_key] += frag_count
                    sample_flength_lines_dict[sample_data_id][cur_status] += '\t'.join(
                        [str(int(accumulated_lengths[length])) for length in range(MIN_FLENGTH, MAX_FLENGTH + 1, 1)])
                    sample_flength_lines_dict[sample_data_id][cur_status] += '\n'  # each sample has 4 lines
                # output accumulated statistics (append if existing at this point)
                with open(output_length_table_path, 'at') as f_tab:
                    if bam_idx == 0:
                        f_tab.write(hdr_line_lengths)
                        f_tab.flush()
                    f_tab.writelines([sample_flength_lines_dict[sample_data_id][cur_status]
                                      for cur_status in ('original', 'corrected')])
                    f_tab.flush()  # continuously write stats per sample to file
                # print result also to cmd
                print(f"Fragment length data for sample '{sample_id}' written to stats TSV file "
                      f"'{output_length_table_path.name}':\n" + '----' * 10 + '\n' + dumps(sample_flength_lines_dict,
                                                                                           sort_keys=True, indent=3))
                for worker in fragment_workers:
                    worker.join()
                for worker in fragment_workers:
                    worker.close()
    # combine statistics tables - per GCparagon preset
    for preset in PRESETS:
        preset_output_table = SCRIPT_PARENT_PATH / \
                              '03_genome-wide_correction_fidelity/ALL-STATISTICS_allSamples_GC-percentage_' \
                              f'perFragmentSequence_{timestamp_str}_preset{preset}.tsv'
        preset_output_table.parent.mkdir(parents=True, exist_ok=True)
        with open(preset_output_table, 'wt') as f_out:
            with open(preset_tables['Griffin'], 'rt') as f_griff:
                f_out.write(f_griff.read())  # simply copy
            # add content form GCparagon preset result
            with open(preset_tables['Griffin'][preset], 'rt') as f_gcp:
                _ = f_gcp.readline()  # skip header
                f_out.writelines(f_gcp.readlines())


if __name__ == '__main__':
    sys.exit(main())


# NEW OUTPUT (no proper-paired problem - concordant fragment length distributions; Multiprocessing problem eradicated;
#             Tagging shift bug solved)

# /home/benjamin/mambaforge-pypy3/envs/GCparagon/bin/python /mnt/NVMeScratch/PycharmProjects/GCparagon_dev/revision/NARGAB/Griffin_correction_validation/genomewide_GC_content_per_fragment_corrected/03_create_corrected_fragment_gc_distributions.py
# INFO - using analysis definitions for GC bias correction algorithm 'GCparagon'
# INFO - processing BAM file 'C01.GCtagged.bam'
# | 2023-08-14 15:23:54,048 - GCPARAGON_validation - INFO - inferred GC-tag: 'GC'
# | 2023-08-14 15:39:44,449 - GCPARAGON_validation - INFO - received stats for sample data 'C01_GCparagon'
# | 2023-08-14 15:39:44,450 - GCPARAGON_validation - INFO - received stats for sample data 'C01_GCparagon'
# | 2023-08-14 15:39:44,450 - GCPARAGON_validation - INFO - received stats for sample data 'C01_GCparagon'
# | 2023-08-14 15:39:44,450 - GCPARAGON_validation - INFO - received stats for sample data 'C01_GCparagon'
# | 2023-08-14 15:39:44,450 - GCPARAGON_validation - INFO - received stats for sample data 'C01_GCparagon'
# | 2023-08-14 15:39:44,451 - GCPARAGON_validation - INFO - received stats for sample data 'C01_GCparagon'
# | 2023-08-14 15:39:44,451 - GCPARAGON_validation - INFO - received stats for sample data 'C01_GCparagon'
# | 2023-08-14 15:39:44,451 - GCPARAGON_validation - INFO - received stats for sample data 'C01_GCparagon'
# | 2023-08-14 15:39:44,451 - GCPARAGON_validation - INFO - received stats for sample data 'C01_GCparagon'
# | 2023-08-14 15:39:44,452 - GCPARAGON_validation - INFO - received stats for sample data 'C01_GCparagon'
# | 2023-08-14 15:39:44,452 - GCPARAGON_validation - INFO - received stats for sample data 'C01_GCparagon'
# | 2023-08-14 15:39:44,452 - GCPARAGON_validation - INFO - received stats for sample data 'C01_GCparagon'
# | 2023-08-14 15:39:44,452 - GCPARAGON_validation - INFO - received stats for sample data 'C01_GCparagon'
# | 2023-08-14 15:39:44,452 - GCPARAGON_validation - INFO - received stats for sample data 'C01_GCparagon'
# | 2023-08-14 15:39:44,452 - GCPARAGON_validation - INFO - received stats for sample data 'C01_GCparagon'
# | 2023-08-14 15:39:44,453 - GCPARAGON_validation - INFO - received stats for sample data 'C01_GCparagon'
# | 2023-08-14 15:39:44,453 - GCPARAGON_validation - INFO - received stats for sample data 'C01_GCparagon'
# | 2023-08-14 15:39:44,453 - GCPARAGON_validation - INFO - received stats for sample data 'C01_GCparagon'
# | 2023-08-14 15:39:44,453 - GCPARAGON_validation - INFO - received stats for sample data 'C01_GCparagon'
# | 2023-08-14 15:39:44,453 - GCPARAGON_validation - INFO - received stats for sample data 'C01_GCparagon'
# | 2023-08-14 15:39:44,454 - GCPARAGON_validation - INFO - received stats for sample data 'C01_GCparagon'
# | 2023-08-14 15:39:44,454 - GCPARAGON_validation - INFO - received stats for sample data 'C01_GCparagon'
# | 2023-08-14 15:39:44,454 - GCPARAGON_validation - INFO - received stats for sample data 'C01_GCparagon'
# INFO - accumulating statistics for sample 'C01'
# GC content data for sample 'C01' written to stats TSV file 'GCPARAGON-STATISTICS_allSamples_GC-percentage_perFragmentSequence_14-08-2023.tsv':
# ----------------------------------------
# {
#    "C01_GCparagon": {
#       "original": "C01\tGCparagon\toriginal\t13\t9\t4\t7\t4\t5\t6\t6\t13\t22\t26\t95\t195\t465\t1154\t2822\t6825\t15519\t34049\t67800\t130776\t226793\t378735\t609982\t877332\t1318193\t1771959\t2344063\t3139646\t3608668\t4521353\t5173885\t5952248\t6459240\t7990153\t8442539\t9499475\t9717954\t11015514\t10890866\t11046625\t10317936\t9892251\t9496761\t9129367\t9026230\t9144000\t8700617\t8948395\t9255960\t8228757\t9203941\t8478746\t7824526\t7055110\t5765055\t4731788\t3907543\t3376077\t2811085\t2516859\t2133940\t1787700\t1674723\t1416879\t1184850\t1039639\t752891\t650764\t502048\t404336\t304800\t255400\t200545\t159395\t135705\t109548\t98308\t88101\t74819\t63665\t45409\t29651\t16678\t7768\t3446\t1451\t551\t193\t85\t29\t21\t8\t4\t2\t4\t9\t46\t38\t10\t15\n"
#    }
# }
# GC content data for sample 'C01' written to stats TSV file 'GCPARAGON-STATISTICS_allSamples_fragment_length_14-08-2023.tsv':
# ----------------------------------------
# {
#    "C01_GCparagon": {
#       "original": "C01\tGCparagon\toriginal\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t1867\t2144\t2365\t48894\t50509\t50396\t50156\t49329\t51655\t54412\t58807\t63948\t71577\t79086\t79205\t75152\t77848\t78415\t84340\t90678\t99110\t111990\t131494\t160936\t185797\t185997\t164455\t152284\t159417\t176547\t200157\t232611\t280816\t357184\t467278\t576525\t630351\t603590\t594876\t659478\t787943\t977110\t1191522\t1389891\t1517521\t1602981\t1699871\t1801748\t1882318\t1978064\t2125188\t2309836\t2532112\t2801255\t3046383\t3181286\t3264601\t3461698\t3795099\t4254543\t4793676\t5335008\t5810205\t6275901\t6859630\t7724267\t8838348\t9720798\t9968505\t9789621\t9390739\t8909936\t8332369\t7738799\t7213352\t6856205\t6558623\t6269774\t5924184\t5497898\t5013880\t4531909\t4082831\t3676280\t3343139\t3061883\t2847661\t2686066\t2542297\t2396912\t2254845\t2093050\t1924053\t1754485\t1599824\t1465344\t1357468\t1269362\t1197725\t1139952\t1082182\t1024027\t966088\t898305\t828877\t762927\t705538\t654258\t612894\t573948\t543428\t513129\t481969\t449685\t415853\t380974\t347621\t318100\t294986\t274002\t256666\t241832\t226520\t211748\t195859\t179225\t165505\t150051\t138676\t130145\t122849\t116830\t112251\t105299\t98679\t92397\t86950\t82628\t78690\t74254\t70957\t70313\t86183\t62546\t3533\t3530\t4155\t3004\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n"
#    }
# }
# INFO - processing BAM file 'H01.GCtagged.bam'
# | 2023-08-14 15:39:44,549 - GCPARAGON_validation - INFO - inferred GC-tag: 'GC'
# | 2023-08-14 16:03:52,219 - GCPARAGON_validation - INFO - received stats for sample data 'H01_GCparagon'
# | 2023-08-14 16:03:52,220 - GCPARAGON_validation - INFO - received stats for sample data 'H01_GCparagon'
# | 2023-08-14 16:03:52,221 - GCPARAGON_validation - INFO - received stats for sample data 'H01_GCparagon'
# | 2023-08-14 16:03:52,221 - GCPARAGON_validation - INFO - received stats for sample data 'H01_GCparagon'
# | 2023-08-14 16:03:52,222 - GCPARAGON_validation - INFO - received stats for sample data 'H01_GCparagon'
# | 2023-08-14 16:03:52,223 - GCPARAGON_validation - INFO - received stats for sample data 'H01_GCparagon'
# | 2023-08-14 16:03:52,223 - GCPARAGON_validation - INFO - received stats for sample data 'H01_GCparagon'
# | 2023-08-14 16:03:52,224 - GCPARAGON_validation - INFO - received stats for sample data 'H01_GCparagon'
# | 2023-08-14 16:03:52,225 - GCPARAGON_validation - INFO - received stats for sample data 'H01_GCparagon'
# | 2023-08-14 16:03:52,225 - GCPARAGON_validation - INFO - received stats for sample data 'H01_GCparagon'
# | 2023-08-14 16:03:52,226 - GCPARAGON_validation - INFO - received stats for sample data 'H01_GCparagon'
# | 2023-08-14 16:03:52,226 - GCPARAGON_validation - INFO - received stats for sample data 'H01_GCparagon'
# | 2023-08-14 16:03:52,227 - GCPARAGON_validation - INFO - received stats for sample data 'H01_GCparagon'
# | 2023-08-14 16:03:52,227 - GCPARAGON_validation - INFO - received stats for sample data 'H01_GCparagon'
# | 2023-08-14 16:03:52,228 - GCPARAGON_validation - INFO - received stats for sample data 'H01_GCparagon'
# | 2023-08-14 16:03:52,229 - GCPARAGON_validation - INFO - received stats for sample data 'H01_GCparagon'
# | 2023-08-14 16:03:52,229 - GCPARAGON_validation - INFO - received stats for sample data 'H01_GCparagon'
# | 2023-08-14 16:03:52,230 - GCPARAGON_validation - INFO - received stats for sample data 'H01_GCparagon'
# | 2023-08-14 16:03:52,230 - GCPARAGON_validation - INFO - received stats for sample data 'H01_GCparagon'
# | 2023-08-14 16:03:52,231 - GCPARAGON_validation - INFO - received stats for sample data 'H01_GCparagon'
# | 2023-08-14 16:03:52,231 - GCPARAGON_validation - INFO - received stats for sample data 'H01_GCparagon'
# | 2023-08-14 16:03:52,232 - GCPARAGON_validation - INFO - received stats for sample data 'H01_GCparagon'
# | 2023-08-14 16:03:52,233 - GCPARAGON_validation - INFO - received stats for sample data 'H01_GCparagon'
# INFO - accumulating statistics for sample 'H01'
# GC content data for sample 'H01' written to stats TSV file 'GCPARAGON-STATISTICS_allSamples_GC-percentage_perFragmentSequence_14-08-2023.tsv':
# ----------------------------------------
# {
#    "C01_GCparagon": {
#       "original": "C01\tGCparagon\toriginal\t13\t9\t4\t7\t4\t5\t6\t6\t13\t22\t26\t95\t195\t465\t1154\t2822\t6825\t15519\t34049\t67800\t130776\t226793\t378735\t609982\t877332\t1318193\t1771959\t2344063\t3139646\t3608668\t4521353\t5173885\t5952248\t6459240\t7990153\t8442539\t9499475\t9717954\t11015514\t10890866\t11046625\t10317936\t9892251\t9496761\t9129367\t9026230\t9144000\t8700617\t8948395\t9255960\t8228757\t9203941\t8478746\t7824526\t7055110\t5765055\t4731788\t3907543\t3376077\t2811085\t2516859\t2133940\t1787700\t1674723\t1416879\t1184850\t1039639\t752891\t650764\t502048\t404336\t304800\t255400\t200545\t159395\t135705\t109548\t98308\t88101\t74819\t63665\t45409\t29651\t16678\t7768\t3446\t1451\t551\t193\t85\t29\t21\t8\t4\t2\t4\t9\t46\t38\t10\t15\n"
#    },
#    "H01_GCparagon": {
#       "original": "H01\tGCparagon\toriginal\t159\t24\t24\t58\t65\t46\t81\t100\t146\t158\t176\t210\t239\t431\t579\t1020\t2034\t4770\t13049\t34294\t84481\t181248\t358349\t643303\t1044519\t1647826\t2407633\t3335407\t4541816\t5486801\t6831467\t8000470\t9121326\t9954613\t11805217\t12502157\t13703826\t14247685\t15705560\t15621958\t14937738\t13682423\t12623700\t11824436\t11152301\t10702933\t10350384\t9703155\t9312687\t9084252\t7840938\t7957161\t6995128\t6286575\t5456217\t4469988\t3658304\t2999150\t2537921\t2134874\t1843619\t1574387\t1304530\t1197298\t996517\t836803\t720575\t540797\t452625\t352577\t282373\t222345\t187033\t151396\t122940\t103385\t88173\t77233\t71146\t60608\t49576\t37539\t26147\t18104\t11991\t7182\t3983\t2171\t1017\t348\t192\t95\t58\t32\t16\t19\t15\t41\t42\t21\t29\n"
#    }
# }
# GC content data for sample 'H01' written to stats TSV file 'GCPARAGON-STATISTICS_allSamples_fragment_length_14-08-2023.tsv':
# ----------------------------------------
# {
#    "C01_GCparagon": {
#       "original": "C01\tGCparagon\toriginal\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t1867\t2144\t2365\t48894\t50509\t50396\t50156\t49329\t51655\t54412\t58807\t63948\t71577\t79086\t79205\t75152\t77848\t78415\t84340\t90678\t99110\t111990\t131494\t160936\t185797\t185997\t164455\t152284\t159417\t176547\t200157\t232611\t280816\t357184\t467278\t576525\t630351\t603590\t594876\t659478\t787943\t977110\t1191522\t1389891\t1517521\t1602981\t1699871\t1801748\t1882318\t1978064\t2125188\t2309836\t2532112\t2801255\t3046383\t3181286\t3264601\t3461698\t3795099\t4254543\t4793676\t5335008\t5810205\t6275901\t6859630\t7724267\t8838348\t9720798\t9968505\t9789621\t9390739\t8909936\t8332369\t7738799\t7213352\t6856205\t6558623\t6269774\t5924184\t5497898\t5013880\t4531909\t4082831\t3676280\t3343139\t3061883\t2847661\t2686066\t2542297\t2396912\t2254845\t2093050\t1924053\t1754485\t1599824\t1465344\t1357468\t1269362\t1197725\t1139952\t1082182\t1024027\t966088\t898305\t828877\t762927\t705538\t654258\t612894\t573948\t543428\t513129\t481969\t449685\t415853\t380974\t347621\t318100\t294986\t274002\t256666\t241832\t226520\t211748\t195859\t179225\t165505\t150051\t138676\t130145\t122849\t116830\t112251\t105299\t98679\t92397\t86950\t82628\t78690\t74254\t70957\t70313\t86183\t62546\t3533\t3530\t4155\t3004\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n"
#    },
#    "H01_GCparagon": {
#       "original": "H01\tGCparagon\toriginal\t2248\t2187\t2299\t2333\t2451\t2282\t2639\t2341\t2710\t2471\t2691\t2633\t2598\t2645\t2933\t2878\t3071\t2992\t3193\t3260\t3414\t3319\t3376\t3455\t3576\t3576\t3766\t3603\t3777\t3889\t4112\t4287\t4359\t4593\t4863\t4985\t4677\t4771\t5138\t5365\t5583\t5375\t5570\t5490\t5577\t5665\t5728\t5758\t6146\t6526\t6734\t6492\t6371\t6364\t6504\t6700\t6692\t7139\t7901\t9080\t10093\t9970\t9038\t8197\t8318\t8409\t8753\t9085\t9872\t11245\t13817\t15511\t16058\t14073\t12432\t12136\t12515\t13600\t15404\t18847\t21335\t22476\t22862\t22525\t22198\t21615\t22113\t24112\t26455\t30716\t36569\t41516\t42239\t38887\t40627\t40467\t42656\t45649\t51678\t60441\t73406\t91959\t106589\t106657\t90759\t78355\t76837\t82109\t90478\t105041\t125179\t157274\t206018\t256950\t277998\t252059\t231661\t243498\t286802\t363715\t461399\t552723\t605429\t637467\t674002\t716464\t744568\t778365\t836794\t918405\t1017646\t1146888\t1275198\t1334118\t1341879\t1389006\t1496337\t1680292\t1920938\t2194642\t2454778\t2672539\t2877438\t3183576\t3625418\t4087360\t4422956\t4656763\t4782399\t4856486\t4851282\t4807448\t4807275\t4942476\t5093511\t5251747\t5337861\t5280452\t5102121\t4866249\t4618672\t4369500\t4180146\t4025809\t3927404\t3883428\t3845222\t3784112\t3692156\t3535681\t3342052\t3125287\t2920964\t2747632\t2605169\t2503837\t2426528\t2359716\t2293181\t2210181\t2112386\t1991794\t1862812\t1737729\t1626855\t1533994\t1457871\t1393637\t1335565\t1279422\t1218407\t1153368\t1079865\t1001081\t928893\t869346\t815363\t769623\t731875\t698264\t663041\t627571\t589083\t548183\t511134\t471688\t440766\t414394\t393017\t373734\t355349\t336264\t316090\t296478\t277860\t259035\t242199\t228845\t218479\t210086\t203239\t195220\t187835\t180011\t170479\t162129\t154852\t148793\t143493\t139841\t137963\t136076\t133757\t131976\t127443\t124357\t120992\t118724\t116128\t113638\t113318\t113488\t113675\t114670\t113689\t112041\t110301\t108996\t107304\t108678\t110784\t113078\t116760\t118593\t119548\t119092\t118354\t118597\t120223\t123815\t129695\t138278\t146412\t151249\t154679\t157071\t161598\t165072\t173701\t181123\t191404\t201028\t211349\t223063\t232806\t242765\t253742\t265979\t279709\t293941\t309404\t321588\t336015\t352417\t372779\t397405\t420224\t444088\t467697\t485427\t505228\t524222\t543567\t564713\t589524\t612308\t634552\t655677\t672748\t694940\t711865\t728318\t744416\t758162\t772643\t784611\t795876\t805058\t819367\t831780\t852025\t867335\t880942\t887473\t888269\t888466\t884705\t877096\t881300\t882430\t887501\t891398\t898619\t900144\t902724\t895990\t887260\t881698\t877316\t870437\t864424\t862531\t858614\t853324\t847597\t840222\t831752\t825770\t816849\t807525\t798329\t786683\t776025\t764000\t752481\t738732\t728761\t718588\t709180\t697138\t684448\t671695\t654386\t639044\t623541\t605318\t592527\t578927\t564368\t551573\t536755\t522702\t505732\t490669\t473949\t456729\t441981\t425590\t412417\t398670\t385557\t372354\t359945\t346197\t332294\t316672\t304260\t290143\t278031\t267162\t257747\t248636\t237829\t227758\t218108\t206725\t196995\t188045\t179304\t171430\t164078\t158122\t151188\t144489\t138024\t131575\t125666\t120224\t113458\t108530\t105325\t101079\t96791\t93239\t89881\t86450\t82261\t79483\t76468\t73995\t71597\t69977\t68341\t66473\t65106\t63944\t61598\t60183\t58885\t58340\t57076\t56619\t56799\t56685\t56796\t57454\t57115\t57058\t56681\t57013\t58027\t58505\t59890\t60956\t62578\t63246\t63581\t64345\t64471\t65462\t67109\t68155\t70516\t72750\t74605\t75366\t75961\t76225\t77057\t77890\t79868\t81988\t83977\t86209\t87356\t88468\t88947\t89754\t90497\t92055\t94197\t94534\t95761\t96910\t97402\t98476\t99062\t99304\t99907\t100373\t101130\t100807\t100553\t101492\t101201\t102293\t102769\t103688\t104060\t105120\t106411\t106220\t106440\t106574\t107745\t108274\t109072\t109621\t109802\t112115\t111600\t112407\t113210\t113899\t114510\t115296\t116538\t116327\t118038\t118151\t119011\t119518\t120991\t121617\t122817\t122354\t122671\t123541\t123483\t124461\t124495\t124963\t124705\t125193\t125787\t125525\t124858\t125477\t125201\t124888\t124115\t124020\t123469\t123380\t122560\t122557\t121727\t119976\t119488\t119077\t117488\t117011\t114790\t115572\t114063\t112615\t111721\t110446\t108192\t107864\t106438\t104534\t102734\t101588\t100171\t98522\t97239\t95200\t93832\t92476\t90371\t88265\t86599\t85110\t82945\t81412\t79751\t78396\t76392\t74704\t73285\t70835\t70021\t67743\t66192\t64566\t63084\t61585\t59335\t58381\t56461\t55381\t53862\t52173\t50884\t49323\t48254\t47120\t45524\t44636\t42983\t41661\t40652\t39481\t38181\t37262\t36640\t35745\t34634\t33722\t32660\t32032\t31658\t30489\t29984\t29279\t28498\t27986\t27509\t26843\t26715\t25668\t25362\t25071\t24552\t24315\t23915\t23769\t23521\t23253\t22765\t22181\t22459\t21996\t22080\t21860\t21859\t21608\t21325\t21397\t20883\t20941\t20704\t20714\t20779\t20547\t20467\t20550\t20165\t20056\t20309\t19940\t19876\t19619\t19978\t19838\t19939\t19750\t19830\t19690\t19636\t19502\t19162\t19392\t19237\t19481\t19534\t19598\t19699\t19620\t19555\t19555\t19574\t19610\t19848\t19741\t20005\t19974\t20067\t20169\t20404\t20310\t20279\t20290\t20415\t20341\t20665\t20386\t20514\t20833\t20869\t20690\t20862\t21196\t21444\t21257\t21506\t21265\t21181\t21747\t21908\t21956\t21868\t21966\t22535\t22278\t22334\t22575\t22795\t22827\t22478\t23034\t22833\t22996\t22897\t22905\t22907\t23142\t22985\t23291\t22822\t23220\t23241\t23170\t23137\t23006\t23012\t22924\t22979\t22823\t22784\t22805\t22706\t22341\t22553\t22336\t22075\t21945\t22223\t21719\t21975\t21614\t21474\t21153\t21034\t20962\t20817\t20298\t20090\t20283\t19980\t19593\t19366\t19129\t18639\t18495\t18305\t18152\t17789\t17563\t17280\t17094\t16638\t16461\t16446\t16039\t15424\t15457\t15261\t14886\t14386\t13864\t13911\t13635\t12783\t12578\t12353\t12019\t11028\t10731\t10858\t10343\t9160\t9031\t9073\t8592\t7115\t6955\t6977\t6428\t4909\t4855\t4931\t4588\t2874\t2760\t2828\t2741\t1224\t1122\t1229\n"
#    }
# }
# INFO - processing BAM file 'B01.GCtagged.bam'
# | 2023-08-14 16:03:52,321 - GCPARAGON_validation - INFO - inferred GC-tag: 'GC'
# | 2023-08-14 16:18:55,044 - GCPARAGON_validation - INFO - received stats for sample data 'B01_GCparagon'
# | 2023-08-14 16:18:55,045 - GCPARAGON_validation - INFO - received stats for sample data 'B01_GCparagon'
# | 2023-08-14 16:18:55,046 - GCPARAGON_validation - INFO - received stats for sample data 'B01_GCparagon'
# | 2023-08-14 16:18:55,047 - GCPARAGON_validation - INFO - received stats for sample data 'B01_GCparagon'
# | 2023-08-14 16:18:55,047 - GCPARAGON_validation - INFO - received stats for sample data 'B01_GCparagon'
# | 2023-08-14 16:18:55,048 - GCPARAGON_validation - INFO - received stats for sample data 'B01_GCparagon'
# | 2023-08-14 16:18:55,048 - GCPARAGON_validation - INFO - received stats for sample data 'B01_GCparagon'
# | 2023-08-14 16:18:55,049 - GCPARAGON_validation - INFO - received stats for sample data 'B01_GCparagon'
# | 2023-08-14 16:18:55,049 - GCPARAGON_validation - INFO - received stats for sample data 'B01_GCparagon'
# | 2023-08-14 16:18:55,050 - GCPARAGON_validation - INFO - received stats for sample data 'B01_GCparagon'
# | 2023-08-14 16:18:55,051 - GCPARAGON_validation - INFO - received stats for sample data 'B01_GCparagon'
# | 2023-08-14 16:18:55,051 - GCPARAGON_validation - INFO - received stats for sample data 'B01_GCparagon'
# | 2023-08-14 16:18:55,052 - GCPARAGON_validation - INFO - received stats for sample data 'B01_GCparagon'
# | 2023-08-14 16:18:55,052 - GCPARAGON_validation - INFO - received stats for sample data 'B01_GCparagon'
# | 2023-08-14 16:18:55,053 - GCPARAGON_validation - INFO - received stats for sample data 'B01_GCparagon'
# | 2023-08-14 16:18:55,053 - GCPARAGON_validation - INFO - received stats for sample data 'B01_GCparagon'
# | 2023-08-14 16:18:55,054 - GCPARAGON_validation - INFO - received stats for sample data 'B01_GCparagon'
# | 2023-08-14 16:18:55,054 - GCPARAGON_validation - INFO - received stats for sample data 'B01_GCparagon'
# | 2023-08-14 16:18:55,055 - GCPARAGON_validation - INFO - received stats for sample data 'B01_GCparagon'
# | 2023-08-14 16:18:55,055 - GCPARAGON_validation - INFO - received stats for sample data 'B01_GCparagon'
# | 2023-08-14 16:18:55,056 - GCPARAGON_validation - INFO - received stats for sample data 'B01_GCparagon'
# | 2023-08-14 16:18:55,056 - GCPARAGON_validation - INFO - received stats for sample data 'B01_GCparagon'
# | 2023-08-14 16:18:55,057 - GCPARAGON_validation - INFO - received stats for sample data 'B01_GCparagon'
# INFO - accumulating statistics for sample 'B01'
# GC content data for sample 'B01' written to stats TSV file 'GCPARAGON-STATISTICS_allSamples_GC-percentage_perFragmentSequence_14-08-2023.tsv':
# ----------------------------------------
# {
#    "B01_GCparagon": {
#       "original": "B01\tGCparagon\toriginal\t285\t99\t113\t107\t200\t231\t403\t679\t1158\t1992\t3309\t5295\t7550\t12824\t18916\t30773\t49839\t82322\t146033\t277603\t482295\t775842\t1138328\t1639861\t2180734\t3123772\t3898125\t4920979\t6201150\t6781688\t7964465\t8624369\t9312134\t9630747\t10996264\t10982276\t11447967\t10953809\t11546315\t10750214\t10333426\t9244136\t8453315\t7729137\t6982718\t6570434\t6254600\t5645380\t5468283\t5170309\t4285208\t4255176\t3580004\t3051652\t2557137\t1989014\t1545779\t1218040\t992148\t782193\t653437\t521637\t412361\t363177\t286226\t225310\t183847\t126761\t102964\t75464\t57497\t42233\t33631\t25680\t19854\t15801\t12559\t10856\t9110\t7543\t6157\t4447\t2969\t2091\t1303\t733\t347\t181\t71\t34\t10\t11\t5\t6\t3\t6\t19\t68\t77\t18\t57\n"
#    },
#    "C01_GCparagon": {
#       "original": "C01\tGCparagon\toriginal\t13\t9\t4\t7\t4\t5\t6\t6\t13\t22\t26\t95\t195\t465\t1154\t2822\t6825\t15519\t34049\t67800\t130776\t226793\t378735\t609982\t877332\t1318193\t1771959\t2344063\t3139646\t3608668\t4521353\t5173885\t5952248\t6459240\t7990153\t8442539\t9499475\t9717954\t11015514\t10890866\t11046625\t10317936\t9892251\t9496761\t9129367\t9026230\t9144000\t8700617\t8948395\t9255960\t8228757\t9203941\t8478746\t7824526\t7055110\t5765055\t4731788\t3907543\t3376077\t2811085\t2516859\t2133940\t1787700\t1674723\t1416879\t1184850\t1039639\t752891\t650764\t502048\t404336\t304800\t255400\t200545\t159395\t135705\t109548\t98308\t88101\t74819\t63665\t45409\t29651\t16678\t7768\t3446\t1451\t551\t193\t85\t29\t21\t8\t4\t2\t4\t9\t46\t38\t10\t15\n"
#    },
#    "H01_GCparagon": {
#       "original": "H01\tGCparagon\toriginal\t159\t24\t24\t58\t65\t46\t81\t100\t146\t158\t176\t210\t239\t431\t579\t1020\t2034\t4770\t13049\t34294\t84481\t181248\t358349\t643303\t1044519\t1647826\t2407633\t3335407\t4541816\t5486801\t6831467\t8000470\t9121326\t9954613\t11805217\t12502157\t13703826\t14247685\t15705560\t15621958\t14937738\t13682423\t12623700\t11824436\t11152301\t10702933\t10350384\t9703155\t9312687\t9084252\t7840938\t7957161\t6995128\t6286575\t5456217\t4469988\t3658304\t2999150\t2537921\t2134874\t1843619\t1574387\t1304530\t1197298\t996517\t836803\t720575\t540797\t452625\t352577\t282373\t222345\t187033\t151396\t122940\t103385\t88173\t77233\t71146\t60608\t49576\t37539\t26147\t18104\t11991\t7182\t3983\t2171\t1017\t348\t192\t95\t58\t32\t16\t19\t15\t41\t42\t21\t29\n"
#    }
# }
# GC content data for sample 'B01' written to stats TSV file 'GCPARAGON-STATISTICS_allSamples_fragment_length_14-08-2023.tsv':
# ----------------------------------------
# {
#    "B01_GCparagon": {
#       "original": "B01\tGCparagon\toriginal\t110\t150\t156\t175\t138\t168\t174\t152\t184\t157\t324\t377\t391\t413\t391\t440\t525\t519\t540\t633\t703\t704\t802\t818\t965\t1026\t1119\t1153\t1393\t1617\t1781\t1923\t1938\t2098\t2233\t2436\t2742\t3030\t3403\t3789\t4249\t4332\t4311\t4445\t4831\t5021\t5113\t5481\t6207\t7122\t7631\t7568\t7662\t7881\t8224\t8786\t9340\t10320\t11761\t14087\t16302\t16898\t15311\t13960\t14398\t15086\t16885\t18129\t20703\t24968\t31264\t37132\t37893\t34336\t31733\t32200\t34866\t39299\t43507\t49077\t54641\t60958\t65665\t64995\t65762\t67046\t68885\t74155\t79704\t88240\t98448\t108466\t111694\t107002\t109593\t109786\t118462\t123199\t138856\t151927\t172066\t207821\t233616\t234952\t215285\t202351\t211324\t234732\t265937\t299671\t355072\t431938\t547293\t657602\t711260\t690806\t692105\t764757\t896906\t1072275\t1242869\t1394885\t1512112\t1611652\t1697449\t1779743\t1849220\t1950328\t2110778\t2297261\t2507543\t2741787\t2978871\t3137803\t3252302\t3428487\t3696657\t4071495\t4513504\t4977229\t5411473\t5844077\t6359029\t7097879\t8000146\t8648117\t8759910\t8457466\t7973811\t7441228\t6841500\t6216190\t5660024\t5222047\t4830343\t4481088\t4107167\t3701677\t3306898\t2932863\t2609705\t2329616\t2095975\t1898324\t1738426\t1612115\t1500468\t1397531\t1298232\t1194812\t1096925\t1003611\t918857\t844759\t783290\t732009\t690171\t654723\t623473\t590505\t556505\t519702\t483639\t450666\t420350\t394546\t372218\t353725\t335032\t316878\t299772\t280210\t260082\t241191\t220843\t203815\t188299\t176508\t165517\t157469\t147500\t139251\t127764\t118193\t108584\t99205\t92290\t87175\t83953\t85610\t80024\t71463\t66945\t65637\t59099\t56226\t53390\t52285\t50637\t50520\t50873\t50395\t50136\t48231\t46391\t44429\t43568\t43015\t42857\t43072\t43887\t45121\t45957\t45427\t44603\t43139\t42480\t41375\t42274\t42999\t45389\t48517\t50673\t52978\t52592\t52290\t51016\t50986\t52482\t55563\t60631\t66893\t72277\t76763\t78422\t79880\t80778\t83204\t87822\t94190\t101644\t109133\t115033\t120040\t122315\t125693\t128622\t134427\t142155\t151633\t161357\t167495\t173775\t178536\t180648\t183964\t189626\t196144\t203732\t208537\t211322\t210319\t209068\t207786\t209126\t210330\t209899\t207686\t205085\t201574\t197857\t193037\t188524\t183127\t178105\t172487\t167195\t161394\t155226\t149449\t144733\t140063\t135020\t130384\t125793\t122176\t116105\t111829\t108482\t104650\t101095\t99192\t96555\t93466\t90458\t87804\t85110\t81991\t80129\t77560\t76317\t73942\t72028\t71070\t69183\t67644\t66348\t64775\t63488\t61846\t59651\t58509\t56585\t55185\t53538\t51957\t50989\t49606\t48670\t47197\t45850\t43917\t43023\t41240\t39615\t38515\t37549\t35621\t34873\t33966\t32668\t31597\t30038\t28338\t27679\t26334\t24797\t24006\t22847\t21666\t21061\t20156\t19141\t18110\t17406\t16644\t15687\t14773\t13903\t13628\t12764\t12433\t11656\t11204\t10429\t9891\t9434\t8934\t8596\t8020\t7583\t7470\t7196\t6921\t6641\t6285\t5994\t5730\t5446\t5382\t5109\t4917\t4933\t4803\t4362\t4322\t4313\t4187\t4084\t3945\t3887\t3841\t3701\t3718\t3726\t3492\t3455\t3619\t3387\t3418\t3288\t3350\t3264\t3194\t3287\t3210\t3262\t3168\t3112\t3164\t3042\t3159\t3032\t3150\t3011\t2905\t2974\t2907\t3010\t2889\t2858\t2828\t2806\t2849\t2690\t2674\t2703\t2595\t2646\t2567\t2684\t2546\t2586\t2497\t2508\t2521\t2466\t2402\t2346\t2327\t2285\t2407\t2368\t2286\t2342\t2327\t2285\t2217\t2287\t2251\t2202\t2200\t2205\t2102\t2131\t2089\t2107\t2060\t2045\t2014\t2042\t1969\t1956\t1941\t1854\t1810\t1867\t1850\t1869\t1875\t1793\t1789\t1880\t1748\t1702\t1821\t1727\t1662\t1599\t1610\t1634\t1532\t1662\t1639\t1522\t1546\t1537\t1461\t1475\t1635\t1535\t1435\t1396\t1441\t1443\t1467\t1465\t1381\t1474\t1301\t1396\t1299\t1352\t1285\t1232\t1252\t1220\t1270\t1271\t1225\t1202\t1217\t1197\t1035\t1144\t1112\t1076\t1076\t1117\t1089\t1032\t1066\t989\t992\t988\t981\t951\t913\t886\t906\t899\t873\t867\t822\t794\t789\t810\t753\t680\t738\t716\t718\t680\t701\t679\t598\t619\t629\t582\t564\t544\t563\t549\t480\t486\t481\t484\t490\t434\t440\t420\t415\t412\t403\t397\t389\t372\t360\t333\t348\t315\t307\t302\t274\t287\t282\t286\t256\t276\t254\t233\t245\t268\t243\t246\t251\t195\t214\t225\t215\t201\t171\t208\t214\t192\t170\t178\t171\t188\t162\t168\t209\t168\t169\t172\t183\t165\t173\t177\t162\t141\t173\t167\t167\t157\t126\t156\t153\t130\t162\t146\t139\t127\t140\t144\t154\t138\t121\t127\t143\t153\t127\t136\t154\t150\t116\t137\t123\t99\t123\t130\t129\t141\t107\t134\t130\t112\t113\t112\t127\t107\t96\t100\t135\t120\t98\t129\t126\t97\t93\t99\t89\t112\t96\t130\t102\t98\t113\t110\t114\t106\t105\t96\t116\t91\t119\t114\t71\t107\t104\t94\t92\t89\t93\t103\t91\t92\t84\t98\t86\t75\t86\t83\t85\t97\t91\t74\t100\t91\t72\t78\t82\t88\t79\t79\t75\t72\t73\t88\t81\t128\t92\t93\t88\t67\t52\t59\t79\t62\t65\t71\t85\t66\t59\t65\t77\t76\t116\t63\t54\t76\t59\t51\t58\t55\t63\t56\t58\t56\t51\t46\t45\t52\t67\t57\t48\t43\t40\t51\t34\t45\t39\t49\t42\t48\t42\t42\t35\t57\t34\t34\t34\t33\t37\t44\t42\t37\t29\t35\t34\t29\n"
#    },
#    "C01_GCparagon": {
#       "original": "C01\tGCparagon\toriginal\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t1867\t2144\t2365\t48894\t50509\t50396\t50156\t49329\t51655\t54412\t58807\t63948\t71577\t79086\t79205\t75152\t77848\t78415\t84340\t90678\t99110\t111990\t131494\t160936\t185797\t185997\t164455\t152284\t159417\t176547\t200157\t232611\t280816\t357184\t467278\t576525\t630351\t603590\t594876\t659478\t787943\t977110\t1191522\t1389891\t1517521\t1602981\t1699871\t1801748\t1882318\t1978064\t2125188\t2309836\t2532112\t2801255\t3046383\t3181286\t3264601\t3461698\t3795099\t4254543\t4793676\t5335008\t5810205\t6275901\t6859630\t7724267\t8838348\t9720798\t9968505\t9789621\t9390739\t8909936\t8332369\t7738799\t7213352\t6856205\t6558623\t6269774\t5924184\t5497898\t5013880\t4531909\t4082831\t3676280\t3343139\t3061883\t2847661\t2686066\t2542297\t2396912\t2254845\t2093050\t1924053\t1754485\t1599824\t1465344\t1357468\t1269362\t1197725\t1139952\t1082182\t1024027\t966088\t898305\t828877\t762927\t705538\t654258\t612894\t573948\t543428\t513129\t481969\t449685\t415853\t380974\t347621\t318100\t294986\t274002\t256666\t241832\t226520\t211748\t195859\t179225\t165505\t150051\t138676\t130145\t122849\t116830\t112251\t105299\t98679\t92397\t86950\t82628\t78690\t74254\t70957\t70313\t86183\t62546\t3533\t3530\t4155\t3004\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n"
#    },
#    "H01_GCparagon": {
#       "original": "H01\tGCparagon\toriginal\t2248\t2187\t2299\t2333\t2451\t2282\t2639\t2341\t2710\t2471\t2691\t2633\t2598\t2645\t2933\t2878\t3071\t2992\t3193\t3260\t3414\t3319\t3376\t3455\t3576\t3576\t3766\t3603\t3777\t3889\t4112\t4287\t4359\t4593\t4863\t4985\t4677\t4771\t5138\t5365\t5583\t5375\t5570\t5490\t5577\t5665\t5728\t5758\t6146\t6526\t6734\t6492\t6371\t6364\t6504\t6700\t6692\t7139\t7901\t9080\t10093\t9970\t9038\t8197\t8318\t8409\t8753\t9085\t9872\t11245\t13817\t15511\t16058\t14073\t12432\t12136\t12515\t13600\t15404\t18847\t21335\t22476\t22862\t22525\t22198\t21615\t22113\t24112\t26455\t30716\t36569\t41516\t42239\t38887\t40627\t40467\t42656\t45649\t51678\t60441\t73406\t91959\t106589\t106657\t90759\t78355\t76837\t82109\t90478\t105041\t125179\t157274\t206018\t256950\t277998\t252059\t231661\t243498\t286802\t363715\t461399\t552723\t605429\t637467\t674002\t716464\t744568\t778365\t836794\t918405\t1017646\t1146888\t1275198\t1334118\t1341879\t1389006\t1496337\t1680292\t1920938\t2194642\t2454778\t2672539\t2877438\t3183576\t3625418\t4087360\t4422956\t4656763\t4782399\t4856486\t4851282\t4807448\t4807275\t4942476\t5093511\t5251747\t5337861\t5280452\t5102121\t4866249\t4618672\t4369500\t4180146\t4025809\t3927404\t3883428\t3845222\t3784112\t3692156\t3535681\t3342052\t3125287\t2920964\t2747632\t2605169\t2503837\t2426528\t2359716\t2293181\t2210181\t2112386\t1991794\t1862812\t1737729\t1626855\t1533994\t1457871\t1393637\t1335565\t1279422\t1218407\t1153368\t1079865\t1001081\t928893\t869346\t815363\t769623\t731875\t698264\t663041\t627571\t589083\t548183\t511134\t471688\t440766\t414394\t393017\t373734\t355349\t336264\t316090\t296478\t277860\t259035\t242199\t228845\t218479\t210086\t203239\t195220\t187835\t180011\t170479\t162129\t154852\t148793\t143493\t139841\t137963\t136076\t133757\t131976\t127443\t124357\t120992\t118724\t116128\t113638\t113318\t113488\t113675\t114670\t113689\t112041\t110301\t108996\t107304\t108678\t110784\t113078\t116760\t118593\t119548\t119092\t118354\t118597\t120223\t123815\t129695\t138278\t146412\t151249\t154679\t157071\t161598\t165072\t173701\t181123\t191404\t201028\t211349\t223063\t232806\t242765\t253742\t265979\t279709\t293941\t309404\t321588\t336015\t352417\t372779\t397405\t420224\t444088\t467697\t485427\t505228\t524222\t543567\t564713\t589524\t612308\t634552\t655677\t672748\t694940\t711865\t728318\t744416\t758162\t772643\t784611\t795876\t805058\t819367\t831780\t852025\t867335\t880942\t887473\t888269\t888466\t884705\t877096\t881300\t882430\t887501\t891398\t898619\t900144\t902724\t895990\t887260\t881698\t877316\t870437\t864424\t862531\t858614\t853324\t847597\t840222\t831752\t825770\t816849\t807525\t798329\t786683\t776025\t764000\t752481\t738732\t728761\t718588\t709180\t697138\t684448\t671695\t654386\t639044\t623541\t605318\t592527\t578927\t564368\t551573\t536755\t522702\t505732\t490669\t473949\t456729\t441981\t425590\t412417\t398670\t385557\t372354\t359945\t346197\t332294\t316672\t304260\t290143\t278031\t267162\t257747\t248636\t237829\t227758\t218108\t206725\t196995\t188045\t179304\t171430\t164078\t158122\t151188\t144489\t138024\t131575\t125666\t120224\t113458\t108530\t105325\t101079\t96791\t93239\t89881\t86450\t82261\t79483\t76468\t73995\t71597\t69977\t68341\t66473\t65106\t63944\t61598\t60183\t58885\t58340\t57076\t56619\t56799\t56685\t56796\t57454\t57115\t57058\t56681\t57013\t58027\t58505\t59890\t60956\t62578\t63246\t63581\t64345\t64471\t65462\t67109\t68155\t70516\t72750\t74605\t75366\t75961\t76225\t77057\t77890\t79868\t81988\t83977\t86209\t87356\t88468\t88947\t89754\t90497\t92055\t94197\t94534\t95761\t96910\t97402\t98476\t99062\t99304\t99907\t100373\t101130\t100807\t100553\t101492\t101201\t102293\t102769\t103688\t104060\t105120\t106411\t106220\t106440\t106574\t107745\t108274\t109072\t109621\t109802\t112115\t111600\t112407\t113210\t113899\t114510\t115296\t116538\t116327\t118038\t118151\t119011\t119518\t120991\t121617\t122817\t122354\t122671\t123541\t123483\t124461\t124495\t124963\t124705\t125193\t125787\t125525\t124858\t125477\t125201\t124888\t124115\t124020\t123469\t123380\t122560\t122557\t121727\t119976\t119488\t119077\t117488\t117011\t114790\t115572\t114063\t112615\t111721\t110446\t108192\t107864\t106438\t104534\t102734\t101588\t100171\t98522\t97239\t95200\t93832\t92476\t90371\t88265\t86599\t85110\t82945\t81412\t79751\t78396\t76392\t74704\t73285\t70835\t70021\t67743\t66192\t64566\t63084\t61585\t59335\t58381\t56461\t55381\t53862\t52173\t50884\t49323\t48254\t47120\t45524\t44636\t42983\t41661\t40652\t39481\t38181\t37262\t36640\t35745\t34634\t33722\t32660\t32032\t31658\t30489\t29984\t29279\t28498\t27986\t27509\t26843\t26715\t25668\t25362\t25071\t24552\t24315\t23915\t23769\t23521\t23253\t22765\t22181\t22459\t21996\t22080\t21860\t21859\t21608\t21325\t21397\t20883\t20941\t20704\t20714\t20779\t20547\t20467\t20550\t20165\t20056\t20309\t19940\t19876\t19619\t19978\t19838\t19939\t19750\t19830\t19690\t19636\t19502\t19162\t19392\t19237\t19481\t19534\t19598\t19699\t19620\t19555\t19555\t19574\t19610\t19848\t19741\t20005\t19974\t20067\t20169\t20404\t20310\t20279\t20290\t20415\t20341\t20665\t20386\t20514\t20833\t20869\t20690\t20862\t21196\t21444\t21257\t21506\t21265\t21181\t21747\t21908\t21956\t21868\t21966\t22535\t22278\t22334\t22575\t22795\t22827\t22478\t23034\t22833\t22996\t22897\t22905\t22907\t23142\t22985\t23291\t22822\t23220\t23241\t23170\t23137\t23006\t23012\t22924\t22979\t22823\t22784\t22805\t22706\t22341\t22553\t22336\t22075\t21945\t22223\t21719\t21975\t21614\t21474\t21153\t21034\t20962\t20817\t20298\t20090\t20283\t19980\t19593\t19366\t19129\t18639\t18495\t18305\t18152\t17789\t17563\t17280\t17094\t16638\t16461\t16446\t16039\t15424\t15457\t15261\t14886\t14386\t13864\t13911\t13635\t12783\t12578\t12353\t12019\t11028\t10731\t10858\t10343\t9160\t9031\t9073\t8592\t7115\t6955\t6977\t6428\t4909\t4855\t4931\t4588\t2874\t2760\t2828\t2741\t1224\t1122\t1229\n"
#    }
# }
# INFO - processing BAM file 'P01.GCtagged.bam'
# | 2023-08-14 16:18:55,150 - GCPARAGON_validation - INFO - inferred GC-tag: 'GC'
# | 2023-08-14 16:53:00,079 - GCPARAGON_validation - INFO - received stats for sample data 'P01_GCparagon'
# | 2023-08-14 16:53:00,080 - GCPARAGON_validation - INFO - received stats for sample data 'P01_GCparagon'
# | 2023-08-14 16:53:00,081 - GCPARAGON_validation - INFO - received stats for sample data 'P01_GCparagon'
# | 2023-08-14 16:53:00,082 - GCPARAGON_validation - INFO - received stats for sample data 'P01_GCparagon'
# | 2023-08-14 16:53:00,082 - GCPARAGON_validation - INFO - received stats for sample data 'P01_GCparagon'
# | 2023-08-14 16:53:00,083 - GCPARAGON_validation - INFO - received stats for sample data 'P01_GCparagon'
# | 2023-08-14 16:53:00,083 - GCPARAGON_validation - INFO - received stats for sample data 'P01_GCparagon'
# | 2023-08-14 16:53:00,084 - GCPARAGON_validation - INFO - received stats for sample data 'P01_GCparagon'
# | 2023-08-14 16:53:00,084 - GCPARAGON_validation - INFO - received stats for sample data 'P01_GCparagon'
# | 2023-08-14 16:53:00,084 - GCPARAGON_validation - INFO - received stats for sample data 'P01_GCparagon'
# | 2023-08-14 16:53:00,085 - GCPARAGON_validation - INFO - received stats for sample data 'P01_GCparagon'
# | 2023-08-14 16:53:00,085 - GCPARAGON_validation - INFO - received stats for sample data 'P01_GCparagon'
# | 2023-08-14 16:53:00,086 - GCPARAGON_validation - INFO - received stats for sample data 'P01_GCparagon'
# | 2023-08-14 16:53:00,086 - GCPARAGON_validation - INFO - received stats for sample data 'P01_GCparagon'
# | 2023-08-14 16:53:00,087 - GCPARAGON_validation - INFO - received stats for sample data 'P01_GCparagon'
# | 2023-08-14 16:53:00,087 - GCPARAGON_validation - INFO - received stats for sample data 'P01_GCparagon'
# | 2023-08-14 16:53:00,088 - GCPARAGON_validation - INFO - received stats for sample data 'P01_GCparagon'
# | 2023-08-14 16:53:00,088 - GCPARAGON_validation - INFO - received stats for sample data 'P01_GCparagon'
# | 2023-08-14 16:53:00,088 - GCPARAGON_validation - INFO - received stats for sample data 'P01_GCparagon'
# | 2023-08-14 16:53:00,089 - GCPARAGON_validation - INFO - received stats for sample data 'P01_GCparagon'
# | 2023-08-14 16:53:00,089 - GCPARAGON_validation - INFO - received stats for sample data 'P01_GCparagon'
# | 2023-08-14 16:53:00,089 - GCPARAGON_validation - INFO - received stats for sample data 'P01_GCparagon'
# | 2023-08-14 16:53:00,090 - GCPARAGON_validation - INFO - received stats for sample data 'P01_GCparagon'
# INFO - accumulating statistics for sample 'P01'
# GC content data for sample 'P01' written to stats TSV file 'GCPARAGON-STATISTICS_allSamples_GC-percentage_perFragmentSequence_14-08-2023.tsv':
# ----------------------------------------
# {
#    "B01_GCparagon": {
#       "original": "B01\tGCparagon\toriginal\t285\t99\t113\t107\t200\t231\t403\t679\t1158\t1992\t3309\t5295\t7550\t12824\t18916\t30773\t49839\t82322\t146033\t277603\t482295\t775842\t1138328\t1639861\t2180734\t3123772\t3898125\t4920979\t6201150\t6781688\t7964465\t8624369\t9312134\t9630747\t10996264\t10982276\t11447967\t10953809\t11546315\t10750214\t10333426\t9244136\t8453315\t7729137\t6982718\t6570434\t6254600\t5645380\t5468283\t5170309\t4285208\t4255176\t3580004\t3051652\t2557137\t1989014\t1545779\t1218040\t992148\t782193\t653437\t521637\t412361\t363177\t286226\t225310\t183847\t126761\t102964\t75464\t57497\t42233\t33631\t25680\t19854\t15801\t12559\t10856\t9110\t7543\t6157\t4447\t2969\t2091\t1303\t733\t347\t181\t71\t34\t10\t11\t5\t6\t3\t6\t19\t68\t77\t18\t57\n"
#    },
#    "C01_GCparagon": {
#       "original": "C01\tGCparagon\toriginal\t13\t9\t4\t7\t4\t5\t6\t6\t13\t22\t26\t95\t195\t465\t1154\t2822\t6825\t15519\t34049\t67800\t130776\t226793\t378735\t609982\t877332\t1318193\t1771959\t2344063\t3139646\t3608668\t4521353\t5173885\t5952248\t6459240\t7990153\t8442539\t9499475\t9717954\t11015514\t10890866\t11046625\t10317936\t9892251\t9496761\t9129367\t9026230\t9144000\t8700617\t8948395\t9255960\t8228757\t9203941\t8478746\t7824526\t7055110\t5765055\t4731788\t3907543\t3376077\t2811085\t2516859\t2133940\t1787700\t1674723\t1416879\t1184850\t1039639\t752891\t650764\t502048\t404336\t304800\t255400\t200545\t159395\t135705\t109548\t98308\t88101\t74819\t63665\t45409\t29651\t16678\t7768\t3446\t1451\t551\t193\t85\t29\t21\t8\t4\t2\t4\t9\t46\t38\t10\t15\n"
#    },
#    "H01_GCparagon": {
#       "original": "H01\tGCparagon\toriginal\t159\t24\t24\t58\t65\t46\t81\t100\t146\t158\t176\t210\t239\t431\t579\t1020\t2034\t4770\t13049\t34294\t84481\t181248\t358349\t643303\t1044519\t1647826\t2407633\t3335407\t4541816\t5486801\t6831467\t8000470\t9121326\t9954613\t11805217\t12502157\t13703826\t14247685\t15705560\t15621958\t14937738\t13682423\t12623700\t11824436\t11152301\t10702933\t10350384\t9703155\t9312687\t9084252\t7840938\t7957161\t6995128\t6286575\t5456217\t4469988\t3658304\t2999150\t2537921\t2134874\t1843619\t1574387\t1304530\t1197298\t996517\t836803\t720575\t540797\t452625\t352577\t282373\t222345\t187033\t151396\t122940\t103385\t88173\t77233\t71146\t60608\t49576\t37539\t26147\t18104\t11991\t7182\t3983\t2171\t1017\t348\t192\t95\t58\t32\t16\t19\t15\t41\t42\t21\t29\n"
#    },
#    "P01_GCparagon": {
#       "original": "P01\tGCparagon\toriginal\t757\t65\t98\t246\t258\t245\t389\t468\t648\t656\t946\t1001\t1246\t2045\t2723\t3996\t6171\t10365\t19013\t35355\t68204\t127497\t224121\t366532\t577393\t926922\t1361333\t1963231\t2774366\t3509221\t4556137\t5605698\t6701224\t7751334\t9638561\t10826181\t12428248\t13500794\t15612302\t16296681\t16701661\t16499107\t16198112\t16210017\t16033548\t16339220\t16577409\t16340501\t16444248\t16991440\t15286938\t16444389\t14982383\t14053702\t12598497\t10750865\t9048643\t7743061\t6747162\t5852843\t5180650\t4544603\t3831440\t3606313\t3017878\t2571299\t2212780\t1678219\t1406024\t1097453\t874778\t691469\t573649\t464632\t380240\t317127\t264247\t230544\t209891\t183902\t156128\t117401\t78163\t48251\t28724\t16204\t8676\t4232\t1841\t730\t393\t234\t123\t96\t61\t72\t206\t333\t209\t104\t151\n"
#    }
# }
# GC content data for sample 'P01' written to stats TSV file 'GCPARAGON-STATISTICS_allSamples_fragment_length_14-08-2023.tsv':
# ----------------------------------------
# {
#    "B01_GCparagon": {
#       "original": "B01\tGCparagon\toriginal\t110\t150\t156\t175\t138\t168\t174\t152\t184\t157\t324\t377\t391\t413\t391\t440\t525\t519\t540\t633\t703\t704\t802\t818\t965\t1026\t1119\t1153\t1393\t1617\t1781\t1923\t1938\t2098\t2233\t2436\t2742\t3030\t3403\t3789\t4249\t4332\t4311\t4445\t4831\t5021\t5113\t5481\t6207\t7122\t7631\t7568\t7662\t7881\t8224\t8786\t9340\t10320\t11761\t14087\t16302\t16898\t15311\t13960\t14398\t15086\t16885\t18129\t20703\t24968\t31264\t37132\t37893\t34336\t31733\t32200\t34866\t39299\t43507\t49077\t54641\t60958\t65665\t64995\t65762\t67046\t68885\t74155\t79704\t88240\t98448\t108466\t111694\t107002\t109593\t109786\t118462\t123199\t138856\t151927\t172066\t207821\t233616\t234952\t215285\t202351\t211324\t234732\t265937\t299671\t355072\t431938\t547293\t657602\t711260\t690806\t692105\t764757\t896906\t1072275\t1242869\t1394885\t1512112\t1611652\t1697449\t1779743\t1849220\t1950328\t2110778\t2297261\t2507543\t2741787\t2978871\t3137803\t3252302\t3428487\t3696657\t4071495\t4513504\t4977229\t5411473\t5844077\t6359029\t7097879\t8000146\t8648117\t8759910\t8457466\t7973811\t7441228\t6841500\t6216190\t5660024\t5222047\t4830343\t4481088\t4107167\t3701677\t3306898\t2932863\t2609705\t2329616\t2095975\t1898324\t1738426\t1612115\t1500468\t1397531\t1298232\t1194812\t1096925\t1003611\t918857\t844759\t783290\t732009\t690171\t654723\t623473\t590505\t556505\t519702\t483639\t450666\t420350\t394546\t372218\t353725\t335032\t316878\t299772\t280210\t260082\t241191\t220843\t203815\t188299\t176508\t165517\t157469\t147500\t139251\t127764\t118193\t108584\t99205\t92290\t87175\t83953\t85610\t80024\t71463\t66945\t65637\t59099\t56226\t53390\t52285\t50637\t50520\t50873\t50395\t50136\t48231\t46391\t44429\t43568\t43015\t42857\t43072\t43887\t45121\t45957\t45427\t44603\t43139\t42480\t41375\t42274\t42999\t45389\t48517\t50673\t52978\t52592\t52290\t51016\t50986\t52482\t55563\t60631\t66893\t72277\t76763\t78422\t79880\t80778\t83204\t87822\t94190\t101644\t109133\t115033\t120040\t122315\t125693\t128622\t134427\t142155\t151633\t161357\t167495\t173775\t178536\t180648\t183964\t189626\t196144\t203732\t208537\t211322\t210319\t209068\t207786\t209126\t210330\t209899\t207686\t205085\t201574\t197857\t193037\t188524\t183127\t178105\t172487\t167195\t161394\t155226\t149449\t144733\t140063\t135020\t130384\t125793\t122176\t116105\t111829\t108482\t104650\t101095\t99192\t96555\t93466\t90458\t87804\t85110\t81991\t80129\t77560\t76317\t73942\t72028\t71070\t69183\t67644\t66348\t64775\t63488\t61846\t59651\t58509\t56585\t55185\t53538\t51957\t50989\t49606\t48670\t47197\t45850\t43917\t43023\t41240\t39615\t38515\t37549\t35621\t34873\t33966\t32668\t31597\t30038\t28338\t27679\t26334\t24797\t24006\t22847\t21666\t21061\t20156\t19141\t18110\t17406\t16644\t15687\t14773\t13903\t13628\t12764\t12433\t11656\t11204\t10429\t9891\t9434\t8934\t8596\t8020\t7583\t7470\t7196\t6921\t6641\t6285\t5994\t5730\t5446\t5382\t5109\t4917\t4933\t4803\t4362\t4322\t4313\t4187\t4084\t3945\t3887\t3841\t3701\t3718\t3726\t3492\t3455\t3619\t3387\t3418\t3288\t3350\t3264\t3194\t3287\t3210\t3262\t3168\t3112\t3164\t3042\t3159\t3032\t3150\t3011\t2905\t2974\t2907\t3010\t2889\t2858\t2828\t2806\t2849\t2690\t2674\t2703\t2595\t2646\t2567\t2684\t2546\t2586\t2497\t2508\t2521\t2466\t2402\t2346\t2327\t2285\t2407\t2368\t2286\t2342\t2327\t2285\t2217\t2287\t2251\t2202\t2200\t2205\t2102\t2131\t2089\t2107\t2060\t2045\t2014\t2042\t1969\t1956\t1941\t1854\t1810\t1867\t1850\t1869\t1875\t1793\t1789\t1880\t1748\t1702\t1821\t1727\t1662\t1599\t1610\t1634\t1532\t1662\t1639\t1522\t1546\t1537\t1461\t1475\t1635\t1535\t1435\t1396\t1441\t1443\t1467\t1465\t1381\t1474\t1301\t1396\t1299\t1352\t1285\t1232\t1252\t1220\t1270\t1271\t1225\t1202\t1217\t1197\t1035\t1144\t1112\t1076\t1076\t1117\t1089\t1032\t1066\t989\t992\t988\t981\t951\t913\t886\t906\t899\t873\t867\t822\t794\t789\t810\t753\t680\t738\t716\t718\t680\t701\t679\t598\t619\t629\t582\t564\t544\t563\t549\t480\t486\t481\t484\t490\t434\t440\t420\t415\t412\t403\t397\t389\t372\t360\t333\t348\t315\t307\t302\t274\t287\t282\t286\t256\t276\t254\t233\t245\t268\t243\t246\t251\t195\t214\t225\t215\t201\t171\t208\t214\t192\t170\t178\t171\t188\t162\t168\t209\t168\t169\t172\t183\t165\t173\t177\t162\t141\t173\t167\t167\t157\t126\t156\t153\t130\t162\t146\t139\t127\t140\t144\t154\t138\t121\t127\t143\t153\t127\t136\t154\t150\t116\t137\t123\t99\t123\t130\t129\t141\t107\t134\t130\t112\t113\t112\t127\t107\t96\t100\t135\t120\t98\t129\t126\t97\t93\t99\t89\t112\t96\t130\t102\t98\t113\t110\t114\t106\t105\t96\t116\t91\t119\t114\t71\t107\t104\t94\t92\t89\t93\t103\t91\t92\t84\t98\t86\t75\t86\t83\t85\t97\t91\t74\t100\t91\t72\t78\t82\t88\t79\t79\t75\t72\t73\t88\t81\t128\t92\t93\t88\t67\t52\t59\t79\t62\t65\t71\t85\t66\t59\t65\t77\t76\t116\t63\t54\t76\t59\t51\t58\t55\t63\t56\t58\t56\t51\t46\t45\t52\t67\t57\t48\t43\t40\t51\t34\t45\t39\t49\t42\t48\t42\t42\t35\t57\t34\t34\t34\t33\t37\t44\t42\t37\t29\t35\t34\t29\n"
#    },
#    "C01_GCparagon": {
#       "original": "C01\tGCparagon\toriginal\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t1867\t2144\t2365\t48894\t50509\t50396\t50156\t49329\t51655\t54412\t58807\t63948\t71577\t79086\t79205\t75152\t77848\t78415\t84340\t90678\t99110\t111990\t131494\t160936\t185797\t185997\t164455\t152284\t159417\t176547\t200157\t232611\t280816\t357184\t467278\t576525\t630351\t603590\t594876\t659478\t787943\t977110\t1191522\t1389891\t1517521\t1602981\t1699871\t1801748\t1882318\t1978064\t2125188\t2309836\t2532112\t2801255\t3046383\t3181286\t3264601\t3461698\t3795099\t4254543\t4793676\t5335008\t5810205\t6275901\t6859630\t7724267\t8838348\t9720798\t9968505\t9789621\t9390739\t8909936\t8332369\t7738799\t7213352\t6856205\t6558623\t6269774\t5924184\t5497898\t5013880\t4531909\t4082831\t3676280\t3343139\t3061883\t2847661\t2686066\t2542297\t2396912\t2254845\t2093050\t1924053\t1754485\t1599824\t1465344\t1357468\t1269362\t1197725\t1139952\t1082182\t1024027\t966088\t898305\t828877\t762927\t705538\t654258\t612894\t573948\t543428\t513129\t481969\t449685\t415853\t380974\t347621\t318100\t294986\t274002\t256666\t241832\t226520\t211748\t195859\t179225\t165505\t150051\t138676\t130145\t122849\t116830\t112251\t105299\t98679\t92397\t86950\t82628\t78690\t74254\t70957\t70313\t86183\t62546\t3533\t3530\t4155\t3004\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n"
#    },
#    "H01_GCparagon": {
#       "original": "H01\tGCparagon\toriginal\t2248\t2187\t2299\t2333\t2451\t2282\t2639\t2341\t2710\t2471\t2691\t2633\t2598\t2645\t2933\t2878\t3071\t2992\t3193\t3260\t3414\t3319\t3376\t3455\t3576\t3576\t3766\t3603\t3777\t3889\t4112\t4287\t4359\t4593\t4863\t4985\t4677\t4771\t5138\t5365\t5583\t5375\t5570\t5490\t5577\t5665\t5728\t5758\t6146\t6526\t6734\t6492\t6371\t6364\t6504\t6700\t6692\t7139\t7901\t9080\t10093\t9970\t9038\t8197\t8318\t8409\t8753\t9085\t9872\t11245\t13817\t15511\t16058\t14073\t12432\t12136\t12515\t13600\t15404\t18847\t21335\t22476\t22862\t22525\t22198\t21615\t22113\t24112\t26455\t30716\t36569\t41516\t42239\t38887\t40627\t40467\t42656\t45649\t51678\t60441\t73406\t91959\t106589\t106657\t90759\t78355\t76837\t82109\t90478\t105041\t125179\t157274\t206018\t256950\t277998\t252059\t231661\t243498\t286802\t363715\t461399\t552723\t605429\t637467\t674002\t716464\t744568\t778365\t836794\t918405\t1017646\t1146888\t1275198\t1334118\t1341879\t1389006\t1496337\t1680292\t1920938\t2194642\t2454778\t2672539\t2877438\t3183576\t3625418\t4087360\t4422956\t4656763\t4782399\t4856486\t4851282\t4807448\t4807275\t4942476\t5093511\t5251747\t5337861\t5280452\t5102121\t4866249\t4618672\t4369500\t4180146\t4025809\t3927404\t3883428\t3845222\t3784112\t3692156\t3535681\t3342052\t3125287\t2920964\t2747632\t2605169\t2503837\t2426528\t2359716\t2293181\t2210181\t2112386\t1991794\t1862812\t1737729\t1626855\t1533994\t1457871\t1393637\t1335565\t1279422\t1218407\t1153368\t1079865\t1001081\t928893\t869346\t815363\t769623\t731875\t698264\t663041\t627571\t589083\t548183\t511134\t471688\t440766\t414394\t393017\t373734\t355349\t336264\t316090\t296478\t277860\t259035\t242199\t228845\t218479\t210086\t203239\t195220\t187835\t180011\t170479\t162129\t154852\t148793\t143493\t139841\t137963\t136076\t133757\t131976\t127443\t124357\t120992\t118724\t116128\t113638\t113318\t113488\t113675\t114670\t113689\t112041\t110301\t108996\t107304\t108678\t110784\t113078\t116760\t118593\t119548\t119092\t118354\t118597\t120223\t123815\t129695\t138278\t146412\t151249\t154679\t157071\t161598\t165072\t173701\t181123\t191404\t201028\t211349\t223063\t232806\t242765\t253742\t265979\t279709\t293941\t309404\t321588\t336015\t352417\t372779\t397405\t420224\t444088\t467697\t485427\t505228\t524222\t543567\t564713\t589524\t612308\t634552\t655677\t672748\t694940\t711865\t728318\t744416\t758162\t772643\t784611\t795876\t805058\t819367\t831780\t852025\t867335\t880942\t887473\t888269\t888466\t884705\t877096\t881300\t882430\t887501\t891398\t898619\t900144\t902724\t895990\t887260\t881698\t877316\t870437\t864424\t862531\t858614\t853324\t847597\t840222\t831752\t825770\t816849\t807525\t798329\t786683\t776025\t764000\t752481\t738732\t728761\t718588\t709180\t697138\t684448\t671695\t654386\t639044\t623541\t605318\t592527\t578927\t564368\t551573\t536755\t522702\t505732\t490669\t473949\t456729\t441981\t425590\t412417\t398670\t385557\t372354\t359945\t346197\t332294\t316672\t304260\t290143\t278031\t267162\t257747\t248636\t237829\t227758\t218108\t206725\t196995\t188045\t179304\t171430\t164078\t158122\t151188\t144489\t138024\t131575\t125666\t120224\t113458\t108530\t105325\t101079\t96791\t93239\t89881\t86450\t82261\t79483\t76468\t73995\t71597\t69977\t68341\t66473\t65106\t63944\t61598\t60183\t58885\t58340\t57076\t56619\t56799\t56685\t56796\t57454\t57115\t57058\t56681\t57013\t58027\t58505\t59890\t60956\t62578\t63246\t63581\t64345\t64471\t65462\t67109\t68155\t70516\t72750\t74605\t75366\t75961\t76225\t77057\t77890\t79868\t81988\t83977\t86209\t87356\t88468\t88947\t89754\t90497\t92055\t94197\t94534\t95761\t96910\t97402\t98476\t99062\t99304\t99907\t100373\t101130\t100807\t100553\t101492\t101201\t102293\t102769\t103688\t104060\t105120\t106411\t106220\t106440\t106574\t107745\t108274\t109072\t109621\t109802\t112115\t111600\t112407\t113210\t113899\t114510\t115296\t116538\t116327\t118038\t118151\t119011\t119518\t120991\t121617\t122817\t122354\t122671\t123541\t123483\t124461\t124495\t124963\t124705\t125193\t125787\t125525\t124858\t125477\t125201\t124888\t124115\t124020\t123469\t123380\t122560\t122557\t121727\t119976\t119488\t119077\t117488\t117011\t114790\t115572\t114063\t112615\t111721\t110446\t108192\t107864\t106438\t104534\t102734\t101588\t100171\t98522\t97239\t95200\t93832\t92476\t90371\t88265\t86599\t85110\t82945\t81412\t79751\t78396\t76392\t74704\t73285\t70835\t70021\t67743\t66192\t64566\t63084\t61585\t59335\t58381\t56461\t55381\t53862\t52173\t50884\t49323\t48254\t47120\t45524\t44636\t42983\t41661\t40652\t39481\t38181\t37262\t36640\t35745\t34634\t33722\t32660\t32032\t31658\t30489\t29984\t29279\t28498\t27986\t27509\t26843\t26715\t25668\t25362\t25071\t24552\t24315\t23915\t23769\t23521\t23253\t22765\t22181\t22459\t21996\t22080\t21860\t21859\t21608\t21325\t21397\t20883\t20941\t20704\t20714\t20779\t20547\t20467\t20550\t20165\t20056\t20309\t19940\t19876\t19619\t19978\t19838\t19939\t19750\t19830\t19690\t19636\t19502\t19162\t19392\t19237\t19481\t19534\t19598\t19699\t19620\t19555\t19555\t19574\t19610\t19848\t19741\t20005\t19974\t20067\t20169\t20404\t20310\t20279\t20290\t20415\t20341\t20665\t20386\t20514\t20833\t20869\t20690\t20862\t21196\t21444\t21257\t21506\t21265\t21181\t21747\t21908\t21956\t21868\t21966\t22535\t22278\t22334\t22575\t22795\t22827\t22478\t23034\t22833\t22996\t22897\t22905\t22907\t23142\t22985\t23291\t22822\t23220\t23241\t23170\t23137\t23006\t23012\t22924\t22979\t22823\t22784\t22805\t22706\t22341\t22553\t22336\t22075\t21945\t22223\t21719\t21975\t21614\t21474\t21153\t21034\t20962\t20817\t20298\t20090\t20283\t19980\t19593\t19366\t19129\t18639\t18495\t18305\t18152\t17789\t17563\t17280\t17094\t16638\t16461\t16446\t16039\t15424\t15457\t15261\t14886\t14386\t13864\t13911\t13635\t12783\t12578\t12353\t12019\t11028\t10731\t10858\t10343\t9160\t9031\t9073\t8592\t7115\t6955\t6977\t6428\t4909\t4855\t4931\t4588\t2874\t2760\t2828\t2741\t1224\t1122\t1229\n"
#    },
#    "P01_GCparagon": {
#       "original": "P01\tGCparagon\toriginal\t4616\t4611\t4932\t4758\t5275\t4896\t5674\t4753\t5502\t4872\t5227\t4903\t4913\t4963\t5231\t5376\t5524\t5328\t5748\t6023\t5950\t5562\t5818\t6131\t6444\t6301\t6379\t6485\t6988\t7518\t7926\t7920\t7934\t8354\t8938\t9187\t9816\t9957\t10158\t10793\t11097\t10930\t11051\t11673\t12116\t12372\t12513\t12840\t13567\t14440\t15082\t14746\t14715\t15285\t16331\t16585\t19179\t19971\t21742\t24518\t26654\t26379\t24321\t23490\t24333\t25312\t26647\t27649\t29777\t33440\t37728\t40874\t41766\t40245\t38762\t39273\t40444\t43300\t47060\t51001\t54624\t58387\t59619\t62149\t63105\t63916\t67324\t71524\t75660\t81421\t89024\t96137\t98035\t96203\t96230\t100230\t107750\t118216\t128906\t144236\t166818\t197560\t224882\t231608\t220554\t219550\t235562\t264367\t296473\t336895\t395115\t478777\t589621\t704460\t764838\t757818\t772141\t844528\t972234\t1143144\t1316490\t1461171\t1563616\t1652536\t1765593\t1885461\t1994006\t2120418\t2284371\t2455364\t2627683\t2800482\t2953796\t3052280\t3115820\t3285803\t3573453\t3960406\t4399827\t4821312\t5172721\t5495299\t5912906\t6626884\t7651868\t8560028\t8963120\t8982245\t8799136\t8528212\t8181381\t7847097\t7622234\t7608614\t7654589\t7676997\t7555935\t7234572\t6757682\t6220587\t5710959\t5251366\t4881984\t4610424\t4407448\t4276371\t4148452\t4008393\t3826308\t3589350\t3320169\t3049833\t2805059\t2603013\t2441975\t2326682\t2234493\t2155052\t2072700\t1979600\t1873266\t1751086\t1629089\t1519229\t1419160\t1339595\t1274326\t1218887\t1170453\t1124635\t1074634\t1019183\t959375\t899557\t840443\t794858\t753861\t723093\t700802\t684687\t666227\t646200\t620846\t594690\t569761\t544808\t523249\t510247\t503166\t502057\t499050\t496805\t490078\t479297\t468598\t458381\t449811\t444623\t442310\t443815\t447623\t450165\t449329\t446716\t444208\t443757\t443314\t440999\t438498\t436986\t440423\t448336\t454925\t460903\t464768\t465839\t467737\t465172\t463790\t461941\t463858\t472314\t485415\t495686\t502231\t504204\t503581\t500068\t499580\t504732\t516462\t531030\t552327\t567141\t577945\t588563\t594738\t604545\t618042\t636460\t655265\t674434\t690305\t700166\t706636\t711283\t726401\t746741\t779665\t809425\t839989\t864504\t882214\t900339\t914891\t934329\t965456\t1003425\t1051217\t1093601\t1115958\t1120950\t1127005\t1137423\t1156596\t1188535\t1216851\t1243624\t1265741\t1271587\t1270408\t1265547\t1253415\t1245037\t1240525\t1233686\t1216865\t1201411\t1183832\t1165280\t1148669\t1128107\t1105917\t1078701\t1055401\t1024485\t997456\t965823\t937598\t917388\t899619\t886311\t869127\t848475\t826625\t804795\t775001\t746236\t722286\t700663\t684369\t670515\t658409\t644885\t631070\t614723\t598081\t581634\t566453\t551589\t537428\t523699\t510976\t501527\t489422\t478484\t466577\t455244\t445187\t433169\t421343\t410533\t399185\t387158\t377289\t366103\t357870\t348871\t341343\t331690\t323576\t314257\t303178\t294287\t284307\t277209\t266680\t260114\t253599\t246373\t239821\t234026\t226169\t219663\t212321\t205721\t199008\t193241\t187702\t182120\t177998\t173377\t169454\t164623\t159576\t156184\t150198\t145949\t142215\t137993\t135474\t132915\t130632\t127597\t124987\t122858\t120208\t117394\t114812\t112764\t110732\t109447\t108551\t106805\t105424\t103954\t102751\t101308\t100025\t99319\t99013\t98126\t97950\t97491\t97100\t96036\t95685\t95151\t94733\t94299\t94791\t94497\t94880\t95240\t94992\t93362\t94361\t94572\t94252\t93359\t93091\t93224\t93093\t94366\t93663\t94500\t94240\t93442\t93345\t93338\t92535\t91871\t93217\t93269\t92885\t92178\t91603\t90370\t88921\t87799\t87269\t87344\t87460\t84341\t84699\t84406\t82102\t78336\t77133\t77174\t76477\t72506\t73818\t74471\t72060\t65915\t65999\t66156\t64599\t57716\t57429\t58330\t55843\t47837\t46993\t47681\t46306\t38412\t38132\t38952\t36934\t31229\t30576\t30886\t28984\t23774\t23944\t24245\t23326\t19557\t19189\t19515\t18776\t15777\t15265\t15871\t15179\t11957\t11856\t12483\t11670\t7208\t7227\t8109\t6981\t3257\t3197\t3531\t3003\t1082\t1099\t1215\t1006\t388\t359\t411\t370\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n"
#    }
# }
