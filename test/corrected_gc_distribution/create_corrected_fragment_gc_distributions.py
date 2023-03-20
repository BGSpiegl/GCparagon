#!/usr/bin/env python3
# tested with GCparagon_oy3.10 conda environment
import os
import sys
import gzip
import pysam
import pathlib
import logging
import numpy as np
from json import dumps
import multiprocessing as mp
from collections import deque
from twobitreader import TwoBitFile
from collections import defaultdict
from time import localtime, strftime
from typing import Union, Optional, Tuple, List

target_alignment_number = None  # 40*10**6  # possible problem that we base our estimate only on chr1
# -> process entire BAMs! # other choice tested: 2*10**7  # twenty million and 40*10**6

SOURCE_CODE_ROOT_PATH = pathlib.Path(__file__).parent.parent.parent
SOURCE_CODE_ROOT_DIR = str(SOURCE_CODE_ROOT_PATH)

hg38_2bit_ref_genome = SOURCE_CODE_ROOT_PATH / '2bit_reference/hg38.analysisSet.2bit'
# above is valid if the EXECUTE_reference_download.sh is used; REPLACE THE PATH ABOVE OTHERWISE!

output_path = SOURCE_CODE_ROOT_PATH / 'test/corrected_gc_distribution'
output_path.mkdir(parents=True, exist_ok=True)

input_path = SOURCE_CODE_ROOT_PATH / 'preset_computation'
bam_correction_matrices_tuples = (
    (input_path / '1/P01/P01.GCtagged.bam',
     input_path / '1/P01/P01_gc_weights_6simsMean.2IQRoutliersRemoved.5IgaussSmoothed.txt.gz'),
    (input_path / '2/P244_7/P244_7.GCtagged.bam',
     input_path / '2/P244_7/P244_7_gc_weights_4simsMean.2IQRoutliersRemoved.'
                  '2IgaussSmoothed.txt.gz'),
    (input_path / 'Preset3/P244_7/P244_7.GCtagged.bam',
     input_path / 'Preset3/P244_7/P244_7_gc_weights_4simsMean.2IQRoutliersRemoved.'
                  '2IgaussSmoothed.txt.gz'),
    (input_path / '1/B01/B01.GCtagged.RGpresent.bam',
     input_path / '1/B01/B58_3_gc_weights_6simsMean.2IQRoutliersRemoved.'
                  '5IgaussSmoothed.txt.gz'),
    (input_path / '1/C01/C01.GCtagged.bam',
     input_path / '1/C01/C219_5_gc_weights_6simsMean.2IQRoutliersRemoved.'
                  '5IgaussSmoothed.txt.gz'),
    (input_path / '1/H01/H01.GCtagged.bam',
     input_path / '1/H01/NPH_011_gc_weights_6simsMean.2IQRoutliersRemoved.'
                  '5IgaussSmoothed.txt.gz'),
    (input_path / '1/P01/P01.GCtagged.bam',
     input_path / '1/P01/P244_6_gc_weights_6simsMean.2IQRoutliersRemoved.'
                  '5IgaussSmoothed.txt.gz'),
    (input_path / '2/B01/B01.GCtagged.RGpresent.bam',
     input_path / '2/B01/B58_3_gc_weights_4simsMean.2IQRoutliersRemoved.'
                  '2IgaussSmoothed.txt.gz'),
    (input_path / '2/C01/C01.GCtagged.bam',
     input_path / '2/C01/C219_5_gc_weights_4simsMean.2IQRoutliersRemoved.'
                  '2IgaussSmoothed.txt.gz'),
    (input_path / '2/H01/H01.GCtagged.bam',
     input_path / '2/H01/NPH_011_gc_weights_4simsMean.2IQRoutliersRemoved.'
                  '2IgaussSmoothed.txt.gz'),
    (input_path / '2/P01/P01.GCtagged.bam',
     input_path / '2/P01/P244_6_gc_weights_4simsMean.2IQRoutliersRemoved.'
                  '2IgaussSmoothed.txt.gz'),
    (input_path / '3/B01/B01.GCtagged.RGpresent.bam',
     input_path / '3/B01/B58_3_gc_weights_4simsMean.2IQRoutliersRemoved.'
                  '2IgaussSmoothed.txt.gz'),
    (input_path / '3/C01/C01.GCtagged.bam',
     input_path / '3/C01/C219_5_gc_weights_4simsMean.2IQRoutliersRemoved.'
                  '2IgaussSmoothed.txt.gz'),
    (input_path / '3/H01/H01.GCtagged.bam',
     input_path / '3/H01/NPH_011_gc_weights_4simsMean.2IQRoutliersRemoved.'
                  '2IgaussSmoothed.txt.gz'),
    (input_path / '3/P01/P01.GCtagged.bam',
     input_path / '3/P01/P244_6_gc_weights_4simsMean.2IQRoutliersRemoved.'
                  '2IgaussSmoothed.txt.gz'))


def gib_cmd_logger() -> logging.Logger:
    any_logger = logging.getLogger(name='GCparagon_dev')
    any_logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter('| %(asctime)s - %(name)s - %(levelname)s - %(message)s')
    stdout_handler = logging.StreamHandler(stream=sys.stdout)
    stdout_handler.setLevel(logging.DEBUG)
    stdout_handler.setStream(sys.stdout)
    stdout_handler.setFormatter(formatter)
    stdout_handler.set_name('zeiln_hendla')
    any_logger.addHandler(stdout_handler)
    return any_logger


LOGGER = gib_cmd_logger()  # formatted cmdline output


def log(message: str, log_level: int, i_log_with: logging.Logger, flush=True, close_handlers=False):
    """

    :param message:
    :param log_level:
    :param i_log_with:
    :param flush:
    :param close_handlers:
    :return:
    """
    match log_level:
        case logging.NOTSET:
            i_log_with.info(message)
        case logging.DEBUG:
            i_log_with.debug(message)
        case logging.INFO:
            i_log_with.info(message)
        case logging.WARNING:
            i_log_with.warning(message)
        case logging.ERROR:
            i_log_with.error(message)
        case logging.CRITICAL:
            i_log_with.critical(message)
    if flush and isinstance(i_log_with, logging.Logger):  # do nothing if i_log_with is undefined
        for hdlr in i_log_with.handlers:
            hdlr.flush()
    if close_handlers:
        for hdlr in i_log_with.handlers:
            hdlr.flush()
            hdlr.close()


def reduce_matrix(matrix_to_trim: np.array, trim_dimensions_exclusively_containing: List[float],
                  border_elements: Optional[int]) -> Tuple[np.array, Tuple[range, range]]:
    """
    Remove non-informative weights so that numpy matrix is smaller (faster access?)
    :param border_elements:
    :param matrix_to_trim:
    :param trim_dimensions_exclusively_containing:
    :return:
    """
    if border_elements is None:  # should not occur but casting according to logic for safety
        border_elements = 0
    # trim leading lines
    delete_initial_rows = 0
    for line_idx in range(matrix_to_trim.shape[0]):
        if all(map(lambda v: v in trim_dimensions_exclusively_containing, matrix_to_trim[line_idx, :])):
            delete_initial_rows += 1
        else:
            break
    if delete_initial_rows > border_elements:  # don't remove lines otherwise
        delete_initial_rows -= border_elements
        matrix_to_trim = matrix_to_trim[delete_initial_rows:, :]
    else:
        delete_initial_rows = 0
    # trim tailing lines
    delete_tailing_rows = 0
    for inverse_line_idx in range(1, matrix_to_trim.shape[0] + 1):
        if all(map(lambda v: v in trim_dimensions_exclusively_containing, matrix_to_trim[-inverse_line_idx, :])):
            delete_tailing_rows += 1
        else:
            break
    if delete_tailing_rows > border_elements:  # don't remove lines otherwise
        delete_tailing_rows -= border_elements
        matrix_to_trim = matrix_to_trim[:-delete_tailing_rows, :]
    else:
        delete_tailing_rows = 0
    # trim leading columns
    delete_initial_columns = 0
    for column_idx in range(matrix_to_trim.shape[1]):
        if all(map(lambda v: v in trim_dimensions_exclusively_containing, matrix_to_trim[:, column_idx])):
            delete_initial_columns += 1
        else:
            break
    if delete_initial_columns > border_elements:  # don't remove columns otherwise
        delete_initial_columns -= border_elements
        matrix_to_trim = matrix_to_trim[:, delete_initial_columns:]
    else:
        delete_initial_columns = 0
    # trim tailing columns
    delete_tailing_columns = 0
    for column_idx in range(1, matrix_to_trim.shape[1] + 1):
        if all(map(lambda v: v in trim_dimensions_exclusively_containing, matrix_to_trim[:, -column_idx])):
            delete_tailing_columns += 1
        else:
            break
    if delete_tailing_columns > border_elements:  # don't remove columns otherwise
        delete_tailing_columns -= border_elements
        matrix_to_trim = matrix_to_trim[:, :-delete_tailing_columns]
    else:
        delete_tailing_columns = 0
    return matrix_to_trim, (range(delete_initial_rows, delete_tailing_rows),
                            range(delete_initial_columns, delete_tailing_columns))


def load_txt_to_matrix_with_meta(filename: Union[str, pathlib.Path], loading_logger: Optional[logging.Logger],
                                 to_dtype=np.float64) -> Tuple[np.array, range]:
    """

    :param loading_logger:
    :param filename:
    :param to_dtype:
    :return:
    """
    if loading_logger is not None:
        log(message=f"Loading statistic matrix from {filename}", log_level=logging.INFO, i_log_with=loading_logger)
    else:
        print(f"Loading statistic matrix from {filename}")
    statistic_matrix = np.loadtxt(filename, delimiter='|', skiprows=0, dtype=to_dtype)
    with gzip.open(str(filename), 'rt') as f_mat:
        hdr = f_mat.readline()
    elements = hdr.split('# rows representing fragment lengths (')[1].split()  # split on whitespace + trim empty fields
    fragment_lengths = range(int(elements[0]), int(elements[3]) + 1)
    return statistic_matrix, fragment_lengths


def gc_stats_counter(bam_path: Union[pathlib.Path, str], tag: str, correction_weights: np.array,
                     sender_connection: mp.Pipe, weight_starts: Tuple[int, int],
                     target_number_alns: Optional[int], sample_data_id: str, max_no_tags=100000):  # ASSERTS pe-reads!
    ref_genome_handle = TwoBitFile(hg38_2bit_ref_genome)
    gc_counter_pre = defaultdict(int)  # max. 101 entries
    gc_counter_corr = defaultdict(float)  # sum of GC-tags
    none_type_bases = 0
    tag_not_present = 0
    flen_offset, gc_offset = weight_starts  # not implemented!
    with pysam.AlignmentFile(bam_path, mode='rb', threads=2) as f_aln:  # process all alignments passing filters
        alignment_iter = filter(lambda p: p.is_proper_pair and 20 < p.template_length < 800 and not p.is_supplementary
                                and not p.is_secondary, f_aln.fetch(multiple_iterators=True, until_eof=True))
        chrom_handles = {}
        held_chrom_handles = deque()
        if target_number_alns is None:
            for aln in alignment_iter:
                frag_len = aln.template_length  # use this filter here because we want
                try:  # to get chromosome handle
                    cur_chrom_handle = chrom_handles[aln.reference_name]
                except KeyError:
                    chrom_handles[aln.reference_name] = ref_genome_handle[aln.reference_name]
                    held_chrom_handles.append(aln.reference_name)
                    cur_chrom_handle = chrom_handles[aln.reference_name]
                    if len(held_chrom_handles) == 3:
                        del chrom_handles[held_chrom_handles.popleft()]
                try:
                    if frag_len <= aln.query_alignment_length:  # template_length filtered to be >0
                        bases = aln.query_alignment_sequence
                    else:
                        start_pos = min(aln.reference_start, aln.next_reference_start)
                        end_pos = start_pos + frag_len  # is filtered to be positive
                        bases = cur_chrom_handle[start_pos:end_pos].upper()  # coords 0-based, end-open
                    gc_pc = int(round(100. * (bases.count('G') + bases.count('C')) / len(bases), ndigits=0))
                    gc_counter_pre[gc_pc] += 1
                    gc_counter_corr[gc_pc] += aln.get_tag(tag)  # returned value is cast into an appropriate python type
                except AttributeError:  # 'NoneType' object has no attribute 'upper'
                    none_type_bases += 1  # no sequence present -> faulty entry?
                except KeyError:
                    tag_not_present += 1
                    if tag_not_present >= max_no_tags:
                        log(message=f"too many (>={max_no_tags:,}) alignments encountered without identifiable "
                                    "GC-tags. Terminating..", log_level=logging.CRITICAL, i_log_with=LOGGER)
                        sender_connection.send((False, None, None, None, None))
                        sender_connection.close()
                        return
        else:
            alns_procd = 0
            aln_scores = []
            for aln in alignment_iter:
                frag_len = aln.template_length  # use this filter here because we want
                try:  # to get chromosome handle
                    cur_chrom_handle = chrom_handles[aln.reference_name]
                except KeyError:
                    chrom_handles[aln.reference_name] = ref_genome_handle[aln.reference_name]
                    held_chrom_handles.append(aln.reference_name)
                    cur_chrom_handle = chrom_handles[aln.reference_name]
                    if len(held_chrom_handles) == 3:
                        del chrom_handles[held_chrom_handles.popleft()]
                try:
                    if frag_len <= aln.query_alignment_length:  # template_length filtered to be >0
                        bases = aln.query_alignment_sequence
                    else:
                        start_pos = min(aln.reference_start, aln.next_reference_start)
                        end_pos = start_pos + frag_len  # is filtered to be positive
                        bases = cur_chrom_handle[start_pos:end_pos].upper()  # no off-by-1-error here; all 0-based, end-open
                    gc_pc = int(round(100. * (bases.count('G') + bases.count('C')) / len(bases), ndigits=0))
                    if gc_pc == 50:
                        aln_scores.append(aln.mapQ)
                        if len(aln_scores) == 50000:
                            print(f"mean alignment score: {np.mean(aln_scores)}, max: {max(aln_scores)}, min: {min(aln_scores)}")
                            sys.exit(1)
                    gc_counter_pre[gc_pc] += 1
                    gc_counter_corr[gc_pc] += aln.get_tag(tag)  # returned value is cast into an appropriate python type
                except AttributeError:  # 'NoneType' object has no attribute 'upper'
                    none_type_bases += 1  # no sequence present -> faulty entry?
                except KeyError:
                    tag_not_present += 1
                    if tag_not_present >= max_no_tags:
                        log(message=f"too many (>={max_no_tags:,}) alignments encountered without identifiable GC-tags. "
                                    f"Terminating..", log_level=logging.CRITICAL, i_log_with=LOGGER)
                        sender_connection.send((False, None, None, None, None))
                        sender_connection.close()
                        return
                alns_procd += 1
                if alns_procd >= target_number_alns:
                    log(message=f"target number of alignments (= {target_number_alns:,}) reached.",
                        log_level=logging.INFO, i_log_with=LOGGER)
                    break
                elif alns_procd % 10**6 == 0:
                    log(message=f"{sample_data_id}: {alns_procd:,} of {target_number_alns:,} alignments processed..",
                        log_level=logging.INFO, i_log_with=LOGGER)
    sender_connection.send((True, gc_counter_pre, gc_counter_corr, none_type_bases, tag_not_present, sample_data_id))
    sender_connection.close()


def infer_gc_tag(preferred: str, bam_path: Union[pathlib.Path, str]):
    options = ('YC',)
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
            for opt in options:
                if aln.has_tag(opt):
                    tag_counts[opt] += 1
                    continue
            tag_counts['unknown'] += 1
    return sorted(list(tag_counts.items(), key=lambda x: x[1], reverse=True))[0][0]


def main() -> Optional[int]:
    # iterate over tuples
    percent_part = '\t'.join([str(gc_pc) for gc_pc in range(0, 101, 1)])
    hdr_line = f'sample\tpreset\tstatus\t{percent_part}\n'
    timestamp_str = strftime('%d-%m-%Y', localtime())
    output_table_path = output_path / f"STATISTICS_allSamples_GC-percentage_perFragmentSequence" \
        f"{('_' + str(target_alignment_number // 10**6) + 'Mreads') if target_alignment_number is not None else ''}_" \
                                      f"{timestamp_str}.tsv"
    fragment_workers = []
    receivers = []
    sample_lines_dict = {}
    with open(output_table_path, 'wt') as f_tab:
        f_tab.write(hdr_line)
        f_tab.flush()
        for bam_path, correction_weights_path in bam_correction_matrices_tuples:
            inferred_gc_tag = infer_gc_tag(preferred='GC', bam_path=bam_path)
            if inferred_gc_tag == 'unknown':
                log(message=f"could not infer GC-tag. Code will likely fail with error.", i_log_with=LOGGER,
                    log_level=logging.WARNING)
            else:
                log(message=f"inferred GC-tag: '{inferred_gc_tag}'", i_log_with=LOGGER, log_level=logging.INFO)
            sample_id_path = bam_path.parent
            sample_id = str(sample_id_path.stem)
            preset = str(sample_id_path.parent.stem).lower().strip('preset')
            sample_data_id = f'{sample_id}_{preset}'
            sample_lines_dict[sample_data_id] = {'original': f'{sample_id}\t{preset}\toriginal\t',
                                                 'corrected': f'{sample_id}\t{preset}\tcorrected\t'}
            correction_matrix, orig_frag_len_range = load_txt_to_matrix_with_meta(
                filename=correction_weights_path, loading_logger=LOGGER)
            complete_mask_focused, (deleted_rows, deleted_columns) = reduce_matrix(
                matrix_to_trim=correction_matrix, trim_dimensions_exclusively_containing=[0., 1.],
                border_elements=None)
            weights_start_offsets = (orig_frag_len_range.start + deleted_rows.start, deleted_columns.start)
            # create GC stats using multiprocessing (for first in pair )
            receiver, sender = mp.Pipe(duplex=False)
            receivers.append(receiver)
            fragment_workers.append(mp.Process(target=gc_stats_counter,
                                               kwargs={'bam_path': bam_path,
                                                       'correction_weights': complete_mask_focused,
                                                       'tag': inferred_gc_tag,
                                                       'sender_connection': sender,
                                                       'weight_starts': weights_start_offsets,
                                                       'target_number_alns': target_alignment_number,
                                                       'sample_data_id': sample_data_id}))
        for read_worker in fragment_workers:
            read_worker.start()
        for rcv in receivers:
            # vars/vals used in worker: True, gc_counter_pre, gc_counter_corr, none_type_bases, tag_not_present
            success, stats_pre, stats_corr, none_bases, no_tag, sample_data_id = rcv.recv()
            # report progress
            log(message=f"received stats for sample data '{sample_data_id}'",
                log_level=logging.INFO, i_log_with=LOGGER, flush=True)
            if none_bases:
                log(message=f"WARNING: unexpected behavior - there were {none_bases:,} None-base alignments in "
                            f"sample data '{sample_data_id}! NOT EXPECTED",
                    log_level=logging.WARNING, i_log_with=LOGGER, flush=True)
            if no_tag:
                log(message=f"WARNING: unexpected behavior - there were {no_tag:,} alignments with no identifiable "
                            f"GC-Tags (GC or YC) in file of sample '{sample_data_id}!",
                    log_level=logging.WARNING, i_log_with=LOGGER, flush=True)
            if not success:
                log(message=f'NOT ENOUGH TAGGED ALNS. TERMINATING..', log_level=logging.CRITICAL, i_log_with=LOGGER,
                    close_handlers=True, flush=True)
                return 1
            # record output statistics
            fragment_stats = {'original': stats_pre, 'corrected': stats_corr}
            for cur_stat in ('original', 'corrected'):
                cur_read_stats = fragment_stats[cur_stat]
                sample_lines_dict[sample_data_id][cur_stat] += '\t'.join([str(int(cur_read_stats[gc_pc]))
                                                                          for gc_pc in range(0, 101, 1)])
                sample_lines_dict[sample_data_id][cur_stat] += '\n'  # each sample has 4 lines -> 2 reads x 2 stati
            f_tab.writelines([sample_lines_dict[sample_data_id][cur_stat]
                              for cur_stat in ('original', 'corrected')])
            f_tab.flush()  # continuously write stats per sample to file
        for worker in fragment_workers:
            worker.join()
        for worker in fragment_workers:
            worker.close()
    # print dully to cmd
    print(f"Sample data written to stats TSV file '{sample_lines_dict}':\n" + '----'*10 + '\n' +
          dumps(sample_lines_dict, sort_keys=True, indent=3))


if __name__ == '__main__':
    sys.exit(main())
