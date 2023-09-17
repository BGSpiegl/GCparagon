#!/usr/bin/env python3
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

target_aln_num = None
# -> process entire BAMs!

########################################################################################################################
# TODO: SET YOUR GCparagon SOURCE DIRECTORY AND TEMPORARY DIRECTORY HERE !!!!!!!!
MY_GRIFFIN_BAMS_OUTPUT_DIRECTORY = Path('your_absolute_path_to_an_output_directory_with_huge_space_from_last_step')
MY_TEMPORARY_DATA_DIR = Path('your_absolute_path_to_a_temporary_directory_with_huge_space')
# TODO: SET YOUR GCparagon SOURCE DIRECTORY AND TEMPORARY DIRECTORY HERE !!!!!!!!
########################################################################################################################
MY_GRIFFIN_BAMS_OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
MY_TEMPORARY_DATA_DIR.mkdir(parents=True, exist_ok=True)

SCRIPT_PARENT_PATH = Path(__file__).parent
SOURCE_CODE_ROOT_PATH = SCRIPT_PARENT_PATH.parent.parent.parent.parent
SOURCE_CODE_ROOT_DIR = str(SOURCE_CODE_ROOT_PATH)

OUTPUT_DIR = SOURCE_CODE_ROOT_PATH / 'validation/genomewide_GC_content_per_fragment_corrected'
# TODO: carry out preset computation first! Use 'driver_scripts/drv_compute_presets.sh'

min_unclipped_aln_fraction = DEFAULT_MIN_UNCLIPPED_ALN_FRACTION  # directly import from main module

hg38_2bit_ref_genome = SOURCE_CODE_ROOT_PATH / 'src/GCparagon/2bit_reference/hg38.analysisSet.2bit'
# TODO: above is valid if the EXECUTE_reference_download.sh is used; REPLACE THE PATH ABOVE OTHERWISE!

assert hg38_2bit_ref_genome.is_file()

MIN_FLENGTH = 20
MAX_FLENGTH = 800

# ---------------------- SELECT CORRECTION ALGORITHM ----------------------
CORRECTION_ALGORITHM = 'Griffin'  # must be from ('Griffin', 'GCparagon')
print(f"INFO - using analysis definitions for GC bias correction algorithm '{CORRECTION_ALGORITHM}'")
# -------------------------------------------------------------------------
# --------------- DEFINE CORRECT preset, if GCparagon correction fidelity should be computed --------------
PRESET = 2
GCPARAGON_PRESET_RESULTS_PATH = SOURCE_CODE_ROOT_PATH / f'preset_computation/preset{PRESET}'
# ---------------------------------------------------------------------------------------------------------
output_path_gcparagon = OUTPUT_DIR / f'original_and_corrected_gc_distribution-GCparagon-preset{PRESET}'
output_path_griffin = OUTPUT_DIR / 'original_and_corrected_gc_distribution-Griffin'
output_paths = {'GCparagon': output_path_gcparagon,
                'Griffin': output_path_griffin}
# TODO: SAMPLE STATS TXT FILES SHOULD BE COMBINED AFTERWARDS IN ONE FILE!
# Griffin correction weights-tagged BAM files
input_path_griffin = MY_GRIFFIN_BAMS_OUTPUT_DIRECTORY
tagged_bams_griffin = (input_path_griffin / 'C01/C01.GCtagged.bam',
                       input_path_griffin / 'H01/H01.GCtagged.bam',
                       input_path_griffin / 'B01/B01.GCtagged.bam',
                       input_path_griffin / 'P01/P01.GCtagged.bam')

# GCparagon - PRESET 2
# -------------------------------------------------------------------------
input_path_gcparagon = GCPARAGON_PRESET_RESULTS_PATH
tagged_bams_gcparagon = (input_path_gcparagon / 'C01/C01.GCtagged.bam',  # took ~ 16 mins on NVMe SSD
                         input_path_gcparagon / 'H01/H01.GCtagged.bam',
                         input_path_gcparagon / 'B01/B01.GCtagged.bam',
                         input_path_gcparagon / 'P01/P01.GCtagged.bam',)
tagged_bam_lists = {'Griffin': tagged_bams_griffin,
                    'GCparagon': tagged_bams_gcparagon}
correction_weight_tags = {'Griffin': 'GG',
                          'GCparagon': 'GC'}  # HIGHLY IMPORTANT !!! use 'GC' for GCparagon and 'GG' for Griffin


# select correct definition set
USE_TAG = correction_weight_tags[CORRECTION_ALGORITHM]
output_path = output_paths[CORRECTION_ALGORITHM]
output_path.mkdir(parents=True, exist_ok=True)
tagged_bams = tagged_bam_lists[CORRECTION_ALGORITHM]


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
                if aln.template_length <= aln.query_alignment_length:  # template_length >0; max. 25% clipped
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
            for opt in options:
                if aln.has_tag(opt):
                    tag_counts[opt] += 1
                    continue
            tag_counts['unknown'] += 1
    return sorted(list(tag_counts.items()), key=lambda x: x[1], reverse=True)[0][0]  # returns tag with highest count


def main() -> Optional[int]:
    # iterate over tuples
    sample_lines_dict = {}
    sample_flength_lines_dict = {}
    percent_part = '\t'.join([str(gc_pc) for gc_pc in range(0, 101, 1)])
    hdr_line = f'sample\talgorithm\tstatus\t{percent_part}\n'
    length_part = '\t'.join([str(length) for length in range(MIN_FLENGTH, MAX_FLENGTH + 1, 1)])
    hdr_line_lengths = f'sample\talgorithm\tstatus\t{length_part}\n'
    timestamp_str = strftime('%d-%m-%Y', localtime())
    output_table_path = \
        Path(output_path /
             (f"{CORRECTION_ALGORITHM.upper()}-STATISTICS_allSamples_GC-percentage_perFragmentSequence"
              f"{('_' + str(target_aln_num // 10 ** 6) + 'Mreads') if target_aln_num is not None else ''}_"
              f"{timestamp_str}.tsv"))
    if output_table_path.is_file():
        output_table_path.unlink()  # delete if exists; will be opened as append otherwise!
    output_length_table_path = \
        Path(output_path /
             (f"{CORRECTION_ALGORITHM.upper()}-STATISTICS_allSamples_fragment_length"
              f"{('_' + str(target_aln_num // 10 ** 6) + 'Mreads') if target_aln_num is not None else ''}_"
              f"{timestamp_str}.tsv"))
    if output_length_table_path.is_file():
        output_length_table_path.unlink()  # delete if exists; will be opened as append otherwise!
    for bam_idx, bam_path in enumerate(tagged_bams):
        print(f"INFO - processing BAM file '{bam_path.name}'")
        std_chroms, std_chrom_lengths = get_standard_chromosomes_from_bam(
            bam_path=str(bam_path), remove_scaffolds=['chrY', 'chrM', 'chrMT', 'chrEBV'])
        fragment_workers = []
        receivers = []
        # GC tag for GCparagon tagged files; additional GG for Griffin results-tagged file
        inferred_gc_tag = infer_gc_tag(preferred=USE_TAG, bam_path=bam_path)
        if inferred_gc_tag == 'unknown':
            log(message=f"could not infer GC-tag. Code will likely fail with error.", logger_name=LOGGER,
                log_level=logging.WARNING)
        else:
            log(message=f"inferred GC-tag: '{inferred_gc_tag}'", logger_name=LOGGER, log_level=logging.INFO)
        sample_id_path = bam_path.parent
        sample_id = str(sample_id_path.stem.split('.')[0])  # just 'P01' etc. not including '.GCtagged'
        sample_data_id = f'{sample_id}_{CORRECTION_ALGORITHM}'
        sample_lines_dict[sample_data_id] = {'original': f'{sample_id}\t{CORRECTION_ALGORITHM}\toriginal\t',
                                             'corrected': f'{sample_id}\t{CORRECTION_ALGORITHM}\tcorrected\t'}
        sample_flength_lines_dict[sample_data_id] = {'original': f'{sample_id}\t{CORRECTION_ALGORITHM}\toriginal\t',
                                                     'corrected': f'{sample_id}\t{CORRECTION_ALGORITHM}\tcorrected\t'}
        # create GC stats using multiprocessing (for first in pair )
        fragment_stats = {'original': [], 'corrected': []}
        fragment_lengths = {'original': [], 'corrected': []}
        for scaffold_to_process in std_chroms:
            receiver, sender = mp.Pipe(duplex=False)
            receivers.append(receiver)
            fragment_workers.append(mp.Process(target=scaffold_gc_stats_counter,
                                               kwargs={'bam_path': bam_path,
                                                       'tag': inferred_gc_tag,
                                                       'sender_connection': sender,
                                                       'scaffold_to_process': scaffold_to_process,
                                                       'sample_data_id': sample_data_id}))
        for read_worker in fragment_workers:
            read_worker.start()
        for rcv in receivers:
            # vars/vals used in worker: True, gc_counter_pre, gc_counter_corr, none_type_bases, tag_not_present
            success, flength_dist_pre, flength_dist_post, \
                stats_pre, stats_corr, none_bases, no_tag, sample_data_id = rcv.recv()
            # report progress
            log(message=f"received stats for sample data '{sample_data_id}'",
                log_level=logging.INFO, logger_name=LOGGER, flush=True)
            if not success:
                log(message=f'NOT ENOUGH TAGGED ALNS. TERMINATING..', log_level=logging.CRITICAL, logger_name=LOGGER,
                    close_handlers=True, flush=True)
                for read_worker in fragment_workers:
                    read_worker.kill()
                return 1
            if none_bases:
                log(message=f"WARNING: unexpected behavior - there were {none_bases:,} None-base alignments in "
                            f"sample data '{sample_data_id}! NOT EXPECTED",
                    log_level=logging.WARNING, logger_name=LOGGER, flush=True)
            if no_tag:
                log(message=f"WARNING: unexpected behavior - there were {no_tag:,} alignments with no identifiable "
                            f"GC-Tags (GC or YC) in file of sample '{sample_data_id}!",
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
            sample_lines_dict[sample_data_id][cur_status] += '\n'  # each sample has 4 lines -> 2 reads x 2 stats
        # output accumulated statistics (append if existing at this point))
        with open(output_table_path, 'at') as f_tab:
            if bam_idx == 0:
                f_tab.write(hdr_line)
                f_tab.flush()
            f_tab.writelines([sample_lines_dict[sample_data_id][cur_status]
                              for cur_status in ('original', 'corrected')])
            f_tab.flush()  # continuously write stats per sample to file
        # print result also to cmd
        print(f"GC content data for sample '{sample_id}' written to stats TSV file '{output_table_path.name}':\n" +
              '----' * 10 + '\n' + dumps(sample_lines_dict, sort_keys=True, indent=3))
        # accumulate fragment length counts
        sample_flength_lines_dict[sample_data_id] = {'original': f'{sample_id}\t{CORRECTION_ALGORITHM}\toriginal\t',
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


if __name__ == '__main__':
    sys.exit(main())


# INFO - using analysis definitions for GC bias correction algorithm 'GCparagon'
# INFO - processing BAM file 'C01.GCtagged.bam'
# | 2023-08-14 15:23:54,048 - GCPARAGON_validation - INFO - inferred GC-tag: 'GC'
# | 2023-08-14 15:39:44,449 - GCPARAGON_validation - INFO - received stats for sample data 'C01_GCparagon'
# ...
# INFO - accumulating statistics for sample 'P01'
# GC content data for sample 'P01' written to stats TSV file 'GCPARAGON-STATISTICS_allSamples_GC-percentage_perFragmentSequence_14-08-2023.tsv':
# ----------------------------------------
# {
#    "B01_GCparagon": {
#       "corrected": "B01\tGCparagon\tcorrected\t285\t99\t113\t107\t200\t231\t403\t679\t1195\t2055\t4349\t8531\t13706\t20058\t24974\t36335\t52867\t78703\t131033\t239787\t407670\t648284\t932099\t1322762\t1741795\t2480958\t3097413\t3921774\t4962588\t5486897\t6501807\t7124909\t7795988\t8183937\t9508776\t9680478\t10290337\t10036503\t10769066\t10186880\t9958197\t9049314\t8384156\t7767629\t7105342\t6777791\t6581607\t6074184\t6073236\t5990025\t5190896\t5408990\t4759248\t4217204\t3671032\t2946661\t2364960\t1921646\t1611882\t1305416\t1125465\t924341\t750488\t685779\t565795\t467704\t406097\t299562\t266359\t215309\t181450\t142084\t117960\t91048\t68469\t54272\t40942\t34929\t28385\t21808\t17222\t9850\t4764\t2419\t1319\t727\t346\t180\t69\t34\t10\t10\t5\t4\t0\t3\t17\t65\t76\t1\t4\n",
#       "original": "B01\tGCparagon\toriginal\t285\t99\t113\t107\t200\t231\t403\t679\t1158\t1992\t3309\t5295\t7550\t12824\t18916\t30773\t49839\t82322\t146033\t277603\t482295\t775842\t1138328\t1639861\t2180734\t3123772\t3898125\t4920979\t6201150\t6781688\t7964465\t8624369\t9312134\t9630747\t10996264\t10982276\t11447967\t10953809\t11546315\t10750214\t10333426\t9244136\t8453315\t7729137\t6982718\t6570434\t6254600\t5645380\t5468283\t5170309\t4285208\t4255176\t3580004\t3051652\t2557137\t1989014\t1545779\t1218040\t992148\t782193\t653437\t521637\t412361\t363177\t286226\t225310\t183847\t126761\t102964\t75464\t57497\t42233\t33631\t25680\t19854\t15801\t12559\t10856\t9110\t7543\t6157\t4447\t2969\t2091\t1303\t733\t347\t181\t71\t34\t10\t11\t5\t6\t3\t6\t19\t68\t77\t18\t57\n"
#    },
#    "C01_GCparagon": {
#       "corrected": "C01\tGCparagon\tcorrected\t13\t9\t4\t7\t4\t5\t6\t6\t13\t22\t26\t95\t195\t465\t1231\t5601\t25483\t66443\t123243\t200363\t338709\t536550\t835845\t1270260\t1742683\t2505163\t3213149\t4065163\t5209870\t5750782\t6907746\t7601141\t8395724\t8772348\t10429152\t10608355\t11497980\t11358975\t12434493\t11916980\t11742821\t10666371\t9947666\t9278553\t8674307\t8349700\t8238621\t7645850\t7660728\t7744545\t6692441\t7312500\t6560320\t5917586\t5249847\t4257242\t3480877\t2852493\t2441392\t2010488\t1774224\t1483155\t1230715\t1141447\t953263\t796706\t705432\t518579\t458270\t368161\t309340\t242055\t209580\t169654\t139285\t121985\t101472\t93967\t87140\t76854\t68152\t54038\t42206\t29980\t15005\t4265\t1467\t546\t189\t84\t29\t18\t8\t2\t2\t3\t9\t44\t37\t1\t1\n",
#       "original": "C01\tGCparagon\toriginal\t13\t9\t4\t7\t4\t5\t6\t6\t13\t22\t26\t95\t195\t465\t1154\t2822\t6825\t15519\t34049\t67800\t130776\t226793\t378735\t609982\t877332\t1318193\t1771959\t2344063\t3139646\t3608668\t4521353\t5173885\t5952248\t6459240\t7990153\t8442539\t9499475\t9717954\t11015514\t10890866\t11046625\t10317936\t9892251\t9496761\t9129367\t9026230\t9144000\t8700617\t8948395\t9255960\t8228757\t9203941\t8478746\t7824526\t7055110\t5765055\t4731788\t3907543\t3376077\t2811085\t2516859\t2133940\t1787700\t1674723\t1416879\t1184850\t1039639\t752891\t650764\t502048\t404336\t304800\t255400\t200545\t159395\t135705\t109548\t98308\t88101\t74819\t63665\t45409\t29651\t16678\t7768\t3446\t1451\t551\t193\t85\t29\t21\t8\t4\t2\t4\t9\t46\t38\t10\t15\n"
#    },
#    "H01_GCparagon": {
#       "corrected": "H01\tGCparagon\tcorrected\t159\t24\t24\t58\t65\t46\t81\t100\t146\t158\t176\t210\t239\t431\t579\t1020\t2048\t6133\t28944\t92197\t228869\t428010\t728077\t1148705\t1675405\t2427690\t3322324\t4368934\t5677213\t6601152\t7954292\t9035160\t10049544\t10721752\t12443949\t12957313\t13965143\t14290391\t15557653\t15322753\t14527307\t13232596\t12158905\t11308231\t10623154\t10157563\t9789986\t9153910\t8756628\t8548352\t7359474\t7476485\t6540265\t5845024\t5049671\t4116158\t3359441\t2761072\t2333855\t1952575\t1679510\t1434987\t1187794\t1087035\t906719\t766455\t667913\t512829\t444532\t360873\t303231\t247145\t212246\t175062\t144449\t124209\t108181\t96090\t90595\t79346\t66001\t50525\t36312\t25237\t15774\t8131\t4010\t2136\t1000\t340\t183\t90\t56\t27\t14\t18\t13\t34\t36\t1\t8\n",
#       "original": "H01\tGCparagon\toriginal\t159\t24\t24\t58\t65\t46\t81\t100\t146\t158\t176\t210\t239\t431\t579\t1020\t2034\t4770\t13049\t34294\t84481\t181248\t358349\t643303\t1044519\t1647826\t2407633\t3335407\t4541816\t5486801\t6831467\t8000470\t9121326\t9954613\t11805217\t12502157\t13703826\t14247685\t15705560\t15621958\t14937738\t13682423\t12623700\t11824436\t11152301\t10702933\t10350384\t9703155\t9312687\t9084252\t7840938\t7957161\t6995128\t6286575\t5456217\t4469988\t3658304\t2999150\t2537921\t2134874\t1843619\t1574387\t1304530\t1197298\t996517\t836803\t720575\t540797\t452625\t352577\t282373\t222345\t187033\t151396\t122940\t103385\t88173\t77233\t71146\t60608\t49576\t37539\t26147\t18104\t11991\t7182\t3983\t2171\t1017\t348\t192\t95\t58\t32\t16\t19\t15\t41\t42\t21\t29\n"
#    },
#    "P01_GCparagon": {
#       "corrected": "P01\tGCparagon\tcorrected\t757\t65\t98\t246\t258\t245\t389\t468\t648\t656\t946\t1001\t1247\t2045\t2723\t4054\t7697\t21502\t68283\t162683\t328207\t588474\t999172\t1595551\t2332626\t3390073\t4537217\t5969237\t7719250\t8930976\t10674196\t12130793\t13379156\t14366701\t16633982\t17398082\t18669700\t18980885\t20635062\t20220384\t19519712\t18228822\t16934963\t16087316\t15171346\t14753778\t14361462\t13606524\t13232637\t13228767\t11572266\t12072640\t10685839\t9731506\t8451127\t6999737\t5730289\t4775392\t4066233\t3446118\t2987176\t2579471\t2140519\t1993095\t1655710\t1413987\t1229912\t938357\t808292\t652761\t540377\t437344\t374081\t313582\t267875\t237386\t209631\t192558\t182341\t164314\t142031\t109108\t77743\t53195\t35021\t19468\t9026\t4137\t1769\t682\t376\t223\t115\t89\t57\t58\t182\t307\t166\t34\t30\n",
#       "original": "P01\tGCparagon\toriginal\t757\t65\t98\t246\t258\t245\t389\t468\t648\t656\t946\t1001\t1246\t2045\t2723\t3996\t6171\t10365\t19013\t35355\t68204\t127497\t224121\t366532\t577393\t926922\t1361333\t1963231\t2774366\t3509221\t4556137\t5605698\t6701224\t7751334\t9638561\t10826181\t12428248\t13500794\t15612302\t16296681\t16701661\t16499107\t16198112\t16210017\t16033548\t16339220\t16577409\t16340501\t16444248\t16991440\t15286938\t16444389\t14982383\t14053702\t12598497\t10750865\t9048643\t7743061\t6747162\t5852843\t5180650\t4544603\t3831440\t3606313\t3017878\t2571299\t2212780\t1678219\t1406024\t1097453\t874778\t691469\t573649\t464632\t380240\t317127\t264247\t230544\t209891\t183902\t156128\t117401\t78163\t48251\t28724\t16204\t8676\t4232\t1841\t730\t393\t234\t123\t96\t61\t72\t206\t333\t209\t104\t151\n"
#    }
# }
