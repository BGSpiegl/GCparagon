#!/usr/bin/env python3

from typing import List, Tuple, Dict, Union
from os.path import join as pth_join, sep as pth_sep
from os import makedirs as os_mkdirs
from shutil import which as sh_which
from sys import exit as sys_exit
from pybedtools import BedTool
from natsort import humansorted
from pathlib import Path

CODE_ROOT_PATH = Path(__file__).parent.parent.parent

bedtools_path = sh_which('bedtools')

CHUNK_SIZE = 1000000  # 1 Mb
SHIFT_N_TIMES = 4
CHUNK_POSITION_OFFSETS = tuple([CHUNK_SIZE // SHIFT_N_TIMES * shift_idx
                                for shift_idx in range(SHIFT_N_TIMES)])

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
GENOME_BUILD = 'hg19'
# PATH DEFINITIONS
exclusion_list = CODE_ROOT_PATH / 'accessory_files/hg19_GCcorrection_ExclusionList.sorted.merged.bed'
genome_file_path = CODE_ROOT_PATH / 'accessory_files/hg19.genome_file.tsv'
output_path = CODE_ROOT_PATH / f'accessory_files/genomic_interval_preselection-shifted{SHIFT_N_TIMES}x_hg19'
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if not bedtools_path:
    print("ERROR - bedtools path not found. Cannot benchmark intersection times.")
    sys_exit(1)


def read_bed_file(bed_path: str) -> List[Tuple[str, int, int]]:
    with open(bed_path, 'rt') as f_bed:
        return [(chrom, int(start), int(stop))
                for chrom, start, stop, *_ in filter(lambda x: x != '', [bed_line.strip().split()
                                                                         for bed_line in f_bed.readlines()])]


def get_stdchrom_intervals(genome_file: Union[str, Path], interval_size=CHUNK_SIZE, offset=0) -> list:
    with open(genome_file, 'rt') as f_gen:
        whole_genome_regions = [(chrom, offset, int(stop))
                                for chrom, stop, *_ in filter(lambda x: x != '',
                                                              [genome_line.strip().split()
                                                               for genome_line in f_gen.readlines()])]
    sorted_std_chroms = humansorted([f'chr{i}' for i in range(1, 23, 1)] + ['chrx', 'chry'])
    std_regs = tuple(filter(lambda c: c is not None,
                            [(chrom, strt, stop) if chrom.lower() in sorted_std_chroms else None
                             for chrom, strt, stop in whole_genome_regions]))
    chroms_sorted = humansorted(list(set([chrom for chrom, _stt, _stp in std_regs])))
    std_intervals = []
    for sorted_chrom in chroms_sorted:
        for chrom, _strt, stop in std_regs:
            if chrom == sorted_chrom:
                n_splits = (stop - offset) // interval_size
                cum_size = offset
                for split_idx in range(0, n_splits, 1):
                    std_intervals.append((chrom, cum_size, cum_size + interval_size))
                    cum_size += interval_size
                std_intervals.append((chrom, cum_size, stop))
    return std_intervals


def get_overlap_statistics(overlaps_bed_path: str) -> Dict[Tuple[str, int, int], Dict[str, int]]:
    with open(overlaps_bed_path, 'rt') as f_overlaps:
        overlap_content = [line.strip().split('\t') for line in f_overlaps.readlines()]
    bad_interval_bases = {}
    for interval_cont, interval_start, interval_stop, _1, _2, _3, *rest in overlap_content:
        overlapping_bases = rest[-1]
        try:
            bad_interval_bases[(interval_cont, int(interval_start), int(interval_stop))].append(int(overlapping_bases))
        except KeyError:
            bad_interval_bases.update({(interval_cont, int(interval_start), int(interval_stop)):
                                       [int(overlapping_bases)]})
    # postprocess: compute sum of overlapping_bases per interval
    for interval in bad_interval_bases.keys():
        bad_interval_bases[interval] = {'bad_regions': len(bad_interval_bases[interval]),
                                        'bad_bases': sum(bad_interval_bases[interval])}
    return {'bad_intervals': bad_interval_bases}  # bad_interval_bases might be empty in cases of a bad regions BED file


def convert_interval_tuple_str(interval_str=None, interval_tuple=None) -> Union[Tuple[str, int, int], str]:
    if all([p is None for p in (interval_str, interval_tuple)]):
        raise AttributeError("at least one input parameter must be not None")
    if interval_str is not None and not isinstance(interval_str, str):
        raise AttributeError("provided interval_str was not a str!")
    if interval_tuple is not None and not isinstance(interval_str, tuple):
        raise AttributeError("provided interval_tuple was not a tuple!")
    if interval_str:
        interval_tuple = (interval_str.split('-')[0],
                       int(interval_str.split('-')[1].split(':')[0]),
                       int(interval_str.split(':')[1]))
        return interval_tuple
    interval_str = f"{interval_tuple[0]}:{interval_tuple[1]}-{interval_tuple[2]}"
    return interval_str


if __name__ == '__main__':
    # read original exclusion_list
    exclusion_list_content = read_bed_file(bed_path=exclusion_list)
    # find exclusion_list intervals larger than 1kb (or multiples of this)
    multiples = (0, 1, 2, 3, 4, 5)  # 0 represents all exclusion listed regions
    num_multiples = {}.fromkeys(multiples)
    sub_exclusion_lists = {}.fromkeys(multiples)
    for k in sub_exclusion_lists.keys():
        sub_exclusion_lists[k] = []
    mult_exclusion_list_beds = {}.fromkeys(multiples)
    sorted_std_chroms = humansorted([f'chr{i}' for i in range(1, 23, 1)] + ['chrx', 'chry'])
    for sorted_chrom in sorted_std_chroms:
        for chrom, start, stop in exclusion_list_content:  # [ DONE ]
            if chrom == sorted_chrom:
                sub_exclusion_lists[0].append(f"{chrom}\t{start}\t{stop}\n")  # always append to basic list
                reg_size = stop - start  # [ DONE. ]
                for mult in multiples:
                    if reg_size >= 1000 * mult:  # filter for regions
                        sub_exclusion_lists[mult].append(f"{chrom}\t{start}\t{stop}\n")
                        try:
                            num_multiples[mult] += 1
                        except TypeError:
                            num_multiples[mult] = 1
    for k in num_multiples.keys():
        if num_multiples[k] is None:
            num_multiples[k] = 0
    for mult in multiples:  # [ DONE ]
        print(f"number of exclusion listed regions larger or equal than {mult}kb: {num_multiples[mult]:,}")
    # write filtered exclusion lists
    for mult, sublist_lines in sub_exclusion_lists.items():
        mult_exclusion_list_beds[mult] = exclusion_list  # THIS was used -> original exclusion list <-> mult=0
    # intersect genomic intervals and exclusion listed regions
    output_path.mkdir(parents=True, exist_ok=True)
    # split up the genome using different position offsets
    for interval_offset in CHUNK_POSITION_OFFSETS:
        print(f"NOW: computing bad regions overlap for intervals with offset {interval_offset} ..")
        if interval_offset == 0:
            current_output_dir = output_path
        else:
            current_output_dir = output_path / f"{interval_offset//1000}kbp_intervalOffset"
            current_output_dir.mkdir(parents=True, exist_ok=True)
        # interval up the genome
        std_intervals = get_stdchrom_intervals(genome_file=genome_file_path, interval_size=CHUNK_SIZE, offset=interval_offset)
        # write whole genome BED file
        whole_genome_bed_path = pth_join(current_output_dir, f'{GENOME_BUILD}_{CHUNK_SIZE//1000}kbp_intervals.bed')
        with open(whole_genome_bed_path, 'wt') as f_whole_g:
            f_whole_g.writelines([f"{chrom}\t{start}\t{stop}\n" for chrom, start, stop in std_intervals])
        # compute overlap between each interval in whole genome BED and the merged, sorted exclusion list instance
        overlap_statistics = {}.fromkeys(multiples)
        for mult in multiples:
            out_overlap_bed = pth_join(current_output_dir,
                                       f'{GENOME_BUILD}_{CHUNK_SIZE//1000}kbpStdChunks_gte{mult}kb_EML-overlap.bed')
            overlap_statistics[mult] = {'overlaps_bed_path': out_overlap_bed}
            intersect_bedtool = BedTool(mult_exclusion_list_beds[mult])
            wg_bedtool = BedTool(whole_genome_bed_path)
            wg_bedtool.intersect(intersect_bedtool, wo=True).saveas(out_overlap_bed)
        for mult in multiples:
            overlap_statistics[mult].update(get_overlap_statistics(overlaps_bed_path=overlap_statistics[mult][
                'overlaps_bed_path']))
        # compute all bad intervals -> create complete table with columns being multiples
        unique_overlapping_intervals = set()
        for mult in multiples:
            unique_overlapping_intervals.update(overlap_statistics[mult]['bad_intervals'].keys())
        sorted_overlapping_intervals = humansorted(list(unique_overlapping_intervals),
                                                   key=lambda x: [x[0], x[1]])
        # write total stats table
        output_overlaps_table = pth_join(current_output_dir,
                                         f'{CHUNK_SIZE//1000}kbp_intervals_bad_regions_overlap_ELRminSizes.tsv')
        sorted_multiples_list = sorted(list(multiples))
        stats_out_lines = ['\t' + '\t\t'.join([f"{m}kb_ELR_min_size" for m in sorted_multiples_list]) +
                           '\t\tTOTAL\t\n' + 'interval\t' +
                           '\t'.join(['regions\tbases'for _m in sorted_multiples_list]) + '\tBASES\tREGIONS\n']
        stats_out_lines.extend([f'{b_interval[0]}:{b_interval[1]:,}-{b_interval[2]:,}\t' +
                                '\t'.join(
                                    ['0\t0'
                                     if overlap_statistics[mult]['bad_intervals'].get(b_interval) is None
                                     else (str(overlap_statistics[mult]['bad_intervals'][b_interval]['bad_regions']) +
                                           '\t' +
                                           str(overlap_statistics[mult]['bad_intervals'][b_interval]['bad_bases']))
                                     for mult in sorted_multiples_list]) + '\t' +
                                f"{sum([0 if overlap_statistics[mult]['bad_intervals'].get(b_interval) is None else overlap_statistics[mult]['bad_intervals'][b_interval]['bad_regions'] for mult in sorted_multiples_list]):,}" + '\t' +
                                f"{sum([0 if overlap_statistics[mult]['bad_intervals'].get(b_interval) is None else overlap_statistics[mult]['bad_intervals'][b_interval]['bad_bases'] for mult in sorted_multiples_list]):,}" + '\n'
                                for b_interval in sorted_overlapping_intervals])
        with open(output_overlaps_table, 'wt') as f_total_table:
            f_total_table.writelines(stats_out_lines)
        # sort list of intervals according to the lowest number of overlapping bases
        with open(output_overlaps_table, 'rt') as f_total_table:
            header_lines = (f_total_table.readline().strip().split('\t'),
                            f_total_table.readline().strip().split('\t'))
            all_stats_content = [line.strip().split('\t') for line in f_total_table.readlines()]
        total_stats_column = header_lines[0].index('TOTAL')
        bases_column = 2  # OR for all thresholds: header_lines[1].index('BASES')
        regions_column = 1  # OR for all thresholds: header_lines[1].index('REGIONS')
        contents_and_overlap = [(interval_content[0],
                                 int(interval_content[regions_column]),
                                 int(interval_content[bases_column]))
                                for interval_content in all_stats_content]
        intervals_per_overlapping_bases = sorted(contents_and_overlap, key=lambda x: [x[2], x[1]])
        intervals_per_overlapping_regions = sorted(contents_and_overlap, key=lambda x: [x[1], x[2]])
        output_base_overlaps_table = pth_join(output_path,
                                              f'sorted-{CHUNK_SIZE//1000}kb-intervals_per_overlapping_bases.tsv')
        with open(output_base_overlaps_table, 'wt') as f_sorted_bases:
            f_sorted_bases.writelines([f"{chrom}\t{regions:,}\t{bases:,}\n"
                                       for chrom, regions, bases in intervals_per_overlapping_bases])
