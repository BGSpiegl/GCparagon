#!/usr/bin/python3

from typing import List, Tuple, Dict, Union
from os.path import join as pth_join, sep as pth_sep
from os import makedirs as os_mkdirs
from shutil import which as sh_which
from sys import exit as sys_exit
from pybedtools import BedTool
from natsort import humansorted
from pathlib import Path

SCRIPT_ROOT_PATH = Path(__file__).parent
SCRIPT_ROOT_DIR = str(SCRIPT_ROOT_PATH)

blacklist = '/home/benjamin/GitHub_Local_Clones/GCparagon_public/accessory_files/hg38_GCcorrection_blacklist.merged.sorted.bed'  # 46.2 MB

# '/home/benjamin/GitHub_Local_Clones/GCparagon_dev/test/BAM_intersection/'\
# 'hg38_GCcorrection_blacklist.merged.bed'

# or: '/home/benjamin/GitHub_Local_Clones/GCparagon_dev/test/BAM_intersection/' \
#      'hg38_GCcorrection_blacklist.merged.bed'
simulate_bed = '/home/benjamin/GitHub_Local_Clones/GCparagon_dev/test/BAM_intersection/' \
               'NPH_004.simulate.bed'  # more entries (10x)
fetch_bed = '/home/benjamin/GitHub_Local_Clones/GCparagon_dev/test/BAM_intersection/' \
            'NPH_004.fetch.bed'  # fewer entries
large_bam_path = '/media/benjamin/Analyses/GC_correction_benchmarks/benchmarking_samples/NPH_004.bam'  # 64.4 GB
small_bam_path = '/media/benjamin/Analyses/GC_correction_benchmarks/benchmarking_samples/B58_3.bam'  # 28.5 GB
genome_file_path = '/media/benjamin/Analyses/GC_correction_benchmarks/benchmarking_samples/NPH_004.genome_file.tsv'
output_path = '/test/BAM_intersection'

bedtools_path = sh_which('bedtools')

CHUNK_SIZE = 1000000  # 1 Mb
CHUNK_POSITION_OFFSETS = (0, CHUNK_SIZE // 4, CHUNK_SIZE // 2, CHUNK_SIZE // 4 * 3)


if not bedtools_path:
    print("ERROR - bedtools path not found. Cannot benchmark intersection times.")
    sys_exit(1)


def read_bed_file(bed_path: str) -> List[Tuple[str, int, int]]:
    with open(bed_path, 'rt') as f_bed:
        return [(chrom, int(start), int(stop))
                for chrom, start, stop, *_ in filter(lambda x: x != '', [bed_line.strip().split()
                                                                         for bed_line in f_bed.readlines()])]


def get_stdchrom_chunks(genome_file: str, chunk_size=CHUNK_SIZE, offset=0) -> list:
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
    std_chunks = []
    for sorted_chrom in chroms_sorted:
        for chrom, _strt, stop in std_regs:
            if chrom == sorted_chrom:
                n_splits = (stop - offset) // chunk_size
                cum_size = offset
                for split_idx in range(0, n_splits, 1):
                    std_chunks.append((chrom, cum_size, cum_size + chunk_size))
                    cum_size += chunk_size
                std_chunks.append((chrom, cum_size, stop))
    return std_chunks


def get_overlap_statistics(overlaps_bed_path: str) -> Dict[Tuple[str, int, int], Dict[str, int]]:
    with open(overlaps_bed_path, 'rt') as f_overlaps:
        overlap_content = [line.strip().split('\t') for line in f_overlaps.readlines()]
    bad_chunk_bases = {}
    for chunk_cont, chunk_start, chunk_stop, _1, _2, _3, *rest in overlap_content:
        overlapping_bases = rest[-1]
        try:
            bad_chunk_bases[(chunk_cont, int(chunk_start), int(chunk_stop))].append(int(overlapping_bases))
        except KeyError:
            bad_chunk_bases.update({(chunk_cont, int(chunk_start), int(chunk_stop)): [int(overlapping_bases)]})
    # postprocess: compute sum of overlapping_bases per chunk
    for chunk in bad_chunk_bases.keys():
        bad_chunk_bases[chunk] = {'bad_regions': len(bad_chunk_bases[chunk]),
                                  'bad_bases': sum(bad_chunk_bases[chunk])}
    return {'bad_chunks': bad_chunk_bases}  # bad_chunk_bases might be empty in cases of a bad bad regions BED file


def convert_chunk_tuple_str(chunk_str=None, chunk_tuple=None) -> Union[Tuple[str, int, int], str]:
    if all([p is None for p in (chunk_str, chunk_tuple)]):
        raise AttributeError("at least one input parameter must be not None")
    if chunk_str is not None and not isinstance(chunk_str, str):
        raise AttributeError("provided chunk_str was not a str!")
    if chunk_tuple is not None and not isinstance(chunk_str, tuple):
        raise AttributeError("provided chunk_tuple was not a tuple!")
    if chunk_str:
        chunk_tuple = (chunk_str.split('-')[0],
                       int(chunk_str.split('-')[1].split(':')[0]),
                       int(chunk_str.split(':')[1]))
        return chunk_tuple
    chunk_str = f"{chunk_tuple[0]}:{chunk_tuple[1]}-{chunk_tuple[2]}"
    return chunk_str


if __name__ == '__main__':
    # read original blacklist
    blacklist_content = read_bed_file(bed_path=blacklist)
    # find blacklist chunks larger than 1kb (or multiples of this)
    multiples = (0, )  # 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 50, 100)  # 0 represents all blacklisted regions
    num_multiples = {}.fromkeys(multiples)
    sub_blacklists = {}.fromkeys(multiples)
    for k in sub_blacklists.keys():
        sub_blacklists[k] = []
    mult_blacklist_beds = {}.fromkeys(multiples)
    sorted_std_chroms = humansorted([f'chr{i}' for i in range(1, 23, 1)] + ['chrx', 'chry'])
    for sorted_chrom in sorted_std_chroms:
        for chrom, start, stop in blacklist_content:  # [ DONE ]
            if chrom == sorted_chrom:
                sub_blacklists[0].append(f"{chrom}\t{start}\t{stop}\n")  # always append to basic list
                reg_size = stop - start  # [ DONE. ]
                for mult in multiples:
                    if reg_size >= 1000 * mult:  # filter for regions
                        sub_blacklists[mult].append(f"{chrom}\t{start}\t{stop}\n")
                        try:
                            num_multiples[mult] += 1
                        except TypeError:
                            num_multiples[mult] = 1
    for k in num_multiples.keys():
        if num_multiples[k] is None:
            num_multiples[k] = 0
    for mult in multiples:  # [ DONE ]
        print(f"number of blacklisted regions larger or equal than {mult}kb: {num_multiples[mult]:,}")
    # write filtered blacklists
    for mult, sublist_lines in sub_blacklists.items():
        if mult == 0:
            mult_blacklist_beds[mult] = blacklist  # THIS was used -> original blacklist <-> mult=0
    # intersect genomic chunks and blacklisted regions
    os_mkdirs(output_path, exist_ok=True)
    # chunk up the genome using different position offsets
    for chunk_offset in CHUNK_POSITION_OFFSETS:
        print(f"NOW: computing bad regions overlap for chunks with offset {chunk_offset} ..")
        if chunk_offset == 0:
            current_output_dir = output_path
        else:
            current_output_dir = f"{output_path.rstrip(pth_sep)}_{chunk_offset//1000}kbp_chunkOffset"
            os_mkdirs(current_output_dir, exist_ok=True)
        # chunk up the genome
        std_chunks = get_stdchrom_chunks(genome_file=genome_file_path, chunk_size=CHUNK_SIZE, offset=chunk_offset)
        # write whole genome BED file
        whole_genome_bed_path = pth_join(current_output_dir, f'hg38_{CHUNK_SIZE//1000}kbpchunks.bed')
        with open(whole_genome_bed_path, 'wt') as f_whole_g:
            f_whole_g.writelines([f"{chrom}\t{start}\t{stop}\n" for chrom, start, stop in std_chunks])
        # compute overlap between each chunk in whole genome BED and the merged, sorted blacklist instance
        overlap_statistics = {}.fromkeys(multiples)
        for mult in multiples:
            out_overlap_bed = pth_join(current_output_dir, f'hg38_{CHUNK_SIZE}bpStdChunks_gte{mult}kb_BL-overlap.bed')
            overlap_statistics[mult] = {'overlaps_bed_path': out_overlap_bed}
            intersect_bedtool = BedTool(mult_blacklist_beds[mult])
            wg_bedtool = BedTool(whole_genome_bed_path)
            wg_bedtool.intersect(intersect_bedtool, wo=True).saveas(out_overlap_bed)
        for mult in multiples:
            overlap_statistics[mult].update(get_overlap_statistics(overlaps_bed_path=overlap_statistics[mult][
                'overlaps_bed_path']))
        # compute all bad chunks -> create complete table with columns being multiples
        unique_overlapping_chunks = set()
        for mult in multiples:
            unique_overlapping_chunks.update(overlap_statistics[mult]['bad_chunks'].keys())
        sorted_overlapping_chunks = humansorted(list(unique_overlapping_chunks),
                                                key=lambda x: [x[0], x[1]])
        # write total stats table
        output_overlaps_table = pth_join(current_output_dir,
                                         f'{CHUNK_SIZE//1000}kbp_chunks_bad_regions_overlap_BRminSizes.tsv')
        sorted_multiples_list = sorted(list(multiples))
        stats_out_lines = ['\t' + '\t\t'.join([f"{m}kb_BR_min_size" for m in sorted_multiples_list]) + '\t\tTOTAL\t\n' +
                           'chunk\t' + '\t'.join(['regions\tbases' for _m in sorted_multiples_list]) + '\tBASES\tREGIONS\n']
        stats_out_lines.extend([f'{b_chunk[0]}:{b_chunk[1]:,}-{b_chunk[2]:,}\t' +
                                '\t'.join(
                                    ['0\t0'
                                     if overlap_statistics[mult]['bad_chunks'].get(b_chunk) is None
                                     else (str(overlap_statistics[mult]['bad_chunks'][b_chunk]['bad_regions']) + '\t' +
                                           str(overlap_statistics[mult]['bad_chunks'][b_chunk]['bad_bases']))
                                     for mult in sorted_multiples_list]) + '\t' +
                                f"{sum([0 if overlap_statistics[mult]['bad_chunks'].get(b_chunk) is None else overlap_statistics[mult]['bad_chunks'][b_chunk]['bad_regions'] for mult in sorted_multiples_list]):,}" + '\t' +
                                f"{sum([0 if overlap_statistics[mult]['bad_chunks'].get(b_chunk) is None else overlap_statistics[mult]['bad_chunks'][b_chunk]['bad_bases'] for mult in sorted_multiples_list]):,}" + '\n'
                                for b_chunk in sorted_overlapping_chunks])
        with open(output_overlaps_table, 'wt') as f_total_table:
            f_total_table.writelines(stats_out_lines)

        # sort list of chunks according to lowest number of overlapping bases
        with open(output_overlaps_table, 'rt') as f_total_table:
            header_lines = (f_total_table.readline().strip().split('\t'),
                            f_total_table.readline().strip().split('\t'))
            all_stats_content = [line.strip().split('\t') for line in f_total_table.readlines()]
        total_stats_column = header_lines[0].index('TOTAL')
        bases_column = 2  # OR for all thresholds: header_lines[1].index('BASES')
        regions_column = 1  # OR for all thresholds: header_lines[1].index('REGIONS')
        contents_and_overlap = [(chunk_content[0],
                                 int(chunk_content[regions_column]),
                                 int(chunk_content[bases_column]))
                                for chunk_content in all_stats_content]
        chunks_per_overlapping_bases = sorted(contents_and_overlap, key=lambda x: [x[2], x[1]])
        chunks_per_overlapping_regions = sorted(contents_and_overlap, key=lambda x: [x[1], x[2]])
        output_base_overlaps_table = pth_join(output_path,
                                              f'sorted-{CHUNK_SIZE//1000}kb-chunks_per_overlapping_bases.tsv')
        with open(output_base_overlaps_table, 'wt') as f_sorted_bases:
            f_sorted_bases.writelines([f"{chrom}\t{regions:,}\t{bases:,}\n"
                                       for chrom, regions, bases in chunks_per_overlapping_bases])
