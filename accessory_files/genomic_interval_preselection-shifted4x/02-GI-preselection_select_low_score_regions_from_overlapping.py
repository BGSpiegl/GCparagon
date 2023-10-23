#!/usr/bin/env python3
import sys
import pathlib
from re import sub as re_sub
from typing import Tuple, List

# path definitions (relative imports):
CODE_ROOT_PATH = pathlib.Path(__file__).parent.parent.parent
CODE_ROOT_DIR = str(CODE_ROOT_PATH)

# add src path to PYTHONPATH variable (= current parent dir):
if CODE_ROOT_DIR not in sys.path:
    print("| INFO - adding GCparagon_dev software parent directory to Python module search path "
          f"variable: {CODE_ROOT_DIR}")
    sys.path.append(CODE_ROOT_DIR)

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# info - all input regions within 1Mbp of the chromosome ends are removed
genome_file_path = CODE_ROOT_PATH / 'accessory_files/GRCh38.genome_file.tsv'  # <---- required input!
# TODO: put here your shifted regions BED files containing exclusion list overlap
# THE FOLLOWING VARIABLES SHOULD BE IDENTICAL TO THE ONES FROM THE "01-GI-preselection_....py" SCRIPT: -----------------
CHUNK_SIZE = 10**6
SHIFT_N_TIMES = 10  # you might want to select a higher number
search_path = CODE_ROOT_PATH / f'accessory_files/genomic_interval_preselection-shifted{SHIFT_N_TIMES}x'
output_path = search_path
# ----------------------------------------------------------------------------------------------------------------------
all_region_shifts_with_overlaps = list(search_path.glob(
    f"*kbp_intervalOffset/{CHUNK_SIZE//10**3}kbp_intervals_bad_regions_overlap_ELRminSizes.tsv"))
all_region_shifts_with_overlaps.extend(list(search_path.glob(
    f"{CHUNK_SIZE//10**3}kbp_intervals_bad_regions_overlap_ELRminSizes.tsv")))
# ^---- required input!
# GCparagon: overlapping-exclusion-listed-bases, shifted up to 3 times by genomic_interval_size/4
#            -> 4 files: 1x un-shifted + 3x shifted assertion: intervals of1 equal size!
REGION_SIZE = 10**6  # 1Mbp -> change according to input BED files!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

with open(genome_file_path, 'rt') as f_genome:
    content = [line.strip().split('\t') for line in f_genome.readlines()]
chromosome_lengths = {}.fromkeys(map(lambda x: x[0], content))
for chrom, length in content:
    chromosome_lengths[chrom] = int(length)

REF_BUILD = 'hg38'
CHROM_END_DISTANCE = 10**6  # regions within 1 Mbp from each chromosome end ware skipped
REGION_OVERLAP_PERCENTAGE_THRESHOLD = 33  # percent max. bad region overlap!
# i : in total, there were 237 events where a selected region overlapped more than the selected threshold of 33% of
#     region size.
# i : 295 cases of lowest-score-selection were skipped (66% allowed overlap)


def somehow_overlap(r1: range, r2: range) -> bool:
    return r1.start in r2 or r1.stop in r2 or r2.start in r1 or r2.stop in r1


def select_lowest_score_element(elements: List[Tuple[range, int]], score_pos: int) \
        -> Tuple[Tuple[range, int], Tuple[range]]:
    ascending_score_elements = sorted(elements, key=lambda x: x[score_pos], reverse=False)
    # sort ascending, take first item
    return ascending_score_elements[0], tuple(map(lambda x: x[0], ascending_score_elements[1:]))


if __name__ == '__main__':
    # read in all regions into dictionary per chromosome
    regions_overlap_dict = {}
    chroms_in_order = []
    for tsv_f in all_region_shifts_with_overlaps:
        with open(tsv_f, 'rt') as f_vrlp:
            hdr1 = f_vrlp.readline().strip().split('\t')
            hdr2 = f_vrlp.readline().strip().split('\t')
            bad_bases_overlap_column_idx = hdr2.index('bases')  # should be 2
            for lin_cont in f_vrlp.readlines():
                interesting_content = lin_cont.strip().split('\t')[:bad_bases_overlap_column_idx+1]
                chrom_str = interesting_content[0]
                overlap = int(re_sub(',', '', interesting_content[-1]))
                chromosome = chrom_str.split(':')[0]
                if chromosome not in chroms_in_order:
                    chroms_in_order.append(chromosome)
                # check if is within 1 Mbp of chromosome end
                start_coord = int(re_sub(',', '', chrom_str.split(':')[1].split('-')[0]))
                stop_coord = int(re_sub(',', '', chrom_str.split(':')[1].split('-')[1])) - 1
                if somehow_overlap(r1=range(start_coord, stop_coord),
                                   r2=(range(chromosome_lengths[chromosome] - CHROM_END_DISTANCE
                                             if chromosome_lengths[chromosome] - CHROM_END_DISTANCE >= 0 else 0,
                                             chromosome_lengths[chromosome])))\
                        or \
                        somehow_overlap(r1=range(start_coord, stop_coord),
                                        r2=(range(0, chromosome_lengths[chromosome]
                                                  if chromosome_lengths[chromosome] < CHROM_END_DISTANCE else
                                                  CHROM_END_DISTANCE))):
                    continue  # ignore end region overlapping intervals!
                region_range = range(start_coord,
                                     stop_coord)  # end-exclusive
                score = int(re_sub(',', '', interesting_content[bad_bases_overlap_column_idx]))
                try:
                    regions_overlap_dict[chromosome].append((region_range, score))
                except KeyError:
                    regions_overlap_dict.update({chromosome: [(region_range, score)]})
    # NOW: for each chromosome, select non-overlapping regions based on lowest score!
    chosen_regions = {}.fromkeys(regions_overlap_dict)
    skip_these_regions = {}.fromkeys(regions_overlap_dict)
    for c in skip_these_regions.keys():
        skip_these_regions[c] = set()
        chosen_regions[c] = []
    total_events_score_too_high = 0
    for chrom, reg_ranges in regions_overlap_dict.items():
        too_high_score_events = 0
        last_selected_range = range(-2, -1)
        for r_idx, (current_reg_range, region_score) in enumerate(reg_ranges):
            if current_reg_range in skip_these_regions[chrom]:
                continue  # ignore region; was already ignored in another overlap trial
            overlapping_regions = list(filter(lambda x: x is not None,
                                              [(reg, scr) if somehow_overlap(r1=current_reg_range, r2=reg) else None
                                               for reg, scr in reg_ranges]))  # includes current region tuple
            if any([somehow_overlap(r1=ovrlp_reg[0], r2=last_selected_range)
                    for ovrlp_reg in overlapping_regions]):  # check if any also overlap with prev. selected
                overlapping_regions.append(last_selected_range)
            selected_tuple, tuple_not_selected_ranges = select_lowest_score_element(elements=overlapping_regions,
                                                                                    score_pos=1)
            # mark as to skip for future
            skip_these_regions[chrom].update(tuple_not_selected_ranges)
            skip_these_regions[chrom].update((selected_tuple[0],))
            # check for score threshold violation:
            if selected_tuple[1] > \
                    int((selected_tuple[0].stop - selected_tuple[0].start) * REGION_OVERLAP_PERCENTAGE_THRESHOLD / 100):
                too_high_score_events += 1
                print(f"WARNING: tried to select region with lowest exclusion list overlap from among "
                      f"{len(overlapping_regions):,} overlapping regions (overlap with "
                      f"'{chrom}:{current_reg_range.start}-{current_reg_range.stop}'), but smallest score was above "
                      f"max. percentage threshold of {REGION_OVERLAP_PERCENTAGE_THRESHOLD}%. None of the involved "
                      f"regions will be present in the output! Continuing ..")
                continue
            # score valid -> test if not overlaps any previous selection!!
            # check if overlaps with something already chosen -> select next best region and check again
            check_again = overlapping_regions
            while any([somehow_overlap(r1=selected_tuple[0], r2=rg[0]) for rg in chosen_regions[chrom]]):
                check_again.remove(selected_tuple)  # remove the overlapping region
                selected_tuple, tuple_not_selected_ranges = select_lowest_score_element(
                    elements=check_again, score_pos=1)
            if selected_tuple:  # might be none
                # find all regions that overlap with finally selected region and mark as skip
                selection_overlapping_ranges = list(
                    map(lambda y: y[0], filter(lambda x: x is not None,
                                               [(reg, scr) if somehow_overlap(r1=selected_tuple[0], r2=reg) else None
                                                for reg, scr in reg_ranges])))
                last_selected_range = selected_tuple[0]  # add to selected regions
                if selected_tuple not in chosen_regions[chrom]:
                    chosen_regions[chrom].append(selected_tuple)
                    skip_these_regions[chrom].update(selection_overlapping_ranges)
        print(f" i : encountered {too_high_score_events:,} events for {chrom} where the selected region had a score "
              f"above the threshold")
        total_events_score_too_high += too_high_score_events
    # report total stats
    print("----------------------------------------------------------------------------------------------------------\n"
          f" i : in total, there were {total_events_score_too_high:,} events where a selected region overlapped more "
          f"than the selected threshold of {REGION_OVERLAP_PERCENTAGE_THRESHOLD}% of region size.")
    # gib selection stats
    print(f"I am done selecting regions. These are the statistics (note: genome is represented {SHIFT_N_TIMES} times "
          f"(i.e. original intervals {SHIFT_N_TIMES-1}x shifted):")
    for chrom in chroms_in_order:
        print(f" i : {chrom} selected regions: {len(chosen_regions[chrom]):,}\n"
              f"     {' '*len(chrom)}  skipped regions: {len(skip_these_regions[chrom]):,}")
    # write output:
    output_file = output_path / f'{REF_BUILD}_minimalExclusionListOverlap_{CHROM_END_DISTANCE//10**6}Mbp_intervals_' \
        f'{REGION_OVERLAP_PERCENTAGE_THRESHOLD}pcOverlapLimited.bed'
    reg_buffer = []
    n_lines_written = 0
    with open(output_file, 'wt') as f_chosen_regs:
        for chrom in chroms_in_order:
            for reg_range, score in sorted(chosen_regions[chrom], key=lambda x: x[0].start, reverse=False):
                reg_buffer.append(f"{chrom}\t{reg_range.start}\t{reg_range.stop+1}\t{score}\n")  # add back 1 to end
                if len(reg_buffer) >= 500:
                    f_chosen_regs.writelines(reg_buffer)
                    n_lines_written += len(reg_buffer)
                    reg_buffer = []
            if reg_buffer:
                f_chosen_regs.writelines(reg_buffer)
                n_lines_written += len(reg_buffer)
                reg_buffer = []
    print(f"done writing {n_lines_written:,} lines (= {REGION_SIZE//10**6}Mbp regions) to BED file "
          f"'{output_file}'")


# RESULT:
# ----------------------------------------------------------------------------------------------------------
#  i : in total, there were 262 events where a selected region overlapped more than the selected threshold of 20% of
#      region size.
# I am done selecting regions. These are the statistics (note: genome is represented 4 times, 4 shifts):
#  i : chr1 selected regions: 135
#               skipped regions: 1,100
#  i : chr2 selected regions: 137
#               skipped regions: 1,069
#  i : chr3 selected regions: 116
#               skipped regions: 878
#  i : chr4 selected regions: 115
#               skipped regions: 843
#  i : chr5 selected regions: 106
#               skipped regions: 803
#  i : chr6 selected regions: 100
#               skipped regions: 750
#  i : chr7 selected regions: 91
#               skipped regions: 704
#  i : chr8 selected regions: 86
#               skipped regions: 640
#  i : chr9 selected regions: 67
#               skipped regions: 598
#  i : chr10 selected regions: 78
#               skipped regions: 583
#  i : chr11 selected regions: 80
#               skipped regions: 591
#  i : chr12 selected regions: 79
#               skipped regions: 590
#  i : chr13 selected regions: 59
#               skipped regions: 498
#  i : chr14 selected regions: 54
#               skipped regions: 458
#  i : chr15 selected regions: 44
#               skipped regions: 429
#  i : chr16 selected regions: 46
#               skipped regions: 390
#  i : chr17 selected regions: 44
#               skipped regions: 360
#  i : chr18 selected regions: 47
#               skipped regions: 349
#  i : chr19 selected regions: 34
#               skipped regions: 249
#  i : chr20 selected regions: 37
#               skipped regions: 275
#  i : chr21 selected regions: 19
#               skipped regions: 189
#  i : chr22 selected regions: 22
#               skipped regions: 212
#  i : chrX selected regions: 90
#               skipped regions: 690
#  i : chrY selected regions: 9
#               skipped regions: 224
# done writing 1,695 lines (= 1Mbp regions) to BED file '/benchmark_solutions/regions/hg38_top_lowest_BlacklistOverlap_1Mbp_regions_max25pcOverlap.bed'


# 10x shifted regions:
#  i : encountered 38 events for chrY where the selected region had a score above the threshold
# ----------------------------------------------------------------------------------------------------------
#  i : in total, there were 204 events where a selected region overlapped more than the selected threshold of 33% of region size.
# I am done selecting regions. These are the statistics (note: genome is represented 10 times (i.e. original intervals 9x shifted):
#  i : chr1 selected regions: 133
#            skipped regions: 2,585
#  i : chr2 selected regions: 138
#            skipped regions: 2,520
#  i : chr3 selected regions: 114
#            skipped regions: 2,060
#  i : chr4 selected regions: 117
#            skipped regions: 1,986
#  i : chr5 selected regions: 103
#            skipped regions: 1,881
#  i : chr6 selected regions: 98
#            skipped regions: 1,770
#  i : chr7 selected regions: 93
#            skipped regions: 1,644
#  i : chr8 selected regions: 84
#            skipped regions: 1,493
#  i : chr9 selected regions: 65
#            skipped regions: 1,414
#  i : chr10 selected regions: 75
#             skipped regions: 1,377
#  i : chr11 selected regions: 75
#             skipped regions: 1,383
#  i : chr12 selected regions: 78
#             skipped regions: 1,372
#  i : chr13 selected regions: 55
#             skipped regions: 1,020
#  i : chr14 selected regions: 50
#             skipped regions: 1,088
#  i : chr15 selected regions: 44
#             skipped regions: 1,028
#  i : chr16 selected regions: 45
#             skipped regions: 914
#  i : chr17 selected regions: 45
#             skipped regions: 843
#  i : chr18 selected regions: 41
#             skipped regions: 809
#  i : chr19 selected regions: 33
#             skipped regions: 588
#  i : chr20 selected regions: 34
#             skipped regions: 646
#  i : chr21 selected regions: 20
#             skipped regions: 455
#  i : chr22 selected regions: 19
#             skipped regions: 411
#  i : chrX selected regions: 89
#            skipped regions: 1,612
#  i : chrY selected regions: 11
#            skipped regions: 553
# done writing 1,659 lines (= 1Mbp regions) to BED file '/accessory_files/genomic_interval_preselection-shifted4x-shifted4x/hg38_minimalExclusionListOverlap_1Mbp_intervals_33pcOverlapLimited.bed'
#
# Process finished with exit code 0