#!/usr/bin/python3
import sys
import pathlib
from re import sub as re_sub
from typing import Tuple, List

# path definitions (relative imports):
SOURCE_CODE_ROOT_PATH = pathlib.Path(__file__).parent.parent
SOURCE_CODE_ROOT_DIR = str(SOURCE_CODE_ROOT_PATH)

# remove all regions from input that are within 1Mbp of the chromosome ends
genome_file_path = SOURCE_CODE_ROOT_PATH / 'accessory_files/GRCh38.genome_file.tsv'
chunks_out_path = SOURCE_CODE_ROOT_PATH / 'accessory_files'

with open(genome_file_path, 'rt') as f_genome:
    content = [line.strip().split('\t') for line in f_genome.readlines()]
chromosome_lengths = {}.fromkeys(map(lambda x: x[0], content))
for chrom, length in content:
    chromosome_lengths[chrom] = int(length)

REF_BUILD = 'hg38'
CHROM_END_DISTANCE = 10**6  # use 1 Mbp safety distance to chromosome start/ends
REGION_OVERLAP_PERCENTAGE_THRESHOLD = 33  # percent max. bad region overlap!
all_region_shifts_with_overlaps = []  # overlapping-blacklisted-bases, shifted up to 3 times by chunk_size/4
# -> 4 files: 1x un-shifted + 3x shifted assertion: chunks of equal size!


def somehow_overlap(r1: range, r2: range) -> bool:
    return r1.start in r2 or r1.stop in r2 or r2.start in r1 or r2.stop in r1


def select_lowest_score_element(elements: List[Tuple[range, int]], score_pos: int) \
        -> Tuple[Tuple[range, int], Tuple[range]]:
    ascending_score_elements = sorted(elements, key=lambda x: x[score_pos], reverse=False)  # sort ascending, take first
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
                line_content_of_interest = lin_cont.strip().split('\t')[:bad_bases_overlap_column_idx+1]
                chrom_str = line_content_of_interest[0]
                overlap = int(re_sub(',', '', line_content_of_interest[-1]))
                chromosome = chrom_str.split(':')[0]
                if chromosome not in chroms_in_order:
                    chroms_in_order.append(chromosome)
                # check if is within 1 Mbp of chromosome end
                start_coord = int(re_sub(',', '', chrom_str.split(':')[1].split('-')[0]))
                stop_coord = int(re_sub(',', '', chrom_str.split(':')[1].split('-')[1])) - 1
                if somehow_overlap(r1=range(start_coord, stop_coord),  # chunk
                                   r2=(range(chromosome_lengths[chromosome] - CHROM_END_DISTANCE
                                             if chromosome_lengths[chromosome] - CHROM_END_DISTANCE >= 0 else 0,
                                             chromosome_lengths[chromosome]))) \
                        or \
                        somehow_overlap(r1=range(start_coord, stop_coord),
                                        r2=(range(0, CHROM_END_DISTANCE
                                                  if chromosome_lengths[chromosome] < CHROM_END_DISTANCE else
                                                  chromosome_lengths[chromosome]))):
                    continue  # ignore end region overlapping chunks!
                region_range = range(start_coord,
                                     stop_coord)  # end-exclusive
                score = int(re_sub(',', '', line_content_of_interest[bad_bases_overlap_column_idx]))
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
            if any([somehow_overlap(r1=ocrlp_reg[0], r2=last_selected_range)
                    for ocrlp_reg in overlapping_regions]):  # check if any also overlap with prev. selected
                overlapping_regions.append(last_selected_range)
            selected_tuple, tuple_not_selected_ranges = select_lowest_score_element(elements=overlapping_regions,
                                                                                    score_pos=1)
            # mark as to skip for future iterations
            skip_these_regions[chrom].update(tuple_not_selected_ranges)
            skip_these_regions[chrom].update((selected_tuple[0],))
            # check for score threshold violation:
            if selected_tuple[1] > \
                    int((selected_tuple[0].stop - selected_tuple[0].start) * REGION_OVERLAP_PERCENTAGE_THRESHOLD / 100):
                too_high_score_events += 1
                print(f"WARNING: tried to select region with lowest blacklist overlap from among "
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
    print(f"I am done selecting regions. These are the statistics (note: genome is represented 4 times, 4 shifts):")
    for chrom in chroms_in_order:
        print(f" i : {chrom} selected regions: {len(chosen_regions[chrom]):,}\n"
              f"     {' '*len(chrom)}  skipped regions: {len(skip_these_regions[chrom]):,}")
    # write output:
    chunks_out_path.mkdir(parents=True, exist_ok=True)
    output_file = chunks_out_path / f'{REF_BUILD}_minimalBlacklistOverlap_1Mbp_chunks_' \
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
    print(f"done writing {n_lines_written:,} lines (= {CHROM_END_DISTANCE//10**6}Mbp regions) to BED file "
          f"'{output_file}'")
