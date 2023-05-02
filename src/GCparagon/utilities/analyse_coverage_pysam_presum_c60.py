# Built-in/Generic Imports
import logging
import sys
import argparse
import numpy
import scipy
import scipy.stats
import os.path
import multiprocessing
import pandas as pd
from scipy import stats
from pathlib import Path
import glob
import pysam  # WARNING - THIS MUST USE THE ADAPTED PYSAM VERSION FROM https://github.com/Faruk-K/pysam!
# OTHERWISE THE '-gc' flag won't work! (i.e. the pysam.AligmentFile.count_coverage_c60() function won't be available)


# Initialise logger
def init_logger():
    logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s',
                        filemode='w',
                        datefmt='%d-%m-%Y %H:%M:%S',
                        level=logging.INFO,
                        stream=sys.stdout)


# Log input arguments
def log_input_arguments(args):
    logging.info(f"Starting program with following arguments: ")
    for arg_key in vars(args):
        logging.info(f"{arg_key.upper()}: {vars(args)[arg_key]}")
    logging.info("")


# Calculate mean and confidence intervals ##############################################################################
def mean_confidence_interval(data, confidence=0.95):
    a = 1.0 * numpy.array(data)
    n = len(a)
    m, se = numpy.mean(a), scipy.stats.sem(a)
    h = se * scipy.stats.t.ppf((1 + confidence) / 2., n - 1)
    return m, m - h, m + h


# Read in copy-number data #############################################################################################
def readCopyData(cna_file):
    cna_path = Path(cna_file)
    extension = ''.join(cna_path.suffixes)
    LOG2 = open(cna_path, "r")
    _header = LOG2.readline()
    log_list = list()
    for line in LOG2.readlines():
        chrom, start, end, log2 = line.rstrip().split("\t")
        if extension == '.correctedDepth.txt':
            if log2 != 'NA':
                log_list.append(['chr' + chrom, int(start), int(end), float(log2)])
        elif extension == '.segments':
            log_list.append([chrom, int(start), int(end), float(log2)])
        else:
            logging.error("Unknown extension " + extension)
            raise Exception("Unknown extension " + extension)
    LOG2.close()
    return log_list


# Calculate mean coverage in 10000bp upstream sequence #################################################################
def calcLocalMean(position, chrom, cnv_data):
    norm = None
    for region in cnv_data:
        reg_chrom, reg_start, reg_end, reg_log2 = region
        if chrom == reg_chrom and reg_start < position < reg_end:
            norm = numpy.power(2, reg_log2)
            break
    if norm:
        return norm
    else:
        return 1


# Calculate Coverage around position ###################################################################################
def pysam_presum_header(q, thread_number, transcript_list, header_map, args):
    sys.stderr.write("Thread " + str(thread_number) + " started\n")
    sys.stderr.flush()
    position_lists = dict()
    position_ctrs = dict()
    tss_visited = list()
    line_count = 0
    skipped = 0
    processed = 0
    bam = pysam.AlignmentFile(args['bam_file'], 'rb')

    for loc in range(-args['start'], args['end'] + 1):
        position_lists[loc] = 0
        position_ctrs[loc] = 0

    for transcript_split in transcript_list:
        chrom = transcript_split[header_map['chr']]
        pos = int(transcript_split[header_map['peak']])


        if (line_count + 1) % 1000 == 0:
            sys.stderr.write("Thread " + str(thread_number) + "\t" + str(line_count + 1) + " positions analyzed\n")
            sys.stderr.flush()
        line_count += 1

        if len(transcript_split) != 6:
            forward = True
        elif transcript_split[5] == '+':
            forward = True
        else:
            forward = False
        if chrom + "_" + str(pos) in tss_visited:
            skipped += 1
            continue
        if pos - args['start'] < 1:
            print("WARNING pos - start: ", transcript_split)
            skipped += 1
            continue
        if chrom.find("_") != -1:
            print("WARNING find _: ", transcript_split)
            skipped += 1
            continue
        if chrom == "chrY":
            print("WARNING chrY: ", transcript_split)
            skipped += 1
            continue
        processed += 1

        cov_start = pos - args['start']
        cov_end = pos + args['end'] + 1

        coverage_tuple = bam.count_coverage_c60(chrom, cov_start, cov_end, quality_threshold=args['mapping_quality'])

        coverage_list = list()
        for i in range(0, len(coverage_tuple[0])):
            print(coverage_tuple)
            coverage_list.append(
                coverage_tuple[0][i] + coverage_tuple[1][i] + coverage_tuple[2][i] + coverage_tuple[3][i])
        # normalize for copy-number variations
        if args['normalisation_file']:
            cnv_data = readCopyData(args['normalisation_file'])
            normcov = calcLocalMean(pos, chrom, cnv_data)
            if normcov == 0:
                normcov = 0.001
        counter = 0
        for i in range(pos - args['start'], pos + args['end'] + 1):
            coverage = float(coverage_list[counter])  # / args['mean_coverage']
            if coverage > args['coverage_limit']:
                counter += 1
                continue
            if args['normalisation_file']:
                coverage /= normcov
            if forward:
                position_lists[i - pos] += coverage
                position_ctrs[i - pos] += 1
            elif not forward:
                position_lists[-(i - pos)] += coverage
                position_ctrs[-(i - pos)] += 1
            counter += 1
        tss_visited.append(chrom + "_" + str(pos))


    sys.stderr.write("Thread " + str(thread_number) + "\t finished\n")
    sys.stderr.flush()
    q.put((position_lists, position_ctrs))


# Calculate Coverage around position ###################################################################################
def pysam_presum(q, thread_number, transcript_list, args):
    sys.stderr.write("Thread " + str(thread_number) + " started\n")
    sys.stderr.flush()
    position_lists = dict()
    position_ctrs = dict()
    tss_visited = list()
    line_count = 0
    skipped = 0
    processed = 0
    coverage_dict = dict()

    bam = pysam.AlignmentFile(args['bam_file'], 'rb')
    for loc in range(-args['start'], args['end'] + 1):
        position_lists[loc] = 0
        position_ctrs[loc] = 0

    for transcript_split in transcript_list:
        chrom = transcript_split[0]
        start = int(transcript_split[1])
        end = int(transcript_split[2])
        ens_id = transcript_split[3]

        if transcript_split[5] == "+":
            pos = start
        else:
            pos = end

        if (line_count + 1) % 1000 == 0:
            sys.stderr.write("Thread " + str(thread_number) + "\t" + str(line_count + 1) + " positions analyzed\n")
            sys.stderr.flush()
        line_count += 1

        if len(transcript_split) < 6 or transcript_split[5] == '+' or transcript_split[5] == '.':
            forward = True
        elif transcript_split[5] == '-':
            forward = False
        if chrom + "_" + str(pos) in tss_visited:
            skipped += 1
            continue
        if pos - args['start'] < 1:
            print("WARNING pos - start: ", transcript_split)
            skipped += 1
            continue
        if chrom.find("_") != -1:
            print("WARNING find _: ", transcript_split)
            skipped += 1
            continue
        if chrom == "chrY":
            print("WARNING chrY: ", transcript_split)
            skipped += 1
            continue
        processed += 1


        cov_start = pos - args['start']
        cov_end = pos + args['end'] + 1

        coverage_tuple = bam.count_coverage_c60(chrom, cov_start - 220, cov_end + 220,
                                                quality_threshold=args['mapping_quality'], gc_corr=args["gc_corr"])
        #print(chrom, cov_start, cov_end)
        coverage_list = list()
        for i in range(0, len(coverage_tuple[0])):
            coverage_list.append(
                coverage_tuple[0][i])
        # normalize for copy-number variations
        if args['normalisation_file']:
            cnv_data = readCopyData(args['normalisation_file'])
            normcov = calcLocalMean(pos, chrom, cnv_data)
            if normcov == 0:
                normcov = 0.001
        counter = 0

        if forward:
            coverage_dict[ens_id + "\t" + chrom + "_" + str(start)] = coverage_list
        elif not forward:
            coverage_dict[ens_id + "\t" + chrom + "_" + str(end)] = coverage_list[::-1]

        for i in range(pos - args['start'], pos + args['end'] + 1):
            coverage = float(coverage_list[counter])  # / args['mean_coverage']
            if not coverage:
                counter += 1
                continue
            if coverage > args['coverage_limit']:
                counter += 1
                continue
            if args['normalisation_file']:
                coverage /= normcov
            if forward:
                position_lists[i - pos] += coverage
                position_ctrs[i - pos] += 1
            elif not forward:
                position_lists[-(i - pos)] += coverage
                position_ctrs[-(i - pos)] += 1
            counter += 1

        tss_visited.append(chrom + "_" + str(pos))

    sys.stderr.write("Thread " + str(thread_number) + "\t finished\n")
    sys.stderr.flush()
    q.put((position_lists, position_ctrs, coverage_dict))


def fetch_file_paths(input_dir, pattern_file_name):
    logging.info("Fetching file paths")
    input_dir_files = glob.glob(os.path.join(input_dir, pattern_file_name))
    return sorted(input_dir_files)


def main(num_content=None):
    # Init stuff
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-bf', '--bam_file', type=str, required=True, help='Sample BAM file')
    parser.add_argument('-nf', '--normalisation_file', type=str, required=False,
                        help='Normalise by local CN from this file')
    parser.add_argument('-ibe', '--input_dir_beds', type=str, required=True,
                        help='Directory containing genomic regions files. Script expects header line for all extensions'
                             'but \".bed\"')
    parser.add_argument('-alg', '--all_genes', type=str, required=False,
                        help='Reference file for all genes')
    parser.add_argument('-s', '--start', type=int, required=False, default=1000,
                        help='Start analyzing coverage at this point before region of interest (default: %(default)d)')
    parser.add_argument('-e', '--end', type=int, required=False, default=1000,
                        help='Stop analyzing coverage at this point after region of interest (default: %(default)d)')
    parser.add_argument('-mc', '--mean_coverage', type=float, required=False, default=1,
                        help='Mean coverage along the genome (default: %(default)d)')
    parser.add_argument('-cl', '--coverage_limit', type=float, required=False, default=1000,
                        help='discard coverage values above this threshold (default: %(default)d)')
    parser.add_argument('-mr', '--max_regions', type=int, required=False, default=0,
                        help='Use a maximum of this amount of regions (default: %(default)d)')
    parser.add_argument('-mq', '--mapping_quality', type=int, required=False, default=0,
                        help='Only count coverage of reads with this mapping quality (default: %(default)d)')
    parser.add_argument('-c', '--cpus', type=int, required=False, default=1,
                        help='Number of CPUs available')
    parser.add_argument('-o', '--output_dir', type=str, required=True,
                        help='Directory to output sample means')
    parser.add_argument('-on', '--output_name', type=str, required=True,
                        help='name of sample')
    parser.add_argument('-gc', '--gc_corr', type=bool, required=False, default=False)
    args = parser.parse_args()
    init_logger()
    log_input_arguments(args)
    args = vars(args)

    # Main code
    bed_paths = fetch_file_paths(args['input_dir_beds'], '*')
    # bed_paths = [bed_path for bed_path in bed_paths if 'Stromal-A' in bed_path]
    for bed_num, bed_path in enumerate(bed_paths):
        # Collect gene data
        bed_path = Path(bed_path)
        extension = bed_path.suffix
        header_map = dict()  # dictionary of header indices
        try:
            with open(bed_path, "rt") as f_in_bed:
                if extension != '.bed':
                    header = f_in_bed.readline().rstrip().split('\t')
                    # parse header
                    header_map['component'] = header.index('component')
                    header_map['chr'] = header.index('seqname')
                    header_map['peak'] = header.index('summit')
                dhss_content = f_in_bed.readlines()
        except IOError:
            logging.error("Failed to open specified files")
            logging.error(bed_path)
            sys.exit(-1)

        bed_name = Path(bed_path).stem

        ctr = 0
        target_content = list()
        # filter out content from file
        for target in dhss_content:
            target_split = target.rstrip().split('\t')

            chrom = target_split[0]

            if len(target_split) < 3:
                continue
            if chrom.find("_") != -1:
                continue
            if chrom.find("Y") != -1:
                continue

            """
            if not header_map or target_component == bed_name:
                target_content.append(target_split)
                ctr += 1
            # check if component name corresponds to file name
            else:
                sys.stderr.write("Target component different " + str(target_component) + ", " + str(bed_name) + "\n")
                sys.stderr.flush()
            """
            target_content.append(target_split)
            ctr += 1

            if num_content and num_content == ctr:
                break
        logging.info(f"{bed_name} {len(target_content)} ({bed_num + 1}/{len(bed_paths)}) BEDs")
        print(len(target_content))
        if args['cpus'] > len(target_content):
            sys.stderr.write(
                'Reducing number of processes from ' + str(args['cpus']) + ' to ' + str(len(target_content)) + '\n')
            sys.stderr.write("--------------------------------------------------\n")
            sys.stderr.flush()
            args['cpus'] = len(target_content)

        # initialize input data
        thread_inputs = dict()
        thread_coverages = dict()
        thread_coverage_ctrs = dict()
        thread_coverage_dicts = dict()

        for thread_number in range(args['cpus']):
            thread_inputs[thread_number] = list()
            thread_coverages[thread_number] = dict()
            thread_coverage_ctrs[thread_number] = dict()
            thread_coverage_dicts[thread_number] = dict()

        # dispatch input data
        max_targets = len(target_content)
        partition_point = max_targets / args["cpus"]
        logging.info(partition_point)
        for thread_number in range(args['cpus']):
            thread_inputs[thread_number] = target_content[int(thread_number * partition_point):int(
                (thread_number + 1) * partition_point)]
        sys.stderr.write(str(len(thread_inputs[list(thread_inputs)[-1]])) + " positions per thread\n")
        sys.stderr.write("--------------------------------------------------\n")
        sys.stderr.flush()

        # start multiple processes
        processes = dict()
        queues = dict()
        for thread in range(args['cpus']):
            queues[thread] = multiprocessing.Queue()
            if header_map:
                processes[thread] = multiprocessing.Process(target=pysam_presum_header, args=(
                    queues[thread], thread, thread_inputs[thread], header_map, args))
            else:
                processes[thread] = multiprocessing.Process(target=pysam_presum, args=(
                    queues[thread], thread, thread_inputs[thread], args))
            processes[thread].start()
        # wait for processes to finish
        for thread in range(args['cpus']):
            thread_coverages[thread], thread_coverage_ctrs[thread], thread_coverage_dicts[thread] = queues[thread].get()
        for thread in range(args['cpus']):
            processes[thread].join()

        # collect all data
        cov_dict = dict()
        coverages = dict()
        ctrs = dict()
        for loc in range(-args['start'], args['end'] + 1):
            coverages[loc] = 0
            ctrs[loc] = 0
        for thread in range(args['cpus']):
            for loc in range(-args['start'], args['end'] + 1):
                coverages[loc] += thread_coverages[thread][loc]
                ctrs[loc] += thread_coverage_ctrs[thread][loc]
        for loc in range(-args['start'], args['end'] + 1):
            coverages[loc] /= 1

        for thread in range(0, args["cpus"]):
            for i in thread_coverage_dicts[thread].keys():
                cov_dict[i] = thread_coverage_dicts[thread][i]
        identifier = []
        for gene in cov_dict.keys():
            identifier.append(gene.split("\t"))

        data = pd.DataFrame(identifier, columns=["Gene", "TSS"])
        data = pd.concat([data, pd.DataFrame(numpy.array(list(cov_dict.values())))], axis=1)
        data.columns = ["Gene", "TSS"] + [str(i) for i in range(-args["start"] - 220, args["end"] + 220 + 1)]
        data = data.drop([str(i) for i in range(-args["start"] - 220, -args["start"] -1)], axis=1)
        data = data.drop([str(i) for i in range(args["end"] + 1, args["end"] + 220)], axis=1)
        trimmed_mean = stats.trim_mean(data[[str(i) for i in range(-args["start"], args["end"] + 1)]].mean(), 0.1)
        data[[str(i) for i in range(-args["start"], args["end"] + 1)]] = \
            data[[str(i) for i in range(-args["start"], args["end"] + 1)]] / trimmed_mean


        data.to_csv(args["output_dir"] + "/" + args["output_name"] + ".csv", index=False)


        # output data
        sys.stderr.write("--------------------------------------------------\n")
        sys.stderr.write(str(ctrs[0]) + " positions analyzed\n")
        sys.stderr.flush()
        Path(args['output_dir']).mkdir(parents=True, exist_ok=True)
        with open(Path(args['output_dir'], bed_name + '.txt'), 'w+') as f_out_component:
            print("Position\tMean Cov\tPositions analyzed", file=f_out_component)
            for loc in range(-args['start'], args['end'] + 1):
                print(str(loc) + "\t" + str(coverages[loc]) + "\t" + str(ctrs[loc]), file=f_out_component)


if __name__ == "__main__":
    main()
