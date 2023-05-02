#!/usr/bin/env python3
from argparse import ArgumentParser
from sys import exit as sys_exit
from os import makedirs as os_mkdirs
from os.path import dirname as pth_dirname, basename as pth_basename, join as pth_join


def inject_cmdline_args():
    cmd_parser = ArgumentParser()
    cmd_parser.add_argument('-g', '--ensembl-id-genes', dest='ensembl_ids_path', required=True,
                            help='Path to list of genes represented as single ENSEMBL ID per line [ REQUIRED ]')
    cmd_parser.add_argument('-r', '--hg38-mane-genes-bed', dest='reference_bed_file_hg38', required=True,
                            help='Path to GRCh38 MANE transcript gene body BED regions file with content in the form '
                                 'of example "chr1	33306765	33321098	ENSG00000184389	.	-". Must contain '
                                 'strand! [ REQUIRED ]')
    cmd_parser.add_argument('-o', '--output-dir', dest='out_dir',
                            help="Optional: an output path can be defined. Will be created if it does not exist. If "
                                 "none is defined, the directory of the ENSEMBL gene ID list will be used.")
    cmd_parser.add_argument('-f', '--first-base-only', action='store_true', dest='reduce_to_first_base',
                            help="Flag: if set, the output BED file will contain only the first base of the gene, "
                                 "viewed in 5'->3' direction of the strand on which the gene is located.")
    return cmd_parser.parse_args().__dict__


if __name__ == '__main__':
    ensembl_ids_path = reference_bed_file_hg38 = out_dir = ''
    reduce_to_first_base = False
    # load commandline parameters
    locals().update(inject_cmdline_args())
    # load gene IDs which the user wants to be output
    with open(ensembl_ids_path, 'rt') as f_ids:
        ensembl_ids = {}.fromkeys([tr.strip() for tr in f_ids.readlines()])
    # load reference
    with open(reference_bed_file_hg38, 'rt') as f_ref:
        if reduce_to_first_base:
            ref_lines = [cont.split('\t') for cont in f_ref.readlines()]
        else:
            ref_lines = f_ref.readlines()
    # assemble output
    out_lines = []
    try:
        if reduce_to_first_base:
            for chrom, strt, stp, gid, dt, strnd, *_grbg in ref_lines:
                if gid in ensembl_ids:
                    out_lines.append('\t'.join([chrom, str(int(stp)-1), stp, gid, dt, strnd.rstrip()]) + '\n'
                                     if strnd.rstrip() == '-' else
                                     '\t'.join([chrom, strt, str(int(strt)+1), gid, dt, strnd.rstrip()]) + '\n')
        else:
            for gene_cont in ref_lines:
                if gene_cont.split('\t')[3] in ensembl_ids:
                    out_lines.append(gene_cont)
    except ValueError:  # if not enough or too many values to unpack
        print()
        sys_exit(1)
    # handle output path
    if not out_dir:
        out_dir = pth_dirname(ensembl_ids_path)
    os_mkdirs(out_dir, exist_ok=True)
    # write output file
    output_bed = pth_join(out_dir, pth_basename('.'.join(ensembl_ids_path.split('.')[:-1]) + '.bed'))
    with open(output_bed, 'wt') as f_bed_out:
        f_bed_out.writelines(out_lines)
    # prompt user feedback
    print(f" i :  successfully wrote {len(out_lines):,} gene entries to: '{output_bed}'. Exiting..")
