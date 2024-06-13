import os, csv, re
import pandas as pd
import numpy as np
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter
from collections import defaultdict
import h5py


def rm_exist_file(file_name):

    if os.path.exists(file_name):
        os.remove(file_name)


def write_row_to_file(file_name, row, delimiter=','):
    with open(file_name, 'a+', newline="", encoding='utf-8') as f:
        writer = csv.writer(f, delimiter=delimiter)
        writer.writerow(row)


def extract_reads(fast5_folder, resquiggle_file, summary_file):

    all_reads = os.listdir(fast5_folder)
    for read_file in all_reads:
        fast5_file = os.path.join(fast5_folder, read_file)
        read_name = read_file.split(".")[0]
        with h5py.File(fast5_file, 'r') as fast5_data:
            GenomeCorrected_object = fast5_data['Analyses']['RawGenomeCorrected_000']
            keys = [i for i in GenomeCorrected_object]
            if keys != ['BaseCalled_template']:
                print(f"{read_name} tombo resquiggle fails.")
                continue
            events_object = GenomeCorrected_object['BaseCalled_template']
            keys = [i for i in events_object]
            if keys != ['Alignment', 'Events']:
                print(f"{read_name} tombo resquiggle fails.")
                continue

            # read all attrs
            lower_lim = events_object.attrs['lower_lim']
            upper_lim = events_object.attrs['upper_lim']
            norm_type = events_object.attrs['norm_type']
            outlier_threshold = events_object.attrs['outlier_threshold']
            rna = events_object.attrs['rna']
            scale = events_object.attrs['scale']
            shift = events_object.attrs['shift']
            signal_match_score = events_object.attrs['signal_match_score']
            status = events_object.attrs['status']

            clipped_bases_end = events_object['Alignment'].attrs['clipped_bases_end']
            clipped_bases_start = events_object['Alignment'].attrs['clipped_bases_start']
            mapped_chrom = events_object['Alignment'].attrs['mapped_chrom']
            mapped_strand = events_object['Alignment'].attrs['mapped_strand']
            mapped_start = events_object['Alignment'].attrs['mapped_start']
            mapped_end = events_object['Alignment'].attrs['mapped_end']
            num_deletions = events_object['Alignment'].attrs['num_deletions']
            num_insertions = events_object['Alignment'].attrs['num_insertions']
            num_matches = events_object['Alignment'].attrs['num_matches']
            num_mismatches = events_object['Alignment'].attrs['num_mismatches']

            read_start_rel_to_raw = events_object['Events'].attrs['read_start_rel_to_raw']
            read_summary = [read_name, read_start_rel_to_raw, mapped_chrom, mapped_strand, mapped_start, mapped_end,
                            num_matches, num_mismatches, num_insertions, num_deletions, clipped_bases_start,
                            clipped_bases_end,  rna, scale, shift, signal_match_score, status, lower_lim, upper_lim,
                            norm_type, outlier_threshold]
            write_row_to_file(summary_file, read_summary, '\t')

            # read resquiggle res
            events = events_object['Events'][:]
            df = pd.DataFrame(events)
            df['base'] = df['base'].apply(lambda x: x.decode("utf-8"))
            df['read_name'] = read_name
            df = df[['read_name', 'norm_mean', 'norm_stdev', 'start', 'length', 'base']]
            df.to_csv(resquiggle_file, index=False, header=False, mode='a', sep='\t')


def arg_parser():
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter, add_help=False)
    parser.add_argument('--root_dir',
                        default='/scratch/cs/nanopore/chengg1/benchmark_datasets/lambda_phage/3_tombo/VER5940')
    parser.add_argument('--folder',
                        default='0')
    return parser


if __name__ == '__main__':

    args = arg_parser().parse_args()

    save_resquiggle_file = os.path.join(args.root_dir, "tombo_resquiggle.txt")
    save_summary_file = os.path.join(args.root_dir, "tombo_summary.txt")
    rm_exist_file(save_resquiggle_file)
    rm_exist_file(save_summary_file)
    reads_fast5_folder = os.path.join(args.root_dir, args.folder)

    resquiggle_file_header = ['read_name', 'norm_mean', 'norm_stdev', 'start', 'length', 'base']
    summary_file_header = ['read_name', 'read_start_rel_to_raw', 'mapped_chrom', 'mapped_strand', 'mapped_start',
                           'mapped_end', 'num_matches', 'num_mismatches', 'num_insertions', 'num_deletions',
                           'clipped_bases_start',  'clipped_bases_end',  'rna', 'scale', 'shift', 'signal_match_score',
                           'status', 'lower_lim', 'upper_lim',  'norm_type', 'outlier_threshold']

    write_row_to_file(save_resquiggle_file, resquiggle_file_header, '\t')
    write_row_to_file(save_summary_file, summary_file_header, '\t')

    extract_reads(reads_fast5_folder, save_resquiggle_file, save_summary_file)
