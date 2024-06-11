import os, csv, re
import pandas as pd
import numpy as np
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter
from collections import defaultdict
from scipy.stats import norm

base_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}


def complement(seq):
    _seq = ""
    for base in seq:
        _seq += base
    return _seq


def read_model_kmer(file_name):
    """
    Read model kmer and return a dict.
    :param file_name:
    :return:
    """
    model_kmer_pd = pd.read_csv(file_name)
    model_kmer_pd['model_mean'] = model_kmer_pd['model_mean'].apply(float)
    model_kmer_pd['model_stdv'] = model_kmer_pd['model_stdv'].apply(float)
    model_kmer_pd = model_kmer_pd.to_dict('records')
    model_kmer_dict = dict()
    for item in model_kmer_pd:
        model_kmer_dict[item['model_kmer']] = item
    return model_kmer_dict


def rm_exist_file(file_name):

    if os.path.exists(file_name):
        os.remove(file_name)


def get_file_lines(file_name):
    count = 0
    file = open(file_name, 'r', encoding='utf-8')
    while 1:
        buffer = file.read(8 * 1024 * 1024)  # 可大概设置
        if not buffer:
            break
        count += buffer.count('\n')
    file.close()
    return count


def write_row_to_file(file_name, row, delimiter=','):
    with open(file_name, 'a+', newline="", encoding='utf-8') as f:
        writer = csv.writer(f, delimiter=delimiter)
        writer.writerow(row)



def arg_parser():
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter, add_help=False)
    parser.add_argument('--file_name',
                        default='/scratch/cs/infantbiome/chengg1/datasets/hek293t_bc_benchmark')
    parser.add_argument('--base_type', default='dna')
    return parser


def cal_pdf(mu, kmer):
    mu = float(mu)
    try:
        pdf = norm.logpdf(mu, loc=model_kmer[kmer]['model_mean'],
                          scale=model_kmer[kmer]['model_stdv'])
    except:
        pdf = None

    return pdf


if __name__ == '__main__':
    """
    cat eventalign_combined.txt | awk -F "\t" '{print $7,$8,$10}' > nanopolish_test.txt
    """

    args = arg_parser().parse_args()
    base_type = args.base_type
    print(f" ========== {base_type} ========== ")
    model_kmer = read_model_kmer(base_type + "_model_kmer.csv")

    nanopolish_df = pd.read_csv(args.file_name, sep=' ')
    if base_type == 'dna':
        nanopolish_df = nanopolish_df[nanopolish_df['model_kmer'] != 'NNNNNN']
    if base_type == 'rna':
        nanopolish_df = nanopolish_df[nanopolish_df['model_kmer'] != 'NNNNN']

    print("Avg. stdv ", nanopolish_df.loc[:, 'event_stdv'].mean())
    nanopolish_df['pdf'] = nanopolish_df.apply(lambda x: cal_pdf(x['event_level_mean'],
                                                                 x['model_kmer']), axis=1)
    print("  ***   ")
    nanopolish_df = nanopolish_df.dropna()
    print("Avg. logpdf ", nanopolish_df.loc[:, 'pdf'].mean())













