import os, csv, re
import pandas as pd
import numpy as np
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter
from collections import defaultdict
from scipy.stats import norm

base_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}


def generate_kmer_dict(file_name):
    data = pd.read_csv(file_name)
    data['index'] = data['index'].apply(lambda x: int(x))
    kmer2idx_dict = data.set_index('kmer')['index'].to_dict()
    idx2kmer_dict = data.set_index('index')['kmer'].to_dict()
    return kmer2idx_dict, idx2kmer_dict


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






def cal_pdf(mu, kmer):
    mu = float(mu)
    return norm.logpdf(mu, loc=model_kmer[kmer]['model_mean'],
                        scale=model_kmer[kmer]['model_stdv'])


def cal_stdv():
    mu = []
    with open("sigma.csv", 'r') as f:
        for line in f:
            line = line.strip().split(",")
            line = [float(r) for r in line]
            mu.append(np.mean(np.array(line)))
    print(np.mean(np.array(mu)))


def arg_parser():
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter, add_help=False)
    parser.add_argument('--file_name', default='eventalign_v1.txt')
    parser.add_argument('--base_type', default='rna')
    return parser


if __name__ == '__main__':
    """
    cat file | awk -F "\t" '{print $7,$8,$10}' > nanopolish_test.txt
    """

    args = arg_parser().parse_args()
    base_type = args.base_type
    print(f" ========== {base_type} ========== ")

    # kmer2idx_dict, idx2kmer_dict = generate_kmer_dict("model_idx_kmer.csv")
    model_kmer = read_model_kmer(base_type + "_model_kmer.csv")

    df = pd.read_csv(args.file_name, header=None, sep='\t')
    df.columns = ['read_index', 'read_name', 'contig', 'model_kmer', 'position',
                  'event_level_mean', 'event_stdv', 'start_idx', 'end_idx', 'event_len']
    df = df[['event_level_mean', 'event_stdv', 'model_kmer']]
    df = df[df['event_level_mean'] > 0]
    df = df[df['model_kmer'] != '-']

    print("Avg. stdv ", df.loc[:, 'event_stdv'].mean())
    df = df[0: 80000]
    df['pdf'] = df.apply(lambda x: cal_pdf(x['event_level_mean'], x['model_kmer']), axis=1)
    print("  ***   ")
    df = df.dropna()

    print("Avg. logpdf ", df.loc[:, 'pdf'].mean())









