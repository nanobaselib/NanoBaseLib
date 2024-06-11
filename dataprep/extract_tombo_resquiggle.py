import os, csv, re
import pandas as pd
import numpy as np
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter
from collections import defaultdict
from scipy import stats
from scipy.stats import norm
import h5py


base_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

POLYA_STANDARD_MU = 108.9
POLYA_STANDARD_SIGMA = 1.67


def smooth(y, box_pts):
    box = np.ones(box_pts) / box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth


def percentile(arr, q=5):
    """
    Remove < 5% and > 95%.
    :param arr:
    :param q:
    :return:
    """
    low_percen = np.percentile(arr, q)
    hign_percen = np.percentile(arr, 100 - q)
    filter_arr = []
    for x in arr:
        if low_percen <= x <= hign_percen:
            filter_arr.append(x)

    return np.array(filter_arr)



def read_transcript_fa(transcript_fasta_file_name):
    fasta = open(transcript_fasta_file_name, "r")
    entries, separate_by_pipe="", False
    for ln in fasta:
        entries += ln
    entries = entries.split(">")
    dict = {}
    for entry in entries:
        entry = entry.split("\n")
        if len(entry[0].split()) > 0:
            id = entry[0].split(' ')[0]
            seq = "".join(entry[1:])
            dict[id] = [seq]
    return dict


def convert_raw_signal_to_pA_value(raw_signal, offset, range, digitisation):
    """
    Transform the raw signal to pico-ampere current values.
      𝑠𝑖𝑔𝑛𝑎𝑙_𝑖𝑛_𝑝𝑖𝑐𝑜_𝑎𝑚𝑝𝑒𝑟𝑒 = (𝑟𝑎𝑤_𝑠𝑖𝑔𝑛𝑎𝑙_𝑣𝑎𝑙𝑢𝑒 + 𝑜𝑓𝑓𝑠𝑒𝑡) ∗ 𝑟𝑎𝑛𝑔𝑒 / 𝑑𝑖𝑔𝑖𝑡𝑖𝑠𝑎𝑡𝑖𝑜𝑛
    :param raw_signal:
    :param offset:
    :param range:
    :param digitisation:
    :return: signal_in_pico_ampere
    """
    return (raw_signal + offset) * range / digitisation


def complement(seq):
    _seq = ""
    for base in seq:
        _seq += base
    return _seq


def reverse_complement(seq):
    return complement(seq)[::-1]


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
    if kmer == "NNNNNN":
        return None
    mu = float(mu)
    try:
        pdf = norm.logpdf(mu, loc=model_kmer[kmer]['model_mean'], scale=model_kmer[kmer]['model_stdv'])
    except:
        pdf = None

    return pdf


def extract_reads(fast5_folder, tombo_file):

    mean = []
    pdf = []

    all_reads = os.listdir(fast5_folder)
    for read_file in all_reads:
        fast5_file = os.path.join(fast5_folder, read_file)
        read_name = read_file.split(".")[0]
        if read_name not in polya_start_dict: continue
        with h5py.File(fast5_file, 'r') as fast5_data:
            events_object = fast5_data['Analyses']['RawGenomeCorrected_000']['BaseCalled_template']
            keys = [i for i in events_object]
            if keys != ['Alignment', 'Events']:
                continue
            print(read_file)
            # get pa current data
            channel_id = fast5_data['UniqueGlobalKey']['channel_id']
            digitisation = channel_id.attrs['digitisation']
            offset = channel_id.attrs['offset']
            range = channel_id.attrs['range']
            Reads = fast5_data['Raw']['Reads']
            for key in Reads:
                raw_signal = Reads[key]['Signal'][:]
                signal_in_pico_ampere = convert_raw_signal_to_pA_value(raw_signal, offset, range, digitisation)

                signal_pts = 50
                polyA_start = polya_start_dict[read_name]
                polyA_end = polya_end_dict[read_name]

                polyA_len = polyA_end - polyA_start
                polyA_sig = signal_in_pico_ampere[polyA_start: polyA_end]

                polyA_extract_start = int(polyA_len * 0.25)
                polyA_extract_end = int(polyA_len * 0.75)
                polyA_extract_sig = polyA_sig[polyA_extract_start: polyA_extract_end]

                polya_extract_smooth_sig = smooth(polyA_sig, signal_pts)[polyA_extract_start: polyA_extract_end]

                polyA_mu = np.mean(polya_extract_smooth_sig)
                polya_res_arr = polyA_extract_sig - polya_extract_smooth_sig
                polyA_sigma = np.std(polya_res_arr)

                ployA_sigma_percentile = np.std(percentile(polya_res_arr))
                normalized_signal = ((signal_in_pico_ampere - polyA_mu) / ployA_sigma_percentile) * POLYA_STANDARD_SIGMA + \
                                        POLYA_STANDARD_MU
                normalized_signal = np.around(np.array(normalized_signal), decimals=3)

            # read tombo res
            mapped_chrom = events_object['Alignment'].attrs['mapped_chrom']
            mapped_end = events_object['Alignment'].attrs['mapped_end']
            mapped_start = events_object['Alignment'].attrs['mapped_start']
            mapped_strand = events_object['Alignment'].attrs['mapped_strand']

            read_start_rel_to_raw = events_object['Events'].attrs['read_start_rel_to_raw']
            events = events_object['Events'][:]
            df = pd.DataFrame(events)
            # df = df[['norm_mean', 'start', 'length', 'base']]
            df['base'] = df['base'].apply(lambda x: x.decode("utf-8"))
            signal_len = df.iloc[-1]['start'] + df.iloc[-1]['length']

            # read_start_rel_to_raw corresponding the raw signal
            # raw_signal = signal_in_pico_ampere[read_start_rel_to_raw: read_start_rel_to_raw + signal_len]
            # raw_signal = raw_signal[::-1]
            #
            # read_start_rel_to_raw corresponding the reversed raw signal
            normalized_signal = normalized_signal[::-1]
            normalized_signal = normalized_signal[read_start_rel_to_raw: read_start_rel_to_raw + signal_len]

            ref_seq = "".join(list(df['base']))

            def get_mean(start, len, signal):
                _signal = signal[start: start + len]
                return np.mean(_signal)

            def get_stdv(start, len, signal):
                _signal = signal[start: start + len]
                return np.std(_signal)

            def get_kmer(i):
                shift = 2
                kmer = ref_seq[i-shift: i+5-shift]
                if base_type == "rna":
                    if len(kmer) != 5:
                        return 'NNNNN'
                    else:
                        return kmer

            df['mean'] = df.apply(lambda x: get_mean(x.start, x.length, normalized_signal), axis=1)
            df['stdv'] = df.apply(lambda x: get_stdv(x.start, x.length, normalized_signal), axis=1)
            df['ref_kmer'] = df.index
            df['ref_kmer'] = df['ref_kmer'].apply(lambda x: get_kmer(x))
            df['model_kmer'] = df['ref_kmer'].apply(lambda x: x[::-1])

            df['pdf'] = df.apply(lambda x: cal_pdf(x['mean'], x['model_kmer']), axis=1)

            df = df.dropna()
            df = df[df['pdf'] > -8]

            mean.append(df.loc[:, 'stdv'].mean())
            pdf.append(df.loc[:, 'pdf'].mean())
            # print("Avg. stdv ", df.loc[:, 'stdv'].mean())
            # print("Avg. logpdf ", df.loc[:, 'pdf'].mean())

        if len(mean) == 100:

            print("\n   **** ")
            print("Avg. stdv ", np.mean(mean))
            print("Avg. logpdf ", np.mean(pdf))
            break



def arg_parser():
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter, add_help=False)
    parser.add_argument('--fast5_folder',
                        default='/scratch/cs/infantbiome/chengg1/datasets/hek293t/tombo_test_wt1/0')
    parser.add_argument('--ploya_file',
                        default='/scratch/cs/infantbiome/chengg1/datasets/hek293t/3_tailfindr/'
                                'HEK293T-WT-rep1/HEK293T-WT-rep1_tailfindr_pass.csv')
    parser.add_argument('--save_folder',
                        default='/scratch/cs/nanopore/chengg1/benchmark_datasets/NA12878/3_tombo/FAB42828')
    parser.add_argument('--base_type', default='rna')
    return parser


if __name__ == '__main__':

    args = arg_parser().parse_args()
    base_type = args.base_type
    print(f" ========== {base_type} ========== ")

    # step 1: read model kmer
    model_kmer = read_model_kmer(base_type + "_model_kmer.csv")
    polya_df = pd.read_csv(args.ploya_file).dropna()
    print(polya_df)
    polya_df['tail_start'] = polya_df['tail_start'].apply(int)
    polya_df['tail_end'] = polya_df['tail_end'].apply(int)
    polya_start_dict = polya_df.set_index("read_id")['tail_start'].to_dict()
    polya_end_dict = polya_df.set_index("read_id")['tail_end'].to_dict()

    # step 2:
    save_tombo_file = os.path.join(args.save_folder, "tombo.txt")
    rm_exist_file(save_tombo_file)

    # step 3:
    extract_reads(args.fast5_folder, save_tombo_file)

    # nanopolish_df = pd.read_csv(args.file_name, sep=' ')
    # if base_type == 'dna':
    #     nanopolish_df = nanopolish_df[nanopolish_df['model_kmer'] != 'NNNNNN']
    # if base_type == 'rna':
    #     nanopolish_df = nanopolish_df[nanopolish_df['model_kmer'] != 'NNNNN']
    #
    # print("Avg. stdv ", nanopolish_df.loc[:, 'event_stdv'].mean())
    # nanopolish_df['pdf'] = nanopolish_df.apply(lambda x: cal_pdf(x['event_level_mean'],
    #                                                              x['model_kmer']), axis=1)
    # print("  ***   ")
    # nanopolish_df = nanopolish_df.dropna()
    # print("Avg. logpdf ", nanopolish_df.loc[:, 'pdf'].mean())













