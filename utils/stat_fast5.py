import os
import csv
import numpy as np
import pandas as pd
from tqdm import tqdm
import argparse
import h5py
import sys
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


def read_polya_results_tsv(polya_results_file_name):
    """
    Read PolyA file
    :param polya_results_file_name:
    :return:
    """
    polya_results_pd = pd.read_csv(polya_results_file_name, sep='\t').to_dict('records')
    polya_results_dict = dict()
    for item in polya_results_pd:
        polya_results_dict[item['readname']] = item
    return polya_results_dict


def read_event_summary_file(eventalign_summary_file, trans2idx_dict):
    """
   ['order', 'contig', 'read_index', 'strand',
    'trans_position_start', 'trans_position_end', 'start_idx', 'end_idx',
    'read_name', 'fast5_path', 'ref']
    """

    data = pd.read_csv(eventalign_summary_file)
    data['read_index'] = data['read_index'].apply(lambda x: int(x))
    data['contig'] = data['contig'].apply(lambda x: trans2idx_dict[x])
    idx2fast5_dict = data.set_index('read_index')['fast5_path'].to_dict()
    idx2readName_dict = data.set_index('read_index')['read_name'].to_dict()
    idx2fullInfo_arr = data.to_dict("records")
    idx2fullInfo_dict = dict()
    for key in idx2fullInfo_arr:
        idx2fullInfo_dict[key['read_index']] = key
    return idx2fast5_dict, idx2readName_dict, idx2fullInfo_dict


def read_filter_file(filter_file_name):
    """
    ["read_name", "ref_idx", "flag", "reference_start", "reference_end"]
    """
    data = pd.read_csv(filter_file_name)
    data['ref_idx'] = data['ref_idx'].apply(lambda x: str(x))
    data['key'] = data['read_name'] + "#" + data['ref_idx']
    return list(data['key'])


def generate_trans_dict(file_name):
    data = pd.read_csv(file_name)
    data['index'] = data['index'].apply(lambda x: int(x))
    trans2idx_dict = data.set_index('transcript')['index'].to_dict()
    idx2trans_dict = data.set_index('index')['transcript'].to_dict()
    return trans2idx_dict, idx2trans_dict


def write_to_file(file_name, data):
    with open(file_name, 'a+', newline="") as f:
        writer = csv.writer(f)
        for row in data:
            writer.writerow(row)


def write_event_combine_file(summary_dict, read_index, contig, read_data, save_dir):

    curent_read_name = summary_dict[int(read_index)]["read_name"]
    _read_file_name = os.path.join(save_dir, read_index + '#' + curent_read_name + "#" + contig + '.txt')
    if os.path.exists(_read_file_name):
        return
        # os.remove(_read_file_name)
    with open(_read_file_name, 'w', newline="") as f:
        writer = csv.writer(f, delimiter='\t')
        for _data in read_data:
            writer.writerow(_data)


def split_eventalign(eventalign_file_name, save_dir, filter_arr, idx2fast5_dict, idx2readname_dict):
    """
    eventalign_file:
               0         1        2      3       4          5           6
         ['read_idx', 'contig', 'pos', 'kmer', 'mean', 'start_idx', 'end_idx']
    """
    print('* Start to processing ... ')
    num_lines = sum(1 for _ in open(eventalign_file_name, 'r'))

    with open(eventalign_file_name, 'r') as f:
        for i, line in enumerate(tqdm(f, total=num_lines)):
            line = line.strip().replace('\n', '').replace('\r', '').split('\t')
            if i == 0:
                last_read_idx = int(line[0])
                last_contig = int(line[1])
                last_pos = int(line[2])
                last_kmer = int(line[3])
                last_mean = float(line[4])
                last_start_idx = int(line[5])
                last_end_idx = int(line[6])
                curr_read_data = [[last_read_idx, last_contig, last_pos, last_kmer, last_mean, last_start_idx,
                                   last_end_idx]]
            else:
                if last_read_idx != int(line[0]):
                    fast5_name = idx2fast5_dict[last_read_idx]
                    read_name = idx2readname_dict[last_read_idx]
                    save_file_name = os.path.join(save_dir, fast5_name + '.csv')
                    if read_name + "#" + str(last_contig) not in filter_arr:
                        write_to_file(save_file_name, curr_read_data)
                    curr_read_data = list()

                last_read_idx = int(line[0])
                last_contig = int(line[1])
                last_pos = int(line[2])
                last_kmer = int(line[3])
                last_mean = float(line[4])
                last_start_idx = int(line[5])
                last_end_idx = int(line[6])
                curr_read_data.append([last_read_idx, last_contig, last_pos, last_kmer, last_mean, last_start_idx,
                                       last_end_idx])

            # the last line
            if i == num_lines - 1:
                fast5_name = idx2fast5_dict[last_read_idx]
                read_name = idx2readname_dict[last_read_idx]
                save_file_name = os.path.join(save_dir, fast5_name + '.csv')
                if read_name + "#" + str(last_contig) not in filter_arr:
                    write_to_file(save_file_name, curr_read_data)


def generate_event_data_dict(tmp_event_file_name, idx2fullInfo_dict):
    """
    :param idx2fullInfo_dict:  ['order', 'contig', 'read_index', 'strand',
                                'trans_position_start', 'trans_position_end', 'start_idx', 'end_idx',
                                'read_name', 'fast5_path', 'ref']
    :return:
    """
    tmp_event_dict = dict()
    tmp_event_data = pd.read_csv(tmp_event_file_name, header=None)
    tmp_event_data.columns = ['read_idx', 'contig', 'pos', 'kmer', 'mean', 'start_idx', 'end_idx']
    tmp_event_data = tmp_event_data.groupby(['read_idx', 'contig'])
    for site, group in tmp_event_data:
        _read_index = site[0]
        _read_info_dict = idx2fullInfo_dict[_read_index]
        _read_name = _read_info_dict['read_name']
        group = group[['pos', 'kmer', 'mean', 'start_idx', 'end_idx']]
        _read_info_dict['event'] = group.values.tolist()
        tmp_event_dict[_read_name] = _read_info_dict
    return tmp_event_dict


# def check_fast5(fast5_dir):
#     all_files = os.listdir(fast5_dir)
#     all_files = ['FAK01563_771a38df9d9b678a485923dc79a6a5d5d59c115c_17.fast5']
#     for fast5_file_name in all_files:
#         print(10 * "-", fast5_file_name, 10 * "-")
#         fast5_file_full_name = os.path.join(fast5_dir, fast5_file_name)
#         with h5py.File(fast5_file_full_name, 'r+') as fast5_data:
#             for read_key in fast5_data:
#                 if "NanopolishEvent" in fast5_data[read_key]:
#                     print("event : ", read_key)
#                     del fast5_data[read_key]["NanopolishEvent"]
#
#                 if "Normalized" in fast5_data[read_key]:
#                     print("normal : ", read_key)
#                     del fast5_data[read_key]["Normalized"]


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


def normalize_fast5_with_polyA(signal_in_pico_ampere, polyA_start, polyA_end):

    signal_pts = 50

    polyA_len = polyA_end - polyA_start
    polyA_sig = signal_in_pico_ampere[polyA_start: polyA_end]

    polyA_extract_start = int(polyA_len * 0.25)
    polyA_extract_end = int(polyA_len * 0.75)
    polyA_extract_sig = polyA_sig[polyA_extract_start: polyA_extract_end]

    polya_extract_smooth_sig = smooth(polyA_sig, signal_pts)[polyA_extract_start: polyA_extract_end]

    polyA_mu = np.mean(polya_extract_smooth_sig)
    polya_res_arr = polyA_extract_sig - polya_extract_smooth_sig
    polyA_sigma = np.std(polya_res_arr)

    if polyA_sigma < 3.0:
        ployA_sigma_percentile = np.std(percentile(polya_res_arr))
        normalized_signal = ((signal_in_pico_ampere - polyA_mu) / ployA_sigma_percentile) * POLYA_STANDARD_SIGMA + \
                           POLYA_STANDARD_MU
        normalized_signal = np.around(np.array(normalized_signal), decimals=3)
        return normalized_signal, polyA_mu, ployA_sigma_percentile
    else:
        return None, None, None


def write_and_normalize_fast5(fast5_file, split_event_file,
                              idx2fullInfo_dict, idx2trans_dict, polya_results_dict):

    tmp_event_dict = generate_event_data_dict(split_event_file, idx2fullInfo_dict)
    print(" * generate event dict over !")

    with h5py.File(fast5_file, 'r+') as fast5_data:
        for read_key in fast5_data:
            read_name = fast5_data[read_key]["Raw"].attrs['read_id'].decode("utf-8")
            if read_name in tmp_event_dict:
                event_data = fast5_data[read_key].create_group("NanopolishEvent")
                event_data['Events'] = tmp_event_dict[read_name]['event']
                event_data['Reference'] = tmp_event_dict[read_name]['ref']
                event_data.attrs.create("read_index", tmp_event_dict[read_name]['read_index'], shape=None, dtype=None)
                event_data.attrs.create("contig_index", tmp_event_dict[read_name]['contig'], shape=None, dtype=None)
                # event_data.attrs.create("contig_name", idx2trans_dict[int(tmp_event_dict[read_name]['contig'])],
                #                         shape=None, dtype="S")
                event_data.attrs.create("order", tmp_event_dict[read_name]['order'], shape=None, dtype=None)
                # event_data.attrs.create("strand", tmp_event_dict[read_name]['strand'], shape=None, dtype="S")
                event_data.attrs.create("trans_position_start", tmp_event_dict[read_name]['trans_position_start'],
                                        shape=None, dtype=None)
                event_data.attrs.create("trans_position_end", tmp_event_dict[read_name]['trans_position_end'],
                                        shape=None, dtype=None)
                event_data.attrs.create("start_idx", tmp_event_dict[read_name]['start_idx'], shape=None, dtype=None)
                event_data.attrs.create("end_idx", tmp_event_dict[read_name]['end_idx'], shape=None, dtype=None)

            if read_name in polya_results_dict:
                polyA_start = int(polya_results_dict[read_name]['polya_start'])
                polyA_end = int(polya_results_dict[read_name]['transcript_start'])
                raw_signal = fast5_data[read_key]["Raw"]['Signal'][:]

                # read channel attributes parameters
                channel_id = fast5_data[read_key]['channel_id']
                digitisation = channel_id.attrs['digitisation']
                offset = channel_id.attrs['offset']
                range = channel_id.attrs['range']

                # transform the raw signal to pico-ampere current values
                signal_in_pico_ampere = convert_raw_signal_to_pA_value(raw_signal, offset, range, digitisation)

                # normalized pico-ampere current values
                normalized_signal, polyA_mu, ployA_sigma_percentile = \
                    normalize_fast5_with_polyA(signal_in_pico_ampere, polyA_start, polyA_end)

                if normalized_signal is not None:
                    norm_data = fast5_data[read_key].create_group("Normalized")
                    norm_data['Signal'] = normalized_signal
                    norm_data.attrs.create("polyA_start", polyA_start, shape=None, dtype=None)
                    norm_data.attrs.create("polyA_end", polyA_end, shape=None, dtype=None)
                    norm_data.attrs.create("polyA_mu", polyA_mu, shape=None, dtype=None)
                    norm_data.attrs.create("polyA_sigma", ployA_sigma_percentile, shape=None, dtype=None)


def count_reads_sum(file_name):

    n = 0
    with h5py.File(file_name, 'r') as fast5_data:
        for read_key in fast5_data:
            n += 1
    return  n



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Reads info statistics.')
    parser.add_argument('--root-dir', type=str, default='/scratch/cs/nanopore/chengg1/base_calling/dataset/dinopore_ivt/'
                                                        'synthetic/0_original/02_single_fast5')
    parser.add_argument('--type', type=str,
                        default='pureI')
    args = parser.parse_args()
    root_dir = args.root_dir
    all_folders = os.listdir(root_dir)

    num_fast5 = 0
    num_reads = 0

    for folder in all_folders:
        if args.type in folder:
            folder_dir = os.path.join(root_dir, folder, "multi")
            files = os.listdir(folder_dir)
            for file in files:
                if ".fast5" in file:
                    num_fast5 += 1
                    num_reads += count_reads_sum(os.path.join(folder_dir, file))

    print(f"{args.type}  fast5 num: {num_fast5},  reads = {num_reads}")
