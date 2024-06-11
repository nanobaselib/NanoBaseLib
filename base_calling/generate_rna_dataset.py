import os, csv, sys, random, time, warnings
import h5py
import pandas as pd
import numpy as np
import multiprocessing
from io import StringIO
from shutil import copyfile
from typing import List, Tuple, Dict, Union
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter


class Consumer(multiprocessing.Process):
    """ For parallelisation """

    def __init__(self, task_queue, task_function, locks=None, result_queue=None):
        multiprocessing.Process.__init__(self)
        self.task_queue = task_queue
        self.locks = locks
        self.task_function = task_function
        self.result_queue = result_queue

    def run(self):
        proc_name = self.name
        while True:
            next_task_args = self.task_queue.get()
            if next_task_args is None:
                self.task_queue.task_done()
                break
            result = self.task_function(*next_task_args, self.locks)
            self.task_queue.task_done()
            if self.result_queue is not None:
                self.result_queue.put(result)


def end_queue(task_queue,n_processes):
    for _ in range(n_processes):
        task_queue.put(None)
    return task_queue


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


def read_summary(summary_file_name):
    """
    Read the eventalign summary.txt and return a dict.

    [read_index,  read_name,  fast5_path,  model_name,   strand, num_events,  num_steps,  num_skips,    num_stays,
     total_duration,  shift,  scale,  drift,  var]
    """
    eventalign_summary_df = pd.read_csv(summary_file_name, sep='\t')
    eventalign_summary_df = eventalign_summary_df[["read_index",  "read_name",  "fast5_path"]]
    eventalign_summary_df["read_index"] = eventalign_summary_df["read_index"].apply(int)
    return eventalign_summary_df


def read_transcript_fa(transcript_fasta_file_name):
    fasta = open(transcript_fasta_file_name, "r")
    entries, separate_by_pipe = "", False
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


def generate(eventalign_result: pd.DataFrame, min_seq_len: int, max_seq_len: int, shift: int,
             out_paths: Dict, sample_reads_dict: Dict, locks: Dict):
    '''
    Function to index the position of a specific read features within eventalign.txt

    Args:
        eventalign_result (pd.DataFrame): A pd.DataFrame object containing a portion of eventalign.txt to be indexed
        os_start (int): An index position within eventalign.txt that corresponds to the start of the eventalign_result
                        portion within eventalign.txt file
        out_paths (Dict): A dictionary containing filepath for all the output files produced by the index function
        locks (Dict): A lock object from multiprocessing library that ensures only one process write to the output file
                      at any given time
    Returns:

    reads_header = ["contig", "read_index", "contig_st",  "contig_en", "signal_st", "signal_en", "order", "first_kmer"]
    '''
    BASE_GAP_MAX_NUM = 5
    keys = ['contig', 'read_index']
    eventalign_result = eventalign_result.groupby(keys)
    with locks['reads'], open(out_paths['reads'], 'a', encoding='utf-8') as f_index:
        for keys, group in eventalign_result:
            group = group.reset_index()
            if len(group) < 2:
                continue
            contig_pos_arr = group['position'].to_numpy()
            contig_pos_prev = contig_pos_arr[0:-1]
            contig_pos_post = contig_pos_arr[1:]
            if np.max(contig_pos_post - contig_pos_prev) > BASE_GAP_MAX_NUM:
                continue
            contig, read_index = keys
            if contig in sample_reads_dict and read_index in sample_reads_dict[contig]:
                if len(sample_reads_dict[contig]) == 1:
                    contig_st = group['position'].min() + shift
                    contig_en = group['position'].max() + shift
                    signal_st = group['start_idx'].min()
                    signal_en = group['end_idx'].max()
                    if group.loc[0, 'start_idx'] > group.loc[1, 'start_idx']:
                        order = 1
                    else:
                        order = 0
                    first_kmer = group.loc[0, 'reference_kmer']
                    f_index.write(
                        '%s,%d,%d,%d,%d,%d,%d,%s\n' % (
                            contig, read_index, contig_st, contig_en, signal_st, signal_en, order, first_kmer))
                else:
                    df_st = random.randint(0, int(len(group)/2))
                    df_en = min(df_st + random.randint(min_seq_len, max_seq_len), len(group))
                    group = group[df_st: df_en].reset_index()
                    contig_st = group['position'].min() + shift
                    contig_en = group['position'].max() + shift
                    signal_st = group['start_idx'].min()
                    signal_en = group['end_idx'].max()
                    if group.loc[0, 'start_idx'] > group.loc[1, 'start_idx']:
                        order = 1
                    else:
                        order = 0
                    first_kmer = group.loc[0, 'reference_kmer']
                    f_index.write(
                        '%s,%d,%d,%d,%d,%d,%d,%s\n' % (
                            contig, read_index, contig_st, contig_en, signal_st, signal_en, order, first_kmer))


def parallel_generate(eventalign_combine_filepath: str, save_reads_filepath: str,  sample_reads_dict: Dict,
                      n_processes: int, chunk_size: int, min_seq_len: int, max_seq_len: int, shift: int):

    # check files
    if not os.path.exists(eventalign_combine_filepath):
        raise FileExistsError("Combined eventalign file does not exist!")
        sys.exit()
    rm_exist_file(save_reads_filepath)

    # Create output paths and locks.
    out_paths, locks = dict(), dict()
    for out_filetype in ['reads']:
        out_paths[out_filetype] = save_reads_filepath
        locks[out_filetype] = multiprocessing.Lock()

    # write the header
    write_row_to_file(out_paths['reads'], reads_header)

    # Create communication queues.
    task_queue = multiprocessing.JoinableQueue(maxsize=n_processes * 2)

    # Create and start consumers.p
    consumers = [Consumer(task_queue=task_queue, task_function=generate, locks=locks) for i in range(n_processes)]
    for process in consumers:
        process.start()

    ## Load tasks into task_queue. A task is eventalign information of one read.
    eventalign_file = open(eventalign_combine_filepath, 'r', encoding='utf-8')
    pos_start = len(eventalign_file.readline()) #remove header
    chunk_split = None
    for chunk in pd.read_csv(eventalign_combine_filepath, chunksize=chunk_size, sep='\t'):
        chunk_complete = chunk[chunk['read_index'] != chunk.iloc[-1]['read_index']]
        chunk_concat = pd.concat([chunk_split, chunk_complete])
        chunk_concat_size = len(chunk_concat.index)
        ## read the file at where it left off because the file is opened once ##
        lines = [len(eventalign_file.readline()) for i in range(chunk_concat_size)]
        chunk_concat.loc[:, 'line_length'] = np.array(lines)
        task_queue.put((chunk_concat, min_seq_len, max_seq_len, shift, out_paths, sample_reads_dict))
        pos_start += sum(lines)
        chunk_split = chunk[chunk['read_index'] == chunk.iloc[-1]['read_index']].copy()

    ## the loop above leaves off w/o adding the last read_index to eventalign.index
    chunk_split_size = len(chunk_split.index)
    lines = [len(eventalign_file.readline()) for i in range(chunk_split_size)]
    chunk_split.loc[:, 'line_length'] = np.array(lines)
    task_queue.put((chunk_split, min_seq_len, max_seq_len, shift, out_paths, sample_reads_dict))

    # Put the stop task into task_queue.
    task_queue = end_queue(task_queue,n_processes)

    # Wait for all of the tasks to finish.
    task_queue.join()


def generate_samples_from_index(index_file, min_reads_per_trans=10):
    """
    :param index_file:
    :param min_reads_per_trans:
    :return: dict
    """
    if not os.path.exists(index_file):
        raise FileExistsError("Index file does not exist!")
        sys.exit()

    df_eventalign_index = pd.read_csv(index_file)[['transcript_id', 'read_index']].groupby(['transcript_id'])
    df_eventalign_index = df_eventalign_index['read_index'].apply(lambda x: list(x)[0: min_reads_per_trans])\
        .reset_index()
    sample_reads_dict = df_eventalign_index.set_index('transcript_id')['read_index'].to_dict()

    sum_reads = 0
    for tid in sample_reads_dict:
        sum_reads += len(sample_reads_dict[tid])
    return sum_reads, sample_reads_dict


def generate_full_reads_info(reads_filepath: str, summary_filepath: str, filter_filepath: str):

    def get_value(x, dict):
        if x in dict:
            return 1
        else:
            return 0

    summary_df = read_summary(summary_filepath)
    readId2name_dict = summary_df.set_index('read_index')['read_name'].to_dict()
    readId2fast5_dict = summary_df.set_index('read_index')['fast5_path'].to_dict()

    reads_df = pd.read_csv(reads_filepath)
    reads_df["read_index"] = reads_df["read_index"].apply(int)
    reads_df['read_name'] = reads_df["read_index"].apply(lambda x: readId2name_dict[x])
    reads_df['fast5_path'] = reads_df["read_index"].apply(lambda x: readId2fast5_dict[x])

    if filter_filepath is not None:
        filter_df = pd.read_csv(filter_filepath)
        filter_dict = filter_df.set_index('read_name')['flag'].to_dict()
        reads_df['select'] = reads_df['read_name'].apply(lambda x: get_value(x, filter_dict))
        reads_df = reads_df[reads_df['select'] == 1]
        reads_df = reads_df.drop(columns=['select'])

    reads_df.to_csv(reads_filepath, header=True, index=False)
    return len(reads_df)


def generate_fast5(reference_filepath: str, reads_filepath: str, fast5_folderpath: str, save_folderpath: str):

    for sub_folder in ['fast5', 'label']:
        if not os.path.exists(os.path.join(save_folderpath, sub_folder)):
            os.mkdir(os.path.join(save_folderpath, sub_folder))

    save_label_file = os.path.join(save_folderpath, "label", "label.csv")
    if os.path.exists(save_label_file):
        os.remove(save_label_file)
    write_row_to_file(save_label_file, label_header)

    ref_dict = read_transcript_fa(reference_filepath)
    print('Read reference done!')

    reads_df = pd.read_csv(reads_filepath)
    reads_df['fast5_path'] = reads_df['fast5_path'].apply(lambda x: x.split("/")[-1])
    reads_df = reads_df.groupby(['fast5_path'])
    for fast5, group in reads_df:
        print(f"* start to process {fast5}")
        source_fast5 = os.path.join(fast5_folderpath, fast5)
        target_fast5 = os.path.join(save_folderpath, 'fast5', fast5)
        copyfile(source_fast5, target_fast5)

        group = group.reset_index().to_dict('records')
        group_dict = dict()
        for item in group:
            group_dict[item['read_name']] = item

        with h5py.File(target_fast5, 'r+') as fast5_data:
            for read_key in fast5_data:
                _read_name = fast5_data[read_key]["Raw"].attrs['read_id'].decode("utf-8")
                if _read_name not in group_dict:
                    del fast5_data[read_key]
                else:
                    _reads_info_dict = group_dict[_read_name]
                    _contig = _reads_info_dict['contig']
                    _read_index = _reads_info_dict['read_index']
                    _contig_st = _reads_info_dict['contig_st']
                    _contig_en = _reads_info_dict['contig_en']
                    _signal_st = _reads_info_dict['signal_st']
                    _signal_en = _reads_info_dict['signal_en']
                    _order = _reads_info_dict['order']
                    _first_kmer = _reads_info_dict['first_kmer']
                    _ref = ref_dict[_contig][0][_contig_st: _contig_en + 1]

                    assert _ref[0:3] == _first_kmer[2:]
                    assert _order == 1
                    _ref = _ref[::-1]

                    value = np.array(fast5_data[read_key]["Raw"]['Signal'][:][_signal_st: _signal_en])
                    del fast5_data[read_key]['Raw/Signal']
                    fast5_data[read_key]['Raw'].attrs['duration'] = _signal_en - _signal_st
                    fast5_data[read_key]['Raw']['Signal'] = value

                    write_row_to_file(save_label_file, [_read_name, _read_index, fast5, _signal_en - _signal_st,
                                                        _contig, _contig_st, _contig_en, len(_ref), _ref])


def arg_parser():
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter,
        add_help=False
    )

    # Required arguments
    parser.add_argument('--input_file',
                        help='input eventalign filepath, the output from nanopolish.',
                        required=True)
    # Optional arguments
    parser.add_argument('--min_reads_per_trans', default=20, type=int)
    parser.add_argument('--min_seq_len', default=100, type=int)
    parser.add_argument('--max_seq_len', default=800, type=int)
    parser.add_argument('--shift', default=2, type=int)
    parser.add_argument('--n_processes', help='number of processes to run.',
                        default=1, type=int)
    parser.add_argument('--chunk_size', help='number of lines from nanopolish eventalign.txt for processing.',
                        default=1000000, type=int)

    parser.add_argument('--reference_file', type=str, default=None)
    parser.add_argument('--summary_file', type=str, default=None)
    parser.add_argument('--filter_file', type=str, default=None)
    parser.add_argument('--fast5_folder', type=str, default=None)
    parser.add_argument('--save_folder', type=str, default=None)
    return parser


def get_reads_num(fast5):
    i = 0
    with h5py.File(fast5, 'r') as fast5_data:
        for read_key in fast5_data:
            i += 1
    return i


if __name__ == '__main__':
    """
    Generate index file for eventalign and combine it.
       Input file:   ../path/to/eventalign/eventalign.txt
       Output file:  ../path/to/eventalign/eventalign.index
                     ../path/to/eventalign/eventalign_combined.txt
    """
    warnings.filterwarnings("ignore")
    eventalign_header = ["contig", "position", "reference_kmer", "read_index", "strand", "event_index",
                         "event_level_mean", "event_stdv", "event_length", "model_kmer", "model_mean", "model_stdv",
                         "standardized_level", "start_idx", "end_idx"]

    reads_header = ["contig", "read_index", "contig_st",  "contig_en", "signal_st", "signal_en", "order", "first_kmer"]
    label_header = ['read_name', 'read_index', 'fast5_path', 'signal_len', 'contig', 'contig_st', 'contig_en',
                    'reference_len', 'reference']

    args = arg_parser().parse_args()
    print(f" ========== rna ========== ")
    args.index_file = args.input_file.replace("txt", "") + "index"
    args.combine_file = args.input_file.replace(".txt", "") + "_combined.txt"
    args.reads_file = args.input_file.replace(".txt", "") + "_reads.txt"

    # step 1: sample reads
    sum_reads, sample_reads_dict = generate_samples_from_index(args.index_file, args.min_reads_per_trans)
    print(f"Sample {sum_reads} reads from {len(sample_reads_dict)} transcripts!")

    # step2: generate sample reads info
    index_start = time.time()
    parallel_generate(args.combine_file, args.reads_file, sample_reads_dict, args.n_processes, args.chunk_size,
                      args.min_seq_len, args.max_seq_len, args.shift)
    index_end = time.time()
    index_run = np.around((index_end - index_start) / 3600, 3)
    print(f"Generate read info done! Running time {index_run} hours.")

    # step3: generate read name and fast5 info
    selected_reads_num = generate_full_reads_info(args.reads_file, args.summary_file, args.filter_file)
    print(f"Generate read full info done! Finally select {selected_reads_num} reads.")

    # step 4: generate base calling input and output
    generate_fast5(args.reference_file, args.reads_file, args.fast5_folder, args.save_folder)
    print("Base calling dataset done!")
