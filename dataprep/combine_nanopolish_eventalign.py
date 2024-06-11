import os, csv
import sys
import pandas as pd
import numpy as np
import multiprocessing
import warnings
import time
from io import StringIO
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


def index(eventalign_result: pd.DataFrame, pos_start: int, out_paths: Dict, locks: Dict):
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
        None
    '''
    eventalign_result = eventalign_result.set_index(['contig', 'read_index'])
    pos_end = pos_start
    with locks['index'], open(out_paths['index'], 'a', encoding='utf-8') as f_index:
        for _index in list(dict.fromkeys(eventalign_result.index)):
            transcript_id, read_index = _index
            pos_end += eventalign_result.loc[_index]['line_length'].sum()
            f_index.write('%s,%d,%d,%d\n' % (transcript_id, read_index, pos_start, pos_end))
            pos_start = pos_end


def parallel_index(eventalign_filepath: str, index_filepath: str,  n_processes: int, chunk_size: int):
    '''
    Function to index every read within eventalign.txt file for faster access later

    Args:
        eventalign_filepath (str): String filepath to the eventalign.txt file
        chunk_size (int): Chunksize argument for pd.read_csv function
        out_dir (str):  String filepath to the output directory of the indexing function
        n_processes (int): Number of processes used for indexing

    Returns:
        None
    '''

    # if exits, remove the saving file
    rm_exist_file(index_filepath)

    # Create output paths and locks.
    out_paths, locks = dict(), dict()
    for out_filetype in ['index']:
        out_paths[out_filetype] = index_filepath
        locks[out_filetype] = multiprocessing.Lock()

    with open(out_paths['index'], 'w', encoding='utf-8') as f:
        f.write('transcript_id,read_index,pos_start,pos_end\n') # header

    # Create communication queues.
    task_queue = multiprocessing.JoinableQueue(maxsize=n_processes * 2)

    # Create and start consumers.p
    consumers = [Consumer(task_queue=task_queue, task_function=index, locks=locks) for i in range(n_processes)]
    for process in consumers:
        process.start()

    ## Load tasks into task_queue. A task is eventalign information of one read.
    eventalign_file = open(eventalign_filepath, 'r', encoding='utf-8')
    pos_start = len(eventalign_file.readline()) #remove header
    chunk_split = None
    index_features = ['contig', 'read_index', 'line_length']
    for chunk in pd.read_csv(eventalign_filepath, chunksize=chunk_size, sep='\t'):
        chunk_complete = chunk[chunk['read_index'] != chunk.iloc[-1]['read_index']]
        chunk_concat = pd.concat([chunk_split, chunk_complete])
        chunk_concat_size = len(chunk_concat.index)
        ## read the file at where it left off because the file is opened once ##
        lines = [len(eventalign_file.readline()) for i in range(chunk_concat_size)]
        chunk_concat.loc[:, 'line_length'] = np.array(lines)
        task_queue.put((chunk_concat[index_features], pos_start, out_paths))
        pos_start += sum(lines)
        chunk_split = chunk[chunk['read_index'] == chunk.iloc[-1]['read_index']].copy()

    ## the loop above leaves off w/o adding the last read_index to eventalign.index
    chunk_split_size = len(chunk_split.index)
    lines = [len(eventalign_file.readline()) for i in range(chunk_split_size)]
    chunk_split.loc[:, 'line_length'] = np.array(lines)
    task_queue.put((chunk_split[index_features], pos_start, out_paths))

    # Put the stop task into task_queue.
    task_queue = end_queue(task_queue,n_processes)

    # Wait for all of the tasks to finish.
    task_queue.join()


def combine_dataframe(events_str: str):

    if base_type == "dna":
        map_fail_kmer = 'NNNNNN'
    if base_type == "rna":
        map_fail_kmer = 'NNNNN'

    f_string = StringIO(events_str)
    eventalign_result = pd.read_csv(f_string, delimiter='\t', names=header)
    f_string.close()
    eventalign_result = eventalign_result[eventalign_result['model_kmer'] != map_fail_kmer]

    if len(eventalign_result) != 0:
        keys = ['read_index', 'contig', 'position', 'reference_kmer']  # for groupby
        eventalign_result.loc[:, 'length'] = pd.to_numeric(eventalign_result['end_idx']) - \
                                             pd.to_numeric(eventalign_result['start_idx'])
        eventalign_result.loc[:, 'sum_norm_mean'] = pd.to_numeric(eventalign_result['event_level_mean']) \
                                                    * eventalign_result['length']
        eventalign_result.loc[:, 'sum_norm_std'] = pd.to_numeric(eventalign_result['event_stdv']) \
                                                   * eventalign_result['length']
        eventalign_result.loc[:, 'sum_dwell_time'] = pd.to_numeric(eventalign_result['event_length']) \
                                                     * eventalign_result['length']

        eventalign_result = eventalign_result.groupby(keys)
        sum_norm_mean = eventalign_result['sum_norm_mean'].sum()
        sum_norm_std = eventalign_result["sum_norm_std"].sum()
        sum_dwell_time = eventalign_result["sum_dwell_time"].sum()

        strand = eventalign_result['strand'].apply(lambda x: list(x)[0])
        model_kmer = eventalign_result['model_kmer'].apply(lambda x: list(x)[0])
        event_idx = eventalign_result['event_index'].min()
        start_idx = eventalign_result['start_idx'].min()
        end_idx = eventalign_result['end_idx'].max()
        total_length = eventalign_result['length'].sum()
        eventalign_result = pd.concat([start_idx, end_idx], axis=1)

        eventalign_result['strand'] = strand
        eventalign_result['event_index'] = event_idx
        eventalign_result['event_level_mean'] = (sum_norm_mean / total_length).round(2)
        eventalign_result["event_stdv"] = (sum_norm_std / total_length).round(3)
        eventalign_result["event_length"] = (sum_dwell_time / total_length).round(5)
        eventalign_result['model_kmer'] = model_kmer

        eventalign_result = eventalign_result.reset_index()
        eventalign_result = eventalign_result.sort_values(by=['read_index', 'event_index'])

        eventalign_result['model_mean'] = 0.0
        eventalign_result['model_stdv'] = 0.0
        eventalign_result['standardized_level'] = 0.0
        eventalign_result = eventalign_result[header]
        return eventalign_result

    return pd.DataFrame([])


def combine(eventalign_result: pd.DataFrame, out_paths: Dict, locks: Dict):
    '''
    Function to index the position of a specific read features within eventalign.txt

    Args:
        eventalign_result (pd.DataFrame): A pd.DataFrame object containing a portion of eventalign.txt to be indexed
        pos_start (int): An index position within eventalign.txt that corresponds to the start of the eventalign_result
                         portion within eventalign.txt file
        out_paths (Dict): A dictionary containing filepath for all the output files produced by the index function
        locks (Dict): A lock object from multiprocessing library that ensures only one process write to the output file
                      at any given time
    Returns:
        None
    '''
    locks['combined'].acquire()
    eventalign_result.to_csv(out_paths['combined'], sep='\t', index=False, mode='a', header=False)
    locks['combined'].release()


def parallel_combine(input_filepath: str, index_filepath: str, combine_filepath: str, n_processes: int):

    # check index file and remove existing combine file
    if not os.path.exists(index_filepath):
        raise FileExistsError("Index file does not exist!")
        sys.exit()
    rm_exist_file(combine_filepath)

    # Create output paths and locks.
    out_paths, locks = dict(), dict()
    for out_filetype in ['combined']:
        out_paths[out_filetype] = combine_filepath
        locks[out_filetype] = multiprocessing.Lock()

    # write the header
    write_row_to_file(out_paths['combined'], header, '\t')

    # Create communication queues.
    task_queue = multiprocessing.JoinableQueue(maxsize=n_processes * 2)

    # Create and start consumers.p
    consumers = [Consumer(task_queue=task_queue, task_function=combine, locks=locks) for i in range(n_processes)]
    for process in consumers:
        process.start()

    # read index file
    df_eventalign_index = pd.read_csv(index_filepath)
    tx_ids = df_eventalign_index['transcript_id'].values.tolist()
    tx_ids = list(dict.fromkeys(tx_ids))
    df_eventalign_index = df_eventalign_index.set_index('transcript_id')

    # start to combine
    with open(input_filepath, 'r', encoding='utf-8') as eventalign_result:
        for tx_id in tx_ids:
            for _, row in df_eventalign_index.loc[[tx_id]].iterrows():
                read_index, pos_start, pos_end = row['read_index'], row['pos_start'], row['pos_end']
                eventalign_result.seek(pos_start, 0)
                events_str = eventalign_result.read(pos_end-pos_start)
                combined_df = combine_dataframe(events_str)
                task_queue.put((combined_df, out_paths))

    # Put the stop task into task_queue.
    task_queue = end_queue(task_queue, n_processes)

    # Wait for all of the tasks to finish.
    task_queue.join()


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
    parser.add_argument('--base_type', default='rna', type=str)
    parser.add_argument('--output_suffix',
                        help='output combined eventalign filename suffix.',
                        default='_combined', type=str)
    parser.add_argument('--n_processes',
                        help='number of processes to run.',
                        default=1, type=int)
    parser.add_argument('--chunk_size',
                        help='number of lines from nanopolish eventalign.txt for processing.',
                        default=1000000, type=int)
    return parser


if __name__ == '__main__':
    """
    Generate index file for eventalign and combine it.
       Input file:   ../path/to/eventalign/eventalign.txt
       Output file:  ../path/to/eventalign/eventalign.index
                     ../path/to/eventalign/eventalign_combined.txt
    """
    warnings.filterwarnings("ignore")
    header = ["contig", "position", "reference_kmer", "read_index", "strand", "event_index", "event_level_mean",
              "event_stdv", "event_length", "model_kmer", "model_mean", "model_stdv", "standardized_level",
              "start_idx", "end_idx"]

    args = arg_parser().parse_args()
    base_type = args.base_type
    print(f" ========== {base_type} ========== ")
    args.index_file = args.input_file.replace("txt", "") + "index"
    args.combine_file = args.input_file.replace(".txt", "") + "_combined.txt"

    print(f"Input file : {args.input_file}")

    print(f"Start to generate index file : {args.index_file}")
    index_start = time.time()
    parallel_index(args.input_file, args.index_file, args.n_processes, args.chunk_size)
    index_end = time.time()
    index_run = np.around((index_end - index_start) / 3600, 3)
    print(f"Index file done! Running time {index_run} hours.")

    print(f"Start to generate combined file : {args.combine_file}")
    combine_start = time.time()
    parallel_combine(args.input_file, args.index_file, args.combine_file, args.n_processes)
    combine_end = time.time()
    combine_run = np.around((combine_end - combine_start) / 3600, 3)
    print(f"Combing done! Running time {combine_run} hours.")
