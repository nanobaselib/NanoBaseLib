import os, csv, re
import parasail
import pandas as pd
import numpy as np
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter
from collections import defaultdict


split_cigar = re.compile(r"(?P<len>\d+)(?P<op>\D+)")
base_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}


def complement(seq):
    _seq = ""
    for base in seq:
        _seq += base
    return _seq


def parasail_to_sam(result, seq):
    """
    Extract reference start and sam compatible cigar string.

    :param result: parasail alignment result.
    :param seq: query sequence.

    :returns: reference start coordinate, cigar string.
    """
    cigstr = result.cigar.decode.decode()
    first = re.search(split_cigar, cigstr)

    first_count, first_op = first.groups()
    prefix = first.group()
    rstart = result.cigar.beg_ref
    cliplen = result.cigar.beg_query

    clip = '' if cliplen == 0 else '{}S'.format(cliplen)
    if first_op == 'I':
        pre = '{}S'.format(int(first_count) + cliplen)
    elif first_op == 'D':
        pre = clip
        rstart = int(first_count)
    else:
        pre = '{}{}'.format(clip, prefix)

    mid = cigstr[len(prefix):]
    end_clip = len(seq) - result.end_query - 1
    suf = '{}S'.format(end_clip) if end_clip > 0 else ''
    new_cigstr = ''.join((pre, mid, suf))
    return rstart, new_cigstr


def accuracy(ref, seq, balanced=False, min_coverage=0.0):
    """
    Calculate the accuracy between `ref` and `seq`
    matches (M), mismatches (X), insertions (I), and deletions (D)
    """
    alignment = parasail.sw_trace_striped_32(seq, ref, 8, 4, parasail.dnafull)
    counts = defaultdict(int)

    q_coverage = len(alignment.traceback.query) / len(seq)
    r_coverage = len(alignment.traceback.ref) / len(ref)


    if r_coverage < min_coverage:
        return []

    _, cigar = parasail_to_sam(alignment, seq)

    for count, op in re.findall(split_cigar, cigar):
        counts[op] += int(count)

    if balanced:
        accuracy = (counts['='] - counts['I']) / (counts['='] + counts['X'] + counts['D'])
    else:
        accuracy = counts['='] / (counts['='] + counts['I'] + counts['X'] + counts['D'])

    alg_len = counts['='] + counts['I'] + counts['X'] + counts['D']
    ref_len = len(ref)

    alg_m = np.around(counts['='] / alg_len * 100, 5)
    alg_i = np.around(counts['I'] / alg_len * 100, 5)
    alg_x = np.around(counts['X'] / alg_len * 100, 5)
    alg_d = np.around(counts['D'] / alg_len * 100, 5)

    ref_m = np.around(counts['='] / ref_len * 100, 5)
    ref_i = np.around(counts['I'] / ref_len * 100, 5)
    ref_x = np.around(counts['X'] / ref_len * 100, 5)
    ref_d = np.around(counts['D'] / ref_len * 100, 5)

    return [alg_m, alg_i, alg_x, alg_d, ref_m, ref_i, ref_x, ref_d]


def print_alignment(ref, seq):
    """
    Print the alignment between `ref` and `seq`
    """
    alignment = parasail.sw_trace_striped_32(seq, ref, 8, 4, parasail.dnafull)

    print('ref len : ', len(alignment.traceback.ref), '   query len : ', len(alignment.traceback.query))
    print(alignment.traceback.ref)
    print(alignment.traceback.comp)
    print(alignment.traceback.query)

    print("  Score=%s" % alignment.score)
    return alignment.score


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


def get_acc_res(save_file, model_dict, gt_df):
    for index, row in gt_df.iterrows():
        read_name = row['read_name']
        reference = row['reference']
        if base_type == 'rna':
            strand = "-"
        else:
            strand = row['strand']
        if base_type == 'rna':
            reference = reference[::-1]
        if read_name in model_dict:
            predict = model_dict[read_name]
            try:
                acc_arr = accuracy(reference, predict)
                write_row_to_file(save_file, [read_name, strand] + acc_arr)
            except:
                print(f"{read_name} error !")


def process_dorado(root_dir, gt_df, model, fomart='sam'):

    if "dorado" in model:
        input_file = os.path.join(root_dir, "output", model, "dorado.sam")
    else:
        input_file = os.path.join(root_dir, "output", model, model + "." + fomart)

    save_file = os.path.join(root_dir, "output", model, "accuracy.csv")
    rm_exist_file(save_file)
    write_row_to_file(save_file, accuracy_file_header)

    dorado_dict = dict()
    with open(input_file, "r") as f:
        for line in f:
            if line[0] == "@": continue
            line = line.split()
            read_name = line[0]
            if read_name not in dorado_dict:
                dorado_dict[read_name] = line[9]
    print(f"{model} data processing done !")

    get_acc_res(save_file, dorado_dict, gt_df)
    print(f"{model} accuracy generation done !")


def process_guppy(root_dir, gt_df, model):

    if "guppy" in model:
        input_file = os.path.join(root_dir, "output", model, "guppy.fastq")
    else:
        input_file = os.path.join(root_dir, "output", model, model + ".fastq")

    save_file = os.path.join(root_dir, "output", model, "accuracy.csv")
    rm_exist_file(save_file)
    write_row_to_file(save_file, accuracy_file_header)

    guppy_dict = dict()
    read_full = None
    with open(input_file, "r") as f:
        for idx, line in enumerate(f):
            if line[0] == "@":
                if read_full:
                    read_full = read_full.split('\n')
                    read_name = read_full[0].split(" ")[0].replace("@", "")
                    pred_seq = read_full[1].replace("U", "T")
                    guppy_dict[read_name] = pred_seq
                read_full = line
            else:
                read_full += line
    print(f"{model} data processing done !")

    j = 0
    for i in guppy_dict:
        print(i)
        j += 1
        if j == 10: break

    get_acc_res(save_file, guppy_dict, gt_df)
    print(f"{model} accuracy generation done !")


def process_rodan(root_dir, gt_df):

    input_file = os.path.join(root_dir, "output", "rodan", "rodan.fasta")
    save_file = os.path.join(root_dir, "output", "rodan", "accuracy.csv")
    rm_exist_file(save_file)
    write_row_to_file(save_file, accuracy_file_header)

    rodan_dict = dict()
    fasta = open(input_file, "r")
    entries = ""
    for ln in fasta:
        entries += ln
    entries = entries.split(">")
    for en in entries:
        if en:
            en = en.split("\n")
            rodan_dict[en[0]] = en[1]
    print("rodan data processing done !")

    get_acc_res(save_file, rodan_dict, gt_df)
    print("rodan accuracy generation done !")


def process_halcyon(root_dir, gt_df):

    input_file = os.path.join(root_dir, "output", "halcyon", "halcyon.fasta")
    save_file = os.path.join(root_dir, "output", "halcyon", "accuracy.csv")
    rm_exist_file(save_file)
    write_row_to_file(save_file, accuracy_file_header)

    rodan_dict = dict()
    fasta = open(input_file, "r")
    entries = ""
    for ln in fasta:
        entries += ln
    entries = entries.split(">")
    for en in entries:
        if en:
            en = en.split("\n")
            read = en[0].split("/")[-1].split(".")[0]
            rodan_dict[read] = en[1]
    print("halcyon data processing done !")

    get_acc_res(save_file, rodan_dict, gt_df)
    print("halcyon accuracy generation done !")


def process_radian(root_dir, gt_df):

    input_file = os.path.join(root_dir, "output", "radian", "radian.fasta")
    save_file = os.path.join(root_dir, "output", "radian", "accuracy.csv")
    rm_exist_file(save_file)
    write_row_to_file(save_file, accuracy_file_header)

    rodan_dict = dict()
    fasta = open(input_file, "r")
    entries = ""
    for ln in fasta:
        entries += ln
    entries = entries.split(">")
    for en in entries:
        if en:
            en = en.split("\n")
            rodan_dict[en[0]] = en[1]
    print("rodan data processing done !")

    get_acc_res(save_file, rodan_dict, gt_df)
    print("rodan accuracy generation done !")


def read_res_acc(file_name, strand=None):
    df = pd.read_csv(file_name)
    if strand:
        df = df[df['strand'] == strand]
    res_arr = []
    acc_cols = accuracy_file_header[2:]
    for col in acc_cols:
        res_arr.append(np.around(df[col].mean(), 4))
    return res_arr


def arg_parser():
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter, add_help=False)
    parser.add_argument('--root_dir',
                        help='input eventalign filepath, the output from nanopolish.',
                        default='/scratch/cs/infantbiome/chengg1/datasets/hek293t_bc_benchmark')
    parser.add_argument('--base_type', default='dna')

    parser.add_argument('--model', default=None)
    parser.add_argument('--plot', default=False)
    parser.add_argument('--strand', default=None)
    return parser


if __name__ == '__main__':

    args = arg_parser().parse_args()
    baselines = ['rodan', 'radian', 'guppy2.3.1',  'guppy4.5.4', 'guppy6.0.1', 'bonito', 'dorado0.5.3', 'dorado0.7.0',
                 'causalcall', 'halcyon']
    accuracy_file_header = ['read_name', 'strand', 'mat_alg_rate', 'ins_alg_rate', 'mis_alg_rate', 'del_alg_rate',
                            'mat_ref_rate', 'ins_ref_rate', 'mis_ref_rate', 'del_ref_rate']

    base_type = args.base_type
    print(f" ========== {base_type} ========== ")
    # data pre
    if args.model:
        assert args.model in baselines
        gt_df = pd.read_csv(os.path.join(args.root_dir, "label", "label.csv"))
        print(f"start to process {args.model} data ... ")
        if 'dorado' in args.model or args.model == 'bonito':
            process_dorado(args.root_dir, gt_df, args.model)
        elif args.model == 'radian':
            process_radian(args.root_dir, gt_df)
        elif 'guppy' in args.model:
            process_guppy(args.root_dir, gt_df, args.model)
        elif args.model == 'rodan':
            process_rodan(args.root_dir, gt_df)
        elif args.model == 'halcyon':
            process_halcyon(args.root_dir, gt_df)
        elif args.model == "causalcall":
            process_guppy(args.root_dir, gt_df, args.model)

    if args.plot:
        for model in baselines:
            model_res = os.path.join(args.root_dir, "output", model, "accuracy.csv")
            if os.path.exists(model_res):
                res = read_res_acc(model_res, args.strand)
                res = [str(i) for i in res]
                print(model.capitalize(), '&', ' & '.join(res), '\\\\')












