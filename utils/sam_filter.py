import csv
import os
import argparse
import pysam
import pandas as pd


def rm_exist_file(file_name):
    if os.path.exists(file_name):
        os.remove(file_name)


def generate_trans_dict(file_name):
    data = pd.read_csv(file_name)
    data['index'] = data['index'].apply(lambda x: int(x))
    trans2idx_dict = data.set_index('transcript')['index'].to_dict()
    idx2trans_dict = data.set_index('index')['transcript'].to_dict()
    return trans2idx_dict, idx2trans_dict


def write_to_csv(file_name, row):
    with open(file_name, 'a+', newline="") as f:
        writer = csv.writer(f)
        writer.writerow(row)


def filter_sam(sam_file, keep_file_name):

    rm_exist_file(keep_file_name)
    write_to_csv(keep_file_name, ["read_name", "contig", "flag", "reference_start", "reference_end"])

    with open(keep_file_name, 'a+', newline="") as f:
        writer = csv.writer(f)
        with pysam.Samfile(sam_file, 'r') as sf:
            for read in sf:
                if read.query_name and read.reference_name:
                    writer.writerow([read.query_name, read.reference_name, read.flag,
                                     read.reference_start, read.reference_end])

    df = pd.read_csv(keep_file_name)
    df['flag'] = df['flag'].apply(lambda x: int(x))
    print(f"* origin num : {len(df)}")

    mapped_df = df[df['flag'] == 0]
    print(f"* mapped num : {len(mapped_df)}")

    mapped_df.to_csv(keep_file_name, index=False)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='.')
    parser.add_argument('--root-dir', type=str, default='benchmark_datasets')
    parser.add_argument('--dataset', type=str, default='ecoli_eligos')
    parser.add_argument('--sample', type=str, default='IVT_m5C')
    args = parser.parse_args()

    sam_file_name = os.path.join(args.root_dir, args.dataset, "4_nanopolish", args.sample,
                                 "reads-ref.sorted.filter.sam")
    keep_file_name = os.path.join(args.root_dir, args.dataset, "4_nanopolish", args.sample,
                                  "reads-ref.sorted.filter.csv")

    print('*', args.dataset, args.sample)
    # filter_sam(sam_file_name, keep_file_name)

    df = pd.read_csv(keep_file_name)
    print(len(df))
    df1 = df[df['contig'] == "C1"]
    print(len(df1), len(df1)/len(df))

    df2 = df[df['contig'] == "C2"]
    print(len(df2), len(df2) / len(df))
