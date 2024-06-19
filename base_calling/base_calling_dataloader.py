import os
import h5py
import pandas as pd
import numpy as np
import torch
from torch.utils.data import Dataset
from torch.utils.data import DataLoader
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter

base2id_dict = {'N': 0, 'A': 1, 'T': 2, 'C': 3, 'G': 4}


class Seg2BaseDataset(Dataset):

    def __init__(self, src_seqs, src_lens, tgt_seqs, tgt_lens):
        super().__init__()
        self.src_seqs = src_seqs
        self.src_lens = src_lens
        self.tgt_seqs = tgt_seqs
        self.tgt_lens = tgt_lens

        assert len(src_seqs) == len(src_lens) == len(tgt_seqs) == len(tgt_lens)

    def __len__(self):
        return len(self.src_seqs)

    def __getitem__(self, index):
        item = {}
        item['src_seq'] = torch.tensor(self.src_seqs[index]).float()
        item['src_len'] = torch.tensor(self.src_lens[index])
        item['tgt_seq'] = torch.tensor(self.tgt_seqs[index]).float()
        item['tgt_len'] = torch.tensor(self.tgt_lens[index])
        return item


def pad_lengths(array, max_len):
    padded = [0.00 for _ in range(max_len)]
    padded[0: len(array)] = array
    return padded, len(array)


def dataset_prep(fast5_folder, label_file, output_directory):
    """
    :param fast5_folder:
    :param label_file:
    :return:  base_len_arr.npy, base_seq_arr.npy, seg_data_arr.npy, seg_len_arr.npy
              saved in output_directory
    """
    label_df = pd.read_csv(label_file)
    label_dict = label_df.set_index(['read_name'])['reference'].to_dict()

    signals = []
    signals_len = []
    labels = []
    labels_len = []

    fast5_files = os.listdir(fast5_folder)
    for fast5_read in fast5_files:
        print(f"start to process {fast5_read} ...")
        fast5_read_name = os.path.join(fast5_folder, fast5_read)

        with h5py.File(fast5_read_name, 'r') as fast5_data:
            for read_key in fast5_data:
                _read_name = fast5_data[read_key]["Raw"].attrs['read_id'].decode("utf-8")
                if _read_name not in label_dict:
                    continue

                _signal = fast5_data[read_key]["Raw"]['Signal'][:]
                _label = label_dict[_read_name]

                _signal_len = len(_signal)
                _label_len = len(_label)

                if _signal_len > max_signal_len or _label_len > max_label_len:
                    continue

                _label = [base2id_dict[i] for i in _label]
                _label, _ = pad_lengths(_label, max_label_len)
                _signal, _ = pad_lengths(_signal, max_signal_len)

                signals.append(_signal)
                signals_len.append(_signal_len)
                labels.append(_label)
                labels_len.append(_label_len)

    signals = np.array(signals, dtype=np.float32)
    signals_len = np.array(signals_len, dtype=np.uint16)
    labels = np.array(labels, dtype=np.uint16)
    labels_len = np.array(labels_len, dtype=np.uint16)

    np.save(os.path.join(output_directory, "signals.npy"), signals)
    np.save(os.path.join(output_directory, "signals_length.npy"), signals_len)
    np.save(os.path.join(output_directory, "reference.npy"), labels)
    np.save(os.path.join(output_directory, "reference_length.npy"), labels_len)

def test_unit(data_dir):

    # load dataset
    src_seqs = np.load(data_dir + '/signals.npy', mmap_mode='r')
    src_lens = np.load(data_dir + '/signals_length.npy', mmap_mode='r').astype(np.int32)
    tgt_seqs = np.load(data_dir + '/reference.npy', mmap_mode='r').astype(np.int32)
    tgt_lens = np.load(data_dir + '/reference_length.npy', mmap_mode='r').astype(np.int32)

    test_num = int(len(src_lens) * 0.3)

    train_dataset = Seg2BaseDataset(src_seqs[0: -test_num], src_lens[0: -test_num],
                                    tgt_seqs[0: -test_num], tgt_lens[0: -test_num])
    train_dataloader = DataLoader(train_dataset, batch_size=16, shuffle=True)

    test_dataset = Seg2BaseDataset(src_seqs[-test_num:], src_lens[-test_num:],
                                   tgt_seqs[-test_num:], tgt_lens[-test_num:])
    test_dataloader = DataLoader(test_dataset, batch_size=16, shuffle=False)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    for batch in train_dataloader:
        src_seq = batch['src_seq'].to(device)
        src_len = batch['src_len'].to(device)
        tgt_seq = batch['tgt_seq'].to(device)
        tgt_len = batch['tgt_len'].to(device)
        print(src_seq[0])
        print(src_len)

        break


def arg_parser():
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter,
        add_help=False
    )

    # Required arguments
    parser.add_argument('--root_dir', required=True, type=str, default="demo_dataset_base_calling")

    # Optional arguments
    parser.add_argument('--fast5_folder', type=str, default="fast5")
    parser.add_argument('--input_folder', type=str, default="input")
    parser.add_argument('--label_file', type=str, default="label/label.csv")
    parser.add_argument('--max_signal_len', type=int, default=50000)
    parser.add_argument('--max_label_len', type=int, default=2000)
    return parser


if __name__ == '__main__':
    """
    Generate input dataset for deep learning model.
    """

    args = arg_parser().parse_args()
    fast5_folder = os.path.join(args.root_dir, args.fast5_folder)
    input_folder = os.path.join(args.root_dir, args.input_folder)
    label_file = os.path.join(args.root_dir, args.label_file)

    max_signal_len = args.max_signal_len
    max_label_len = args.max_label_len

    if not os.path.exists(input_folder):
        os.mkdir(input_folder)

    print("start to generate dataset for base calling.")
    dataset_prep(fast5_folder, label_file, input_folder)
    print("\ndataloader example:")
    test_unit(input_folder)