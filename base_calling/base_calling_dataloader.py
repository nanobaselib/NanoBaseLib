import torch
from torch.utils.data import Dataset


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


def dataset_prep(root_dir, fast5_reads_folder, label_file):
    """
    :param fast5_reads_folder:
    :param label_file:
    :return:  base_len_arr.npy, base_seq_arr.npy, seg_data_arr.npy, seg_len_arr.npy
    """
    pass


def test_unit(exampel_data_dir):

    # load dataset
    src_seqs = np.load(exampel_data_dir + 'seg_data_arr.npy', mmap_mode='r')
    src_lens = np.load(exampel_data_dir + 'seg_len_arr.npy', mmap_mode='r').astype(np.int32)
    tgt_seqs = np.load(exampel_data_dir + 'base_seq_arr.npy', mmap_mode='r')
    tgt_lens = np.load(exampel_data_dir + 'base_len_arr.npy', mmap_mode='r').astype(np.int32)

    test_num = 3000

    train_dataset = Seg2BaseDataset(src_seqs[0: -test_num], src_lens[0: -test_num],
                                    tgt_seqs[0: -test_num], tgt_lens[0: -test_num])
    train_dataloader = DataLoader(train_dataset, batch_size=32, shuffle=True)

    test_dataset = Seg2BaseDataset(src_seqs[-test_num:], src_lens[-test_num:],
                                   tgt_seqs[-test_num:], tgt_lens[-test_num:])
    test_dataloader = DataLoader(test_dataset, batch_size=20, shuffle=False)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    for batch in train_dataloader:
        src_seq = batch['src_seq'].to(device)
        src_len = batch['src_len'].to(device)
        tgt_seq = batch['tgt_seq'].to(device)
        tgt_len = batch['tgt_len'].to(device)
        print(src_seq[0])
        print(src_len)

        break

if __name__ == '__main__':

    dataset_prep(exampel_data_dir, fast5_reads_folder, label_file)
    test_unit()