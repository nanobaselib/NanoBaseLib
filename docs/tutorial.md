## NanoBaseLib Tutorials

### 1. dataprep

#### combine_nanopolish_eventalign

`python dataprep/combine_nanopolish_eventalign.py --input_file eventalign.txt`

Input:

| Parameter         | Description |
| -----------       | ----------- |
| --input_file      | Input eventalign filepath, the output from nanopolish eventalign.  |
| --output_suffix   | Output combined eventalign filename suffix        |
| --n_processes     |  Number of processes to run. |
| --chunk_size      |  Number of lines from nanopolish eventalign.txt for processing. |


Output:

| File                     | Description |
| -----------              | ----------- |
| eventalign_combined.txt  | The combined eventalign file. |
| eventalign.index         | The index file for original eventalign file.    |


#### extract_tombo_resquiggle

`python dataprep/extract_tombo_resquiggle.py`

Input:

| Parameter         | Description |
| -----------       | ----------- |
| --root_dir        | The root dir of tombo resquiggle resluts. |
| --folder          | The folder of  tombo resquiggle resluts. Default: "0".      |


Output:

| File                     | Description |
| -----------              | ----------- |
| tombo_resquiggle.txt     | Tombo resquiggle resluts. |
| tombo_summary.txt        | The arrtibutes of tombo resquiggle resluts.   |

### 2. base_calling

#### generate_dna_dataset / generate_rna_dataset

Input:

| Parameter         | Description |
| -----------       | ----------- |
| --input_file      |  Input combined eventalign filepath.  |
| --reference_file  |  Genome reference_file. |
| --fast5_folder    |  Orignal fast5 folder. |
| --save_folder     |  Save fast5 folder. |

Output:

| File                     | Description |
| -----------              | ----------- |
| Generated fast5 files.   | Generated fast5 files saved in `--save_folder`, which correspond with the ground truth sequence.  |
| label.csv                | The the ground truth sequence file.    |

#### base_calling_dataloader

Input:

| Parameter         | Description |
| -----------       | ----------- |
| --fast5_folder    | Input fast5 folder.  |
| --label_file      | Ground turth (label) file for each reads. Default: "label.csv". |

Output:

| File                     | Description |
| -----------              | ----------- |
| seg_data_arr.npy         |  Current signals.           |
| seg_len_arr.npy          | The length of current signals. |
| base_len_arr.npy         | Ground turth sequence. |
| base_seq_arr.npy         | The length of ground turth sequence.    |

#### base_calling_evaluation

Input:

| Parameter         | Description |
| -----------       | ----------- |
| --fast5_folder    | Input fast5 folder.  |
| --label_file      | Ground turth (label) file for each reads. Default: "label.csv". |

Output: The average accuracy of matches ($M$), mismatches ($X$), insertions ($I$), and deletions ($D$).

### 3. polya_detection

### 4. segment_align


