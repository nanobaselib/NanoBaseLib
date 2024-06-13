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

#### generate_dna_dataset

#### generate_rna_dataset

#### base_calling_dataloader

#### base_calling_evaluation

### 3. polya_detection

### 4. segment_align


