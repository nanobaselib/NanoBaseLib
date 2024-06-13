## NanoBaseLib Tutorials

### 1. dataprep

#### combine_nanopolish_eventalign

`python dataprep/combine_nanopolish_eventalign.py --input_file eventalign.txt`

Input:

| parameter         | Description |
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
