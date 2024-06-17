# Demo

## Example Dataset
```
demo_dataset
    | -- 0_reference
            | -- ref.fa
    | -- 1_raw_signal
            | -- single_fast5
            | -- multi_fast5
            | -- multi_pod5
    | -- 2_base_called
            | -- dorado
            | -- guppy
    | -- 3_tailfindr
            | -- tailfindr_pass.csv
    | -- 4_nanopolish
            | -- eventalign_combined.txt
            | -- polya-pass-only-with-head.tsv
    | -- 5_tombo
            | -- tombo_resquiggle.txt
            | -- tombo_summary.txt
    | -- 6_segpore        
```
## Preprocessing Pipeline

### step 1: Data Standardization

If the format of raw data is single-fast5, convert it into multi-fast5.
```
single_to_multi_fast5 --input_path demo_dataset/1_raw_signal/single_fast5 --save_path demo_dataset/1_raw_signal/multi_fast5 --filename_base demo --batch_size 4000 --recursive; 
```
If Dorado base caller is needed, convert multi-fast5 to pod5.
```
pod5 convert fast5 demo_dataset/1_raw_signal/multi_fast5/*.fast5 --output demo_dataset/1_raw_signal/multi_pod5/ --one-to-one demo_dataset/1_raw_signal/multi_fast5
```
### step 2: Base Calling

## NanoBaseLib Software Package
