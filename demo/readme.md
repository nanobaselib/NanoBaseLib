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

### step 1: Data standardization

```
single_to_multi_fast5 --input_path demo_dataset/1_raw_signal/single_fast5 --save_path demo_dataset/1_raw_signal/multi_fast5 --filename_base demo --batch_size 4000 --recursive; 
```


## NanoBaseLib Software Package
