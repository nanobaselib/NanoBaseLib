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

```
cd demo_dataset
# download the reference (ref.fa) to 0_reference
# download the raw signal to 1_raw_signal
```

### Step 1: Data Standardization

If the format of raw data is single-fast5, convert it into multi-fast5.
```
single_to_multi_fast5 --input_path 1_raw_signal/single_fast5 --save_path 1_raw_signal/multi_fast5 --filename_base demo --batch_size 4000 --recursive; 
```
If Dorado base caller is needed, convert multi-fast5 to pod5.
```
pod5 convert fast5 1_raw_signal/multi_fast5/*.fast5 --output 1_raw_signal/multi_pod5/ --one-to-one 1_raw_signal/multi_fast5
```
### Step 2: Base Calling

```
guppy_basecaller -c rna_r9.4.1_70bps_hac.cfg --num_callers 20 --cpu_threads_per_caller 20 -i 1_raw_signal/multi_fast5 -s 2_base_called/guppy  --fast5_out
```

## NanoBaseLib Software Package
