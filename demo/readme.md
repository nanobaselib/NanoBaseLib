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
                    | -- dorado.bam
            | -- guppy
                    | -- fail
                    | -- pass
                    | -- workspace
            | -- demo_dataset.fastq
            | -- reads-ref.sorted.filter.bam
    | -- 3_tailfindr
            | -- tails.csv
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

```
dorado basecaller rna002_70bps_hac@v3 1_raw_signal/multi_pod5 --estimate-poly-a > 2_base_called/dorado/dorado.bam 
```

### Step 3: Mapping

```
cd 2_base_called
find guppy -type f -name "*.fastq" -exec cat {} + > demo_dataset.fastq
nanopolish index --directory=../1_raw_signal/multi_fast5 demo_dataset.fastq
minimap2 -ax map-ont --MD -t 8 --secondary=no ../0_reference/ref.fa demo_dataset.fastq | samtools sort -o reads-ref.sorted.bam -T reads.tmp
```
![image](https://github.com/nanobaselib/NanoBaseLib/assets/166529164/c92a60fb-8fa7-47de-8456-fe119c020d42)

```
samtools view -b -F 2324 reads-ref.sorted.bam > reads-ref.sorted.filter.bam
samtools index reads-ref.sorted.filter.bam
samtools quickcheck reads-ref.sorted.filter.bam
samtools view -h -o reads-ref.sorted.filter.sam reads-ref.sorted.filter.bam
```

### Step 4: PolyA Detection

```
library(tailfindr)
df <- find_tails(fast5_dir = '2_base_called/guppy/workspace',
                 save_dir = '3_tailfindr',
                 csv_filename = 'tails.csv',
                 num_cores = 20)
```


## NanoBaseLib Software Package
