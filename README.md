
<div align="center">
    <img src=docs/logo_small.svg width=60% />
</div>
<p align="center" ">
 <a href="https://nanobaselib.github.io">NanoBaseLib Website</a> |  <a href="docs/tutorial.md">NanoBaseLib Tutorials</a> |  <a href="https://drive.usercontent.google.com/download?id=1CQ_KkR4_HeTGjhr-zcHznhpUUcUEx1PL&export=download">Download Demo Dataset</a>
</p>

NanoBaseLib contains five modules: `dataprep`, `base_calling`, `ploya_detection`, `segment_align`, and `rna_mod_detection`. The `dataprep` module processes the output of different tools into the same format. The other four modules are designed to benchmark different tasks.

## Demo Dataset
```
NanoBaseLib
    | -- demo_dataset
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
                | -- eventalign.txt
                | -- eventalign_combined.txt
                | -- polya.tsv
                | -- polya-pass-only-with-head.tsv
        | -- 5_tombo
                | -- single_reads
                        | -- tombo_resquiggle.txt
                        | -- tombo_summary.txt
        | -- 6_segpore        
```
## Preprocessing Pipeline

```
git clone https://github.com/nanobaselib/NanoBaseLib.git
cd NanoBaseLib
# download demo dataset
```

### Step 1: Data Standardization

If the format of raw data is single-fast5, convert it into multi-fast5.
```
single_to_multi_fast5 --input_path 1_raw_signal/single_fast5 --save_path 1_raw_signal/multi_fast5 \
                      --filename_base demo --batch_size 4000 --recursive; 
```
If Dorado base caller is needed, convert multi-fast5 to pod5.
```
pod5 convert fast5 1_raw_signal/multi_fast5/*.fast5 --output 1_raw_signal/multi_pod5 \
                   --one-to-one 1_raw_signal/multi_fast5
```
### Step 2: Base Calling

```
guppy_basecaller -c rna_r9.4.1_70bps_hac.cfg --num_callers 20 --cpu_threads_per_caller 20 \
                 -i 1_raw_signal/multi_fast5 -s 2_base_called/guppy  --fast5_out
```

```
dorado basecaller rna002_70bps_hac@v3 1_raw_signal/multi_pod5 \
                  --estimate-poly-a > 2_base_called/dorado/dorado.bam 
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
nanopolish polya --threads=32 --reads=demo_dataset.fastq --bam=reads-ref.sorted.filter.bam \
                 --genome=../0_reference/ref.fa > ../4_nanopolish/polya.tsv
grep -E 'PASS|readname' ../4_nanopolish/polya.tsv > ../4_nanopolish/polya-pass-only-with-head.tsv
```

```
library(tailfindr)
df <- find_tails(fast5_dir = '2_base_called/guppy/workspace',
                 save_dir = '3_tailfindr',
                 csv_filename = 'tails.csv',
                 num_cores = 20)
```

### Step 5: Segmentation and Event Alignmnet
```
nanopolish eventalign --reads demo_dataset.fastq \
--bam reads-ref.sorted.filter.bam \
--genome ../0_reference/ref.fa \
--signal-index \
--scale-events \
--summary ../4_nanopolish/summary.txt \
--threads 32 > ../4_nanopolish/eventalign.txt
```

```
cd ..
mkdir 5_tombo/single_reads
for i in `ls 2_base_called/guppy/workspace`; do multi_to_single_fast5 --input_path 2_base_called/guppy/workspace/${i}  --save_path 5_tombo/single_reads --recursive -t 20; done
tombo resquiggle 5_tombo/single_reads 0_reference/ref.fa --overwrite --processes 20
```
![image](https://github.com/nanobaselib/NanoBaseLib/assets/166529164/de974327-21fd-4fb4-b85e-2e7d819381a8)


## NanoBaseLib Commands

### combine_nanopolish_eventalign
```
python dataprep/combine_nanopolish_eventalign.py \
--input_file demo_dataset/4_nanopolish/eventalign.txt  \
--base_type rna --n_processes 10
```
### extract_tombo_resquiggle
```
python dataprep/extract_tombo_resquiggle.py --root_dir demo_dataset/5_tombo/single_reads --folder 0
```

### generate_rna_dataset
```
mkdir demo_dataset_base_calling
python base_calling/generate_rna_dataset.py \
--input_file demo_dataset/4_nanopolish/eventalign.txt \
--summary_file demo_dataset/4_nanopolish/summary.txt \
--reference_file demo_dataset/0_reference/ref.fa \
--fast5_folder demo_dataset/1_raw_signal/multi_fast5 \
--save_folder demo_dataset_base_calling
```
![image](https://github.com/nanobaselib/NanoBaseLib/assets/166529164/134e290c-8318-47d3-b995-6805ada5abd8)

### base_calling_dataloader
```
python base_calling/base_calling_dataloader.py --root_dir demo_dataset_base_calling
```
