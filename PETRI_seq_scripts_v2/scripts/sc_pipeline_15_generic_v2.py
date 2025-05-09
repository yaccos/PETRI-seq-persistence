## Supplement to "Prokaryotic single-cell RNA sequencing by in situ combinatorial indexing" (doi: 10.1038/s41564-020-0729-6)
## Written by Sydney Blattman
## Tavazoie Lab, Columbia University
## Last updated January 2025

import time
import os
import sys
import rpy2.robjects as robjects

start = time.time()
script_dir = sys.argv[0].split(sys.argv[0].split('/')[-1])[0]
i = 1
sample = sys.argv[i][0:sys.argv[i].find('_S')]
if sys.argv[i].find('_S') == -1:
    raise ValueError('Must include _S after sample')
n_lanes = int(sys.argv[2])

print(f'Preprocessing {sample}')

cutadapt_command = (
    f'seq {n_lanes} | time parallel --bar -j4 cutadapt -m 75 -o results/{sample}/{sample}_QF_L00{{}}_paired_min75_001.fastq '
    f'results/{sample}/{sample}_QF_L00{{}}_p.assembled.fastq'
)
os.system(cutadapt_command)

# Split paired reads back into two files of read 1 (58 bases) and read 2 (remaining sequence - reverse comp)
split_r1_command = (
    f'seq {n_lanes} | time parallel --bar -j4 cutadapt -l 58 -o results/{sample}/{sample}_QF_L00{{}}_R1_paired.fastq '
    f'results/{sample}/{sample}_QF_L00{{}}_paired_min75_001.fastq'
)
os.system(split_r1_command)

split_r2_command = (
    f'seq {n_lanes} | time parallel --bar -j4 cutadapt -u 58 -o results/{sample}/{sample}_QF_L00{{}}_preR2_paired.fastq '
    f'results/{sample}/{sample}_QF_L00{{}}_paired_min75_001.fastq'
)
os.system(split_r2_command)
reverse_comp_command = (
    f'seq {n_lanes} | time parallel --bar -j4 seqkit seq -r -p -t DNA '
    f'results/{sample}/{sample}_QF_L00{{}}_preR2_paired.fastq -o results/{sample}/{sample}_QF_L00{{}}_R2_paired.fastq'
)
os.system(reverse_comp_command)

# Merge back the reads that overlapped by pear and those that didn't to get clean non-overlapping reads 1 and 2

for l in range(1, n_lanes + 1):
    merge_r1_command = (
        f'cat results/{sample}/{sample}_QF_L00{l}_p.unassembled.forward.fastq '
        f'results/{sample}/{sample}_QF_L00{l}_R1_paired.fastq > results/{sample}/{sample}_QF_merged_L00{l}_R1.fastq'
    )
    os.system(merge_r1_command)
    merge_r2_command = (
        f'cat results/{sample}/{sample}_QF_L00{l}_p.unassembled.reverse.fastq '
        f'results/{sample}/{sample}_QF_L00{l}_R2_paired.fastq > results/{sample}/{sample}_QF_merged_L00{l}_R2.fastq'
    )
    os.system(merge_r2_command)


robjects.globalenv['n_lanes'] = n_lanes
robjects.globalenv['script_dir'] = script_dir
robjects.globalenv['sample'] = sample

robjects.r.source(f'{script_dir}/demultiplexer.R')

end = time.time()
print(f"Time elapsed during Python pipeline: {end - start}")
