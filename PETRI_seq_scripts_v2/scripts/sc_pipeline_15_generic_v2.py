## Supplement to "Prokaryotic single-cell RNA sequencing by in situ combinatorial indexing" (doi: 10.1038/s41564-020-0729-6)
## Written by Sydney Blattman
## Tavazoie Lab, Columbia University
## Last updated January 2025

import time
import os
import sys
import preprocess_2_sc15 as preprocess
import rpy2.robjects as robjects
import os.path

start = time.time()
script_dir = sys.argv[0].split(sys.argv[0].split('/')[-1])[0]
i = 1
sample = sys.argv[i][0:sys.argv[i].find('_S')]
if sys.argv[i].find('_S') == -1:
    print('error: must include _S after sample')
    exit()
n_lanes = int(sys.argv[2])

# Clean up old files
cleanup_commands = [
    f'rm -r {sample}/*QF_* 2> /dev/null',
    f'rm -r {sample}_logs/sc_pipeline_15 2> /dev/null',
]

for cmd in cleanup_commands:
    os.system(cmd)

os.system(f'mkdir {sample}_logs')
os.system(f'mkdir {sample}_logs/sc_pipeline_15')

print(f'Preprocessing {sample}')

# Run Fastqc on all lanes
os.system(f'mkdir {sample}_logs/sc_pipeline_15/fastqc')
fastqc_command = (
    f'ls {sample}/*_001.fastq.gz | time parallel --bar --results '
    f'{sample}_logs/sc_pipeline_15/fastqc -j8 fastqc {{}}'
)
os.system(fastqc_command)
print('Fastqc done')

# Trim low quality reads
trim_command = (
    f'seq {n_lanes} | time parallel --bar -j4 cutadapt -q 10,10 --minimum-length 55:14 '
    f'--max-n 3 --pair-filter=any -o {sample}/{sample}_QF_L00{{}}_R1_001.fastq '
    f'-p {sample}/{sample}_QF_L00{{}}_R2_001.fastq {sample}/{sys.argv[i]}_L00{{}}_R1_001.fastq.gz '
    f'{sample}/{sys.argv[i]}_L00{{}}_R2_001.fastq.gz > {sample}_logs/sc_pipeline_15/QF.log'
)
os.system(trim_command)

print('Quality Trim Done')

# Use pear to match read 1 and read 2; for those that overlap, remove reads less than 75bp
pear_command = (
    f'seq {n_lanes} | time parallel --bar -j5 pear -f {sample}/{sample}_QF_L00{{}}_R1_001.fastq '
    f'-r {sample}/{sample}_QF_L00{{}}_R2_001.fastq -o {sample}/{sample}_QF_L00{{}}_p -v 8 -p 0.001 -n 0'
)
os.system(pear_command)

cutadapt_command = (
    f'seq {n_lanes} | time parallel --bar -j4 cutadapt -m 75 -o {sample}/{sample}_QF_L00{{}}_paired_min75_001.fastq '
    f'{sample}/{sample}_QF_L00{{}}_p.assembled.fastq'
)
os.system(cutadapt_command)

# Split paired reads back into two files of read 1 (58 bases) and read 2 (remaining sequence - reverse comp)
split_r1_command = (
    f'seq {n_lanes} | time parallel --bar -j4 cutadapt -l 58 -o {sample}/{sample}_QF_L00{{}}_R1_paired.fastq '
    f'{sample}/{sample}_QF_L00{{}}_paired_min75_001.fastq'
)
os.system(split_r1_command)

split_r2_command = (
    f'seq {n_lanes} | time parallel --bar -j4 cutadapt -u 58 -o {sample}/{sample}_QF_L00{{}}_preR2_paired.fastq '
    f'{sample}/{sample}_QF_L00{{}}_paired_min75_001.fastq'
)
os.system(split_r2_command)
reverse_comp_command = (
    f'seq {n_lanes} | time parallel --bar -j4 seqkit seq -r -p -t DNA '
    f'{sample}/{sample}_QF_L00{{}}_preR2_paired.fastq -o {sample}/{sample}_QF_L00{{}}_R2_paired.fastq'
)
os.system(reverse_comp_command)

# Merge back the reads that overlapped by pear and those that didn't to get clean non-overlapping reads 1 and 2

for l in range(1, n_lanes + 1):
    merge_r1_command = (
        f'cat {sample}/{sample}_QF_L00{l}_p.unassembled.forward.fastq '
        f'{sample}/{sample}_QF_L00{l}_R1_paired.fastq > {sample}/{sample}_QF_merged_L00{l}_R1.fastq'
    )
    os.system(merge_r1_command)
    merge_r2_command = (
        f'cat {sample}/{sample}_QF_L00{l}_p.unassembled.reverse.fastq '
        f'{sample}/{sample}_QF_L00{l}_R2_paired.fastq > {sample}/{sample}_QF_merged_L00{l}_R2.fastq'
    )
    os.system(merge_r2_command)

# Remove intermediate files
remove_intermediate_commands = [
    f'rm -r {sample}/{sample}*p.unassembled.forward.fastq',
    f'rm -r {sample}/{sample}*p.unassembled.reverse.fastq',
    f'rm -r {sample}/{sample}*_paired_min75_001.fastq',
    f'rm -r {sample}/{sample}*preR2_paired.fastq',
    f'rm -r {sample}/{sample}*_p.assembled.fastq',
    f'rm -r {sample}/{sample}*_p.discarded.fastq',
    f'rm -r {sample}/{sample}*_R1_paired.fastq',
    f'rm -r {sample}/{sample}*_R2_paired.fastq'
]

for cmd in remove_intermediate_commands:
    os.system(cmd)

robjects.globalenv['n_lanes'] = n_lanes
robjects.globalenv['script_dir'] = script_dir
robjects.globalenv['sample'] = sample

robjects.r.source(f'{script_dir}/demultiplexer.R')

end = time.time()
print(f"Time elapsed during Python pipeline: {end - start}")
