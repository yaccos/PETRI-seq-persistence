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
    f'rm -r {sample}_bc3 2> /dev/null',
    f'rm -r {sample}_bc2 2> /dev/null',
    f'rm -r {sample}_bc1 2> /dev/null',
    f'rm -r {sample}_bc1_table.txt 2> /dev/null',
    f'rm -r {sample}_bc2_table.txt 2> /dev/null',
    f'rm -r {sample}_bc3_table.txt 2> /dev/null',
    f'rm -r {sample}_logs/sc_pipeline_15 2> /dev/null',
    f'rm -r {sample}_bc1_cumulative_frequency_table.txt 2> /dev/null',
    f'rm -r {sample}_bc2_cumulative_frequency_table.txt 2> /dev/null',
    f'rm -r {sample}_bc3_cumulative_frequency_table.txt 2> /dev/null'
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

robjects.r.source(f'{script_dir}/demultiplexer.R')

# Trim low quality reads
trim_command = (
    f'seq {n_lanes} | time parallel --bar -j4 cutadapt -q 10,10 --minimum-length 55:14 '
    f'--max-n 3 --pair-filter=any -o {sample}/{sample}_QF_L00{{}}_R1_001.fastq.gz '
    f'-p {sample}/{sample}_QF_L00{{}}_R2_001.fastq.gz {sample}/{sys.argv[i]}_L00{{}}_R1_001.fastq.gz '
    f'{sample}/{sys.argv[i]}_L00{{}}_R2_001.fastq.gz > {sample}_logs/sc_pipeline_15/QF.log'
)
os.system(trim_command)

print('Quality Trim Done')

# Use pear to match read 1 and read 2; for those that overlap, remove reads less than 75bp
pear_command = (
    f'seq {n_lanes} | time parallel --bar -j5 pear -f {sample}/{sample}_QF_L00{{}}_R1_001.fastq.gz '
    f'-r {sample}/{sample}_QF_L00{{}}_R2_001.fastq.gz -o {sample}/{sample}_QF_L00{{}}_p -v 8 -p 0.001 -n 0'
)
os.system(pear_command)

cutadapt_command = (
    f'seq {n_lanes} | time parallel --bar -j4 cutadapt -m 75 -o {sample}/{sample}_QF_L00{{}}_paired_min75_001.fastq.gz '
    f'{sample}/{sample}_QF_L00{{}}_p.assembled.fastq'
)
os.system(cutadapt_command)

# Split paired reads back into two files of read 1 (58 bases) and read 2 (remaining sequence - reverse comp)
split_r1_command = (
    f'seq {n_lanes} | time parallel --bar -j4 cutadapt -l 58 -o {sample}/{sample}_QF_L00{{}}_R1_paired.fastq.gz '
    f'{sample}/{sample}_QF_L00{{}}_paired_min75_001.fastq.gz'
)
os.system(split_r1_command)

split_r2_command = (
    f'seq {n_lanes} | time parallel --bar -j4 cutadapt -u 58 -o {sample}/{sample}_QF_L00{{}}_preR2_paired.fastq.gz '
    f'{sample}/{sample}_QF_L00{{}}_paired_min75_001.fastq.gz'
)
os.system(split_r2_command)
reverse_comp_command = (
    f'seq {n_lanes} | time parallel --bar -j4 seqkit seq -r -p -t DNA '
    f'{sample}/{sample}_QF_L00{{}}_preR2_paired.fastq.gz -o {sample}/{sample}_QF_L00{{}}_R2_paired.fastq.gz'
)
os.system(reverse_comp_command)

# Merge back the reads that overlapped by pear and those that didn't to get clean non-overlapping reads 1 and 2
gzip_forward_command = (
    f'seq {n_lanes} | time parallel --bar -j4 gzip {sample}/{sample}_QF_L00{{}}_p.unassembled.forward.fastq'
)
os.system(gzip_forward_command)

gzip_reverse_command = (
    f'seq {n_lanes} | time parallel --bar -j4 gzip {sample}/{sample}_QF_L00{{}}_p.unassembled.reverse.fastq'
)
os.system(gzip_reverse_command)
for l in range(1, n_lanes + 1):
    merge_r1_command = (
        f'gunzip -c {sample}/{sample}_QF_L00{l}_p.unassembled.forward.fastq.gz '
        f'{sample}/{sample}_QF_L00{l}_R1_paired.fastq.gz | gzip > {sample}/{sample}_QF_merged_L00{l}_R1.fastq.gz'
    )
    os.system(merge_r1_command)
    merge_r2_command = (
        f'gunzip -c {sample}/{sample}_QF_L00{l}_p.unassembled.reverse.fastq.gz '
        f'{sample}/{sample}_QF_L00{l}_R2_paired.fastq.gz | gzip > {sample}/{sample}_QF_merged_L00{l}_R2.fastq.gz'
    )
    os.system(merge_r2_command)

# Extract UMI
umi_command = (
    f'seq {n_lanes} | time parallel --bar -j5 umi_tools extract --stdin={sample}/{sample}_QF_merged_L00{{}}_R1.fastq.gz '
    f'--bc-pattern=NNNNNNN --read2-in={sample}/{sample}_QF_merged_L00{{}}_R2.fastq.gz '
    f'--log={sample}_logs/sc_pipeline_15/UMI_extract.log --stdout {sample}/{sample}_QF_UMI_L00{{}}_R1_001.fastq.gz '
    f'--read2-out={sample}/{sample}_QF_UMI_L00{{}}_R2_001.fastq.gz'
)
os.system(umi_command)

# Remove intermediate files
remove_intermediate_commands = [
    f'rm -r {sample}/{sample}*p.unassembled.forward.fastq.gz',
    f'rm -r {sample}/{sample}*p.unassembled.reverse.fastq.gz',
    f'rm -r {sample}/{sample}*_paired_min75_001.fastq.gz',
    f'rm -r {sample}/{sample}*preR2_paired.fastq.gz',
    f'rm -r {sample}/{sample}*_p.assembled.fastq',
    f'rm -r {sample}/{sample}*_p.discarded.fastq',
    f'rm -r {sample}/{sample}*_R1_paired.fastq.gz',
    f'rm -r {sample}/{sample}*_R2_paired.fastq.gz'
]

for cmd in remove_intermediate_commands:
    os.system(cmd)

# Demultiplex by bc3
os.system(f'mkdir {sample}_bc3')
if os.path.exists(f'{sample}_logs/sc_pipeline_15/bc3.log'):
    os.system(f'rm {sample}_logs/sc_pipeline_15/bc3.log')

# Note: bc3 adapter sequences are 23 nt long. The error rate is set to 0.05, but this creates ambiguities
# causing trouble for demultiplexing.
bc3_command = (
    f'seq {n_lanes} | time parallel --bar -j5 cutadapt -g file:{script_dir}sc_barcodes_v2/BC3_anchored.fa '
    f'-e 0.05 --overlap 21 '
    f'--untrimmed-output {sample}_bc3/{sample}_no_bc3_L00{{}}_R1_001.fastq.gz '
    f'--untrimmed-paired-output {sample}_bc3/{sample}_no_bc3_L00{{}}_R2_001.fastq.gz '
    f'-o {sample}_bc3/{sample}_{{name}}x_L00{{}}_R1_001.fastq.gz -p {sample}_bc3/{sample}_{{name}}x_L00{{}}_R2_001.fastq.gz '
    f'{sample}/{sample}_QF_UMI_L00{{}}_R1_001.fastq.gz {sample}/{sample}_QF_UMI_L00{{}}_R2_001.fastq.gz >> '
    f'{sample}_logs/sc_pipeline_15/bc3.log'
)
os.system(bc3_command)
print('script dir:', script_dir)

os.system(f'cd {sample}_bc3 && python {script_dir}merge_lanes_mac_compatible.py')
print('bc3 done')


# Demultiplex by bc2
os.system(f'mkdir {sample}_bc2')
if os.path.exists(f'{sample}_logs/sc_pipeline_15/bc2.log'):
    os.system(f'rm {sample}_logs/sc_pipeline_15/bc2.log')
bc3_list = ''
for i in range(1, 97):
    if os.path.exists(f'{sample}_bc3/{sample}_bc3_{i}x_R1_all_lanes.fastq.gz'):
        if bc3_list == '':
            bc3_list = str(i)
        else:
            bc3_list = f'{bc3_list}\n{i}'
bc2_command = (
    f'echo "{bc3_list}" | time parallel --bar -j12 cutadapt -g file:{script_dir}sc_barcodes_v2/BC2_anchored.fa '
    f'-e 0.05 --overlap 20 --untrimmed-output '
    f'{sample}_bc2/{sample}_bc1_{{}}_R1_no_bc2.fastq.gz '
    f'--untrimmed-paired-output {sample}_bc2/{sample}_bc3_{{}}_R2_no_bc2.fastq.gz '
    f'-o {sample}_bc2/{sample}_R1_{{name}}_bc3_{{}}.fastq.gz -p {sample}_bc2/{sample}_R2_{{name}}_bc3_{{}}.fastq.gz '
    f'{sample}_bc3/{sample}_bc3_{{}}x_R1_all_lanes.fastq.gz {sample}_bc3/{sample}_bc3_{{}}x_R2_all_lanes.fastq.gz >> '
    f'{sample}_logs/sc_pipeline_15/bc2.log'
)
os.system(bc2_command)
print('bc2 done')

# Checkpoint to be sure all bc3 files were demultiplexed
n_R1 = len([name for name in os.listdir(f'{sample}_bc3') if 'R1' in name and 'no_bc3' not in name])
n_R2 = len([name for name in os.listdir(f'{sample}_bc3') if 'R2' in name and 'no_bc3' not in name])
bc2_log = open(f'{sample}_logs/sc_pipeline_15/bc2.log', 'r').read()
expected_n = bc2_log.count("Summary") + bc2_log.count("No reads processed!")
if (n_R1 != expected_n) | (n_R2 != expected_n):
    print('ERROR: total demultiplexed bc2 files do not match expected input from bc3. Maybe process was disrupted?')
    quit()

os.system(f'rm -r {sample}_bc3')

# Demultiplex by bc1
os.system(f'mkdir {sample}_bc1')
if os.path.exists(f'{sample}_logs/sc_pipeline_15/bc1.log'):
    os.system(f'rm {sample}_logs/sc_pipeline_15/bc1.log')
for bc3 in range(1, 97):
    bc2_list = ''
    for bc2 in range(1, 97):
        if os.path.exists(f'{sample}_bc2/{sample}_R1_bc2_{bc2}_bc3_{bc3}.fastq.gz'):
            if bc2_list == '':
                bc2_list = str(bc2)
            else:
                bc2_list = f'{bc2_list}\n{bc2}'
    print(bc3)
    if bc2_list != '':
        bc1_command = (
            f'echo "{bc2_list}" | time parallel --bar -j12 cutadapt -g file:{script_dir}sc_barcodes_v2/BC1_5p_anchor_v2.fa '
            f'-e 0.2 --no-indels --overlap 7 --discard-untrimmed --action=none '
            f'-o {sample}_bc1/{sample}_R1_{{name}}_bc2_{{}}_bc3_{bc3}.fastq.gz '
            f'-p {sample}_bc1/{sample}_R2_{{name}}_bc2_{{}}_bc3_{bc3}.fastq.gz '
            f'{sample}_bc2/{sample}_R1_bc2_{{}}_bc3_{bc3}.fastq.gz {sample}_bc2/{sample}_R2_bc2_{{}}_bc3_{bc3}.fastq.gz >> '
            f'{sample}_logs/sc_pipeline_15/bc1.log'
        )
        os.system(bc1_command)

# Checkpoint to be sure all bc2 files were demultiplexed
n_R1 = len([name for name in os.listdir(f'{sample}_bc2') if 'R1' in name and 'no_bc2' not in name])
n_R2 = len([name for name in os.listdir(f'{sample}_bc2') if 'R2' in name and 'no_bc2' not in name])
bc_1_log = open(f'{sample}_logs/sc_pipeline_15/bc1.log', 'r').read()
expected_n = bc_1_log.count("Summary") + bc_1_log.count("No reads processed!")
if (n_R1 != expected_n) | (n_R2 != expected_n):
    print('ERROR: total demultiplexed bc1 files do not match expected input from bc2. Maybe process was disrupted?')
    quit()

os.system(f'rm -r {sample}_bc2')
if os.path.exists(f'{sample}_bc1_cumulative_frequency_table.txt'):
    os.system(f'rm {sample}_bc1_cumulative_frequency_table.txt')

preprocess.log_to_table(sample, 'bc1')
preprocess.log_plot(sample, 'bc1', 0)
preprocess.freq_plot(sample, 'bc1', 0)
os.system(f'rm {sample}_bc1_table.txt')
print('bc1 done')

end = time.time()
print(f"Time elapsed during Python pipeline: {end - start}")
