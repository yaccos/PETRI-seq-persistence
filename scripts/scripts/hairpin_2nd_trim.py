import sys
import os
import multiprocessing
import pandas as pd
print(sys.argv[0])
script_dir = sys.argv[0].split('hairpin_2nd_trim.py')[0]
print(script_dir)
def trim_R2(sample,ID,name,bc1_to_seq):
    R2_file_name = ID + '_R2' + name[name.find('_bc1_'):len(name)]
    R2_file_name = R2_file_name.replace(' ' , '')
    R2_short_name = R2_file_name.replace(ID+'_','')
    barcodes = R2_file_name.split('_bc')
    bc1 = barcodes[1].split('_')[1]
    print('2nd trim:' + R2_file_name)
    bc2 = barcodes[2].split('_')[1]
    bc3 = barcodes[3].split('_')[1]
    bc1_sequence = bc1_to_seq[bc1]
    bc2_sequence = bc2_to_seq[bc2]
    bc3_sequence = bc3_to_seq[bc3]
    full_seq = f"{bc3_sequence}GGTCCTTGGCTTCGC{bc2_sequence}CCTCCTACGCCAGA{bc1_sequence}"
    cutadapt_command = (
        f"cutadapt -b {full_seq} -O 7 -e 0.2 --minimum-length 16 "
        f"-o {sample}_2nd_trim/{R2_file_name}_2trim.fastq "
        f"{sample}_R2_trimmed/{R2_file_name}_R2_trimmed.fastq "
        f">> {sample}_logs/2nd_trim/2nd_trim.log"
    )
    os.system(cutadapt_command)

bc1_list = open(script_dir + '/sc_barcodes_v2/BC1.fa')
bc1_to_seq = {}
for line in bc1_list:
    if '>' in line:
        barcode = line.split('bc1_')[1].split('\n')[0]
    else:
        seq = line.split('N$')[0]
        bc1_to_seq[barcode] = seq
bc2_list = open(script_dir + '/sc_barcodes_v2/BC2_anchored.fa')
bc2_to_seq = {}
for line in bc2_list:
    if '>' in line:
        barcode = line.split('bc2_')[1].split('\n')[0]
    else:
        seq = line.split('CCTCCTACGCCAGA')[0].replace('^','')
        bc2_to_seq[barcode] = seq
bc3_list = open(script_dir + '/sc_barcodes_v2/BC3_anchored.fa')
bc3_to_seq = {}
for line in bc3_list:
    if '>' in line:
        barcode = line.split('bc3_')[1].split('\n')[0]
    else:
        seq = line.split('GGTCCTTGGCTTCGC')[0].replace('^','')
        bc3_to_seq[barcode] = seq

sample = sys.argv[1]
if len(sys.argv) > 2:
    ID = sys.argv[2]
else:
    ID = sample
table = open(sample + '_selected_cumulative_frequency_table.txt')
os.system(f'mkdir {sample}_logs/2nd_trim')
bc1_list = open(script_dir + '/sc_barcodes_v2/BC1.fa')

os.system(f'mkdir {sample}_2nd_trim')
for line in table:
    name = line.split('\t')[1]
    trim_R2(sample,ID,name,bc1_to_seq)
