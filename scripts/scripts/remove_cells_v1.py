import pandas as pd
import sys
import os

sample = sys.argv[1]
threshold = sys.argv[2]

table = open(f'{sample}_bc1_cumulative_frequency_table.txt')
output = open(f'{sample}_selected_cumulative_frequency_table.txt', 'w')
next(table)

for line in table:
    name = line.split('\t')[1]
    R1_file_name = f"{sample}_bc1/{sample}_R1{name[name.find('_bc1_'):len(name)]}.fastq".replace(' ', '')
    R2_file_name = f"{sample}_bc1/{sample}_R2{name[name.find('_bc1_'):len(name)]}.fastq".replace(' ', '')
    if int(line.split('\t')[5]) < int(threshold):
        output.write(line)
    # We remove all fastq files that are not selected even if they are empty
    else:
        if os.path.exists(R1_file_name):
            os.unlink(R1_file_name)
        if os.path.exists(R2_file_name):
            os.unlink(R2_file_name)
os.system(f'rm {sample}_bc1/*_unknown*')
os.system(f'mv {sample}_bc1/ {sample}_selected_cells')
output.close()
table.close()
