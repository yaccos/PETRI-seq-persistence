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
    R1_file_name = f"{sample}_bc1/{sample}_R1{name[name.find('_bc1_'):len(name)]}.fastq.gz".replace(' ', '')
    R2_file_name = f"{sample}_bc1/{sample}_R2{name[name.find('_bc1_'):len(name)]}.fastq.gz".replace(' ', '')
    if int(line.split('\t')[5]) < int(threshold):
        output.write(line)
    elif int(line.split('\t')[2]) > 0:
        os.system(f'rm {R1_file_name}')
        os.system(f'rm {R2_file_name}')
os.system(f'rm {sample}_bc1/*_unknown*')
os.system(f'mv {sample}_bc1/ {sample}_selected_cells')
output.close()
table.close()
