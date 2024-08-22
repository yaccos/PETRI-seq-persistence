## Supplement to "Prokaryotic single-cell RNA sequencing by in situ combinatorial indexing" (doi: 10.1038/s41564-020-0729-6)
## Written by Sydney Blattman
## Tavazoie Lab, Columbia University
## Last updated July 2024 

import time
import os
import sys
import preprocess_2_sc15 as preprocess
import os.path
import os, errno

start = time.time()
script_dir = sys.argv[0].split(sys.argv[0].split('/')[-1])[0]
i=1
sample = sys.argv[i][0:sys.argv[i].find('_S')]
if(sys.argv[i].find('_S') == -1):
    print('error: must include _S after after sample')
    exit()
n_lanes = int(sys.argv[2])

## Clean up old files
os.system('rm -r ' + sample + '/*QF_* 2> rm.log')
os.system('rm -r ' + sample + '_bc3 2> rm.log')
os.system('rm -r ' + sample + '_bc2 2> rm.log')
os.system('rm -r ' + sample + '_bc1 2> rm.log')
os.system('rm -r ' + sample + '_bc1_table.txt 2> rm.log')
os.system('rm -r ' + sample + '_bc2_table.txt 2> rm.log')
os.system('rm -r ' + sample + '_bc3_table.txt 2> rm.log')
os.system('rm -r ' + sample + '_logs/sc_pipeline_15 2> rm.log')
os.system('rm -r ' + sample + '_bc1_cumulative_frequency_table.txt 2> rm.log')
os.system('rm -r ' + sample + '_bc2_cumulative_frequency_table.txt 2> rm.log')
os.system('rm -r ' + sample + '_bc3_cumulative_frequency_table.txt 2> rm.log')
if(os.path.exists('rm.log')):
    os.system('rm rm.log')

os.system('mkdir ' + sample + '_logs')
os.system('mkdir ' + sample + '_logs/sc_pipeline_15')

print('Preprocessing ' + sample)

## Run Fastqc on all lanes
os.system('mkdir ' + sample + '_logs/sc_pipeline_15/fastqc')
os.system('ls ' + sample + '/*_001.fastq.gz | time parallel --bar --results ' + sample + '_logs/sc_pipeline_15/fastqc -j8 fastqc {}')
print('Fastqc done')

## trim low quality reads
os.system('seq ' + str(n_lanes) + ' | time parallel --bar -j4 cutadapt -q 10,10 --minimum-length 55:14 --max-n 3 --pair-filter=any -o ' + sample + '/' + sample + '_QF_L00{}_R1_001.fastq.gz -p ' + sample + '/' + sample + '_QF_L00{}_R2_001.fastq.gz ' + sample + '/' + sys.argv[i] + '_L00{}_R1_001.fastq.gz ' + sample + '/' + sys.argv[i] + '_L00{}_R2_001.fastq.gz > ' + sample + '_logs/sc_pipeline_15/QF.log')
print('Quality Trim Done')

## use pear to match read 1 and read 2; for those that overlap, remove reads less than 75bp (too little cDNA to align) 
os.system('seq ' + str(n_lanes) + ' | time parallel --bar -j5 pear -f ' + sample + '/' + sample + '_QF_L00{}_R1_001.fastq.gz -r ' + sample + '/' + sample + '_QF_L00{}_R2_001.fastq.gz -o ' + sample + '/' + sample + '_QF_L00{}_p -v 8 -p 0.001 -n 0 -k')
os.system('seq ' + str(n_lanes) + ' | time parallel --bar -j4 cutadapt -m 75 -o ' + sample + '/' + sample + '_QF_L00{}_paired_min75_001.fastq.gz ' + sample + '/' + sample + '_QF_L00{}_p.assembled.fastq')

## split paired reads back into two files of read 1 (58 bases) and read 2 (remaining sequence - reverse comp)
os.system('seq ' + str(n_lanes) + ' | time parallel --bar -j4 cutadapt -l 58 -o ' + sample + '/' + sample + '_QF_L00{}_R1_paired.fastq.gz ' + sample + '/' + sample + '_QF_L00{}_paired_min75_001.fastq.gz') 
os.system('seq ' + str(n_lanes) + ' | time parallel --bar -j4 cutadapt -u 58 -o ' + sample + '/' + sample + '_QF_L00{}_preR2_paired.fastq.gz ' + sample + '/' + sample + '_QF_L00{}_paired_min75_001.fastq.gz')
os.system('seq ' + str(n_lanes) + ' | time parallel --bar -j4 seqkit seq -r -p -t DNA ' + sample + '/' + sample + '_QF_L00{}_preR2_paired.fastq.gz -o ' + sample + '/' + sample + '_QF_L00{}_R2_paired.fastq.gz') # reverse comp for R2

## merge back the reads that overlapped by pear and those that didn't to get clean non-overlapping reads 1 and 2
os.system('seq ' + str(n_lanes) + ' | time parallel --bar -j4 gzip ' + sample + '/' + sample + '_QF_L00{}_p.unassembled.forward.fastq')
os.system('seq ' + str(n_lanes) + ' | time parallel --bar -j4 gzip ' + sample + '/' + sample + '_QF_L00{}_p.unassembled.reverse.fastq')
for l in range(1,n_lanes+1):
    os.system('gunzip -c ' + sample + '/' + sample + '_QF_L00' + str(l) + '_p.unassembled.forward.fastq.gz ' + sample + '/' + sample + '_QF_L00' + str(l) + '_R1_paired.fastq.gz | gzip > ' + sample + '/' + sample + '_QF_merged_L00' + str(l) + '_R1.fastq.gz')
    os.system('gunzip -c ' + sample + '/' + sample + '_QF_L00' + str(l) + '_p.unassembled.reverse.fastq.gz ' + sample + '/' + sample + '_QF_L00' + str(l) + '_R2_paired.fastq.gz | gzip > ' + sample + '/' + sample + '_QF_merged_L00' + str(l) + '_R2.fastq.gz')

## extract UMI
os.system('seq ' + str(n_lanes) + ' | time parallel --bar -j5 umi_tools extract --stdin=' + sample + '/' + sample + '_QF_merged_L00{}_R1.fastq.gz --bc-pattern=NNNNNNN --read2-in=' + sample + '/' + sample + '_QF_merged_L00{}_R2.fastq.gz --log=' + sample + '_logs/sc_pipeline_15/UMI_extract.log --stdout ' + sample + '/' + sample + '_QF_UMI_L00{}_R1_001.fastq.gz --read2-out=' + sample + '/' +  sample + '_QF_UMI_L00{}_R2_001.fastq.gz')

os.system('rm -r ' + sample + '/' + sample + '*p.unassembled.forward.fastq.gz')
os.system('rm -r ' + sample + '/' + sample + '*p.unassembled.reverse.fastq.gz')
os.system('rm -r ' + sample + '/' + sample + '*_paired_min75_001.fastq.gz')
os.system('rm -r ' + sample + '/' + sample + '*preR2_paired.fastq.gz')
os.system('rm -r ' + sample + '/' + sample + '*_p.assembled.fastq')
os.system('rm -r ' + sample + '/' + sample + '*_p.discarded.fastq')
os.system('rm -r ' + sample + '/' + sample + '*_R1_paired.fastq.gz') 
os.system('rm -r ' + sample + '/' + sample + '*_R2_paired.fastq.gz')
## demultiplex by bc3
os.system('mkdir ' + sample + '_bc3')
if(os.path.exists(sample + '_logs/sc_pipeline_15/bc3.log')):
    os.system('rm ' + sample + '_logs/sc_pipeline_15/bc3.log')
os.system('seq ' + str(n_lanes) + ' | time parallel --bar -j5 cutadapt -g file:'+script_dir+'sc_barcodes_v2/BC3_anchored.fa -e 0.05 --overlap 21 --untrimmed-output ' + sample + '_bc3/' + sample + '_no_bc3_L00{}_R1_001.fastq.gz  --untrimmed-paired-output ' + sample + '_bc3/' + sample + '_no_bc3_L00{}_R2_001.fastq.gz  -o ' + sample + '_bc3/' + sample + '_{name}x_L00{}_R1_001.fastq.gz -p ' + sample + '_bc3/' + sample + '_{name}x_L00{}_R2_001.fastq.gz ' + sample + '/' + sample + '_QF_UMI_L00{}_R1_001.fastq.gz ' + sample + '/' + sample + '_QF_UMI_L00{}_R2_001.fastq.gz >> ' + sample + '_logs/sc_pipeline_15/bc3.log')
print('script dir: ' + script_dir)

os.system('cd ' + sample + '_bc3 && python ' + script_dir + 'merge_lanes_mac_compatible.py')
print('bc3 done')

## demultiplex by bc2
os.system('mkdir ' + sample + '_bc2')
if(os.path.exists(sample + '_logs/sc_pipeline_15/bc2.log')):
    os.system('rm ' + sample + '_logs/sc_pipeline_15/bc2.log')
bc3_list = ''
for i in range(1,97):
    if(os.path.exists(sample + '_bc3/' + sample + '_bc3_' + str(i) + 'x_R1_all_lanes.fastq.gz')):
        if bc3_list == '':
            bc3_list = str(i)
        else:
            bc3_list = bc3_list + '\n' + str(i)
os.system('echo "' + bc3_list + '" | time parallel --bar -j12 cutadapt -g file:'+script_dir+'sc_barcodes_v2/BC2_anchored.fa -e 0.05 --overlap 20 --untrimmed-output ' + sample + '_bc2/' + sample + '_bc1_{}_R1_no_bc2.fastq.gz  --untrimmed-paired-output ' + sample + '_bc2/' + sample + '_bc3_{}_R2_no_bc2.fastq.gz -o ' + sample + '_bc2/' + sample + '_R1_{name}_bc3_{}.fastq.gz -p ' + sample + '_bc2/' + sample + '_R2_{name}_bc3_{}.fastq.gz ' + sample + '_bc3/' + sample +'_bc3_{}x_R1_all_lanes.fastq.gz ' + sample + '_bc3/' + sample + '_bc3_{}x_R2_all_lanes.fastq.gz >> ' + sample + '_logs/sc_pipeline_15/bc2.log')
print('bc2 done')

#### checkpoint to be sure all bc3 files were demultiplexed
n_R1 = len([name for name in os.listdir(sample + '_bc3') if (('R1' in name) and ('no_bc3' not in name))])
n_R2 = len([name for name in os.listdir(sample + '_bc3') if (('R2' in name) and ('no_bc3' not in name))])
expected_n = open(sample + '_logs/sc_pipeline_15/bc2.log', 'r').read().count("Summary")
if (n_R1 != expected_n) | (n_R2 != expected_n):
    print('ERROR: total demultiplexed bc2 files do not match expected input from bc3. Maybe process was disrupted?') 
    quit()

os.system('rm -r ' + sample + '_bc3')

## demultiplex by bc1
os.system('mkdir ' + sample + '_bc1')
if(os.path.exists(sample + '_logs/sc_pipeline_15/bc1.log')):
    os.system('rm ' + sample + '_logs/sc_pipeline_15/bc1.log')
for bc3 in range(1,97):
    bc2_list = ''
    for bc2 in range(1,97):
        if(os.path.exists(sample + '_bc2/' + sample + '_R1_bc2_' + str(bc2) + '_bc3_' + str(bc3) + '.fastq.gz')):
            if bc2_list == '':
                bc2_list = str(bc2)
            else:
                bc2_list = bc2_list + '\n' + str(bc2)
    print(bc3)
    if bc2_list != '':
        os.system('echo "' + bc2_list + '" | time parallel --bar -j12 cutadapt -g file:'+script_dir+'sc_barcodes_v2/BC1_5p_anchor_v2.fa -e 0.2 --no-indels --overlap 7 --no-trim -o ' + sample + '_bc1/' + sample + '_R1_{name}_bc2_{}_bc3_' + str(bc3) + '.fastq.gz -p ' + sample + '_bc1/' + sample + '_R2_{name}_bc2_{}_bc3_' + str(bc3) + '.fastq.gz ' + sample + '_bc2/' + sample +'_R1_bc2_{}_bc3_' + str(bc3) + '.fastq.gz ' + sample + '_bc2/' + sample + '_R2_bc2_{}_bc3_' + str(bc3) + '.fastq.gz >> ' + sample + '_logs/sc_pipeline_15/bc1.log')

#### checkpoint to be sure all bc2 files were demultiplexed
n_R1 = len([name for name in os.listdir(sample + '_bc2') if (('R1' in name) and ('no_bc2' not in name))])
n_R2 = len([name for name in os.listdir(sample + '_bc2') if (('R2' in name) and ('no_bc2' not in name))])
expected_n = open(sample + '_logs/sc_pipeline_15/bc1.log', 'r').read().count("Summary")
if (n_R1 != expected_n) | (n_R2 != expected_n):
    print('ERROR: total demultiplexed bc1 files do not match expected input from bc2. Maybe process was disrupted?')
    quit()


os.system('rm -r ' + sample + '_bc2')
if(os.path.exists(sample + '_bc1_cumulative_frequency_table.txt')):
    os.system('rm ' + sample + '_bc1_cumulative_frequency_table.txt')


preprocess.log_to_table(sample,'bc1')
preprocess.log_plot(sample,'bc1',0)
preprocess.freq_plot(sample,'bc1',0)
os.system('rm ' + sample + '_bc1_table.txt')
print('bc1 done')

end = time.time()
print(end-start)

