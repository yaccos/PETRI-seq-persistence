## edited May 26 2020 (after upload) to handle different feature label in gff

import sys
import os

def pipeline(sample):
    
    umi_tools_command = (
        f'umi_tools group --per-gene --gene-tag=XT --per-cell --cell-tag=??? -I {sample}_FC/{sample}_no_XT_sorted.bam.featureCounts.bam '
        f'--group-out={sample}_FC_directional_grouped_2/{sample}_UMI_counts.tsv '
        f'--method=directional --output-bam -S {sample}_FC_directional_grouped_2/{sample}_group_FC.bam '
        f'>> {sample}_logs/featureCounts_directional_5/{sample}_umi_group.log'
    )
    os.system(umi_tools_command)

sample = sys.argv[1] 

os.system('mkdir ' + sample + '_FC_directional_grouped_2')
os.system('mkdir ' + sample + '_logs/featureCounts_directional_5')
pipeline(sample)
