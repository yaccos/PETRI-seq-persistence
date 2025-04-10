## edited May 26 2020 (after upload) to handle different feature label in gff

import sam_edit_tools
import sys
import os

def pipeline(sample, gff):
    sam_edit_tools.no_xt_new(f"{sample}/{sample}_bwa_sam.sam", f"{sample}/{sample}_no_XT.sam")
    samtools_view_command = (
        f'samtools view -bS {sample}/{sample}_no_XT.sam | '
        f'samtools sort - > "{sample}/{sample}_no_XT_sorted.bam'
    )
    os.system(samtools_view_command)
    
    samtools_index_command_1 = (
        f'samtools index "{sample}/{sample}_no_XT_sorted.bam'
    )
    os.system(samtools_index_command_1)
    
    featureCounts_command = (
        f"featureCounts -t 'Coding_or_RNA' -g 'name' -s 1 -a {gff} -o {sample}_FC/{sample} "
        f"-R BAM {sample}/{sample}_no_XT_sorted.bam"
    )
    os.system(featureCounts_command)
    
    samtools_index_command_2 = (
        f'samtools index {sample}_FC/{sample}_no_XT_sorted.bam.featureCounts.bam'
    )
    os.system(samtools_index_command_2)
    
    umi_tools_command = (
        f'umi_tools group --per-gene --gene-tag=XT --per-cell --cell-tag=??? -I {sample}_FC/{sample}_no_XT_sorted.bam.featureCounts.bam '
        f'--group-out={sample}_FC_directional_grouped_2/{sample}_UMI_counts.tsv '
        f'--method=directional --output-bam -S {sample}_FC_directional_grouped_2/{sample}_group_FC.bam '
        f'>> {sample}_logs/featureCounts_directional_5/{sample}_umi_group.log'
    )
    os.system(umi_tools_command)

sample = sys.argv[1] 
gff = sys.argv[2]

os.system('mkdir ' + sample + '_FC')
os.system('mkdir ' + sample + '_FC_directional_grouped_2')
os.system('mkdir ' + old_sample + '_logs/featureCounts_directional_5')
pipeline(sample, gff)
