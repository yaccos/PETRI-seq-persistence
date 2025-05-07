dir=$(echo $0 | sed "s/pipeline_v2_generic.sh//")
sample=$1
fasta=$2
gff=$3


bwa aln -n 0.06 ${fasta} results/${sample}/${sample}_2trim.fastq > results/${sample}/${sample}_bwa.sai
bwa samse -n 14 ${fasta} results/${sample}/${sample}_bwa.sai results/${sample}/${sample}_2trim.fastq > results/${sample}/${sample}_bwa.sam
sed "s/XT:/XN:/" results/${sample}/${sample}_bwa.sam > results/${sample}/${sample}_no_XT.sam # Removes the XT tag from the sam file because it reportably interfer with featureCount
samtools view -bS results/${sample}/${sample}_no_XT.sam | samtools sort - > results/${sample}/${sample}_sorted.bam
samtools index results/${sample}/${sample}_sorted.bam
mkdir results/${sample}_FC
featureCounts -t 'Coding_or_RNA' -g 'name' -s 1 -a ${gff} -o results/${sample}_FC/${sample} -R BAM results/${sample}/${sample}_sorted.bam # annotate features from gff
samtools index results/${sample}_FC/${sample}_sorted.bam.featureCounts.bam
Rscript ${dir}/add_cell_barcode.R $sample
mkdir results/${sample}_FC_directional_grouped_2
umi_tools group --per-gene --gene-tag=XT --per-cell --cell-tag=CB --extract-umi-method=tag --umi-tag=BX \
-I results/${sample}_FC/${sample}_sorted.bam.featureCounts_with_celltag.bam \
	--group-out=results/${sample}_FC_directional_grouped_2/${sample}_UMI_counts.tsv \
	--method=directional --output-bam -S results/${sample}_FC_directional_grouped_2/${sample}_group_FC.bam
samtools view results/${sample}_FC_directional_grouped_2/${sample}_group_FC.bam > results/${sample}_FC_directional_grouped_2/${sample}_group_FC.sam
python $dir/sc_sam_processor_11_generic.py 0 ${sample} # generates a single file of collapsed UMIs (output suffix is _filtered_mapped_UMIs.txt})
python $dir/make_matrix_mixed_species.py ${sample}_v11_threshold_0 # make matrix from filtered UMIs
