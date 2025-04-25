## edited May 26 2020 so logs folder is made if not already there

dir=$(echo $0 | sed "s/pipeline_v2_generic.sh//")
sample=$1
fasta=$2
gff=$3

if [ ! -d "${sample}_logs" ]
then
	mkdir ${sample}_logs
fi

num=1
if [ ! -f "${sample}_selected_frequency_table.txt" ]
then
	echo "File ${sample}_selected_frequency_table.txt does not exist. Maybe need to run demultiplexer.R"
	num=0
fi

if [ $num -gt 0 ]
then
	if [ ! -f "${fasta}.bwt" ]
	then
		bwa index ${fasta}
	fi

	bwa aln -n 0.06 ${fasta} ${sample}/${sample}_2trim.fastq > ${sample}/${sample}_bwa.sai
	bwa samse -n 14 ${fasta} ${sample}/${sample}_bwa.sai ${sample}/${sample}_2trim.fastq > ${sample}/${sample}_bwa.sam
	sed "s/XT:/XN:/" ${sample}/${sample}_bwa.sam > ${sample}/${sample}_no_XT.sam # Removes the XT tag from the sam file because it reportably interfer with featureCount
	samtools view -bS ${sample}/${sample}_no_XT.sam | samtools sort - > ${sample}/${sample}_sorted.bam
	samtools index ${sample}/${sample}_sorted.bam
	mkdir ${sample}_FC
	featureCounts -t 'Coding_or_RNA' -g 'name' -s 1 -a ${gff} -o ${sample}_FC/${sample} -R BAM ${sample}/${sample}_sorted.bam # annotate features from gff
	samtools index ${sample}_FC/${sample}_sorted.bam.featureCounts.bam
	Rscript ${dir}/add_cell_barcode.R $sample
	mkdir ${sample}_FC_directional_grouped_2/
	mkdir ${sample}_logs/featureCounts_directional_5
	umi_tools group --per-gene --gene-tag=XT --per-cell --cell-tag=CB --extract-umi-method=tag --umi-tag=BX \
	-I ${sample}_FC/${sample}_sorted.bam.featureCounts_with_celltag.bam \
        --group-out=${sample}_FC_directional_grouped_2/${sample}_UMI_counts.tsv \
        --method=directional --output-bam -S ${sample}_FC_directional_grouped_2/${sample}_group_FC.bam \
        >> ${sample}_logs/featureCounts_directional_5/${sample}_umi_group.log
	samtools view ${sample}_FC_directional_grouped_2/${sample}_group_FC.bam > ${sample}_FC_directional_grouped_2/${sample}_group_FC.sam
	python $dir/sc_sam_processor_11_generic.py 0 ${custom_name} ${sample} # generates a single file of collapsed UMIs (output suffix is _filtered_mapped_UMIs.txt})
	python $dir/make_matrix_mixed_species.py ${custom_name}_v11_threshold_0 # make matrix from filtered UMIs
	## Below are cleanup steps - comment out to see intermediate files
	rm -r ${sample}/${sample}*QF*.fastq
	rm -r ${sample}_R2_trimmed #fastq files after trim_R2_v4
	rm -r ${sample}_2nd_trim
	rm -r ${sample}_bwa_sam # sam files after align_v4
	rm -r ${custom_name}_no_XT # sam and bam files without XT tag (generated during featureCounts_directional_5)
	rm -r ${custom_name}_FC # bam files after feature counts (featureCounts_directional_5)
	rm -r ${custom_name}_FC_directional_grouped_2 # bam files after umi tools group (featureCounts_directional_5)
fi
