## edited May 26 2020 so logs folder is made if not already there

dir=$(echo $0 | sed "s/pipeline_v2_generic.sh//")
sample=$1
n_BCs=$2
fasta=$3
gff=$4
custom_name=$5

if [ ! -d "${sample}_logs" ]
then
	mkdir ${sample}_logs
fi

num=1
if [! -d "${sample}_selected_frequency_table.txt" ]
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
	python $dir/align_HiSeq_after_hairpin_trim_v2.py ${sample} ${fasta} ${sample} # align read 2 to genomes
	python $dir/featureCounts_directional_5.py ${sample} ${custom_name} ${sample} ${gff} # annotate features from gff and identify UMI groups
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
	echo $msg
fi
