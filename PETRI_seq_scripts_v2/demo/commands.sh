script_dir=$1 #full path to scripts folder

python $script_dir/sc_pipeline_15_generic_v2.py random20000_S1 1
$script_dir/pipeline_v2_generic.sh random20000 100 ../scripts/U00096_JE2.fa ../scripts/U00096_JE2_rRNA.gff random20000

