# PETRI-seq-persistence
Code for "Identification and genetic dissection of convergent persister cell states"

# Processing PETRI-seq data
-	See PETRI_seq_scripts_v2, which is updated from https://tavazoielab.c2b2.columbia.edu/PETRI-seq/. This software is designed to be broadly useful for processing PETRI-seq data by different users. The folder contains another README with software requirements and instructions for the demo.

# Notebooks to generate all figures for  “Identification and genetic dissection of convergent persister cell states”. System requirements are listed in each notebook.
o	Many notebooks and additional scripts require processed data from GEO (GSE229435). Download supplementary files and put in ‘source_data/from_GEO’.
o	Note that within source_data, only ‘included’ directory contains files. The others are filled as scripts are run or from GEO.

Contents:
-	Figure notebooks\
                - (run first to make intermediate source file)\
 	- figure1_py_v2.ipynb\
                - (figure1_R_v2.ipynb)\
        - figure1_R_v2.ipynb
    o	figure2_py.ipynb
    o	figure2_R_v4.ipynb
    o	figure3_py_v4.ipynb
        	(figure3_R.ipynb, fig3_persister_only.R [first run save_seurat_objects.R])
    o	figure3_R.ipynb
        	(fig3_persister_only.R [first run save_seurat_objects.R])
    o	figure4_py_v7.ipynb
        	(fig3_persister_only.R [first run save_seurat_objects.R])
    o	figure4_R.ipynb
        	(fig3_persister_only.R [first run save_seurat_objects.R])
    o	figure5_v8_py.ipynb
    o	figureED1_R.ipynb
    o	figureED1_py.ipynb
        	(figureED1_R.ipynb, fig1_markers.R)
    o	figureED2_py.ipynb
        	(figED2_markers_vs_clust5.R)
    o	figureED3_py.ipynb
    o	figureED3_R.ipynb
    o	figureED4_R.ipynb
    o	figureED5_py.ipynb
        	(v4_markers_vs_clust1.R, v4_markers_vs_clust4.R)
    o	figureED5_R.ipynb
        	(hipA7_only.R, 6day_only.R, 6day_metG_only.R, ds30_seurat.R, save_seurat_objects.R)
    o	figureED6_R.ipynb
        	(seurat_ds38_min10_CFTonly_v4.R [first run save_seurat_objects.R])
    o	figureED6_py.ipynb
        	(v4_CFT_marker.R, figED2_markers_vs_clust5.R)
    o	figureED7_R.ipynb
        	(proteomics_seurat_v3.R)
    o	figureED7_py.ipynb
        	(fig1_markers.R)
    o	figureED8_py.ipynb
    o	figureED8_R.ipynb
    o	figureED9_py.ipynb
        	(fig3_persister_only.R)
    o	figureED9_R.ipynb
        	(fig3_persister_only.R)
    o	figureED10_R.ipynb
        	(fig3_persister_only.R [first run save_seurat_objects.R])
    o	figureED10_py.ipynb
    o	figureED11_py.ipynb
    o	figureED12_py.ipynb

-	Additional scripts (*.R) to generate needed files for notebooks
-	Source data (‘source_data/included/’)
-	enrichment_and_p_val.py: example script for CRISPRi processing
-	PETRI_seq_scripts_v2
    o	Updated scripts and demo for processing raw PETRI-seq data to generate count matrix; see README.txt in this folder
