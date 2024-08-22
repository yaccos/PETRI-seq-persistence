# PETRI-seq-persistence
Code for "Identification and genetic dissection of convergent persister cell states"

# Processing PETRI-seq data
-	See PETRI_seq_scripts_v2, which is updated from https://tavazoielab.c2b2.columbia.edu/PETRI-seq/. The folder contains another README with software requirements and instructions for the demo.

# Notebooks to generate all figures for  “Identification and genetic dissection of convergent persister cell states”. 
- System requirements are listed in each notebook.
- Many notebooks and additional scripts require processed data from GEO (GSE229435). Download supplementary files and put in ‘source_data/from_GEO’.
- Note that within source_data, only ‘included’ directory contains files. The others are filled as scripts are run or from GEO.

Contents:
-	Figure notebooks
 	- figure1_py_v2.ipynb
 	- figure1_R_v2.ipynb
 	- figure2_py.ipynb
 	- figure2_R_v4.ipynb
 	- figure3_py_v4.ipynb\
                - (first run save_seurat_objects.R, figure3_R.ipynb, fig3_persister_only.R)
 	- figure3_R.ipynb\
                - (first run save_seurat_objects.R, fig3_persister_only.R)
 	- figure4_py_v7.ipynb\
                - (first run save_seurat_objects.R, fig3_persister_only.R)
 	- figure4_R.ipynb\
                - (first run save_seurat_objects.R, fig3_persister_only.R)
 	- figure5_v8_py.ipynb
 	- figureED1_R.ipynb
 	- figureED1_py.ipynb\
                - (first run figureED1_R.ipynb, fig1_markers.R)
 	- figureED2_py.ipynb\
                - (first run figED2_markers_vs_clust5.R)
 	- figureED3_py.ipynb
 	- figureED3_R.ipynb
 	- figureED4_R.ipynb
 	- figureED5_py.ipynb\
                - (first run v4_markers_vs_clust1.R, v4_markers_vs_clust4.R)
 	- figureED5_R.ipynb\
                - (first run save_seurat_objects.R, hipA7_only.R, 6day_only.R, 6day_metG_only.R, ds30_seurat.R)
 	- figureED6_R.ipynb\
                - (first run save_seurat_objects.R, seurat_ds38_min10_CFTonly_v4.R)
 	- figureED6_py.ipynb\
                - (first run v4_CFT_marker.R, figED2_markers_vs_clust5.R)
 	- figureED7_R.ipynb\
                - (first run proteomics_seurat_v3.R)
 	- figureED7_py.ipynb\
                - (first run fig1_markers.R)
 	- figureED8_py.ipynb
 	- figureED8_R.ipynb
 	- figureED9_py.ipynb\
                - (first run fig3_persister_only.R)
 	- figureED9_R.ipynb\
                - (first run fig3_persister_only.R)
 	- figureED10_R.ipynb\
                - (first run save_seurat_objects.R,fig3_persister_only.R)
 	- figureED10_py.ipynb
 	- figureED11_py.ipynb
 	- figureED12_py.ipynb

-	Additional scripts (*.R) to generate needed files for notebooks
-	Source data (‘source_data/included/’)
-	enrichment_and_p_val.py: example script for CRISPRi processing
-	PETRI_seq_scripts_v2
 	 - Updated scripts and demo for processing raw PETRI-seq data to generate count matrix; see README.txt in this folder
