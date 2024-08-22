library(Seurat)

seurat <- readRDS('source_data/from_GEO/seurat_full_with_CFT.rds') ## from GEO
saveRDS(subset(seurat,idents=c('CFT073_1hr_20231016','CFT073_amp_20231016')),'source_data/generated/seurat_CFT.rds') ## save just CFT cells
saveRDS(subset(seurat,invert=TRUE,idents=c('CFT073_1hr_20231016','CFT073_amp_20231016')),'source_data/generated/seurat_full.rds') ## remove the CFT cells
