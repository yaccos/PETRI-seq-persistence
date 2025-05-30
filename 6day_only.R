library(Seurat) #v4.1.1
library(umap) #v0.2.8.0
library(EnvStats) #v2.7.0

seurat_file = 'source_data/from_GEO/ds_seurat_PCA_UMAP_clusters_ds38_min10_v4.rds'
seurat <- readRDS(seurat_file)
seurat@active.ident <- seurat$cell.type
seurat <- subset(seurat,idents=c('d6_1hr_1_20220422','d6_1hr_2_20220422'))

all.genes <- rownames(seurat)
geo.mean <- geoMean(colSums(data.frame(seurat@assays$RNA@counts)))
seurat <- NormalizeData(seurat, normalization.method = "LogNormalize", scale.factor = geo.mean)
seurat <- ScaleData(seurat,features=all.genes)
seurat <- FindVariableFeatures(seurat, selection.method = "vst")
seurat <- RunPCA(seurat, verbose = TRUE,features=all.genes)

seurat <- FindNeighbors(seurat, dims = 1:6)
seurat <- FindClusters(seurat, resolution = 0.3)

seurat <- DietSeurat(
              seurat,
              counts = TRUE,
              data = TRUE,
              scale.data = FALSE,
              assays = 'RNA',
              dimreducs = c('pca'),
              graphs = NULL
            )
saveRDS(seurat,'source_data/generated/day6_only_seurat.rds')
