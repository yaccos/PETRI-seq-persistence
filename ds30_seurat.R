library(Seurat) #v4.1.1
library(umap) #v0.2.8.0
library(EnvStats) #v2.7.0

set.seed(5)
seurat <- readRDS('source_data/from_GEO/ds_seurat_PCA_UMAP_clusters_ds38_min10_v4.rds')
cell.type <- seurat$cell.type

#### downsample and filter ###
counts = as.matrix(x = GetAssayData(object = seurat, assay = "RNA", slot = "counts"))
downsampled = SampleUMI(data = counts,max.umi = 30, upsample = FALSE, verbose = TRUE)
seurat = CreateSeuratObject(counts = downsampled, project = "hvgs")
seurat <- subset(seurat, subset = nCount_RNA > 28 & nCount_RNA < 31)

## save cell type slot ##
seurat$cell.type = cell.type
all.genes <- rownames(seurat)

#### normalize and scale ###
geo.mean <- geoMean(colSums(data.frame(seurat@assays$RNA@counts)))
seurat <- NormalizeData(seurat, normalization.method = "LogNormalize", scale.factor = geo.mean)
seurat <- ScaleData(seurat,features=all.genes)
seurat <- FindVariableFeatures(seurat, selection.method = "vst")
seurat <- RunPCA(seurat, verbose = TRUE,features=all.genes)
seurat <- DietSeurat(
              seurat,
              counts = TRUE,
              data = TRUE,
              scale.data = FALSE,
              assays = 'RNA',
              dimreducs = c('pca'),
              graphs = NULL
            )
seurat <- RunUMAP(seurat, dims = 1:10)
seurat <- FindNeighbors(seurat, dims = 1:10)
seurat <- FindClusters(seurat, resolution = 0.34)

### reorder cluster numbers ##

levels(seurat$seurat_clusters) <- c(levels(seurat$seurat_clusters),10,11,12,13,14,15,16)
seurat$seurat_clusters[seurat$seurat_clusters==0] = 15
seurat$seurat_clusters[seurat$seurat_clusters==1] = 14
seurat$seurat_clusters[seurat$seurat_clusters==2] = 16
seurat$seurat_clusters[seurat$seurat_clusters==3] = 12
seurat$seurat_clusters[seurat$seurat_clusters==4] = 10
seurat$seurat_clusters[seurat$seurat_clusters==5] = 11
seurat$seurat_clusters[seurat$seurat_clusters==6] = 13
seurat$seurat_clusters = droplevels(seurat$seurat_clusters)
seurat$seurat_clusters <- factor(as.numeric(seurat$seurat_clusters)-1) ## revert back to numbers 0-6

seurat$umap@cell.embeddings[,'UMAP_1'] = -seurat$umap@cell.embeddings[,'UMAP_1'] ## flip UMAP1

saveRDS(seurat, file = 'source_data/generated/ds30_seurat_PCA_UMAP_clusters.rds')

