library(Seurat) #v4.1.1
library(umap) #v0.2.8.0
library(EnvStats) #v2.7.0

file_name = 'source_data/generated/onlyCFT_ds_seurat_PCA_UMAP_clusters_v4.rds' 
seurat_full <- readRDS('source_data/generated/seurat_CFT.rds') # from GEO
counts <- GetAssayData(seurat_full, assay = "RNA")
counts <- counts[-(which(rownames(counts) %in% c('rRNA'))),]
seurat_full <- subset(seurat_full, features = rownames(counts))

#### downsample and filter ###
set.seed(30)
seurat <- subset(seurat_full, subset = nCount_RNA > 9 & nCount_RNA < 1000)
cell.type <- seurat$cell.type
counts = as.matrix(x = GetAssayData(object = seurat, assay = "RNA", slot = "counts"))
downsampled = SampleUMI(data = counts,max.umi = 38, upsample = FALSE, verbose = TRUE)
seurat = CreateSeuratObject(counts = downsampled, project = "hvgs")
seurat$cell.type = cell.type
all.genes <- rownames(seurat)


#### normalize and scale ###
geo.mean <- geoMean(colSums(data.frame(seurat@assays$RNA@counts)))
seurat <- NormalizeData(seurat, normalization.method = "LogNormalize", scale.factor = geo.mean)
seurat <- ScaleData(seurat,features=all.genes)
seurat <- FindVariableFeatures(seurat, selection.method = "vst")
seurat <- RunPCA(seurat, verbose = TRUE,features=VariableFeatures(seurat))
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
seurat <- FindClusters(seurat, resolution = 0.31, random.seed=2)

saveRDS(seurat, file = file_name)
