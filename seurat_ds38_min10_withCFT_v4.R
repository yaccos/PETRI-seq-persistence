library(Seurat) #v4.1.1
library(umap) #v0.2.8.0
library(EnvStats) #v2.7.0

file_name = 'source_data/from_GEO/withCFT_ds_seurat_PCA_UMAP_clusters_v4.rds'
seurat_full <- readRDS('source_data/from_GEO/seurat_full_with_CFT.rds')
counts <- GetAssayData(seurat_full, assay = "RNA")
counts <- counts[-(which(rownames(counts) %in% c('rRNA'))),]

CFT_to_MG_df = read.delim('CFT_to_MG_gene_names.txt',sep='\t')
genes = CFT_to_MG_df[(CFT_to_MG_df['MG_name']!='none') & (CFT_to_MG_df['MG_name']!='not_annotated'),'MG_name']
counts <- counts[(which(rownames(counts) %in% genes)),]
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

saveRDS(seurat, file = file_name)

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
saveRDS(seurat, file = paste(file_name,sep=''))

seurat <- RunUMAP(seurat, dims = 1:10)
seurat <- FindNeighbors(seurat, dims = 1:10)
seurat <- FindClusters(seurat, resolution = 0.38, random.seed=2)

seurat$umap@cell.embeddings[,'UMAP_1'] = -seurat$umap@cell.embeddings[,'UMAP_1'] ## flip UMAP1
### reorder cluster numbers ##

seurat$pca@cell.embeddings[,'PC_1'] = -seurat$pca@cell.embeddings[,'PC_1'] ## flip PC1
seurat$pca@feature.loadings[,'PC_1'] = -seurat$pca@feature.loadings[,'PC_1'] ## flip PC1
seurat$pca@cell.embeddings[,'PC_2'] = -seurat$pca@cell.embeddings[,'PC_2'] ## flip PC2
seurat$pca@feature.loadings[,'PC_2'] = -seurat$pca@feature.loadings[,'PC_2'] ## flip PC2


levels(seurat$seurat_clusters) <- c(levels(seurat$seurat_clusters),10,11,12,13,14,15,16)
seurat$seurat_clusters[seurat$seurat_clusters==0] = 15
seurat$seurat_clusters[seurat$seurat_clusters==1] = 14
seurat$seurat_clusters[seurat$seurat_clusters==2] = 16
seurat$seurat_clusters[seurat$seurat_clusters==3] = 12
seurat$seurat_clusters[seurat$seurat_clusters==4] = 11
seurat$seurat_clusters[seurat$seurat_clusters==5] = 13
seurat$seurat_clusters[seurat$seurat_clusters==6] = 10
seurat$seurat_clusters = droplevels(seurat$seurat_clusters)
seurat$seurat_clusters <- factor(as.numeric(seurat$seurat_clusters)-1) ## revert back to numbers 0-6
seurat@active.ident <- seurat$seurat_clusters

saveRDS(seurat, file = file_name)
