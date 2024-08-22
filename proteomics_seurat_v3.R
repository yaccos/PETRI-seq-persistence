# R version 4.1.0 #
# Load packages: Seurat_4.1.1, RColorBrewer_1.1-3, ggplot2_3.3.6, extrafont_0.18,stringr_1.4.0 
packages = c('Seurat','RColorBrewer','ggplot2','extrafont','stringr','EnvStats')

for (p in packages){
    suppressPackageStartupMessages(eval(bquote(library(.(p)))))
}
seurat <- readRDS('source_data/from_GEO/ds_seurat_PCA_UMAP_clusters_ds38_min10_v4.rds')
orig.embeddings = seurat[['pca']]@cell.embeddings
geo.mean <- geoMean(colSums(data.frame(seurat@assays$RNA@counts)))
all.genes <- rownames(seurat)

data_table <- read.delim('source_data/included/SB1354_Fusion_MaxQuant_LFQ_intensity_wt_metG.txt',row.names=1)
seurat_p <- CreateSeuratObject(counts = data_table, project = "hvgs")  
seurat_p <- NormalizeData(seurat_p, normalization.method = "LogNormalize", scale.factor = geo.mean)

seurat_p <- merge(seurat,seurat_p) 

seurat_p$protein <- factor(grepl('LFQ',names(seurat_p$cell.type)))
seurat_p@active.ident <- seurat_p$protein
counts <- GetAssayData(subset(seurat_p,idents=TRUE), assay = "RNA")
counts <- (counts[rowSums(counts)>0,])
seurat_p <- subset(seurat_p, features = rownames(counts))
seurat_p <- ScaleData(seurat_p,features=all.genes)
loadings <- seurat[['pca']]@feature.loadings
reference <- subset(seurat_p,idents=FALSE)
o <- subset(seurat_p,idents=TRUE)
rm(seurat)
rm(seurat_p)

new.names <- c()
for (val in names(o@active.ident))
{
    temp.name = unlist(strsplit(val,'_1'))[1]
    temp.name = unlist(strsplit(temp.name,'_2'))[1]
    temp.name = unlist(strsplit(temp.name,'_3'))[1]
    new.names <- c(new.names,paste(temp.name,sep=''))
}
o$cell.type <- factor(new.names)
o@active.ident <- o$cell.type

pca.features <- rownames(counts[rowSums(counts)>0,])
sub_loadings <- loadings[rownames(loadings) %in% pca.features,]
pca.features <- rownames(sub_loadings)

o.pca.embeddings <- t(GetAssayData(object = o, assay = 'RNA',slot = "scale.data")[pca.features,]) %*% sub_loadings
o[['pca']] <- CreateDimReducObject(embeddings = o.pca.embeddings, key = "PC_", assay = "RNA"  )

ref.pca.embeddings <- t(GetAssayData(object = reference, assay = 'RNA',slot = "scale.data")[pca.features,]) %*% sub_loadings
reference[['pca']] <- CreateDimReducObject(embeddings = ref.pca.embeddings, key = "PC_", assay = "RNA"  )

merged = merge(reference,o,merge.dr = c("pca"))
merged <- RunUMAP(merged, dims = 1:6,return.model=TRUE,n.neighbors=50)
merged <- FindNeighbors(merged, dims = 1:6,k.param=50)
merged <- FindClusters(merged, resolution = 0.32)

seurat <- merged
rm(merged)
seurat$umap@cell.embeddings[,'UMAP_2'] = -seurat$umap@cell.embeddings[,'UMAP_2'] ## flip UMAP2

#### reorder cluster numbers ##
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

saveRDS(seurat,'source_data/generated/v3_proteomics_ds_seurat_PCA_UMAP_clusters_ds38_min10_v4.rds')
