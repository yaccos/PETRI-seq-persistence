library('Seurat') #v4.1.1
library('limma') # 1.5-0
library('pbapply') # 3.50.3
library('edgeR')
library('EnvStats') #v2.7.0
library('stringr')

save_markers <- function(seurat,seurat_full,cells.1,cells.2,filename){
	length(cells.1)
	length(cells.2)
	data.use = seurat$RNA@data
	data.use <- data.use[, c(cells.1, cells.2), drop = FALSE]
	j <- seq_len(length.out = length(x = cells.1))
	pvals <- pbsapply(
	        X = 1:nrow(x = data.use),
	        FUN = function(x) {
	                return(2 * limma::rankSumTestWithCorrelation(index = j, statistics = data.use[x, ]))
	        }
	    )
	colnames(pvals) <- rownames(data.use)
	bulk.1 = rowSums(subset(seurat_full,cells=cells.1)$RNA@counts)
	bulk.2 = rowSums(subset(seurat_full,cells=cells.2)$RNA@counts)
	bulk.df = t(rbind(bulk.1,bulk.2))
	write.table(bulk.df,str_replace(filename,'.txt','_counts.txt'),sep='\t',quote=FALSE)
 	y <- DGEList(counts=bulk.df, genes=rownames(bulk.df), group=c(0,1))
 	countsPerMillion <- cpm(y)
 	countCheck <- countsPerMillion > 5
 	keep <- which(rowSums(countCheck) >= 1)
 	y <- y[keep,]
 	countsPerMillion <- as.data.frame(countsPerMillion)
 	countsPerMillion$operon <- rownames(countsPerMillion)
 	y <- calcNormFactors(y, method="TMM")
 	design <- model.matrix(~0+group, data=y$samples)
 	colnames(design) <- levels(y$samples$group)
 	bcv <- 0.1
 	et <- exactTest(y, pair=c('1','0'),dispersion = bcv^2)
 	result <- cbind(operon=rownames(et$table), et$table ) #add rownames as a column
 	result = merge(result,t(pvals),by.x=0,by.y=0)
 	write.table(result,filename,sep='\t',quote=FALSE)
}
seurat <- readRDS('source_data/generated/onlyCFT_ds_seurat_PCA_UMAP_clusters_v4.rds')
seurat@active.ident <- seurat$cell.type
seurat_full <- readRDS('source_data/from_GEO/seurat_full_with_CFT.rds')
seurat_full <- subset(seurat_full,cells=names(seurat$cell.type))
counts <- GetAssayData(seurat_full, assay = "RNA")
counts <- counts[-(which(rownames(counts) %in% c('rRNA'))),]
seurat_full <- subset(seurat_full, features = rownames(counts))
seurat_full <- subset(seurat_full, subset = nCount_RNA > 9 & nCount_RNA < 1000)
seurat_full$cell.type <- seurat$cell.type
seurat_full$seurat_clusters <- seurat$seurat_clusters

clustid = names(which.max(table(subset(seurat,idents='CFT073_amp_20231016')$seurat_clusters)))
ids = c('CFT073_amp_20231016')
cells.1 = names(subset(seurat,subset=(cell.type %in% ids) & (seurat_clusters==clustid))$cell.type)
ids = c('CFT073_1hr_20231016')
cells.2 = names(subset(seurat,subset=(cell.type %in% ids) & (seurat_clusters==clustid))$cell.type)
save_markers(seurat,seurat_full,cells.1,cells.2,'source_data/generated/v4_CFTamp_vs_CFT1hr.txt')

