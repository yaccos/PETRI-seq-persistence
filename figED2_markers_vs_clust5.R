library('Seurat') #v4.1.1
library('limma') # 1.5-0
library('pbapply') # 3.50.3
library('edgeR')
library('EnvStats') #v2.7.0
library('stringr')
seurat <- readRDS('source_data/from_GEO/ds_seurat_PCA_UMAP_clusters_ds38_min10_v4.rds') # download from GEO; also saved in seurat_ds38_min10_v4.R
seurat_full <- readRDS('source_data/generated/seurat_full.rds')
counts <- GetAssayData(seurat_full, assay = "RNA")
counts <- counts[-(which(rownames(counts) %in% c('rRNA'))),]
seurat_full <- subset(seurat_full, features = rownames(counts))
seurat_full <- subset(seurat_full, subset = nCount_RNA > 9 & nCount_RNA < 1000)
seurat_full$cell.type <- seurat$cell.type
seurat_full$seurat_clusters <- seurat$seurat_clusters

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

cells.2 = sample(names(subset(seurat,idents='5')$cell.type),40000)

cells.1 = names(subset(seurat,subset=(seurat_clusters=='1'))$cell.type)
save_markers(seurat,seurat_full,cells.1,cells.2,'source_data/generated/v4_clust1_vs_5.txt')

ids = c('tet1_20211014','tet5_20211014','tet1hr_20230919')
cells.1 = names(subset(seurat,subset=(cell.type %in% ids))$cell.type)
save_markers(seurat,seurat_full,cells.1,cells.2,'source_data/generated/v4_tet_vs_5.txt')

ids = levels(seurat$cell.type)[(grepl('metG',levels(seurat$cell.type))) & !(grepl('metG3|undil|amp|abx',levels(seurat$cell.type)))]
cells.1 = names(subset(seurat,subset=(cell.type %in% ids) & (seurat_clusters=='2'))$cell.type)
save_markers(seurat,seurat_full,cells.1,cells.2,'source_data/generated/v4_metG_clust2_vs_5.txt')

ids = levels(seurat$cell.type)[(grepl('d6_1hr',levels(seurat$cell.type)))]
cells.1 = names(subset(seurat,subset=(cell.type %in% ids) & (seurat_clusters=='2'))$cell.type)
save_markers(seurat,seurat_full,cells.1,cells.2,'source_data/generated/v4_d6_clust2_vs_5.txt')

ids = levels(seurat$cell.type)[(grepl('h1|h2|h3',levels(seurat$cell.type))) & !(grepl('10min|20min',levels(seurat$cell.type)))]
cells.1 = names(subset(seurat,subset=(cell.type %in% ids) & (seurat_clusters=='2'))$cell.type)
save_markers(seurat,seurat_full,cells.1,cells.2,'source_data/generated/v4_hipA7_clust2_vs_5.txt')

ids = c('metG_amp_20231016','metG_imm_abx_20230919')
cells.1 = names(subset(seurat,subset=(cell.type %in% ids))$cell.type)
save_markers(seurat,seurat_full,cells.1,cells.2,'source_data/generated/v4_metG_amp_vs_5.txt')

ids = c('wt_amp_20231016')
cells.1 = names(subset(seurat,subset=(cell.type %in% ids))$cell.type)
save_markers(seurat,seurat_full,cells.1,cells.2,'source_data/generated/v4_wt_amp_vs_5.txt')

ids = c('d6_amp_3_20220422')
cells.1 = names(subset(seurat,subset=(cell.type %in% ids))$cell.type)
save_markers(seurat,seurat_full,cells.1,cells.2,'source_data/generated/v4_d6_amp_vs_5.txt')
