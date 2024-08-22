# R 4.1.0
library('Seurat') #v4.1.1
library('limma') # 1.5-0
library('pbapply') # 3.50.3
library('edgeR') # 3.36.0
library('EnvStats') #v2.7.0
library('stringr') #1.4.0

seurat <- readRDS('source_data/from_GEO/ds_seurat_PCA_UMAP_clusters_ds38_min10_v4.rds') # download from GEO; also saved in seurat_ds38_min10_v4.R
seurat_full <- readRDS('source_data/generated/seurat_full.rds') # download full from GEO then run save_seurat_objects.R to remove CFT
counts <- GetAssayData(seurat_full, assay = "RNA")
counts <- counts[-(which(rownames(counts) %in% c('rRNA'))),]
seurat_full <- subset(seurat_full, features = rownames(counts))
seurat_full <- subset(seurat_full, subset = nCount_RNA > 9 & nCount_RNA < 1000)
seurat_full$cell.type <- seurat$cell.type
seurat_full$seurat_clusters <- seurat$seurat_clusters

to_include = c('t1_2hr_WT_20210722','t1_20min_WT_20210722','t1_undil_WT_20210722',
               't2_20min_WT_20210722','t2_undil_WT_20210722','t3_undil_WT_20210722',
               't4_2hr_WT_20210722','t4_1hr_WT_20210722','t4_20min_WT_20210722','t4_undil_WT_20210722',
               't1_2hr_metG_20210623','t1_20min_metG_20210623','t1_undil_metG_20210623',
               't2_2hr_metG_20210623','t2_20min_metG_20210623','t2_undil_metG_20210623','t3_undil_metG_20210623',
               't4_2hr_metG_20210623','t4_20min_metG_20210623','t4_undil_metG_20210623',
               'WT_1hr_20220331','WT_30min_20220331','WT_10min_20220331','WT_3min_20220331','WT_stat_20220331',
               'metG6_20201018_1','metG5_20201018_1','metG4_20201018_1','metG3_20201018_1','metG2_20200926','metG1_20200926','WT_10min_20220210')

seurat@active.ident <- seurat$cell.type
seurat <- subset(seurat,idents=to_include)
seurat@active.ident <- seurat$seurat_clusters

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

for (id1 in c('2','3')){
	for (id2 in c('1','2','3','4','5')){
		if (id1!=id2){
			cells.1 = names(subset(seurat,idents=id1)$cell.type)
			cells.2 = names(subset(seurat,idents=id2)$cell.type)
			if (id2=='5'){
				cells.2 = na.omit(names(subset(seurat,idents=id2)$cell.type)[1:50000])
			}
			filename = paste('source_data/generated/v2_clust',id1,'_vs_',id2,'_markers_for_fig1.txt',sep='')
			save_markers(seurat,seurat_full,cells.1,cells.2,filename)
		}
	}
}

