# R 4.1.0
library(Seurat) #v4.1.1
library(EnvStats) #v2.7.0
library(limma) # 1.5-0
library(pbapply) # 3.50.3
library(edgeR)
library(stringr)

seurat <- readRDS('source_data/from_GEO/ds_seurat_PCA_UMAP_clusters_ds38_min10_v4.rds')
seurat@active.ident <- seurat$seurat_clusters
seurat <- subset(seurat,idents=2)

seurat@active.ident <- seurat$cell.type
WT_d6 <- c('d6_1hr_1_20220422','d6_1hr_2_20220422')
hipA7 <- c('h1_42min_20211217','h2_42min_20211217','h3_50min_20211217')
metG <- c('metG_52min_20211217','metG1_20200926','metG2_20200926','metG5_20201018_1','metG6_20201018_1','t4_2hr_metG_20210623','metG_pre_20231016')
tet <- c('tet1_20211014','tet5_20211014','tet1hr_20230919')
to_include = c(hipA7,WT_d6,tet,metG)
seurat <- subset(seurat,idents=to_include)

seurat$cell.group <- seurat$cell.type
levels(seurat$cell.group) <- c(levels(seurat$cell.group),'hipA7','WT_d6','tet','metG')
seurat$cell.group[seurat$cell.group %in% WT_d6] <- 'WT_d6'
seurat$cell.group[seurat$cell.group %in% hipA7] <- 'hipA7'
seurat$cell.group[seurat$cell.group %in% metG] <- 'metG'
seurat$cell.group[seurat$cell.group %in% tet] <- 'tet'
seurat$cell.group <- droplevels(seurat$cell.group)
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

saveRDS(seurat, file = paste('source_data/generated/fig3_persister_only.rds',sep=''))

seurat_full <- readRDS('source_data/generated/seurat_full.rds')
## find and save marker genes
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

seurat@active.ident <- seurat$cell.group
for (id1 in c('WT_d6','hipA7','metG')){
	id2 = 'tet'
	cells.2 = names(subset(seurat,idents=id2)$cell.type)
        cells.1 = names(subset(seurat,idents=id1)$cell.type)
	save_markers(seurat,seurat_full,cells.1,cells.2,paste('source_data/generated/',id1,'_vs_',id2,'_markers_for_fig3.txt',sep=''))
}


