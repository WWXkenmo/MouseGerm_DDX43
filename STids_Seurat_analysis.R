##########################################################
##Integrate three population
SeuratMerge <- function(obj1, obj2){
   require("Seurat")
   cat('integrate assay...')
   assay1 <- GetAssayData(obj1)
   assay2 <- GetAssayData(obj2)
   
   cat('\nintegrate pca...')
   pca1 <- obj1[['pca']]@cell.embeddings
   pca2 <- obj2[['pca']]@cell.embeddings
   pca <- rbind(pca1,pca2)
   
   cat('\nintegrate meta.data...')
   meta.data1 <- obj1@meta.data
   meta.data2 <- obj2@meta.data
   meta.data <- rbind(meta.data1,meta.data2)
   
   cat('\nDone')
   assay <- cbind(assay1,assay2)
   assay <- CreateSeuratObject(assay)
   assay@meta.data <- meta.data
   assay[['pca']] <- CreateDimReducObject(pca)
   assay
}

SeuratAssay <- function(obj1, obj2,obj3,assay){
   require("Seurat")
   cat('integrate assay...')
   assay1 <- GetAssayData(obj1,assay=assay)
   assay2 <- GetAssayData(obj2,assay=assay)
   assay3 <- GetAssayData(obj3,assay=assay)
   
   geneID <- intersect(rownames(assay1),rownames(assay2))
   geneID <- intersect(geneID,rownames(assay3))
   assay <- cbind(assay1[geneID,],assay2[geneID,],assay3[geneID,])
   assay
}


wilcoxTest <- function(x,label){
    fc <- mean(x[label==1])-mean(x[label==0])
    pvalue <- wilcox.test(x[label==1],x[label==0])$p.value
	line <- c(fc,pvalue)
	line
}

KivsWt <- function(mat, batch,ncore=6){
    SPG_Ki_rep1 <- mat[,batch=="MouseSperm366"]
    SPG_Ki_rep2 <- mat[,batch=="MouseSperm368"]
    SPG_Wt_rep1 <- mat[,batch=="MouseSperm499"]
    SPG_Wt_rep2 <- mat[,batch=="MouseSperm500"]
    
    cat("Processing 1th")
	cat("\nRunning the first group...")
	test1 <- cbind(SPG_Ki_rep1, SPG_Wt_rep1)
	label <- c(rep(1,ncol(SPG_Ki_rep1)),rep(0,ncol(SPG_Wt_rep1)))
	library(parallel)
	cat(paste("\nUsing ",ncore," cores...",sep=""))
    cl <- makeCluster(ncore)
    clusterExport(cl, 'label')
    clusterExport(cl,'test1')
	clusterExport(cl,'wilcoxTest')
    dt <- parApply(cl,X=test1,1,wilcoxTest,label=label)
	stopCluster(cl)
	cat("\nDone")
	dt1 <- t(dt)

	cat("\nRunning the second group...")
	test1 <- cbind(SPG_Ki_rep1, SPG_Wt_rep2)
	label <- c(rep(1,ncol(SPG_Ki_rep1)),rep(0,ncol(SPG_Wt_rep2)))
	library(parallel)
	cat(paste("\nUsing ",ncore," cores...",sep=""))
    cl <- makeCluster(ncore)
    clusterExport(cl, 'label')
    clusterExport(cl,'test1')
	clusterExport(cl,'wilcoxTest')
    dt <- parApply(cl,X=test1,1,wilcoxTest,label=label)
	stopCluster(cl)
	cat("\nDone")
	dt2 <- t(dt)
	
	cat("\nRunning the third group...")
	test1 <- cbind(SPG_Ki_rep2, SPG_Wt_rep1)
	label <- c(rep(1,ncol(SPG_Ki_rep2)),rep(0,ncol(SPG_Wt_rep1)))
	library(parallel)
	cat(paste("\nUsing ",ncore," cores...",sep=""))
    cl <- makeCluster(ncore)
    clusterExport(cl, 'label')
    clusterExport(cl,'test1')
	clusterExport(cl,'wilcoxTest')
    dt <- parApply(cl,X=test1,1,wilcoxTest,label=label)
	stopCluster(cl)
	cat("\nDone")
	dt3 <- t(dt)
	
	
	cat("\nRunning the fourth group...")
	test1 <- cbind(SPG_Ki_rep2, SPG_Wt_rep2)
	label <- c(rep(1,ncol(SPG_Ki_rep2)),rep(0,ncol(SPG_Wt_rep2)))
	library(parallel)
	cat(paste("\nUsing ",ncore," cores...",sep=""))
    cl <- makeCluster(ncore)
    clusterExport(cl, 'label')
    clusterExport(cl,'test1')
	clusterExport(cl,'wilcoxTest')
    dt <- parApply(cl,X=test1,1,wilcoxTest,label=label)
	stopCluster(cl)
	cat("\nDone")
	dt4 <- t(dt)
	
	fc_1 <- dt1[,1]
	fc_2 <- dt2[,1]
	fc_3 <- dt3[,1]
	fc_4 <- dt4[,1]
	pvalue_1 <- dt1[,2]
	pvalue_2 <- dt2[,2]
	pvalue_3 <- dt3[,2]
	pvalue_4 <- dt4[,2]
    cat("\nProcessing 2th")
    fc <- sign(fc_1)+sign(fc_2)+sign(fc_3)+sign(fc_4)
    pvalue <- 1-(1-pvalue_1)*(1-pvalue_2)*(1-pvalue_3)*(1-pvalue_4)
    
    #initial filtering
    fc <- which(abs(fc)==4)
    pvalue <- which(pvalue<0.5)
    gene <- rownames(mat)[intersect(fc,pvalue)]
    cat("\nProcessing OverAll")
    
	cat("\nRunning the final group...")
	test1 <- cbind(SPG_Ki_rep1,SPG_Ki_rep2,SPG_Wt_rep1, SPG_Wt_rep2)[gene,]
	label <- c(rep(1,ncol(cbind(SPG_Ki_rep1,SPG_Ki_rep2))),rep(0,ncol(cbind(SPG_Wt_rep1, SPG_Wt_rep2))))
	library(parallel)
	cat(paste("\nUsing ",ncore," cores...",sep=""))
    cl <- makeCluster(ncore)
    clusterExport(cl, 'label')
    clusterExport(cl,'test1')
	clusterExport(cl,'wilcoxTest')
    dt <- parApply(cl,X=test1,1,wilcoxTest,label=label)
	stopCluster(cl)
	cat("\nDone")
	dt_final <- t(dt)
	pvalue <- dt_final[,2]
	
    fdr <- p.adjust(pvalue,method="BH")
    
    res <- cbind(dt_final,fdr)
    colnames(res) <- c("Diff","P-value","FDR")
    rownames(res) <- gene
    res
}

KivsWt_all <- function(mat, batch,ncore=6){
    SPG_Ki_rep1 <- mat[,batch=="MouseSperm366"]
    SPG_Ki_rep2 <- mat[,batch=="MouseSperm368"]
    SPG_Wt_rep1 <- mat[,batch=="MouseSperm499"]
    SPG_Wt_rep2 <- mat[,batch=="MouseSperm500"]
    
    cat("Processing 1th")
	cat("\nRunning the first group...")
	test1 <- cbind(SPG_Ki_rep1, SPG_Wt_rep1)
	label <- c(rep(1,ncol(SPG_Ki_rep1)),rep(0,ncol(SPG_Wt_rep1)))
	library(parallel)
	cat(paste("\nUsing ",ncore," cores...",sep=""))
    cl <- makeCluster(ncore)
    clusterExport(cl, 'label')
    clusterExport(cl,'test1')
	clusterExport(cl,'wilcoxTest')
    dt <- parApply(cl,X=test1,1,wilcoxTest,label=label)
	stopCluster(cl)
	cat("\nDone")
	dt1 <- t(dt)

	cat("\nRunning the second group...")
	test1 <- cbind(SPG_Ki_rep1, SPG_Wt_rep2)
	label <- c(rep(1,ncol(SPG_Ki_rep1)),rep(0,ncol(SPG_Wt_rep2)))
	library(parallel)
	cat(paste("\nUsing ",ncore," cores...",sep=""))
    cl <- makeCluster(ncore)
    clusterExport(cl, 'label')
    clusterExport(cl,'test1')
	clusterExport(cl,'wilcoxTest')
    dt <- parApply(cl,X=test1,1,wilcoxTest,label=label)
	stopCluster(cl)
	cat("\nDone")
	dt2 <- t(dt)
	
	cat("\nRunning the third group...")
	test1 <- cbind(SPG_Ki_rep2, SPG_Wt_rep1)
	label <- c(rep(1,ncol(SPG_Ki_rep2)),rep(0,ncol(SPG_Wt_rep1)))
	library(parallel)
	cat(paste("\nUsing ",ncore," cores...",sep=""))
    cl <- makeCluster(ncore)
    clusterExport(cl, 'label')
    clusterExport(cl,'test1')
	clusterExport(cl,'wilcoxTest')
    dt <- parApply(cl,X=test1,1,wilcoxTest,label=label)
	stopCluster(cl)
	cat("\nDone")
	dt3 <- t(dt)
	
	
	cat("\nRunning the fourth group...")
	test1 <- cbind(SPG_Ki_rep2, SPG_Wt_rep2)
	label <- c(rep(1,ncol(SPG_Ki_rep2)),rep(0,ncol(SPG_Wt_rep2)))
	library(parallel)
	cat(paste("\nUsing ",ncore," cores...",sep=""))
    cl <- makeCluster(ncore)
    clusterExport(cl, 'label')
    clusterExport(cl,'test1')
	clusterExport(cl,'wilcoxTest')
    dt <- parApply(cl,X=test1,1,wilcoxTest,label=label)
	stopCluster(cl)
	cat("\nDone")
	dt4 <- t(dt)
	
	fc_1 <- dt1[,1]
	fc_2 <- dt2[,1]
	fc_3 <- dt3[,1]
	fc_4 <- dt4[,1]
	pvalue_1 <- dt1[,2]
	pvalue_2 <- dt2[,2]
	pvalue_3 <- dt3[,2]
	pvalue_4 <- dt4[,2]
    cat("\nProcessing 2th")
    fc <- sign(fc_1)+sign(fc_2)+sign(fc_3)+sign(fc_4)
    pvalue <- 1-(1-pvalue_1)*(1-pvalue_2)*(1-pvalue_3)*(1-pvalue_4)
    
    #initial filtering
    fc <- which(abs(fc)>=0)
    pvalue <- which(pvalue<1)
    gene <- rownames(mat)[intersect(fc,pvalue)]
    cat("\nProcessing OverAll")
    
	cat("\nRunning the final group...")
	test1 <- cbind(SPG_Ki_rep1,SPG_Ki_rep2,SPG_Wt_rep1, SPG_Wt_rep2)[gene,]
	label <- c(rep(1,ncol(cbind(SPG_Ki_rep1,SPG_Ki_rep2))),rep(0,ncol(cbind(SPG_Wt_rep1, SPG_Wt_rep2))))
	library(parallel)
	cat(paste("\nUsing ",ncore," cores...",sep=""))
    cl <- makeCluster(ncore)
    clusterExport(cl, 'label')
    clusterExport(cl,'test1')
	clusterExport(cl,'wilcoxTest')
    dt <- parApply(cl,X=test1,1,wilcoxTest,label=label)
	stopCluster(cl)
	cat("\nDone")
	dt_final <- t(dt)
	pvalue <- dt_final[,2]
	
    fdr <- p.adjust(pvalue,method="BH")
    
    res <- cbind(dt_final,fdr)
    colnames(res) <- c("Diff","P-value","FDR")
    rownames(res) <- gene
    res
}


KovsWt <- function(mat, batch,ncore=6){
    SPG_Ko_rep1 <- mat[,batch=="MouseSperm340"]
    SPG_Ko_rep2 <- mat[,batch=="MouseSperm341"]
    SPG_Wt_rep1 <- mat[,batch=="MouseSperm499"]
    SPG_Wt_rep2 <- mat[,batch=="MouseSperm500"]
    
    cat("Processing 1th")
	cat("\nRunning the first group...")
	test1 <- cbind(SPG_Ko_rep1, SPG_Wt_rep1)
	label <- c(rep(1,ncol(SPG_Ko_rep1)),rep(0,ncol(SPG_Wt_rep1)))
	library(parallel)
	cat(paste("\nUsing ",ncore," cores...",sep=""))
    cl <- makeCluster(ncore)
    clusterExport(cl, 'label')
    clusterExport(cl,'test1')
	clusterExport(cl,'wilcoxTest')
    dt <- parApply(cl,X=test1,1,wilcoxTest,label=label)
	stopCluster(cl)
	cat("\nDone")
	dt1 <- t(dt)

	cat("\nRunning the second group...")
	test1 <- cbind(SPG_Ko_rep1, SPG_Wt_rep2)
	label <- c(rep(1,ncol(SPG_Ko_rep1)),rep(0,ncol(SPG_Wt_rep2)))
	library(parallel)
	cat(paste("\nUsing ",ncore," cores...",sep=""))
    cl <- makeCluster(ncore)
    clusterExport(cl, 'label')
    clusterExport(cl,'test1')
	clusterExport(cl,'wilcoxTest')
    dt <- parApply(cl,X=test1,1,wilcoxTest,label=label)
	stopCluster(cl)
	cat("\nDone")
	dt2 <- t(dt)
	
	cat("\nRunning the third group...")
	test1 <- cbind(SPG_Ko_rep2, SPG_Wt_rep1)
	label <- c(rep(1,ncol(SPG_Ko_rep2)),rep(0,ncol(SPG_Wt_rep1)))
	library(parallel)
	cat(paste("\nUsing ",ncore," cores...",sep=""))
    cl <- makeCluster(ncore)
    clusterExport(cl, 'label')
    clusterExport(cl,'test1')
	clusterExport(cl,'wilcoxTest')
    dt <- parApply(cl,X=test1,1,wilcoxTest,label=label)
	stopCluster(cl)
	cat("\nDone")
	dt3 <- t(dt)
	
	
	cat("\nRunning the fourth group...")
	test1 <- cbind(SPG_Ko_rep2, SPG_Wt_rep2)
	label <- c(rep(1,ncol(SPG_Ko_rep2)),rep(0,ncol(SPG_Wt_rep2)))
	library(parallel)
	cat(paste("\nUsing ",ncore," cores...",sep=""))
    cl <- makeCluster(ncore)
    clusterExport(cl, 'label')
    clusterExport(cl,'test1')
	clusterExport(cl,'wilcoxTest')
    dt <- parApply(cl,X=test1,1,wilcoxTest,label=label)
	stopCluster(cl)
	cat("\nDone")
	dt4 <- t(dt)
	
	fc_1 <- dt1[,1]
	fc_2 <- dt2[,1]
	fc_3 <- dt3[,1]
	fc_4 <- dt4[,1]
	pvalue_1 <- dt1[,2]
	pvalue_2 <- dt2[,2]
	pvalue_3 <- dt3[,2]
	pvalue_4 <- dt4[,2]
    cat("\nProcessing 2th")
    fc <- sign(fc_1)+sign(fc_2)+sign(fc_3)+sign(fc_4)
    pvalue <- 1-(1-pvalue_1)*(1-pvalue_2)*(1-pvalue_3)*(1-pvalue_4)
    
    #initial filtering
    fc <- which(abs(fc)>=0)
    pvalue <- which(pvalue<1)
    gene <- rownames(mat)[intersect(fc,pvalue)]
    cat("\nProcessing OverAll")
    
	cat("\nRunning the final group...")
	test1 <- cbind(SPG_Ko_rep1,SPG_Ko_rep2,SPG_Wt_rep1, SPG_Wt_rep2)[gene,]
	label <- c(rep(1,ncol(cbind(SPG_Ko_rep1,SPG_Ko_rep2))),rep(0,ncol(cbind(SPG_Wt_rep1, SPG_Wt_rep2))))
	library(parallel)
	cat(paste("\nUsing ",ncore," cores...",sep=""))
    cl <- makeCluster(ncore)
    clusterExport(cl, 'label')
    clusterExport(cl,'test1')
	clusterExport(cl,'wilcoxTest')
    dt <- parApply(cl,X=test1,1,wilcoxTest,label=label)
	stopCluster(cl)
	cat("\nDone")
	dt_final <- t(dt)
	pvalue <- dt_final[,2]
	
    fdr <- p.adjust(pvalue,method="BH")
    
    res <- cbind(dt_final,fdr)
    colnames(res) <- c("Diff","P-value","FDR")
    rownames(res) <- gene
    res
}

KovsWt_all <- function(mat, batch,ncore=6){
    SPG_Ko_rep1 <- mat[,batch=="MouseSperm340"]
    SPG_Ko_rep2 <- mat[,batch=="MouseSperm341"]
    SPG_Wt_rep1 <- mat[,batch=="MouseSperm499"]
    SPG_Wt_rep2 <- mat[,batch=="MouseSperm500"]
    
    cat("Processing 1th")
	cat("\nRunning the first group...")
	test1 <- cbind(SPG_Ko_rep1, SPG_Wt_rep1)
	label <- c(rep(1,ncol(SPG_Ko_rep1)),rep(0,ncol(SPG_Wt_rep1)))
	library(parallel)
	cat(paste("\nUsing ",ncore," cores...",sep=""))
    cl <- makeCluster(ncore)
    clusterExport(cl, 'label')
    clusterExport(cl,'test1')
	clusterExport(cl,'wilcoxTest')
    dt <- parApply(cl,X=test1,1,wilcoxTest,label=label)
	stopCluster(cl)
	cat("\nDone")
	dt1 <- t(dt)

	cat("\nRunning the second group...")
	test1 <- cbind(SPG_Ko_rep1, SPG_Wt_rep2)
	label <- c(rep(1,ncol(SPG_Ko_rep1)),rep(0,ncol(SPG_Wt_rep2)))
	library(parallel)
	cat(paste("\nUsing ",ncore," cores...",sep=""))
    cl <- makeCluster(ncore)
    clusterExport(cl, 'label')
    clusterExport(cl,'test1')
	clusterExport(cl,'wilcoxTest')
    dt <- parApply(cl,X=test1,1,wilcoxTest,label=label)
	stopCluster(cl)
	cat("\nDone")
	dt2 <- t(dt)
	
	cat("\nRunning the third group...")
	test1 <- cbind(SPG_Ko_rep2, SPG_Wt_rep1)
	label <- c(rep(1,ncol(SPG_Ko_rep2)),rep(0,ncol(SPG_Wt_rep1)))
	library(parallel)
	cat(paste("\nUsing ",ncore," cores...",sep=""))
    cl <- makeCluster(ncore)
    clusterExport(cl, 'label')
    clusterExport(cl,'test1')
	clusterExport(cl,'wilcoxTest')
    dt <- parApply(cl,X=test1,1,wilcoxTest,label=label)
	stopCluster(cl)
	cat("\nDone")
	dt3 <- t(dt)
	
	
	cat("\nRunning the fourth group...")
	test1 <- cbind(SPG_Ko_rep2, SPG_Wt_rep2)
	label <- c(rep(1,ncol(SPG_Ko_rep2)),rep(0,ncol(SPG_Wt_rep2)))
	library(parallel)
	cat(paste("\nUsing ",ncore," cores...",sep=""))
    cl <- makeCluster(ncore)
    clusterExport(cl, 'label')
    clusterExport(cl,'test1')
	clusterExport(cl,'wilcoxTest')
    dt <- parApply(cl,X=test1,1,wilcoxTest,label=label)
	stopCluster(cl)
	cat("\nDone")
	dt4 <- t(dt)
	
	fc_1 <- dt1[,1]
	fc_2 <- dt2[,1]
	fc_3 <- dt3[,1]
	fc_4 <- dt4[,1]
	pvalue_1 <- dt1[,2]
	pvalue_2 <- dt2[,2]
	pvalue_3 <- dt3[,2]
	pvalue_4 <- dt4[,2]
    cat("\nProcessing 2th")
    fc <- sign(fc_1)+sign(fc_2)+sign(fc_3)+sign(fc_4)
    pvalue <- 1-(1-pvalue_1)*(1-pvalue_2)*(1-pvalue_3)*(1-pvalue_4)
    
    #initial filtering
    fc <- which(abs(fc)>=0)
    pvalue <- which(pvalue<1)
    gene <- rownames(mat)[intersect(fc,pvalue)]
    cat("\nProcessing OverAll")
    
	cat("\nRunning the final group...")
	test1 <- cbind(SPG_Ko_rep1,SPG_Ko_rep2,SPG_Wt_rep1, SPG_Wt_rep2)[gene,]
	label <- c(rep(1,ncol(cbind(SPG_Ko_rep1,SPG_Ko_rep2))),rep(0,ncol(cbind(SPG_Wt_rep1, SPG_Wt_rep2))))
	library(parallel)
	cat(paste("\nUsing ",ncore," cores...",sep=""))
    cl <- makeCluster(ncore)
    clusterExport(cl, 'label')
    clusterExport(cl,'test1')
	clusterExport(cl,'wilcoxTest')
    dt <- parApply(cl,X=test1,1,wilcoxTest,label=label)
	stopCluster(cl)
	cat("\nDone")
	dt_final <- t(dt)
	pvalue <- dt_final[,2]
	
    fdr <- p.adjust(pvalue,method="BH")
    
    res <- cbind(dt_final,fdr)
    colnames(res) <- c("Diff","P-value","FDR")
    rownames(res) <- gene
    res
}

SeuratCount <- function(obj1, obj2,obj3,celltype="STids1"){
    require("Seurat")
    cat('integrate count...')
    assay1 <- obj1@assays$RNA@counts[,obj1$celltype %in% celltype]
    assay2 <- obj2@assays$RNA@counts[,obj2$celltype %in% celltype]
    assay3 <- obj3@assays$RNA@counts[,obj3$celltype %in% celltype]
    
    geneID <- intersect(rownames(assay1),rownames(assay2))
    geneID <- intersect(geneID,rownames(assay3))
    assay <- cbind(assay1[geneID,],assay2[geneID,],assay3[geneID,])
    assay
}

KivsWt_all <- function(mat, batch,ncore=6){
    SPG_Ki_rep1 <- mat[,batch=="MouseSperm366"]
    SPG_Ki_rep2 <- mat[,batch=="MouseSperm368"]
    SPG_Wt_rep1 <- mat[,batch=="MouseSperm499"]
    SPG_Wt_rep2 <- mat[,batch=="MouseSperm500"]
    gene <- rownames(mat)
    cat("\nRunning the final group...")
    test1 <- cbind(SPG_Ki_rep1,SPG_Ki_rep2,SPG_Wt_rep1, SPG_Wt_rep2)[gene,]
    label <- c(rep(1,ncol(cbind(SPG_Ki_rep1,SPG_Ki_rep2))),rep(0,ncol(cbind(SPG_Wt_rep1, SPG_Wt_rep2))))
    library(parallel)
    cat(paste("\nUsing ",ncore," cores...",sep=""))
    cl <- makeCluster(ncore)
    clusterExport(cl, 'label')
    clusterExport(cl,'test1')
    clusterExport(cl,'wilcoxTest')
    dt <- parApply(cl,X=test1,1,wilcoxTest,label=label)
    stopCluster(cl)
    cat("\nDone")
    dt_final <- t(dt)
    pvalue <- dt_final[,2]
    
    fdr <- p.adjust(pvalue,method="BH")
    
    res <- cbind(dt_final,fdr)
    colnames(res) <- c("Diff","P-value","FDR")
    rownames(res) <- gene
    res}
################################
MouseGerm_Ki <- MouseGerm_Ki[,MouseGerm_Ki$celltype %in% c("Endothelial","Macrophage","Unknown") == FALSE]
MouseGerm_KO <- MouseGerm_KO[,MouseGerm_KO$celltype %in% c("Endothelial","Macrophage","Unknown") == FALSE]
MouseGerm_WT <- MouseGerm_WT[,MouseGerm_WT$celltype %in% c("Endothelial","Macrophage","Unknown") == FALSE]

MouseGerm.integrate <- SeuratMerge(MouseGerm_Ki, MouseGerm_KO)
MouseGerm.integrate <- SeuratMerge(MouseGerm.integrate, MouseGerm_WT)
##############################
#Select STids2
STids <- MouseGerm.integrate[,MouseGerm.integrate$celltype %in% "STids" == TRUE]
STids <- ScaleData(STids)
STids <- RunPCA(object = STids, verbose = FALSE,features=rownames(STids))
ElbowPlot(STids,ndims=30)
STids <- RunUMAP(object = STids, dims = 1:20)
STids <- FindNeighbors(STids, dims = 1:20,reduction="pca")
STids <- FindClusters(STids, resolution = 0.25,algorithm=2)
DimPlot(STids,reduction="umap",label=1)
##Find All Markers
STids.SCT <- SeuratAssay(MouseGerm_Ki, MouseGerm_KO,MouseGerm_WT,assay="SCT")
STids.SCT <- STids.SCT[,MouseGerm.integrate$celltype %in% "STids" == TRUE]
STids[['SCT']] <- CreateAssayObject(STids.SCT)
stids.markers <- FindAllMarkers(STids, only.pos = TRUE, test.use="bimod",min.diff.pct=0.2)

####identify the DEG of fail differentiation
STids.SCT <- GetAssayData(STids)
STids.SCT <- STids.SCT[,STids$seurat_clusters %in% 2 == TRUE]

STids.SCT.ki <- STids.SCT[,colnames(STids.SCT) %in% colnames(MouseGerm_Ki)]
STids.SCT.ko <- STids.SCT[,colnames(STids.SCT) %in% colnames(MouseGerm_KO)]
STids.SCT.wt <- STids.SCT[,colnames(STids.SCT) %in% colnames(MouseGerm_WT)]
batch.ki <- MouseGerm_Ki$batch[colnames(STids.SCT.ki)]
batch.ko <- MouseGerm_KO$batch[colnames(STids.SCT.ko)]
batch.wt <- MouseGerm_WT$batch[colnames(STids.SCT.wt)]


DEG_test_kivswt <- KivsWt(cbind(STids.SCT.ki,STids.SCT.wt), c(batch.ki,batch.wt))
DEG_test_kovswt <- KovsWt(cbind(STids.SCT.ko,STids.SCT.wt), c(batch.ko,batch.wt))

###############################################################
#possible confounded by sequencing depth, using read count with SCTransform to reproduce analysis
STids.new <- SeuratCount(MouseGerm_Ki, MouseGerm_KO,MouseGerm_WT)
STids.new <- CreateSeuratObject(STids.new)
STids.new <-SCTransform(STids.new,vars.to.regress = c("nFeature_RNA","nCount_RNA"))

STids_splice <- readRDS("/mnt/data2/wangweixu/mm10/spliced_count.rds")
STids_splice <- CreateSeuratObject(STids_splice$spliced_count)
STids_splice <-SCTransform(STids_splice,vars.to.regress = c("nFeature_RNA","nCount_RNA"))

STids.SCT <- as.matrix(GetAssayData(STids_splice)[,STids$seurat_clusters %in% 3 == TRUE])

STids.SCT.ki <- STids.SCT[,colnames(STids.SCT) %in% colnames(MouseGerm_Ki)]
STids.SCT.ko <- STids.SCT[,colnames(STids.SCT) %in% colnames(MouseGerm_KO)]
STids.SCT.wt <- STids.SCT[,colnames(STids.SCT) %in% colnames(MouseGerm_WT)]
batch.ki <- MouseGerm_Ki$batch[colnames(STids.SCT.ki)]
batch.ko <- MouseGerm_KO$batch[colnames(STids.SCT.ko)]
batch.wt <- MouseGerm_WT$batch[colnames(STids.SCT.wt)]

DEG_test_kivswt_3 <- KivsWt_all(cbind(STids.SCT.ki,STids.SCT.wt), c(batch.ki,batch.wt))
#DEG_test_kovswt_3 <- KovsWt_final(cbind(STids.SCT.ko,STids.SCT.wt), c(batch.ko,batch.wt))

DEG_test_kivswt_3 <- DEG_test_kivswt_3[!is.nan(DEG_test_kivswt_3[,2]),]

#DEG_test_kovswt_3 <- as.data.frame(DEG_test_kovswt_3)
DEG_test_kivswt_3 <- as.data.frame(DEG_test_kivswt_3)

#DEG_test_kovswt_3$logFC <- log2(rowMeans(STids.SCT.ko[rownames(DEG_test_kovswt_3),])/rowMeans(STids.SCT.wt[rownames(DEG_test_kovswt_3),]))
DEG_test_kivswt_3$logFC <- log2(rowMeans(STids.SCT.ki[rownames(DEG_test_kivswt_3),])/rowMeans(STids.SCT.wt[rownames(DEG_test_kivswt_3),]))
data <- DEG_test_kivswt_3

data$sig[(data$FDR > 0.05|data$FDR=="NA")|(data$logFC < 0)& data$logFC > -0] <- "no"
data$sig[data$FDR <= 0.05 & data$logFC >= 0] <- "up"
data$sig[data$FDR <= 0.05 & data$logFC <= 0] <- "down"
data <- na.omit(data)
logFC <- (data[,"logFC"])
x_lim <- max(logFC,-logFC)
library(ggplot2)
library(RColorBrewer)

theme_set(theme_bw())
p <- ggplot(data,aes(logFC,-1*log10(FDR),
                     color = sig))+geom_point()+
    xlim(-10,10) +  labs(x="log2(FoldChange)",y="-log10(FDR)")
p <- p + scale_color_manual(values =c("#0072B5","grey","#BC3C28"))+
    geom_hline(yintercept=-log10(0.05),linetype=4)+
    geom_vline(xintercept=c(-0,0),linetype=4)
p <- p +theme(panel.grid =element_blank())+
    theme(axis.line = element_line(size=0))+ylim(0,15)
p <- p  +guides(colour = FALSE)
p_ko_3 <- p +theme(axis.text=element_text(size=10),axis.title=element_text(size=10))





########################################
STids.SCT <- as.matrix(GetAssayData(STids_splice)[,STids$seurat_clusters %in% 3 == TRUE])

#STids.SCT.ki <- STids.SCT[,colnames(STids.SCT) %in% colnames(MouseGerm_Ki)]
STids.SCT.ko <- STids.SCT[,colnames(STids.SCT) %in% colnames(MouseGerm_KO)]
STids.SCT.wt <- STids.SCT[,colnames(STids.SCT) %in% colnames(MouseGerm_WT)]
#batch.ki <- MouseGerm_Ki$batch[colnames(STids.SCT.ki)]
batch.ko <- MouseGerm_KO$batch[colnames(STids.SCT.ko)]
batch.wt <- MouseGerm_WT$batch[colnames(STids.SCT.wt)]

#DEG_test_kivswt_3 <- KivsWt_all(cbind(STids.SCT.ki,STids.SCT.wt), c(batch.ki,batch.wt))
DEG_test_kovswt_3 <- KovsWt_all(cbind(STids.SCT.ko,STids.SCT.wt), c(batch.ko,batch.wt))
DEG_test_kovswt_3 <- as.data.frame(DEG_test_kovswt_3)
#DEG_test_kivswt_3 <- as.data.frame(DEG_test_kivswt_3)
###build vocalno plot
DEG_test_kovswt_3$logFC <- log2(rowMeans(STids.SCT.ko[rownames(DEG_test_kovswt_3),])/rowMeans(STids.SCT.wt[rownames(DEG_test_kovswt_3),]))
data <- DEG_test_kovswt_3

data$sig[(data$FDR > 0.05|data$FDR=="NA")|(data$logFC < 0)& data$logFC > -0] <- "no"
data$sig[data$FDR <= 0.05 & data$logFC >= 0] <- "up"
data$sig[data$FDR <= 0.05 & data$logFC <= 0] <- "down"
data <- na.omit(data)
logFC <- (data[,"logFC"])
x_lim <- max(logFC,-logFC)
library(ggplot2)
library(RColorBrewer)
theme_set(theme_bw())
p <- ggplot(data,aes(logFC,-1*log10(FDR),
                     color = sig))+geom_point()+
    xlim(-10,10) +  labs(x="log2(FoldChange)",y="-log10(FDR)")
p <- p + scale_color_manual(values =c("#0072B5","grey","#BC3C28"))+
    geom_hline(yintercept=-log10(0.05),linetype=4)+
    geom_vline(xintercept=c(-0,0),linetype=4)
p <- p +theme(panel.grid =element_blank())+
    theme(axis.line = element_line(size=0))+ylim(0,15)
p <- p  +guides(colour = FALSE)
p_ko <- p +theme(axis.text=element_text(size=10),axis.title=element_text(size=10)) + 

png("DEG_ko_cluster0.png",res=300,height=1000,width=2000)
p_ko
dev.off()

data <- DEG_test_kivswt_0
logFC <- (data[,"logFC"])
data$sig[(data$FDR > 0.05|data$FDR=="NA")|(data$logFC < 0)& data$logFC > -0] <- "no"
data$sig[data$FDR <= 0.05 & data$logFC >= 0] <- "up"
data$sig[data$FDR <= 0.05 & data$logFC <= 0] <- "down"

x_lim <- max(logFC,-logFC)
library(ggplot2)
library(RColorBrewer)
#pdf(file = "miRNA_volcano.pdf",width=8,height=8)
theme_set(theme_bw())
p <- ggplot(data,aes(logFC,-1*log10(FDR),
                     color = sig))+geom_point()+
    xlim(-10,10) +  labs(x="log2(FoldChange)",y="-log10(FDR)")
p <- p + scale_color_manual(values =c("#0072B5","grey","#BC3C28"))+
    geom_hline(yintercept=-log10(0.05),linetype=4)+
    geom_vline(xintercept=c(-0,0),linetype=4)
p <- p +theme(panel.grid =element_blank())+
    theme(axis.line = element_line(size=0))+ylim(0,15)
p <- p  +guides(colour = FALSE)
p_ki <- p +theme(axis.text=element_text(size=10),axis.title=element_text(size=10))

tiff("DEG_in_cluster_ki_ko.tiff",res=300,height=1000,width=2000)
p_ko+p_ki
dev.off()

########################################################
background <- rownames(STids_splice)
background <- bitr(background,fromType="SYMBOL",toType="ENTREZID",org.Mm.eg.db)

gene_3_ko <- bitr(rownames(DEG_test_kovswt_3)[DEG_test_kovswt_3[,3]<0.05],fromType="SYMBOL",toType="ENTREZID",org.Mm.eg.db)
gene_3_ko <- gene_3_ko[gene_3_ko$SYMBOL %in% Ribosome == FALSE,]
goterm_3 <- enrichGO(names(table(c(gene_3_ko[,2]))), ont="BP",pvalueCutoff = 0.05,qvalueCutoff = 0.05,OrgDb = org.Mm.eg.db,universe=background[,2],readable=TRUE)
goterm_3_filter <- clusterProfiler::gofilter(goterm_3,level = 4)


background <- rownames(STids_splice)
background <- bitr(background,fromType="SYMBOL",toType="ENTREZID",org.Mm.eg.db)

gene_0_ko <- bitr(rownames(DEG_test_kovswt_0)[DEG_test_kovswt_0[,3]<0.05],fromType="SYMBOL",toType="ENTREZID",org.Mm.eg.db)
gene_0_ko <- gene_0_ko[gene_0_ko$SYMBOL %in% Ribosome == FALSE,]
goterm_0 <- enrichGO(names(table(c(gene_0_ko[,2]))), ont="BP",pvalueCutoff = 0.05,qvalueCutoff = 0.05,OrgDb = org.Mm.eg.db,universe=background[,2],readable=TRUE)
goterm_0_filter <- clusterProfiler::gofilter(goterm_0,level = 4)



background <- rownames(STids_splice)
background <- bitr(background,fromType="SYMBOL",toType="ENTREZID",org.Mm.eg.db)

gene_1_ko <- bitr(rownames(DEG_test_kovswt_1)[DEG_test_kovswt_1[,3]<0.05],fromType="SYMBOL",toType="ENTREZID",org.Mm.eg.db)
gene_1_ko <- gene_1_ko[gene_1_ko$SYMBOL %in% Ribosome == FALSE,]
goterm_1 <- enrichGO(names(table(c(gene_1_ko[,2]))), ont="BP",pvalueCutoff = 0.05,qvalueCutoff = 0.05,OrgDb = org.Mm.eg.db,universe=background[,2],readable=TRUE)
goterm_1_filter <- clusterProfiler::gofilter(goterm_1,level = 4)



background <- rownames(STids_splice)
background <- bitr(background,fromType="SYMBOL",toType="ENTREZID",org.Mm.eg.db)

gene_2_ko <- bitr(rownames(DEG_test_kovswt_2)[DEG_test_kovswt_2[,3]<0.05],fromType="SYMBOL",toType="ENTREZID",org.Mm.eg.db)
gene_2_ko <- gene_2_ko[gene_2_ko$SYMBOL %in% Ribosome == FALSE,]
goterm_2 <- enrichGO(names(table(c(gene_2_ko[,2]))), ont="BP",pvalueCutoff = 0.05,qvalueCutoff = 0.05,OrgDb = org.Mm.eg.db,universe=background[,2],readable=TRUE)
goterm_2_filter <- clusterProfiler::gofilter(goterm_2,level = 4)


bk_genes <- bk_genes[,2]
cp = list(subtype_0 = gene_3_ko[,2], subtype_1 = gene_0_ko[,2],subtype_2 = gene_1_ko[,2],subtype_3 = gene_2_ko[,2])  
xx <- compareCluster(cp, fun="enrichGO",OrgDb= "org.Mm.eg.db", pvalueCutoff=0.05,pAdjustMethod = "BH",qvalueCutoff = 0.05,ont="BP",universe=background[,2]) 
xx_filter <- clusterProfiler::gofilter(xx,level = 4)
png("DEG_dotplot.png",width=1600,height=1500,res=2500)
dotplot(xx_filter,showCategory=20,includeAll=TRUE)
dev.off()

########################################################
#plot the heatmap along with velocity pesudotime
metadata_wt <- readRDS("/mnt/data2/wangweixu/mm10/WT_metadata.rds")
geneID <- names(table(c(gene_0_ko_down[,1],gene_3_ko_down[,1])))
geneID <- geneID[geneID %in% Ribosome == FALSE]

##extract the expression profile of SCT
STids.SCT <- GetAssayData(STids_splice)[,STids$batch %in% c("MouseSperm499","MouseSperm500")]
STids.SCT <- STids.SCT[,rownames(metadata_wt)]
STids.SCT <- STids.SCT[geneID,]

##
#Select the genes that associate with pesudotime
cor.gene <- NULL
for(i in 1:nrow(STids.SCT)){
  cor.gene <- c(cor.gene, cor.test(STids.SCT[i,],metadata_wt$velocity_pseudotime,method="sp")$p.value)
}
cor.gene <- p.adjust(cor.gene,method="BH")
STids.SCT <- STids.SCT[cor.gene<0.05,]

##build pheatmap
##identify co-expression module

mat <- as.matrix(STids.SCT)
gsg = goodSamplesGenes(t(mat), verbose = 3);
# Remove the offending genes and samples from the data:
germ_exp = mat[gsg$goodGenes,gsg$goodSamples]

datExpr <- t(germ_exp);rm(germ_exp);gc()
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
## calculate spearman
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5,networkType = "unsigned")

softPower = 2;
adjacency = adjacency(datExpr, power = softPower,type="unsigned");
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);

minModuleSize = 50;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = TRUE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
#tiff("dendrogram_graph.tiff",res=300,height=1000,width=1700)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
#dev.off()

#Select the most representative genes from each clusters and perform go enrichments
kME <- WGCNA::intramodularConnectivity(adjacency,dynamicColors)
kME$module <- dynamicColors

###Perform GO enrichment for each module
geneSets <- bitr(rownames(kME)[kME$module %in% "blue"], fromType="SYMBOL",toType="ENTREZID","org.Mm.eg.db")
goterm_blue <- enrichGO(geneSets[,2], ont="BP",pvalueCutoff = 0.1,qvalueCutoff = 0.1,
            OrgDb = org.Mm.eg.db,readable=TRUE)

geneSets <- bitr(rownames(kME)[kME$module %in% "brown"], fromType="SYMBOL",toType="ENTREZID","org.Mm.eg.db")
goterm_brown <- enrichGO(geneSets[,2], ont="BP",pvalueCutoff = 0.1,qvalueCutoff = 0.1,
            OrgDb = org.Mm.eg.db,readable=TRUE)

geneSets <- bitr(rownames(kME)[kME$module %in% "turquoise"], fromType="SYMBOL",toType="ENTREZID","org.Mm.eg.db")
goterm_turquoise <- enrichGO(geneSets[,2], ont="BP",pvalueCutoff = 0.1,qvalueCutoff = 0.1,
            OrgDb = org.Mm.eg.db,readable=TRUE)

geneSets <- bitr(rownames(kME)[kME$module %in% "yellow"], fromType="SYMBOL",toType="ENTREZID","org.Mm.eg.db")
goterm_yellow <- enrichGO(geneSets[,2], ont="BP",pvalueCutoff = 0.1,qvalueCutoff = 0.1,
            OrgDb = org.Mm.eg.db,readable=TRUE)

####Using top20 genes to plot heatmap
brown_hub <- subset(kME, module %in% "brown" == TRUE)
brown_hub <- brown_hub[order(brown_hub$kWithin,decreasing=TRUE),]
brown_hub <- rownames(brown_hub)[1:20]

blue_hub <- subset(kME, module %in% "blue" == TRUE)
blue_hub <- blue_hub[order(blue_hub$kWithin,decreasing=TRUE),]
blue_hub <- rownames(blue_hub)[1:20]

turquoise_hub <- subset(kME, module %in% "turquoise" == TRUE)
turquoise_hub <- turquoise_hub[order(turquoise_hub$kWithin,decreasing=TRUE),]
turquoise_hub <- rownames(turquoise_hub)[1:20]

yellow_hub <- subset(kME, module %in% "yellow" == TRUE)
yellow_hub <- yellow_hub[order(yellow_hub$kWithin,decreasing=TRUE),]
yellow_hub <- rownames(yellow_hub)[1:20]

heatmap_matrix=as.matrix(STids.SCT[c(brown_hub,blue_hub,turquoise_hub,yellow_hub),order(metadata_wt$velocity_pseudotime,decreasing=FALSE)])

heatmap_matrix=heatmap_matrix[!apply(heatmap_matrix, 1, sd)==0,]
heatmap_matrix=Matrix::t(scale(Matrix::t(heatmap_matrix),center=TRUE))
heatmap_matrix=heatmap_matrix[is.na(row.names(heatmap_matrix)) == FALSE,]
heatmap_matrix[is.nan(heatmap_matrix)] = 0
heatmap_matrix[heatmap_matrix>2.5] = 2.5
heatmap_matrix[heatmap_matrix< -2.5] = -2.5

row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix)))/2)
row_dist[is.na(row_dist)] <- 1

exp_rng <- range(heatmap_matrix) #bks is based on the expression range
bks <- seq(exp_rng[1] - 0.1, exp_rng[2] + 0.1, by=0.1)
if(is.null(hmcols)) {
  hmcols <- blue2green2red(length(bks) - 1)
}

# prin  t(hmcols)
ph_stids <- pheatmap(heatmap_matrix, 
               useRaster = T,
               cluster_cols=FALSE, 
               cluster_rows=TRUE, 
               show_rownames=TRUE, 
               show_colnames=F, 
               #scale="row",
               clustering_distance_rows=row_dist,
               clustering_method = "ward.D2",
               cutree_rows=4,
               silent=TRUE,
               filename=NA,
               breaks=bks,
               color=hmcols
               #color=hmcols#,
               # filename="expression_pseudotime_pheatmap.pdf",
)
annotation_row <- data.frame(Cluster=factor(cutree(ph_stids$tree_row, 4)))
cluster <- rep("empty",nrow(annotation_row))
cluster[annotation_row$Cluster==1] <- "brown"
cluster[annotation_row$Cluster==2] <- "blue"
cluster[annotation_row$Cluster==3] <- "yellow"
cluster[annotation_row$Cluster==4] <- "turquoise"
annotation_row2 <- data.frame(Cluster=factor(cluster))
rownames(annotation_row2) <- rownames(annotation_row)

annotation_col <- data.frame(VelocityTime = metadata_wt[colnames(heatmap_matrix),]$velocity_pseudotime)
rownames(annotation_col) <- colnames(heatmap_matrix)

ann_colors = list(
  Cluster = c(blue = "blue", brown = "brown", yellow = "yellow",turquoise = "turquoise")
)

ph_stids <- pheatmap(heatmap_matrix, 
               useRaster = T,
               cluster_cols=FALSE, 
               cluster_rows=TRUE, 
               show_rownames=TRUE, 
               show_colnames=F, 
               #scale="row",
               clustering_distance_rows=row_dist,
               clustering_method = "ward.D2",
               annotation_row = annotation_row2,
			   annotation_col = annotation_col,
               cutree_rows=4,
               silent=TRUE,
               filename=NA,
               breaks=bks,
               color=hmcols,
			   annotation_colors = ann_colors
               #color=hmcols#,
               # filename="expression_pseudotime_pheatmap.pdf",
)

png("ph_stids.png",res=300,height=2500,width=2500)
ph_stids
dev.off()

###################################
#measure the spliced RNA expression activity
cells_rankings <- AUCell_buildRankings(as.matrix(GetAssayData(STids_splice[,STids$batch %in% c("MouseSperm499","MouseSperm500")])), nCores = 6)
cells_AUC <- AUCell_calcAUC(list(turquoise = rownames(kME)[kME$module %in% "turquoise"]), cells_rankings, aucMaxRank=3000, nCores =6)
exp_turquoise <- t(getAUC(cells_AUC))

exp_turquoise_order <- exp_turquoise[colnames(heatmap_matrix),]
###build loess curve to fit the variation
dat <- cbind(exp_turquoise_order,metadata_wt[colnames(heatmap_matrix),"velocity_pseudotime"])
colnames(dat) <- c("ExpVar","Pesudotime")
dat <- as.data.frame(dat)

model1=loess(ExpVar~Pesudotime,data=dat,span=0.1)
model2=loess(VelocityVar~Pesudotime,data=dat,span=0.1)
dat$ExpVar_smooth <- model1$fitted
dat$VelocityVar_smooth <- model2$fitted

dat_var <- cbind(rep(dat$Pesudotime,1),c(dat$ExpVar_smooth),c(rep("spliced RNA expression",nrow(dat))))
dat_var <- as.data.frame(dat_var)
colnames(dat_var) <- c("Pesudotime","Variation","Curve Type")
dat_var$Pesudotime <- as.numeric(as.matrix(dat_var$Pesudotime))
dat_var$Variation <- as.numeric(as.matrix(dat_var$Variation))

ggplot(data=dat_var, aes(x=Pesudotime, y=Variation, group=`Curve Type`, color=`Curve Type`)) +  geom_line()

###########################################
cells_rankings <- AUCell_buildRankings(as.matrix(GetAssayData(STids_splice[,STids$batch %in% c("MouseSperm341","MouseSperm340")])), nCores = 6)
cells_AUC <- AUCell_calcAUC(list(turquoise = rownames(kME)[kME$module %in% "turquoise"]), cells_rankings, aucMaxRank=3000, nCores =6)
exp_turquoise_ko <- t(getAUC(cells_AUC))

metadata_ko <- readRDS("KO_metadata.rds")
exp_turquoise_ko_order <- exp_turquoise_ko[rownames(metadata_ko),]

dat2 <- cbind(velocity_turquoise_ko_order,exp_turquoise_ko_order,metadata_ko[,"velocity_pseudotime"])
colnames(dat2) <- c("VelocityVar","ExpVar","Pesudotime")
dat2 <- as.data.frame(dat2)

model1=loess(ExpVar~Pesudotime,data=dat2,span=0.1)
dat2$ExpVar_smooth <- model1$fitted

dat_var2 <- cbind(rep(dat2$Pesudotime,1),c(dat2$ExpVar_smooth),c(rep("spliced RNA expression",nrow(dat2))))
dat_var2 <- as.data.frame(dat_var2)
colnames(dat_var2) <- c("Pesudotime","Variation","Curve Type")
dat_var2$Pesudotime <- as.numeric(as.matrix(dat_var2$Pesudotime))
dat_var2$Variation <- as.numeric(as.matrix(dat_var2$Variation))

ggplot(data=dat_var2, aes(x=Pesudotime, y=Variation, group=`Curve Type`, color=`Curve Type`)) +  geom_line()

dat_var2$`Curve Type` <- paste(dat_var2$`Curve Type`,"-KO",sep="")
dat_var_all <- rbind(dat_var,dat_var2)
dat_var_all <- subset(dat_var_all, `Curve Type` %in% c("spliced RNA expression","spliced RNA expression-KO") == TRUE)

ggplot(data=dat_var_all, aes(x=Pesudotime, y=Variation, group=`Curve Type`, color=`Curve Type`)) +  geom_line()

tiff("spliced_RNA_expression_turquoise.tiff",res=300,height=1500,width=2500)
ggplot(data=dat_var_all, aes(x=Pesudotime, y=Variation, group=`Curve Type`, color=`Curve Type`)) +  geom_line(size=3) + theme_classic2()+theme(axis.text=element_text(size=20),axis.title=element_text(size=20))+labs(x="Velocity Pseudotime")
dev.off()

##########################################################################################################
##Running the random walk with restart to check the priorization
mat <- GetAssayData(STids_splice)[,STids$batch %in% c("MouseSperm499","MouseSperm500")]
mat <- mat[,rownames(metadata_wt)]

##build pheatmap
##identify co-expression module

mat <- as.matrix(mat)
gsg = goodSamplesGenes(t(mat), verbose = 3);
# Remove the offending genes and samples from the data:
germ_exp = mat[gsg$goodGenes,gsg$goodSamples]

datExpr <- t(germ_exp);rm(germ_exp);gc()
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
## calculate spearman
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5,networkType = "unsigned")

softPower = 2;
adjacency_all = adjacency(datExpr, power = softPower,type="unsigned");
TOM_all = TOMsimilarity(adjacency_all);

################################################
RWRtrace <- function(start_vec, trP_m = trP, beta=0.5){
   s_new <- start_vec
   s <- as.matrix(sample(c(0,1),length(s_new),replace=TRUE))
   repeat{
   
   if(abs(cor(s_new,s))>0.99999){
     break
   }
   
   s <- s_new
   s_new <- (1-beta) * t(s) %*% trP + beta * t(s)
   s_new <- t(as.matrix(s_new))
   print(abs(cor(s_new,s)))
   }
   
   return(s_new)
}

###Build transition matrix based on TOM
trP <- diag(1/rowSums(TOM_all)) %*% TOM_all
start_vec <- rep(0,nrow(trP)); names(start_vec)<-colnames(trP);start_vec[rownames(TOM_all) %in% rownames(kME)[kME$module %in% "brown"]] <- 1
start_vec <- as.matrix(start_vec)
vec_brown <- RWRtrace(start_vec,trP,beta=0.5)

start_vec <- rep(0,nrow(trP)); names(start_vec)<-colnames(trP);start_vec[rownames(TOM_all) %in% rownames(kME)[kME$module %in% "turquoise"]] <- 1
start_vec <- as.matrix(start_vec)
vec_turquoise <- RWRtrace(start_vec,trP,beta=0.5)

start_vec <- rep(0,nrow(trP)); names(start_vec)<-colnames(trP);start_vec[rownames(TOM_all) %in% rownames(kME)[kME$module %in% "blue"]] <- 1
start_vec <- as.matrix(start_vec)
vec_blue <- RWRtrace(start_vec,trP,beta=0.5)

start_vec <- rep(0,nrow(trP)); names(start_vec)<-colnames(trP);start_vec[rownames(TOM_all) %in% rownames(kME)[kME$module %in% "yellow"]] <- 1
start_vec <- as.matrix(start_vec)
vec_yellow <- RWRtrace(start_vec,trP,beta=0.5)

vec_brown_order <- vec_brown[order(vec_brown,decreasing=TRUE),]


vec_label <- rep(0,nrow(trP)); names(vec_label)<-colnames(trP);vec_label[rownames(TOM_all) %in% rownames(kME)] <- 1
####################################################################
##Build connect graph based on page ranking
##Select the genes with pr score large than 0.1
blueGene <- rownames(vec_blue)[vec_blue[,1]>0.1]
turquoiseGene <- rownames(vec_turquoise)[vec_turquoise[,1]>0.1]
yellowGene <- rownames(vec_yellow)[vec_yellow[,1]>0.1]
brownGene <- rownames(vec_brown)[vec_brown[,1]>0.1]

keyGene <- c(blueGene,turquoiseGene,yellowGene,brownGene)
keyGene <- names(table(keyGene)[table(keyGene)>=3])

###Build graph based on page ranking
graph_mat_id <- c(keyGene)
graph_mat <- matrix(0,nrow=length(graph_mat_id),ncol=4)
rownames(graph_mat) <- graph_mat_id
colnames(graph_mat) <- c("blue","turquoise","yellow","brown")

for(i in keyGene){
  if(sum(blueGene %in% i)==1) graph_mat[i,"blue"] <- vec_blue[i,]
  if(sum(brownGene %in% i)==1) graph_mat[i,"brown"] <- vec_brown[i,]
  if(sum(yellowGene %in% i)==1) graph_mat[i,"yellow"] <- vec_yellow[i,]
  if(sum(turquoiseGene %in% i)==1) graph_mat[i,"turquoise"] <- vec_turquoise[i,]
}

nch.net<- bip_init_network(graph_mat)

coordP <- cbind(rep(3, dim(graph_mat)[1]), seq(1, dim(graph_mat)[1]) + 
                    2)
coordA <- cbind(rep(3, dim(graph_mat)[2]), seq(1, dim(graph_mat)[2]) + 
                    2)
coordA[,2] <- coordA[,2]+3
coordA[1:3,1] <- coordA[1:3,1]-3
coordA[4,1] <- coordA[4,1]+3
coordA[4,2] <- coordA[4,2]-2
coordA[1,2] <- coordA[1,2]-1
coordA[3,2] <- coordA[3,2]+1
mylayout <- as.matrix(rbind(coordP, coordA))

col = c(Module = "grey", Gene = "gold")
pp <- GGally::ggnet2(nch.net, shape = "mode", label = TRUE, color = "mode", 
        palette = "Set1", size = 9, legend.size = 9, mode = mylayout, 
        label.size = 8, layout.par = NULL, layout.exp = 0, 
        size.legend = NA, label.trim = FALSE, edge.lty = "solid", 
        edge.label = NULL, edge.size = bip_edgewt(graph_mat, 
            5), edge.alpha = 0.25) + coord_flip()
##########################################################
#plot brown module
adjacency_brown <- adjacency[kME$module %in% c("grey") == FALSE,kME$module %in% c("grey") == FALSE]
diag(adjacency_brown) <- 0
adjacency_brown[adjacency_brown<0.05] <- 0
adjacency_brown[adjacency_brown!=0] <- 1
adjacency_brown <- adjacency_brown[rowSums(adjacency_brown)>0,rowSums(adjacency_brown)>0]
#net <- graph_from_adjacency_matrix(adjacency_brown,mode="undirected")
net = network(adjacency_brown, directed = FALSE)

geneClass <- as.character(kME[rownames(adjacency_brown),]$module)
geneClass[rownames(adjacency_brown) %in% c("Elfn2","Acrv1")] <- c("Elfn2","Acrv1")

y <- names(table(c(geneClass,"Elfn2")))
col <- names(table(c(kME[rownames(adjacency_brown),]$module,"red")))
names(col) <- y

net %v% "module" = geneClass
ggnet2(net, color = "module", palette = col, alpha = 0.75, size = 4, edge.alpha = 0.5)

############################################################
#barplot of brown module importance
kME_brown <- kME[kME$module %in% "brown",]
kME_brown <- kME_brown[order(kME_brown$kWithin,decreasing=TRUE),]

df <- data.frame(gene=rownames(kME_brown),degree=kME_brown$kWithin)[1:20,]
df <- df[order(df$degree,decreasing=FALSE),]
df$gene <- factor(df$gene,levels=as.character(as.matrix(df$gene)))
p<-ggplot(df, aes(x=gene, y=degree,fill=degree)) +
    geom_bar(stat="identity")+theme_minimal() +scale_fill_continuous(type = "viridis")
tiff("ranking_brown.tiff",res=300,height=1800,width=2300)
#p+mytheme+labs(x="",y="Degree",main="Connectivity in Brown")+coord_flip()
p+labs(x="",y="Degree",main="Connectivity in Brown")+theme_classic2(base_size=18)+coord_flip()
#p+mytheme+labs(x="",y="RWR score",main="RWR distance")+coord_flip()
dev.off()

p_brown <- ((length(vec_brown)-rank(vec_brown,ties.method = "average"))+1)/(length(vec_brown)+1)
p_blue <- ((length(vec_blue)-rank(vec_blue,ties.method = "average"))+1)/(length(vec_blue)+1)
p_yellow <- ((length(vec_yellow)-rank(vec_yellow,ties.method = "average"))+1)/(length(vec_yellow)+1)
p_turquoise <- ((length(vec_turquoise)-rank(vec_turquoise,ties.method = "average"))+1)/(length(vec_turquoise)+1)

p <- (1-p_brown)*(1-p_blue)*(1-p_yellow)*(1-p_turquoise)
names(p) <- rownames(vec_brown)
#p <- p[keyGene]
#p <- -log10(p)
p <- p[order(p,decreasing=TRUE)]

#################
library(ggplot2)
library(ggridges)
d_p <- data.frame(gene=names(p),page_ranking_score=p,module=rep("module",length(p)))
p_page <- ggplot(
    d_p, 
    aes(x = page_ranking_score, y = module, fill = stat(x))
) +
    geom_density_ridges_gradient() +
    labs(title = 'Association to the modules')+xlim(0,1) +scale_fill_continuous(type = "viridis")+labs(x="Association Score",y="No. genes")+theme_classic2(base_size=18)+theme(axis.text.y = element_blank(),axis.ticks.y = element_blank())

tiff("/mnt/data2/wangweixu/MouseGerm_figure/RWR_pdf.tiff",res=300,height=1300,width=2500)
p_page+NoLegend()
dev.off()

p <- vec_brown[names(p)[1:30],]
p <- p[order(p,decreasing=TRUE)]

df_p <- data.frame(gene=names(p),page_ranking_score=p)
df_p <- df_p[order(df_p$page_ranking_score,decreasing=FALSE),]
df_p$gene <- factor(df_p$gene,levels=as.character(as.matrix(df_p$gene)))
df_p$brown_score <- vec_brown[as.character(as.matrix(df_p$gene)),]
df_p <- df_p[order(df_p$brown_score,decreasing=FALSE),]
df_p$gene <- factor(df_p$gene,levels=as.character(as.matrix(df_p$gene)))
p<-ggplot(df_p, aes(x=gene, y=brown_score,fill=brown_score)) +
    geom_bar(stat="identity")+theme_minimal() +scale_fill_continuous(type = "viridis")
p+labs(x="",y="RWR score",main="RWR distance")+theme_classic2(base_size=18)+coord_flip()
#p+mytheme+theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())+labs(x="",y="RWR score")

tiff("/mnt/data2/wangweixu/MouseGerm_figure/RWR_ranking_main.tiff",res=300,height=1800,width=2300)
p+labs(x="",y="RWR score",main="RWR distance")+theme_classic2(base_size=18)+coord_flip()
dev.off()

#################################################################
###build loess curve to fit the variation
expr_rwr <- t(GetAssayData(MouseGerm_WT)[as.character(as.matrix(df_p$gene)),])
dat <- cbind(expr_rwr,metadata_wt[rownames(expr_rwr),"velocity_pseudotime"])
colnames(dat) <- c(colnames(expr_rwr),"Pesudotime")
dat <- as.data.frame(dat)
#rm(expr_rwr);gc()

for(i in colnames(expr_rwr)){
model=loess(dat[,i]~dat[,21],span=0.1)
pre <- predict(model, dat[,21])
dat[,i] <- pre
}

dat_var <- NULL
for(i in colnames(expr_rwr)){
  dat_var <- rbind(dat_var, cbind(dat$Pesudotime,scale(dat[,i],scale=FALSE),rep(i,length(dat[,i]))))
}
dat_var <- as.data.frame(dat_var)
colnames(dat_var) <- c("Pesudotime","Variation","Curve Type")
dat_var$Pesudotime <- as.numeric(as.matrix(dat_var$Pesudotime))
dat_var$Variation <- as.numeric(as.matrix(dat_var$Variation))
dat_var$Elfn2_ind <- rep("Non Elfn2",nrow(dat_var))
dat_var$Elfn2_ind[dat_var$`Curve Type` == "Elfn2"] <- "Elfn2"

tiff("/mnt/data2/wangweixu/MouseGerm_figure/expression_pattern.tiff",res=300,height=1500,width=2300)
ggplot(data=dat_var, aes(x=Pesudotime, y=Variation, group=`Curve Type`, color=`Elfn2_ind`)) +  geom_line(size=1) + theme_classic2(base_size=18)+NoLegend()
dev.off()

#################################################
df_brown <- data.frame(module=rep("brown",length(vec_brown)),order=c(1:length(vec_brown)),page_ranking_score=vec_brown[order(vec_brown,decreasing=TRUE)])
p_brown <-ggplot(data=df_brown,aes(x=order,y=page_ranking_score)) + geom_point(size=1) + xlab("")+ylab("RWR_score") + theme_classic2() + labs(title="Brown")

df_blue <- data.frame(module=rep("blue",length(vec_blue)),order=c(1:length(vec_blue)),page_ranking_score=vec_blue[order(vec_blue,decreasing=TRUE)])
p_blue <-ggplot(data=df_blue,aes(x=order,y=page_ranking_score)) + geom_point(size=1) + xlab("")+ylab("RWR_score") + theme_classic2() + labs(title="blue")

df_turquoise <- data.frame(module=rep("turquoise",length(vec_turquoise)),order=c(1:length(vec_turquoise)),page_ranking_score=vec_turquoise[order(vec_turquoise,decreasing=TRUE)])
p_turquoise <-ggplot(data=df_turquoise,aes(x=order,y=page_ranking_score)) + geom_point(size=1) + xlab("")+ylab("RWR_score") + theme_classic2() + labs(title="turquoise")

df_yellow <- data.frame(module=rep("yellow",length(vec_yellow)),order=c(1:length(vec_yellow)),page_ranking_score=vec_yellow[order(vec_yellow,decreasing=TRUE)])
p_yellow <-ggplot(data=df_yellow,aes(x=order,y=page_ranking_score)) + geom_point(size=1) + xlab("")+ylab("RWR_score") + theme_classic2() + labs(title="yellow")



#################################


STids.SCT <- GetAssayData(STids_splice)[,STids$seurat_clusters %in% 3 == TRUE]

STids.SCT.ki <- STids.SCT[,colnames(STids.SCT) %in% colnames(MouseGerm_Ki)]
STids.SCT.ko <- STids.SCT[,colnames(STids.SCT) %in% colnames(MouseGerm_KO)]
STids.SCT.wt <- STids.SCT[,colnames(STids.SCT) %in% colnames(MouseGerm_WT)]
batch.ki <- MouseGerm_Ki$batch[colnames(STids.SCT.ki)]
batch.ko <- MouseGerm_KO$batch[colnames(STids.SCT.ko)]
batch.wt <- MouseGerm_WT$batch[colnames(STids.SCT.wt)]

DEG_test_kivswt_3 <- KivsWt(cbind(STids.SCT.ki,STids.SCT.wt), c(batch.ki,batch.wt))
DEG_test_kovswt_3 <- KovsWt_all(cbind(STids.SCT.ko,STids.SCT.wt), c(batch.ko,batch.wt))


STids.SCT <- GetAssayData(STids_splice)[,STids$seurat_clusters %in% 1 == TRUE]

STids.SCT.ki <- STids.SCT[,colnames(STids.SCT) %in% colnames(MouseGerm_Ki)]
STids.SCT.ko <- STids.SCT[,colnames(STids.SCT) %in% colnames(MouseGerm_KO)]
STids.SCT.wt <- STids.SCT[,colnames(STids.SCT) %in% colnames(MouseGerm_WT)]
batch.ki <- MouseGerm_Ki$batch[colnames(STids.SCT.ki)]
batch.ko <- MouseGerm_KO$batch[colnames(STids.SCT.ko)]
batch.wt <- MouseGerm_WT$batch[colnames(STids.SCT.wt)]

DEG_test_kivswt_1 <- KivsWt(cbind(STids.SCT.ki,STids.SCT.wt), c(batch.ki,batch.wt))
DEG_test_kovswt_1 <- KovsWt_all(cbind(STids.SCT.ko,STids.SCT.wt), c(batch.ko,batch.wt))

################################################################
#Select sub-column of DEG
DEG_test_kivswt_1 <- DEG_test_kivswt_1[,10:12]
DEG_test_kivswt_2 <- DEG_test_kivswt_2[,10:12]
DEG_test_kivswt_3 <- DEG_test_kivswt_3[,10:12]

DEG_test_kovswt_1 <- DEG_test_kovswt_1[,10:12]
DEG_test_kovswt_2 <- DEG_test_kovswt_2[,10:12]
DEG_test_kovswt_3 <- DEG_test_kovswt_3[,10:12]

#########################################
#test GO enrichment pathway on each stage
background <- rownames(STids.new)
background <- bitr(background,fromType="SYMBOL",toType="ENTREZID",org.Mm.eg.db)
gene_0_ki_down <- bitr(rownames(DEG_test_kivswt_0)[DEG_test_kivswt_0[,3]<0.05&DEG_test_kivswt_0[,1]<0],fromType="SYMBOL",toType="ENTREZID",org.Mm.eg.db)
goterm_0_ki_down <- enrichGO(gene_0_ki_down[,2], ont="BP",pvalueCutoff = 0.05,qvalueCutoff = 0.05,OrgDb = org.Mm.eg.db,universe=background[,2],readable=TRUE)

gene_0_ko_down <- bitr(rownames(DEG_test_kovswt_0)[DEG_test_kovswt_0[,3]<0.05&DEG_test_kovswt_0[,1]<0],fromType="SYMBOL",toType="ENTREZID",org.Mm.eg.db)
goterm_0_ko_down <- enrichGO(gene_0_ko_down[,2], ont="BP",pvalueCutoff = 0.05,qvalueCutoff = 0.05,OrgDb = org.Mm.eg.db,universe=background[,2],readable=TRUE)

gene_2_ki_down <- bitr(rownames(DEG_test_kivswt_2)[DEG_test_kivswt_2[,3]<0.05&DEG_test_kivswt_2[,1]<0],fromType="SYMBOL",toType="ENTREZID",org.Mm.eg.db)
goterm_2_ki_down <- enrichGO(gene_2_ki_down[,2], ont="BP",pvalueCutoff = 0.05,qvalueCutoff = 0.05,OrgDb = org.Mm.eg.db,readable=TRUE)

gene_2_ko_down <- bitr(rownames(DEG_test_kovswt_2)[DEG_test_kovswt_2[,3]<0.05&DEG_test_kovswt_2[,1]<0],fromType="SYMBOL",toType="ENTREZID",org.Mm.eg.db)
goterm_2_ko_down <- enrichGO(gene_2_ko_down[,2], ont="BP",pvalueCutoff = 0.05,qvalueCutoff = 0.05,OrgDb = org.Mm.eg.db,readable=TRUE)

gene_1_ki_down <- bitr(rownames(DEG_test_kivswt_1)[DEG_test_kivswt_1[,3]<0.05&DEG_test_kivswt_1[,1]<0],fromType="SYMBOL",toType="ENTREZID",org.Mm.eg.db)
goterm_1_ki_down <- enrichGO(gene_1_ki_down[,2], ont="BP",pvalueCutoff = 0.05,qvalueCutoff = 0.05,OrgDb = org.Mm.eg.db,readable=TRUE)

gene_1_ko_down <- bitr(rownames(DEG_test_kovswt_1)[DEG_test_kovswt_1[,3]<0.05&DEG_test_kovswt_1[,1]<0],fromType="SYMBOL",toType="ENTREZID",org.Mm.eg.db)
goterm_1_ko_down <- enrichGO(gene_1_ko_down[,2], ont="BP",pvalueCutoff = 0.05,qvalueCutoff = 0.05,OrgDb = org.Mm.eg.db,readable=TRUE)

gene_3_ko_down <- bitr(rownames(DEG_test_kovswt_3)[DEG_test_kovswt_3[,3]<0.05&DEG_test_kovswt_3[,1]<0],fromType="SYMBOL",toType="ENTREZID",org.Mm.eg.db)
goterm_3_ko_down <- enrichGO(gene_3_ko_down[,2], ont="BP",pvalueCutoff = 0.05,qvalueCutoff = 0.05,OrgDb = org.Mm.eg.db,readable=TRUE)

gene_3_ki_down <- bitr(rownames(DEG_test_kivswt_3)[DEG_test_kivswt_3[,3]<0.05&DEG_test_kivswt_3[,1]<0],fromType="SYMBOL",toType="ENTREZID",org.Mm.eg.db)
goterm_3_ki_down <- enrichGO(gene_3_ki_down[,2], ont="BP",pvalueCutoff = 0.05,qvalueCutoff = 0.05,OrgDb = org.Mm.eg.db,readable=TRUE)

#################################################################
##
goterm_0_ki_down
df_barMarker2 <- data.frame(pval = -log10(summary(goterm_0_ki_down)[1:5,]$p.adjust),goterm=summary(goterm_0_ki_down)[1:5,]$Description)
df_barMarker2$goterm_name <- factor(c(summary(goterm_0_ki_down)[1:5,]$Description),levels=c(c(summary(goterm_0_ki_down)[1:5,]$Description)[order(-log10(summary(goterm_0_ki_down)[1:5,]$p.adjust),decreasing=FALSE)]))
p_0_go_ki <-ggplot(data=df_barMarker2, aes(x=goterm_name, y=pval)) +
    geom_bar(stat="identity",fill=rgb(250,157,150,maxColorValue = 255))+coord_flip()+theme_classic()+theme(text = element_text(family = "Arial", size=15))+labs(x="",y="-log10(P-value)")
tiff("Ki_down_0.tiff",res=300,height=800,width=2100)
p_0_go_ki
dev.off()

df_barMarker2 <- data.frame(pval = -log10(summary(goterm_0_ko_down)[1:5,]$p.adjust),goterm=summary(goterm_0_ko_down)[1:5,]$Description)
df_barMarker2$goterm_name <- factor(c(summary(goterm_0_ko_down)[1:5,]$Description),levels=c(c(summary(goterm_0_ko_down)[1:5,]$Description)[order(-log10(summary(goterm_0_ko_down)[1:5,]$p.adjust),decreasing=FALSE)]))
p_0_go_ko <-ggplot(data=df_barMarker2, aes(x=goterm_name, y=pval)) +
    geom_bar(stat="identity",fill=rgb(250,157,150,maxColorValue = 255))+coord_flip()+theme_classic()+theme(text = element_text(family = "Arial", size=15))+labs(x="",y="-log10(P-value)")
tiff("KO_down_0.tiff",res=300,height=800,width=2300)
p_0_go_ko
dev.off()



df_barMarker2 <- data.frame(pval = -log10(summary(goterm_2_ki_down)[1:5,]$p.adjust),goterm=summary(goterm_2_ki_down)[1:5,]$Description)
df_barMarker2$goterm_name <- factor(c(summary(goterm_2_ki_down)[1:5,]$Description),levels=c(c(summary(goterm_2_ki_down)[1:5,]$Description)[order(-log10(summary(goterm_2_ki_down)[1:5,]$p.adjust),decreasing=FALSE)]))
p_2_go_ki <-ggplot(data=df_barMarker2, aes(x=goterm_name, y=pval)) +
    geom_bar(stat="identity",fill=rgb(0,188,194,maxColorValue = 255))+coord_flip()+theme_classic()+theme(text = element_text(family = "Arial", size=15))+labs(x="",y="-log10(P-value)")
tiff("Ki_down_2.tiff",res=300,height=800,width=2100)
p_2_go_ki
dev.off()

df_barMarker2 <- data.frame(pval = -log10(summary(goterm_2_ko_down)[1:5,]$p.adjust),goterm=summary(goterm_2_ko_down)[1:5,]$Description)
df_barMarker2$goterm_name <- factor(c(summary(goterm_2_ko_down)[1:5,]$Description),levels=c(c(summary(goterm_2_ko_down)[1:5,]$Description)[order(-log10(summary(goterm_2_ko_down)[1:5,]$p.adjust),decreasing=FALSE)]))
p_2_go_ko <-ggplot(data=df_barMarker2, aes(x=goterm_name, y=pval)) +
    geom_bar(stat="identity",fill=rgb(0,188,194,maxColorValue = 255))+coord_flip()+theme_classic()+theme(text = element_text(family = "Arial", size=15))+labs(x="",y="-log10(P-value)")
tiff("KO_down_2.tiff",res=300,height=800,width=2300)
p_2_go_ko
dev.off()

                                

df_barMarker2 <- data.frame(pval = -log10(summary(goterm_1_ki_down)[1:5,]$p.adjust),goterm=summary(goterm_1_ki_down)[1:5,]$Description)
df_barMarker2$goterm_name <- factor(c(summary(goterm_1_ki_down)[1:5,]$Description),levels=c(c(summary(goterm_1_ki_down)[1:5,]$Description)[order(-log10(summary(goterm_1_ki_down)[1:5,]$p.adjust),decreasing=FALSE)]))
p_1_go_ki <-ggplot(data=df_barMarker2, aes(x=goterm_name, y=pval)) +
    geom_bar(stat="identity",fill=rgb(184,211,116,maxColorValue = 255))+coord_flip()+theme_classic()+theme(text = element_text(family = "Arial", size=15))+labs(x="",y="-log10(P-value)")
tiff("Ki_down_1.tiff",res=300,height=800,width=2100)
p_1_go_ki
dev.off()

df_barMarker2 <- data.frame(pval = -log10(summary(goterm_1_ko_down)[1:5,]$p.adjust),goterm=summary(goterm_1_ko_down)[1:5,]$Description)
df_barMarker2$goterm_name <- factor(c(summary(goterm_1_ko_down)[1:5,]$Description),levels=c(c(summary(goterm_1_ko_down)[1:5,]$Description)[order(-log10(summary(goterm_1_ko_down)[1:5,]$p.adjust),decreasing=FALSE)]))
p_1_go_ko <-ggplot(data=df_barMarker2, aes(x=goterm_name, y=pval)) +
    geom_bar(stat="identity",fill=rgb(184,211,116,maxColorValue = 255))+coord_flip()+theme_classic()+theme(text = element_text(family = "Arial", size=15))+labs(x="",y="-log10(P-value)")
tiff("KO_down_1.tiff",res=300,height=800,width=2300)
p_1_go_ko
dev.off()
###############################################################
##Test the binding target enrichments
length(intersect(names(table(peak1_annot$Symbol)),gene_0_ki_down$SYMBOL))
DEG_vec <- rep(0,length(rownames(STids_splice)))
DEG_vec[rownames(STids_splice) %in% gene_0_ki_down$SYMBOL] <- 1
target_wt1_vec <- rep(0,length(rownames(STids_splice)))
target_wt1_vec[rownames(STids_splice) %in% names(table(peak1_annot$Symbol))] <- 1
target_wt2_vec <- rep(0,length(rownames(STids_splice)))
target_wt2_vec[rownames(STids_splice) %in% names(table(peak2_annot$Symbol))] <- 1

fisher.test(table(DEG_vec,target_wt1_vec),alternative="greater")$p.value
fisher.test(table(DEG_vec,target_wt2_vec),alternative="greater")$p.value

DEG_vec <- rep(0,length(rownames(STids_splice)))
DEG_vec[rownames(STids_splice) %in% gene_0_ko_down$SYMBOL] <- 1
target_wt1_vec <- rep(0,length(rownames(STids_splice)))
target_wt1_vec[rownames(STids_splice) %in% names(table(peak1_annot$Symbol))] <- 1
target_wt2_vec <- rep(0,length(rownames(STids_splice)))
target_wt2_vec[rownames(STids_splice) %in% names(table(peak2_annot$Symbol))] <- 1

fisher.test(table(DEG_vec,target_wt1_vec),alternative="greater")
fisher.test(table(DEG_vec,target_wt2_vec),alternative="greater")
#######################################################

DEG_vec <- rep(0,length(rownames(STids_splice)))
DEG_vec[rownames(STids_splice) %in% gene_2_ki_down$SYMBOL] <- 1
target_wt1_vec <- rep(0,length(rownames(STids_splice)))
target_wt1_vec[rownames(STids_splice) %in% names(table(peak1_annot$Symbol))] <- 1
target_wt2_vec <- rep(0,length(rownames(STids_splice)))
target_wt2_vec[rownames(STids_splice) %in% names(table(peak2_annot$Symbol))] <- 1

fisher.test(table(DEG_vec,target_wt1_vec),alternative="greater")
fisher.test(table(DEG_vec,target_wt2_vec),alternative="greater")

DEG_vec <- rep(0,length(rownames(STids_splice)))
DEG_vec[rownames(STids_splice) %in% gene_2_ko_down$SYMBOL] <- 1
target_wt1_vec <- rep(0,length(rownames(STids_splice)))
target_wt1_vec[rownames(STids_splice) %in% names(table(peak1_annot$Symbol))] <- 1
target_wt2_vec <- rep(0,length(rownames(STids_splice)))
target_wt2_vec[rownames(STids_splice) %in% names(table(peak2_annot$Symbol))] <- 1

fisher.test(table(DEG_vec,target_wt1_vec),alternative="greater")
fisher.test(table(DEG_vec,target_wt2_vec),alternative="greater")

###############################################
DEG_vec <- rep(0,length(rownames(STids_splice)))
DEG_vec[rownames(STids_splice) %in% gene_1_ki_down$SYMBOL] <- 1
target_wt1_vec <- rep(0,length(rownames(STids_splice)))
target_wt1_vec[rownames(STids_splice) %in% names(table(peak1_annot$Symbol))] <- 1
target_wt2_vec <- rep(0,length(rownames(STids_splice)))
target_wt2_vec[rownames(STids_splice) %in% names(table(peak2_annot$Symbol))] <- 1

fisher.test(table(DEG_vec,target_wt1_vec),alternative="greater")
fisher.test(table(DEG_vec,target_wt2_vec),alternative="greater")


DEG_vec <- rep(0,length(rownames(STids_splice)))
DEG_vec[rownames(STids_splice) %in% gene_1_ko_down$SYMBOL] <- 1
target_wt1_vec <- rep(0,length(rownames(STids_splice)))
target_wt1_vec[rownames(STids_splice) %in% names(table(peak1_annot$Symbol))] <- 1
target_wt2_vec <- rep(0,length(rownames(STids_splice)))
target_wt2_vec[rownames(STids_splice) %in% names(table(peak2_annot$Symbol))] <- 1

fisher.test(table(DEG_vec,target_wt1_vec),alternative="greater")
fisher.test(table(DEG_vec,target_wt2_vec),alternative="greater")

DEG_vec <- rep(0,length(rownames(STids_splice)))
DEG_vec[rownames(STids_splice) %in% gene_3_ki_down$SYMBOL] <- 1
target_wt1_vec <- rep(0,length(rownames(STids_splice)))
target_wt1_vec[rownames(STids_splice) %in% names(table(peak1_annot$Symbol))] <- 1
target_wt2_vec <- rep(0,length(rownames(STids_splice)))
target_wt2_vec[rownames(STids_splice) %in% names(table(peak2_annot$Symbol))] <- 1

fisher.test(table(DEG_vec,target_wt1_vec),alternative="greater")
fisher.test(table(DEG_vec,target_wt2_vec),alternative="greater")



DEG_vec <- rep(0,length(rownames(STids_splice)))
DEG_vec[rownames(STids_splice) %in% gene_3_ko_up$SYMBOL] <- 1
target_wt1_vec <- rep(0,length(rownames(STids_splice)))
target_wt1_vec[rownames(STids_splice) %in% names(table(peak1_annot$Symbol))] <- 1
target_wt2_vec <- rep(0,length(rownames(STids_splice)))
target_wt2_vec[rownames(STids_splice) %in% names(table(peak2_annot$Symbol))] <- 1

fisher.test(table(DEG_vec,target_wt1_vec),alternative="greater")
##0.4131 130 overlapped
fisher.test(table(DEG_vec,target_wt2_vec),alternative="greater")
#########################################################################
##delate ribosome
marker <- stids.markers[stids.markers[,6]==0&stids.markers[,5]<0.05,7]
marker <- marker[marker %in% Ribosome ==FALSE]
gene_0_marker <- bitr(marker,fromType="SYMBOL",toType="ENTREZID",org.Mm.eg.db)
goterm_0_marker <- enrichGO(gene_0_marker[,2], ont="BP",pvalueCutoff = 0.05,qvalueCutoff = 0.05,OrgDb = org.Mm.eg.db,readable=TRUE)

marker <- stids.markers[stids.markers[,6]==1&stids.markers[,5]<0.05,7]
marker <- marker[marker %in% Ribosome ==FALSE]
gene_1_marker <- bitr(marker,fromType="SYMBOL",toType="ENTREZID",org.Mm.eg.db)
goterm_1_marker <- enrichGO(gene_1_marker[,2], ont="BP",pvalueCutoff = 0.05,qvalueCutoff = 0.05,OrgDb = org.Mm.eg.db,readable=TRUE)

marker <- stids.markers[stids.markers[,6]==2&stids.markers[,5]<0.05,7]
marker <- marker[marker %in% Ribosome ==FALSE]
gene_2_marker <- bitr(marker,fromType="SYMBOL",toType="ENTREZID",org.Mm.eg.db)
goterm_2_marker <- enrichGO(gene_2_marker[,2], ont="BP",pvalueCutoff = 0.05,qvalueCutoff = 0.05,OrgDb = org.Mm.eg.db,readable=TRUE)

marker <- stids.markers[stids.markers[,6]==3&stids.markers[,5]<0.05,7]
marker <- marker[marker %in% Ribosome ==FALSE]
gene_3_marker <- bitr(marker,fromType="SYMBOL",toType="ENTREZID",org.Mm.eg.db)
goterm_3_marker <- enrichGO(gene_3_marker[,2], ont="BP",pvalueCutoff = 0.05,qvalueCutoff = 0.05,OrgDb = org.Mm.eg.db,readable=TRUE)

####Select important go term to plot
revigo_0_marker <- summary(goterm_0_marker)$qvalue
revigo_0_marker_GO <- rownames(summary(goterm_0_marker)) 
revigo_0_marker_GO <- cbind(revigo_0_marker_GO, revigo_0_marker)
revigo_0_marker_GO <- revigo_0_marker_GO[order(revigo_0_marker,decreasing=FALSE),]
revigo_0_marker_GO <- as.data.frame(revigo_0_marker_GO)


revigo_1_marker <- summary(goterm_1_marker)$qvalue
revigo_1_marker_GO <- rownames(summary(goterm_1_marker)) 
revigo_1_marker_GO <- cbind(revigo_1_marker_GO, revigo_1_marker)
revigo_1_marker_GO <- revigo_1_marker_GO[order(revigo_1_marker,decreasing=FALSE),]
revigo_1_marker_GO <- as.data.frame(revigo_1_marker_GO)


revigo_2_marker <- summary(goterm_2_marker)$qvalue
revigo_2_marker_GO <- rownames(summary(goterm_2_marker)) 
revigo_2_marker_GO <- cbind(revigo_2_marker_GO, revigo_2_marker)
revigo_2_marker_GO <- revigo_2_marker_GO[order(revigo_2_marker,decreasing=FALSE),]
revigo_2_marker_GO <- as.data.frame(revigo_2_marker_GO)

revigo_3_marker <- summary(goterm_3_marker)$qvalue
revigo_3_marker_GO <- rownames(summary(goterm_3_marker)) 
revigo_3_marker_GO <- cbind(revigo_3_marker_GO, revigo_3_marker)
revigo_3_marker_GO <- revigo_3_marker_GO[order(revigo_3_marker,decreasing=FALSE),]
revigo_3_marker_GO <- as.data.frame(revigo_3_marker_GO)
#################################################################
#Load Go filtering
GO_filter <- read.csv("GO_filter.csv",header=TRUE,stringsAsFactors=FALSE)
GO_all <- rbind(summary(goterm_0_marker)[GO_filter$GO[GO_filter$Group==0],],summary(goterm_3_marker)[GO_filter$GO[GO_filter$Group==3],],summary(goterm_2_marker)[GO_filter$GO[GO_filter$Group==2],],summary(goterm_1_marker)[GO_filter$GO[GO_filter$Group==1],])

df_barMarker2 <- data.frame(pval = -log10(GO_all$p.adjust),goterm=GO_all$Description,group=GO_filter$Group)
df_barMarker2$goterm <- as.character(as.matrix(df_barMarker2$goterm))
df_barMarker2 <- subset(df_barMarker2,group %in% 3)
df_barMarker2$goterm_name <- factor(c(df_barMarker2$goterm),levels=c(c(df_barMarker2$goterm)[order(df_barMarker2$pval,decreasing=FALSE)]))
go_all_3 <-ggplot(data=df_barMarker2, aes(x=goterm_name, y=pval)) +
    geom_bar(stat="identity",fill=rgb(232,202,255,maxColorValue = 255))+coord_flip()+theme_classic()+theme(text = element_text(family = "Arial", size=15))+labs(x="",y="-log10(P-value)")

df_barMarker2 <- data.frame(pval = -log10(GO_all$p.adjust),goterm=GO_all$Description,group=GO_filter$Group)
df_barMarker2$goterm <- as.character(as.matrix(df_barMarker2$goterm))
df_barMarker2 <- subset(df_barMarker2,group %in% 0)
df_barMarker2$goterm_name <- factor(c(df_barMarker2$goterm),levels=c(c(df_barMarker2$goterm)[order(df_barMarker2$pval,decreasing=FALSE)]))
go_all_0 <-ggplot(data=df_barMarker2, aes(x=goterm_name, y=pval)) +
    geom_bar(stat="identity",fill=rgb(250,157,150,maxColorValue = 255))+coord_flip()+theme_classic()+theme(text = element_text(family = "Arial", size=15))+labs(x="",y="-log10(P-value)")
	
df_barMarker2 <- data.frame(pval = -log10(GO_all$p.adjust),goterm=GO_all$Description,group=GO_filter$Group)
df_barMarker2$goterm <- as.character(as.matrix(df_barMarker2$goterm))
df_barMarker2 <- subset(df_barMarker2,group %in% 1)
df_barMarker2$goterm_name <- factor(c(df_barMarker2$goterm),levels=c(c(df_barMarker2$goterm)[order(df_barMarker2$pval,decreasing=FALSE)]))
go_all_1 <-ggplot(data=df_barMarker2, aes(x=goterm_name, y=pval)) +
    geom_bar(stat="identity",fill=rgb(184,211,116,maxColorValue = 255))+coord_flip()+theme_classic()+theme(text = element_text(family = "Arial", size=15))+labs(x="",y="-log10(P-value)")
	
df_barMarker2 <- data.frame(pval = -log10(GO_all$p.adjust),goterm=GO_all$Description,group=GO_filter$Group)
df_barMarker2$goterm <- as.character(as.matrix(df_barMarker2$goterm))
df_barMarker2 <- subset(df_barMarker2,group %in% 2)
df_barMarker2$goterm_name <- factor(c(df_barMarker2$goterm),levels=c(c(df_barMarker2$goterm)[order(df_barMarker2$pval,decreasing=FALSE)]))
go_all_2 <-ggplot(data=df_barMarker2, aes(x=goterm_name, y=pval)) +
    geom_bar(stat="identity",fill=rgb(0,188,194,maxColorValue = 255))+coord_flip()+theme_classic()+theme(text = element_text(family = "Arial", size=15))+labs(x="",y="-log10(P-value)")



#################################################################
#Heatmap of gene expression
STids_kivswt_exp <- cbind(STids.SCT.ki,STids.SCT.wt)[rownames(DEG_test_kivswt),]
STids_kovswt_exp <- cbind(STids.SCT.ko,STids.SCT.wt)[rownames(DEG_test_kovswt),]
STids_kovswt_exp <- STids_kovswt_exp[apply(STids_kovswt_exp,1,var)>0.1,]
STids_kivswt_exp <- STids_kivswt_exp[apply(STids_kivswt_exp,1,var)>0.1,]

###
sce <- SingleCellExperiment(assays = list(logcounts = STids_kovswt_exp))
sce$batch <- c(rep("KO",ncol(STids.SCT.ko)),rep("WT",ncol(STids.SCT.wt)))
p <- scater::plotHeatmap(sce, features = rownames(STids_kovswt_exp), 
                    center = T, zlim = c(-0.6, 0.6), colour_columns_by = "batch", show_colnames = F,show_rownames = F, 
                    cluster_cols = F, fontsize_row = 6, color = colorRampPalette(c("purple", "black", 
                                                                                   "yellow"))(90))
p <- scater::plotHeatmap(sce[p$tree_row$order,], features = rownames(STids_kovswt_exp)[p$tree_row$order], 
                    center = T, zlim = c(-0.6, 0.6),cluster_rows=FALSE, colour_columns_by = "batch", show_colnames = F,show_rownames = F, 
                    cluster_cols = F, fontsize_row = 6, color = colorRampPalette(c("purple", "black", 
                                                                                   "yellow"))(90))																				   
sce <- SingleCellExperiment(assays = list(logcounts = STids_kivswt_exp))
sce$batch <- c(rep("Ki",ncol(STids.SCT.ki)),rep("WT",ncol(STids.SCT.wt)))
p2 <- scater::plotHeatmap(sce, features = rownames(STids_kivswt_exp), 
                    center = T, zlim = c(-0.6, 0.6), colour_columns_by = "batch", show_colnames = F,show_rownames = F, 
                    cluster_cols = F, fontsize_row = 6, color = colorRampPalette(c("purple", "black", 
                                                                                   "yellow"))(90))
p2 <- scater::plotHeatmap(sce[p2$tree_row$order,], features = rownames(STids_kivswt_exp)[p2$tree_row$order], 
                    center = T, zlim = c(-0.6, 0.6), cluster_rows=FALSE,colour_columns_by = "batch", show_colnames = F,show_rownames = F, 
                    cluster_cols = F, fontsize_row = 6, color = colorRampPalette(c("purple", "black", 
   