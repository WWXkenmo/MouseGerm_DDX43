####################################################################
#Build high-quality normal reference
MouseGerm_499 <- readRDS("MouseGerm_499.rds")
MouseGerm_500 <- readRDS("MouseGerm_500.rds")

MouseGerm <- cbind(MouseGerm_500,MouseGerm_499)
MouseGerm <- CreateSeuratObject(counts = MouseGerm, project = "MouseGerm",min.cell=0.001*ncol(MouseGerm))
MouseGerm$batch <- c(rep("MouseSperm500",ncol(MouseGerm_500)),rep("MouseSperm499",ncol(MouseGerm_499)))
rm(MouseGerm_500,MouseGerm_499)
MouseGerm[["percent.mt"]] <- PercentageFeatureSet(MouseGerm, pattern = "^mt-")

VlnPlot(MouseGerm, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
quantile(MouseGerm@meta.data$nCount_RNA,c(0.025,0.975))
quantile(MouseGerm@meta.data$nFeature_RNA,c(0.025,0.95))
plot(MouseGerm@meta.data$nCount_RNA,MouseGerm@meta.data$nFeature_RNA,pch=16,cex=0.7,bty="n")
abline(h=c(571.00,6081.75),v=c(731.00,38407.25),lty=2,lwd=1,col="red")
MouseGerm <- subset(MouseGerm, subset = nFeature_RNA > 571 & nFeature_RNA < 6081.75 & percent.mt < 10)

MouseGerm<-NormalizeData(MouseGerm)
## Scale Data
MouseGerm <- ScaleData(MouseGerm, features = rownames(MouseGerm))

MouseGerm<-CellCycleScoring(MouseGerm,s.features = cc.genes$s.genes,g2m.features = cc.genes$g2m.genes)
MouseGerm <- FindVariableFeatures(MouseGerm, selection.method ="vst", nfeatures =3000)
MouseGerm <- RunPCA(MouseGerm, verbose=FALSE)
DimPlot(MouseGerm, reduction = "pca",dims = c(1,2),group.by = "Phase")

## Using SCT to regress out the co-variable
MouseGerm<-SCTransform(MouseGerm,vars.to.regress = c("S.Score","G2M.Score","percent.mt","nCount_RNA"))
MouseGerm <- RunPCA(object = MouseGerm, verbose = FALSE)
ElbowPlot(MouseGerm,ndims=30)
MouseGerm <- RunUMAP(object = MouseGerm, dims = 1:20)
MouseGerm <- FindNeighbors(MouseGerm, dims = 1:20,reduction="pca")
MouseGerm <- FindClusters(MouseGerm, resolution = 0.25,algorithm=2)
DimPlot(MouseGerm, reduction = 'umap', label=1)

MouseGerm$batch <- rep("Wild Type",ncol(MouseGerm))
Ica.mousegerm <- ICAcomputing(MouseGerm,ICA.type="JADE",seurat.obj = TRUE,RMT=TRUE)
Ica.filter  <- crossBatchMapping(Ica.mousegerm$ica.pooling,cor="spearman")
MouseGerm <- RunICAnet(MouseGerm,Ica.filter$ica.filter,W.top=2.5,aucMaxRank=600,species = 10090)


#####################################################
## load KO dataset     
      
MouseGerm_366 <- readRDS("MouseGerm_366.rds")         
MouseGerm_368 <- readRDS("MouseGerm_368.rds")

MouseGerm_366 <- CreateSeuratObject(counts = MouseGerm_366, project = "MouseGerm1", min.cells = 0)
MouseGerm_368 <- CreateSeuratObject(counts = MouseGerm_368, project = "MouseGerm2", min.cells = 0)

MouseGerm_368[["percent.mt"]] <- PercentageFeatureSet(MouseGerm_368, pattern = "^mt-")
VlnPlot(MouseGerm_368, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
v <- quantile(MouseGerm_368@meta.data$nCount_RNA,c(0.025,0.975))
h <- quantile(MouseGerm_368@meta.data$nFeature_RNA,c(0.025,0.975))
MouseGerm_368 <- subset(MouseGerm_368, subset = nFeature_RNA > 440 & nFeature_RNA < 3745.8 & percent.mt < 5)

MouseGerm_366[["percent.mt"]] <- PercentageFeatureSet(MouseGerm_366, pattern = "^mt-")
VlnPlot(MouseGerm_366, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
v <- quantile(MouseGerm_366@meta.data$nCount_RNA,c(0.025,0.975))
h <- quantile(MouseGerm_366@meta.data$nFeature_RNA,c(0.025,0.975))
plot(MouseGerm_366@meta.data$nCount_RNA,MouseGerm_366@meta.data$nFeature_RNA,pch=16,cex=0.7,bty="n")
abline(h=h,v=v,lty=2,lwd=1,col="red")
MouseGerm_366 <- subset(MouseGerm_366, subset = nFeature_RNA > 1070 & nFeature_RNA < 6289.7 & percent.mt < 5)

####
MouseGerm_Ki <- cbind(GetAssayData(MouseGerm_366),
                   GetAssayData(MouseGerm_368))
MouseGerm_Ki <- CreateSeuratObject(counts = MouseGerm_Ki, project = "MouseGerm")
MouseGerm_Ki <- MouseGerm_Ki[rownames(MouseGerm),]
MouseGerm_Ki$batch <- c(rep("MouseSperm366",ncol(MouseGerm_366)),rep("MouseSperm368",ncol(MouseGerm_368)))
MouseGerm_Ki <- PercentageFeatureSet(MouseGerm_Ki, pattern = "^mt-", col.name = "percent.mt")
MouseGerm_Ki <- SCTransform(MouseGerm_Ki, vars.to.regress = "percent.mt", verbose = FALSE)
rm(MouseGerm_366,MouseGerm_368);gc()

cells_rankings <- AUCell_buildRankings(GetAssayData(MouseGerm_Ki), nCores = 1)
cells_AUC <- AUCell_calcAUC(MouseGerm@misc$IcaNet_geneSets, cells_rankings, aucMaxRank=600, nCores =3)
MouseGerm_Ki[['IcaNet']] <- CreateAssayObject(as.matrix(getAUC(cells_AUC)))
DefaultAssay(MouseGerm_Ki) <- "IcaNet"

##################################################
##################################################################################
###Load KO          
MouseGerm_340 <- readRDS("MouseGerm_340.rds")
MouseGerm_341 <- readRDS("MouseGerm_341.rds")

MouseGerm_340 <- CreateSeuratObject(counts = MouseGerm_340, project = "MouseGerm1", min.cells = 0)
MouseGerm_341 <- CreateSeuratObject(counts = MouseGerm_341, project = "MouseGerm2", min.cells = 0)

MouseGerm_341[["percent.mt"]] <- PercentageFeatureSet(MouseGerm_341, pattern = "^mt-")
VlnPlot(MouseGerm_341, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
v <- quantile(MouseGerm_341@meta.data$nCount_RNA,c(0.025,0.975))
h <- quantile(MouseGerm_341@meta.data$nFeature_RNA,c(0.025,0.975))
plot(MouseGerm_341@meta.data$nCount_RNA,MouseGerm_341@meta.data$nFeature_RNA,pch=16,cex=0.7,bty="n")
abline(h=h,v=v,lty=2,lwd=1,col="red")
MouseGerm_341 <- subset(MouseGerm_341, subset = nFeature_RNA > 696.125 & nFeature_RNA < 5708.875 & percent.mt < 5)

MouseGerm_340[["percent.mt"]] <- PercentageFeatureSet(MouseGerm_340, pattern = "^mt-")
VlnPlot(MouseGerm_340, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
v <- quantile(MouseGerm_340@meta.data$nCount_RNA,c(0.025,0.975))
h <- quantile(MouseGerm_340@meta.data$nFeature_RNA,c(0.025,0.975))
plot(MouseGerm_340@meta.data$nCount_RNA,MouseGerm_340@meta.data$nFeature_RNA,pch=16,cex=0.7,bty="n")
abline(h=h,v=v,lty=2,lwd=1,col="red")
MouseGerm_340 <- subset(MouseGerm_340, subset = nFeature_RNA > 716 & nFeature_RNA < 6229 & percent.mt < 5)


MouseGerm_KO <- cbind(GetAssayData(MouseGerm_340),
                   GetAssayData(MouseGerm_341))
MouseGerm_KO <- CreateSeuratObject(counts = MouseGerm_KO, project = "MouseGerm")
MouseGerm_KO <- MouseGerm_KO[rownames(MouseGerm),]
MouseGerm_KO$batch <- c(rep("MouseSperm340",ncol(MouseGerm_340)),rep("MouseSperm341",ncol(MouseGerm_341)))
MouseGerm_KO <- PercentageFeatureSet(MouseGerm_KO, pattern = "^mt-", col.name = "percent.mt")
MouseGerm_KO <- SCTransform(MouseGerm_KO, vars.to.regress = "percent.mt", verbose = FALSE)
rm(MouseGerm_340,MouseGerm_341);gc()

cells_rankings <- AUCell_buildRankings(GetAssayData(MouseGerm_KO), nCores = 3)
cells_AUC <- AUCell_calcAUC(MouseGerm@misc$IcaNet_geneSets, cells_rankings, aucMaxRank=600, nCores =6)
MouseGerm_KO[['IcaNet']] <- CreateAssayObject(as.matrix(getAUC(cells_AUC)))
DefaultAssay(MouseGerm_KO) <- "IcaNet"

##integration
DefaultAssay(MouseGerm) <- "IcaNet"
DefaultAssay(MouseGerm_Ki) <- "IcaNet"
MouseGerm_list <- list()
MouseGerm_list[[1]] <- MouseGerm_Ki
MouseGerm_list[[2]] <- MouseGerm_KO
MouseGerm_list[[3]] <- MouseGerm
names(MouseGerm_list) <- c("Ki","KO","Normal")
for (i in 1:length(MouseGerm_list)) {
    VariableFeatures(MouseGerm_list[[i]]) <- rownames(MouseGerm_list[[i]])
}
MouseGerm.anchors <- FindIntegrationAnchors(object.list = MouseGerm_list, anchor.features = rownames(MouseGerm), reference = 3)
MouseGerm.anchors_Ki <- FindTransferAnchors(query = MouseGerm_list[[1]], reference=MouseGerm_list[[3]], features = rownames(MouseGerm))
MouseGerm.anchors_KO <- FindTransferAnchors(query = MouseGerm_list[[2]], reference=MouseGerm_list[[3]], features = rownames(MouseGerm))
predictions_Ki <- TransferData(anchorset = MouseGerm.anchors_Ki, refdata = MouseGerm_list[[3]]$celltype, 
    dims = 1:30)
predictions_KO <- TransferData(anchorset = MouseGerm.anchors_KO, refdata = MouseGerm_list[[3]]$celltype, 
    dims = 1:30)
MouseGerm.integrated <- IntegrateData(anchorset = MouseGerm.anchors)

MouseGerm.integrated <- ScaleData(MouseGerm.integrated)
MouseGerm.integrated <- RunPCA(object = MouseGerm.integrated, verbose = FALSE)
ElbowPlot(MouseGerm.integrated,ndims=30)
MouseGerm.integrated <- RunUMAP(object = MouseGerm.integrated, dims = 1:20)
MouseGerm.integrated <- RunTSNE(object = MouseGerm.integrated, dims = 1:20)


moduleExp <- as.matrix(GetAssayData(MouseGerm.integrated))
norm2scaling <- function(x){
   x <- x-mean(x)
   x <- x/norm(as.matrix(x),"2")
   x
}
moduleExp_scale <- apply(moduleExp,2,norm2scaling)   
svd <- rARPACK::svds((moduleExp_scale),k=30,nu=30,nv=30)
svd.reduction <- svd$v %*% diag(sqrt(svd$d))
rownames(svd.reduction) <- colnames(moduleExp_scale) 
MouseGerm.integrated[['IcaNetSVD']] <- CreateDimReducObject(embeddings =svd.reduction,key="svd_")
MouseGerm.integrated <- RunUMAP(MouseGerm.integrated,reduction="IcaNetSVD",dims=1:20,reduction.name = "IcaNetUMAP")

###Seperate as Ki Ko WT
MouseGerm_Ki <- MouseGerm.integrated[,MouseGerm.integrated$batch %in% c("MouseSperm366","MouseSperm368")]
MouseGerm_Ki <- MouseGerm_Ki[,which(predictions_Ki$prediction.score.max > 0)]
MouseGerm_Ki$celltype <- predictions_Ki$predicted.id[which(predictions_Ki$prediction.score.max > 0)]
MouseGerm_KO <- MouseGerm.integrated[,MouseGerm.integrated$batch %in% c("MouseSperm340","MouseSperm341")]
MouseGerm_KO <- MouseGerm_KO[,which(predictions_KO$prediction.score.max > 0)]
MouseGerm_KO$celltype <- predictions_KO$predicted.id[which(predictions_KO$prediction.score.max > 0)]
MouseGerm_WT <- MouseGerm.integrated[,MouseGerm.integrated$batch %in% c("MouseSperm499","MouseSperm500")]

DimPlot(MouseGerm_WT, reduction = 'IcaNetUMAP', label=1, group.by = 'celltype')
DimPlot(MouseGerm_KO, reduction = 'IcaNetUMAP', label=1, group.by = 'celltype')
DimPlot(MouseGerm_Ki, reduction = 'IcaNetUMAP', label=1, group.by = 'celltype')

tiff("reference_UMAP_1213.tiff",width=2300,height=2000,res=300)
DimPlot(MouseGerm_WT, reduction = 'IcaNetUMAP', label=1, group.by = 'celltype')
dev.off()
tiff("KO_UMAP_1213.tiff",width=2300,height=2000,res=300)
DimPlot(MouseGerm_KO, reduction = 'IcaNetUMAP', label=1, group.by = 'celltype')
dev.off()
tiff("Ki_UMAP_1213.tiff",width=2300,height=2000,res=300)
DimPlot(MouseGerm_Ki, reduction = 'IcaNetUMAP', label=1, group.by = 'celltype')
dev.off()
save(MouseGerm_Ki,MouseGerm_KO,MouseGerm_WT,file="MouseGerm_1213.Rdata")


####################################

png("KO_UMAP.png",width=2300,height=2000,res=300)
DimPlot(MouseGerm_overall[,MouseGerm_overall$batch == "Mapped KO"], reduction = 'umap', label=1, group.by = 'celltype2')
dev.off()

png("WT_UMAP.png",width=2300,height=2000,res=300)
DimPlot(MouseGerm_overall[,MouseGerm_overall$batch == "Reference"], reduction = 'umap', label=1, group.by = 'celltype2')
dev.off()

png("Mapping_UMAP.png",width=2300,height=2000,res=300)
DimPlot(MouseGerm_overall, reduction = 'umap', label=1, group.by = 'batch')
dev.off()

png("KI_UMAP.png",width=2300,height=2000,res=300)
DimPlot(MouseGerm_Ki, reduction = 'umap', label=1)
dev.off()

###########################plot the result#########################
MouseGerm_overall <- CreateSeuratObject(cbind(GetAssayData(MouseGerm_WT),GetAssayData(MouseGerm_KO)))
UMAP <- rbind(MouseGerm_WT@reductions$umap@cell.embeddings,MouseGerm_KO@reductions$umap@cell.embeddings)
MouseGerm_overall[['umap']] <- CreateDimReducObject(UMAP)
MouseGerm_overall$celltype <- c(MouseGerm_WT$celltype, MouseGerm_KO$celltype)
celltype <- MouseGerm_overall$celltype
MouseGerm_overall$celltype <- celltype
MouseGerm_overall$batch <- c(rep("Reference",ncol(MouseGerm_WT)),rep("Mapped KO",ncol(MouseGerm_KO)))
MouseGerm_WT$celltype2 <- celltype

tiff("mapping.tiff",res=300,height=1500,width=2300)
DimPlot(MouseGerm_overall,group.by = "batch",reduction="umap")+labs(x="UMAP-1",y="UMAP-2")
dev.off()

png("reference_umap.png",res=300,height=1500,width=2000)
DimPlot(MouseGerm_WT_all,group.by = "celltype2",reduction="umap",label=1,label.size = 4.2)+NoLegend()+labs(x="UMAP-1",y="UMAP-2")+theme(text = element_text(size=15))+FontSize(x.text=15,y.text=15)
dev.off()

Marker <- FindAllMarkers(MouseGerm_WT_all,only.pos=TRUE,min.diff.pct = 0.2,logfc.threshold=0.5)
###Calculate Average Expression Value
FS <- names(table(Marker$gene)[table(Marker$gene)==1])
Marker.f <- subset(Marker,gene %in% FS)
FS_list <- NULL;for(i in c("Unknown","Macrophage","Endothelial","SPG","SCytes1","SCytes2","STids","Elongating1","Elongating2","Elongating3")){
   FS_list <- c(FS_list, rownames(subset(Marker.f,cluster %in% i)))
}
Exp <- AverageExpression(MouseGerm_WT_all,assays = "SCT")
Exp <- Exp$SCT
Exp <- Exp[FS_list,]
Exp <- Exp[,c("Unknown","Macrophage","Endothelial","SPG","SCytes1","SCytes2","STids","Elongating1","Elongating2","Elongating3")]

png("marker_heatmap.png",res=300,height=2500,width=1000)
pheatmap(Exp,cluster_rows = FALSE,cluster_cols = FALSE,scale="row",show_rownames=FALSE)
dev.off()


################################################
## calculate correlation matrix and return
MII <- read.csv("MII_gset.csv")
RS1o2 <- read.csv("RS1o2.csv")
TangSingleCell_Seurat <- readRDS("TangSingleCell_Seurat.rds")
TangSingleCell_Seurat <- NormalizeData(TangSingleCell_Seurat)
Exp <- AverageExpression(MouseGerm_WT,assays = "SCT")
Exp <- Exp$SCT
refExp <- AverageExpression(TangSingleCell_Seurat,assays = "RNA")
refExp <- refExp$RNA
Exp <- Exp[rowMeans(Exp)>0.5,]
refExp <- refExp[rowMeans(refExp)>0.5,]
g <- intersect(rownames(Exp), rownames(refExp))
Exp <- Exp[g,]
refExp <- refExp[g,]

corMat <- cor(log2(Exp+1),log2(refExp+1))
corMat[c("SPG","SCytes1","SCytes2","STids","Elongating1","Elongating2","Elongating3"),] -> corMat
corMat <- corMat[,c("A1","ln","TypeBS","TypeBG2M","G1","ePL","mPL","lPL","L","Z","eP","mP","lP","D","MI","MII","RS1o2","RS3o4","RS5o6","RS7o8")]
#png("cor_heatmap_WT.png",res=300,height=1000,width=2500)
pheatmap::pheatmap(corMat,cluster_rows = FALSE,cluster_cols = FALSE)
#dev.off()


Exp <- AverageExpression(MouseGerm_KO,assays = "SCT")
Exp <- Exp$SCT
refExp <- AverageExpression(TangSingleCell_Seurat,assays = "RNA")
refExp <- refExp$RNA
Exp <- Exp[rowMeans(Exp)>0.5,]
refExp <- refExp[rowMeans(refExp)>0.5,]
g <- intersect(rownames(Exp), rownames(refExp))
Exp <- Exp[g,]
refExp <- refExp[g,]

corMat <- cor(log2(Exp+1),log2(refExp+1))
corMat[c("SPG","SCytes1","SCytes2","STids","Elongating1","Elongating2","Elongating3"),] -> corMat
corMat <- corMat[,c("A1","ln","TypeBS","TypeBG2M","G1","ePL","mPL","lPL","L","Z","eP","mP","lP","D","MI","MII","RS1o2","RS3o4","RS5o6","RS7o8")]
png("cor_heatmap_KO.png",res=300,height=1000,width=2500)
pheatmap::pheatmap(corMat,cluster_rows = FALSE,cluster_cols = FALSE)
dev.off()

cells_rankings <- AUCell_buildRankings(as.matrix(GetAssayData(MouseGerm_KO)), nCores = 6)
cells_AUC <- AUCell_calcAUC(list(MII=intersect(MII$Gene,rownames(MouseGerm_KO)),
                                 RS2=intersect(RS1o2$Gene,rownames(MouseGerm_KO))), cells_rankings, aucMaxRank=3000, nCores =6)
cells_AUC_matrix <- getAUC(cells_AUC)

p1 <- Boxplot(MouseGerm_KO,"RS1o2",idents = rev(c("SCytes2","STids")))
p2 <- Boxplot(MouseGerm_KO,"MII",idents = rev(c("SCytes2","STids")))
png("MII_vs_RS2(KO).png",res=300,height=1500,width=2200)
cowplot::plot_grid(p2,p1,nrow=1)
dev.off()
####################################################################
Exp <- AverageExpression(MouseGerm_Ki,assays = "SCT")
Exp <- Exp$SCT
refExp <- AverageExpression(TangSingleCell_Seurat,assays = "RNA")
refExp <- refExp$RNA
Exp <- Exp[rowMeans(Exp)>0.5,]
refExp <- refExp[rowMeans(refExp)>0.5,]
g <- intersect(rownames(Exp), rownames(refExp))
Exp <- Exp[g,]
refExp <- refExp[g,]

corMat <- cor(log2(Exp+1),log2(refExp+1))
corMat[c("SPG","SCytes1","SCytes2","STids","Elongating1","Elongating2","Elongating3"),] -> corMat
corMat <- corMat[,c("A1","ln","TypeBS","TypeBG2M","G1","ePL","mPL","lPL","L","Z","eP","mP","lP","D","MI","MII","RS1o2","RS3o4","RS5o6","RS7o8")]
png("cor_heatmap_Ki.png",res=300,height=1000,width=2500)
pheatmap::pheatmap(corMat,cluster_rows = FALSE,cluster_cols = FALSE)
dev.off()

cells_rankings <- AUCell_buildRankings(as.matrix(GetAssayData(MouseGerm_Ki)), nCores = 6)
cells_AUC <- AUCell_calcAUC(list(MII=intersect(MII$Gene,rownames(MouseGerm_KO)),
                                 RS2=intersect(RS1o2$Gene,rownames(MouseGerm_KO))), cells_rankings, aucMaxRank=5000, nCores =6)
cells_AUC_matrix <- getAUC(cells_AUC)
p1 <- Boxplot(MouseGerm_Ki,"RS1o2",idents = rev(c("SCytes2","STids")))
p2 <- Boxplot(MouseGerm_Ki,"MII",idents = rev(c("SCytes2","STids")))
png("MII_vs_RS2(Ki).png",res=300,height=1500,width=2200)
cowplot::plot_grid(p2,p1,nrow=1)
dev.off()


Boxplot <- function(object,feature,idents = NULL){
   ## 
   require("ggplot2")
   
   cell_type <- as.character(as.matrix(Idents(object)))
   feature_vec <- object@meta.data[,feature][cell_type %in% idents]
   cell_type <- cell_type[cell_type %in% idents]
   
   datm <- data.frame(feature = feature_vec, celltype = cell_type)
   datm$celltype <- as.factor(datm$celltype)
   datm$celltype <- factor(datm$celltype, levels = idents)
      
   p<-ggplot(datm, aes(x=celltype, y=feature, fill=celltype)) +
    geom_boxplot() + theme_classic() + theme(text = element_text(size=15)) + labs(y = paste(feature," gene set enrichment score",sep=""))
   p
}

Boxplot(MouseGerm_Ki,"RS1o2",idents = rev(c("SCytes2","STids")))

###################################
gene <- intersect(rownames(MouseGerm_WT), rownames(MouseGerm_KO))
PCA <- rbind(MouseGerm_WT@reductions$pca@cell.embeddings,MouseGerm_KO@reductions$pca@cell.embeddings)
UMAP <- rbind(MouseGerm_WT@reductions$umap@cell.embeddings,MouseGerm_KO@reductions$umap@cell.embeddings)
MouseGerm_overall <- CreateSeuratObject(cbind(GetAssayData(MouseGerm_WT)[gene,],GetAssayData(MouseGerm_KO)[gene,]))
MouseGerm_overall$sample <- c(MouseGerm_WT$batch,MouseGerm_KO$batch)
MouseGerm_overall$celltype <- c(as.character(as.matrix(Idents(MouseGerm_WT))),as.character(as.matrix(Idents(MouseGerm_KO))))
MouseGerm_overall$batch <- c(rep("Reference",ncol(MouseGerm_WT)),rep("Mapped KO",ncol(MouseGerm_KO)))


##Create SingleCellExperiment Object
sce <- SingleCellExperiment(assays=list(counts=GetAssayData(MouseGerm_overall)),
    reducedDims=SimpleList(PCA=PCA, UMAP=UMAP))

library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(patchwork)

#Create Milo Object
sce$batch <- MouseGerm_overall$batch
sce$sample <- MouseGerm_overall$sample
sce$celltype <- MouseGerm_overall$celltype2
sce[,colnames(MouseGerm)] -> sce
milo <- Milo(sce)
milo <- buildGraph(milo, k = 30, d = 20, reduced.dim = "PCA")
milo <- makeNhoods(milo, prop = 0.1, k = 30, d=20, refined = TRUE, reduced_dims = "PCA")
p1 <- plotReducedDim(sce, colour_by="batch", dimred = "UMAP") 
p2 <- plotNhoodSizeHist(milo)

milo <- countCells(milo, meta.data = as.data.frame(colData(milo)), sample="sample")

design <- data.frame(colData(milo))[,c("sample", "batch")]

## Convert batch info from integer to factor
design$batch <- as.factor(design$batch) 
design$sample <- as.factor(design$sample) 
design <- dplyr::distinct(design)
rownames(design) <- design$sample

milo <- calcNhoodDistance(milo, d=20, reduced.dim = "PCA")
da_results <- testNhoods(milo, design = ~ batch, design.df = design)

milo <- buildNhoodGraph(milo)

## Plot single-cell UMAP
umap_pl <- plotReducedDim(milo, dimred = "PCA", colour_by="batch", text_by = "celltype", 
                          text_size = 3, point_size=0.5) + theme(text = element_text(size=15)) + 
  guides(fill="none")

## Plot neighbourhood graph
nh_graph_pl <- plotNhoodGraphDA(milo, da_results, layout="PCA",alpha=0.05) + theme(text = element_text(family = "Arial",size=15))

png("milo_pca.png",res=300,height=1800,width=2500)  
umap_pl	 + nh_graph_pl + plot_layout(guides="collect")
dev.off()
########################################################################
PCA <- rbind(MouseGerm_WT@reductions$pca@cell.embeddings,MouseGerm_Ki@reductions$pca@cell.embeddings)
UMAP <- rbind(MouseGerm_WT@reductions$umap@cell.embeddings,MouseGerm_Ki@reductions$umap@cell.embeddings)
MouseGerm_overall <- CreateSeuratObject(cbind(GetAssayData(MouseGerm_WT),GetAssayData(MouseGerm_Ki)))
MouseGerm_overall$sample <- c(MouseGerm_WT$batch,MouseGerm_Ki$batch)
MouseGerm_overall$celltype <- c(as.character(as.matrix(Idents(MouseGerm_WT))),as.character(as.matrix(Idents(MouseGerm_Ki))))
MouseGerm_overall$batch <- c(rep("Reference",ncol(MouseGerm_WT)),rep("Mapped KO",ncol(MouseGerm_Ki)))

##Create SingleCellExperiment Object
sce <- SingleCellExperiment(assays=list(counts=GetAssayData(MouseGerm_overall)),
    reducedDims=SimpleList(PCA=PCA, UMAP=UMAP))
sce$batch <- MouseGerm_overall$batch
sce$sample <- MouseGerm_overall$sample
sce$celltype <- MouseGerm_overall$celltype

milo <- Milo(sce)
milo <- buildGraph(milo, k = 30, d = 20, reduced.dim = "PCA")
milo <- makeNhoods(milo, prop = 0.1, k = 30, d=20, refined = TRUE, reduced_dims = "PCA")
p1 <- plotReducedDim(sce, colour_by="batch", dimred = "UMAP") 
p2 <- plotNhoodSizeHist(milo)

milo <- countCells(milo, meta.data = as.data.frame(colData(milo)), sample="sample")

design <- data.frame(colData(milo))[,c("sample", "batch")]

## Convert batch info from integer to factor
design$batch <- as.factor(design$batch) 
design$sample <- as.factor(design$sample) 
design <- dplyr::distinct(design)
rownames(design) <- design$sample

milo <- calcNhoodDistance(milo, d=20, reduced.dim = "PCA")
da_results <- testNhoods(milo, design = ~ batch, design.df = design)
da_results$SpatialFDR <- da_results$PValue
milo <- buildNhoodGraph(milo)

## Plot single-cell UMAP
umap_pl <- plotReducedDim(milo, dimred = "PCA", colour_by="batch", text_by = "celltype", 
                          text_size = 3, point_size=0.5) + theme(text = element_text(size=15)) + 
  guides(fill="none")

## Plot neighbourhood graph
nh_graph_pl <- plotNhoodGraphDA(milo, da_results, layout="PCA",alpha=0.05) + theme(text = element_text(family = "Arial",size=15))

png("milo_pca(Ki).png",res=300,height=1800,width=2500)  
umap_pl	 + nh_graph_pl + patchwork::plot_layout(guides="collect")
dev.off()

########################################################################
library(SCENT)
library(RSCORE)
library(reticulate)

### calculate WT
PPI_mousegerm <- getPPI_String2(MouseGerm,species=10090)
score <- SCENT::CompCCAT(GetAssayData(MouseGerm),PPI_mousegerm)
MouseGerm$SCENT_score <- score

dat <- MouseGerm@meta.data[,c("batch","celltype","SCENT_score")]
dat$celltype <- factor(dat$celltype, levels=c("SPG","SCytes1","SCytes2","STids","Elongating3","Elongating2","Elongating1"))
p<-ggplot(dat, aes(x=celltype, y=SCENT_score, fill=batch)) +
  geom_boxplot(position=position_dodge(1))
  
### calculate KO
PPI_mousegerm <- getPPI_String2(MouseGerm_KO,species=10090)
score <- SCENT::CompCCAT(GetAssayData(MouseGerm_KO),PPI_mousegerm)
MouseGerm_KO$SCENT_score <- score

dat <- MouseGerm_KO@meta.data[,c("batch","celltype","SCENT_score")]
dat$celltype <- factor(dat$celltype, levels=c("SPG","SCytes1","SCytes2","STids","Elongating3","Elongating2","Elongating1"))
p<-ggplot(dat, aes(x=celltype, y=SCENT_score, fill=batch)) +
  geom_boxplot(position=position_dodge(1))
  
### calculate Ki
PPI_mousegerm <- getPPI_String2(MouseGerm_Ki,species=10090)
score <- SCENT::CompCCAT(GetAssayData(MouseGerm_Ki),PPI_mousegerm)
MouseGerm_Ki$SCENT_score <- score

dat <- MouseGerm_Ki@meta.data[,c("batch","celltype","SCENT_score")]
dat$celltype <- factor(dat$celltype, levels=c("SPG","SCytes1","SCytes2","STids","Elongating3","Elongating2","Elongating1"))
p<-ggplot(dat, aes(x=celltype, y=SCENT_score, fill=batch)) +
  geom_boxplot(position=position_dodge(1)) 


###############################################################################
### Uncover cellular identity through GO enrichment analysis

p1 <- FeaturePlot(MouseGerm_WT_all,"Pdgfra")
p2 <- FeaturePlot(MouseGerm_WT_all,"Tcf21")
p3 <- FeaturePlot(MouseGerm_WT_all,"Ly6a")
p4 <- FeaturePlot(MouseGerm_WT_all,"Thra")
png("unknow_celltype_marker.png",res=300,height=2000,width=2000)
cowplot::plot_grid(p1,p2,p3,p4,nrow=2)
dev.off()

markers <- FindMarkers(MouseGerm_WT_all,ident.1 = "Unknown",only.pos=TRUE,min.diff.pct = 0.2,logfc.threshold=0.5)
background <- rownames(MouseGerm_WT_all)
background <- bitr(background,fromType="SYMBOL",toType="ENTREZID",org.Mm.eg.db)
markers <- markers[markers[,"p_val_adj"]<0.05,]

markers <- bitr(rownames(markers),fromType="SYMBOL",toType="ENTREZID",org.Mm.eg.db)
goterm <- enrichGO(names(table(c(markers[,2]))), ont="BP",pvalueCutoff = 0.05,qvalueCutoff = 0.05,OrgDb = org.Mm.eg.db,universe=background[,2],readable=TRUE)

png("GO_unknown_marker.png",res=300,height=1500,width=1000)
barplot(goterm,showCategory = 20)
dev.off()
























