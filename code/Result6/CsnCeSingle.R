
singlecell_path <- "../singlecell/GSE115978"
library(Seurat)
count_data <- read.csv(paste0(singlecell_path,"/GSE115978_counts.csv"),row.names = 1)
dim(count_data)#[1] 23686  7186
anno_data <- read.csv(paste0(singlecell_path,"/GSE115978_cell.annotations.csv"))
rownames(anno_data) <- anno_data$cells
seuractobject <- CreateSeuratObject(counts = count_data, project = "single" ,min.cells = 3, min.features = 200)
seuractobject <- NormalizeData(seuractobject, normalization.method = "LogNormalize", scale.factor = 10000)
seuractobject <- FindVariableFeatures(seuractobject, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seuractobject)
seuractobject <- ScaleData(seuractobject, features = all.genes)
seuractobject <- RunPCA(seuractobject, features = VariableFeatures(object = seuractobject))
ElbowPlot(seuractobject,ndims = 40)
seuractobject <- FindNeighbors(seuractobject, dims = 1:20)
seuractobject <- FindClusters(seuractobject, resolution = 0.8)#
seuractobject <- RunTSNE(seuractobject, dims = 1:20,check_duplicates=FALSE)
seuractobject <- RunUMAP(seuractobject, dims = 1:20,check_duplicates=FALSE)
DimPlot(seuractobject, reduction = "tsne", label = TRUE)
DimPlot(seuractobject, reduction = "umap", label = TRUE)
seuractobject@meta.data[["celltype"]] <- anno_data[rownames(seuractobject@meta.data),"cell.types"]
anno_data$treatment <- ifelse(anno_data$treatment.group=="post.treatment","T","NT")
seuractobject@meta.data[["treatment"]] <- anno_data[rownames(seuractobject@meta.data),"treatment"]
newIdent <- seuractobject@meta.data[["celltype"]]
names(newIdent) <- rownames(seuractobject@meta.data)
seuractobject@active.ident <- as.factor(newIdent)
tsne_umap_data <- as.data.frame(seuractobject@reductions[["tsne"]]@cell.embeddings)
tsne_umap_data <- cbind(tsne_umap_data,as.data.frame(seuractobject@reductions[["umap"]]@cell.embeddings))
tsne_umap_data$celltype <- seuractobject@meta.data[["celltype"]]
tsne_umap_data$treatment <- seuractobject@meta.data[["treatment"]]
tsne_umap_data$dataset <- "GSE115978"
tsne_umap_data$cellID <- rownames(tsne_umap_data)
write.csv(tsne_umap_data,"../singlecell/GSE115978/tsne_umap_data.csv")
save.image("../singlecell/GSE115978/GSE115978.Rdata")

DimPlot(seuractobject, reduction = "tsne", label = TRUE)
DimPlot(seuractobject, reduction = "umap", label = TRUE)
newIdent <- seuractobject@meta.data[["treatment"]]
names(newIdent) <- rownames(seuractobject@meta.data)
seuractobject@active.ident <- as.factor(newIdent)
DimPlot(seuractobject, reduction = "tsne", label = TRUE)
colsNT <- c("#A4BBB6","#E1CE6F")
plot1 <- DimPlot(seuractobject, reduction = "tsne",label = F, pt.size = 1,cols = colsNT)



################################CSN(cell-specific network)################################################################
source("../CSN.R")
immcluster_path <- "../single_immcluster"

for(d in 1:10){

  library(Seurat)
  exp_data <- GetAssayData(seuractobject,slot="data",assay="RNA")
  cancer_ceRNA <- read.csv("../cancer_exp_ceRNA_data.csv",header = T,row.names = 1)
  if(length(intersect(rownames(exp_data),cancer_ceRNA$lncRNA))==0 | length(intersect(rownames(exp_data),cancer_ceRNA$gene))==0) next
  dim(exp_data)#[1]  534 7186
  exp_data <- exp_data[rownames(exp_data) %in% unique(c(cancer_ceRNA$lncRNA,cancer_ceRNA$gene)),]
  csn_data <- CSN(exp_data)
  cancer_ceRNA <- cancer_ceRNA[,c(2,4)]
  #save.image("../singlecell/GSE115978_example/GSE115978_csn_data.Rdata")
  load("../singlecell/GSE115978_example/GSE115978_csn_data.Rdata")
  library(psych)
  library(tidyr)
  library(tidyverse)
  
  cell_regu_data <- data.frame()
  for(c in 1:length(csn_data)){
    temp <- as.data.frame(csn_data[[c]])
    #temp <- as.data.frame(list_mat[[1]])
    temp <- temp[intersect(rownames(temp),cancer_ceRNA$lncRNA),intersect(colnames(temp),cancer_ceRNA$gene)]
    temp$lncRNA <- rownames(temp)
    temp <- gather(temp, key = "gene",value = "exp",-`lncRNA`) 
    temp$ce <- paste(temp$lncRNA,temp$gene,sep = "_")
    temp_data <- data.frame(t(temp$exp))
    colnames(temp_data) <- temp$ce
    cancer_ceRNA$ce <- paste(cancer_ceRNA$lncRNA,cancer_ceRNA$gene,sep = "_")
    if(length(intersect(colnames(temp_data),unique(cancer_ceRNA$ce)))==0) next
    temp_data <- temp_data[,intersect(colnames(temp_data),unique(cancer_ceRNA$ce))]
    cell_regu_data <- rbind(cell_regu_data,temp_data)
    print(c)
  }
  rownames(cell_regu_data) <- colnames(exp_data)
}

cell_regu_data <- read.csv("../singlecell/GSE115978/cell_regu_data.csv",header = T,row.names = 1)
dim(cell_regu_data)#[1] 7186 1322
tsne_umap_data <- cbind(tsne_umap_data,cell_regu_data)
write.csv(tsne_umap_data,"../singlecell/GSE115978/tsne_umap_data.csv")
