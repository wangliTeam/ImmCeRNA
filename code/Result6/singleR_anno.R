singleR_anno <- function(seuractobject){
  # if (!requireNamespace("BiocManager", quietly = TRUE))
  #    install.packages("BiocManager")
  # options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
  # options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
  # BiocManager::install("SingleR")
  library(Seurat)
  library(SingleR)
  library(celldex)

  # load("../HumanPrimaryCellAtlas_hpca.se_human.RData")
  load("../DatabaseImmuneCellExpressionData_dice.se_human.RData")
  load("../BlueprintEncode_bpe.se_human.RData")
  load("../MonacoImmuneData_mid.se_human.RData")
  load("../ImmGenData_igd.se_human.RData")

  #load("BlueprintEncode_bpe.se_human.RData")
  meta <- seuractobject@meta.data
 ## dice
  pbmc.dice <- SingleR(test = GetAssayData(seuractobject, slot="data"), ref = dice.se, 
                       labels = dice.se$label.main,
                       method = "cluster", clusters = meta$seurat_clusters)
  dice_celltype = data.frame(ClusterID=rownames(pbmc.dice), celltype=pbmc.dice$labels, stringsAsFactors = FALSE)
  rownames(dice_celltype) <- dice_celltype$ClusterID
  seuractobject@meta.data[["dice_celltype"]] <- dice_celltype[seuractobject@meta.data[["seurat_clusters"]],"celltype"]
  ## bpe
  pbmc.bpe <- SingleR(test = GetAssayData(seuractobject, slot="data"), ref = bpe.se, 
                      labels = bpe.se$label.main,
                      method = "cluster", clusters = meta$seurat_clusters)
  bpe_celltype = data.frame(ClusterID=rownames(pbmc.bpe), celltype=pbmc.bpe$labels, stringsAsFactors = FALSE)
  rownames(bpe_celltype) <- bpe_celltype$ClusterID
  seuractobject@meta.data[["bpe_celltype"]] <- bpe_celltype[seuractobject@meta.data[["seurat_clusters"]],"celltype"]
  ## mid
  pbmc.mid <- SingleR(test = GetAssayData(seuractobject, slot="data"), ref = mid.se, 
                      labels = mid.se$label.main,
                      method = "cluster", clusters = meta$seurat_clusters)
  mid_celltype = data.frame(ClusterID=rownames(pbmc.mid), celltype=pbmc.mid$labels, stringsAsFactors = FALSE)
  rownames(mid_celltype) <- mid_celltype$ClusterID
  seuractobject@meta.data[["mid_celltype"]] <- mid_celltype[seuractobject@meta.data[["seurat_clusters"]],"celltype"]
  ## igd
  pbmc.igd <- SingleR(test = GetAssayData(seuractobject, slot="data"), ref = igd.se, 
                      labels = igd.se$label.main,
                      method = "cluster", clusters = meta$seurat_clusters)
  igd_celltype = data.frame(ClusterID=rownames(pbmc.igd), celltype=pbmc.igd$labels, stringsAsFactors = FALSE)
  rownames(igd_celltype) <- igd_celltype$ClusterID
  seuractobject@meta.data[["igd_celltype"]] <- igd_celltype[seuractobject@meta.data[["seurat_clusters"]],"celltype"]
  #汇总
  celltypes <- data.frame(seuractobject@meta.data)
  celltypes$cellID <- rownames(celltypes)
  celltypes <- celltypes[,c("cellID","dice_celltype","bpe_celltype","mid_celltype","igd_celltype")]
  return(celltypes)
}