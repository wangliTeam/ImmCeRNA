#

library(monocle)
packageVersion("monocle")  
#install.packages('SeuratObject')
library(Seurat)
input_path <- "../inputData/impy/single_cell/"
output_path <- "../output/output_Traject/"
dataset <- read.csv(paste0(input_path,"sc_dats.csv"),header = T,row.names = 1)
dataset <- dataset$dats
Imm_celltypes <- c("T.CD8" ,"T.CD4","T.cell","T_cells" ,"NK" ,"NK_cell" ,"B_cell" ,"B.cell" ,"CD4","CD8")
Imm_cell_uni <- data.frame(row.names =  Imm_celltypes, 
                           uni_celtype=c("CD8 T cells" ,"CD4 T cells","T cells","T cells" ,"NK cells" ,"NK cells" ,"B cells" ,"B cells" ,"CD4 T cells","CD8 T cells" ))

for(d in 1:length(dataset)){#
  print(dataset[d])
  sc_data <- load(paste0(input_path,dataset[d],"_scexp.Rdata"))
  temp_celltypes <- intersect(unique(seuractobject@meta.data$celltype),Imm_celltypes)
  if(length(temp_celltypes)==0) next
  print(temp_celltypes)
  if(length(temp_celltypes) > 1){
    for(i in 1:length(temp_celltypes)){
      select_temp_celltypes <- temp_celltypes[i]
      print(dataset[d])
      print(select_temp_celltypes)
      pd <- data.frame(celltype = seuractobject@meta.data$celltype,row.names = colnames(seuractobject))
      pd <- pd[pd$celltype==select_temp_celltypes,, drop = FALSE]
      matrix <- GetAssayData(seuractobject,layer="data",assay="RNA")
      geneid <-gsub('[_]','-',row.names(matrix))
      row.names(matrix)<-geneid
      data <- matrix[,rownames(pd)]
      fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
      #
      pd <- new("AnnotatedDataFrame", data = pd)
      fd <- new("AnnotatedDataFrame", data = fData)
      
      library(monocle)
      library("dplyr",warn.conflicts = FALSE)
      monocle_cds <- newCellDataSet(data,
                                    phenoData = pd,
                                    featureData = fd,
                                    lowerDetectionLimit = 0.5,
                                    expressionFamily = negbinomial.size())
      head(pData(monocle_cds))
      monocle_cds <- estimateSizeFactors(monocle_cds)
      monocle_cds <- estimateDispersions(monocle_cds)
      monocle_cds <- detectGenes(monocle_cds, min_expr = 3 )
      print(head(fData(monocle_cds)))
      
      expressed_genes <- row.names(subset(fData(monocle_cds), num_cells_expressed >= 20))
      diff_test_res <- dispersionTable(monocle_cds)
      ordering_genes <- subset(diff_test_res, mean_expression >= 0.1  & dispersion_empirical >= 1 * dispersion_fit)$gene_id
      #write.table(ordering_genes,paste0(output_path,dataset[i],"_ordering_genes.csv"))
      monocle_cds <-setOrderingFilter(monocle_cds,ordering_genes)
      monocle_cds <- reduceDimension(monocle_cds, max_components = 2, reduction_method = "DDRTree")
      monocle_cds <- orderCells(monocle_cds)
      colnames(pData(monocle_cds))
      # [1] "celltype"            "Size_Factor"         "num_genes_expressed"
      # [4] "Pseudotime"          "State"
      reduced_coords <- as.data.frame(t(reducedDimS(monocle_cds)))
      colnames(reduced_coords) <- c("x","y")
      pseudotime <- pData(monocle_cds)$Pseudotime
      CellTraject_data <- cbind(data.frame(cell_id = rownames(pData(monocle_cds)),pseudotime = pseudotime,row.names=rownames(pData(monocle_cds))),
                                reduced_coords)
      temp_unicelltype <- unique(Imm_cell_uni[select_temp_celltypes,"uni_celtype"])
      save(monocle_cds,file = paste0(output_path,dataset[d],"_",temp_unicelltype,"_monoTraject.Rdata"))
      tsne_umap_data <- read.csv(paste0("../output/output_impy/singlecell/",dataset[d],"_tsne_umap_data.csv"),header=T,row.names=1,check.names = F)
      uni_cellID <- intersect(rownames(tsne_umap_data),rownames(CellTraject_data))
      tsne_umap_data <- tsne_umap_data[uni_cellID,]
      CellTraject_data <- CellTraject_data[uni_cellID,]
      ceCellTraject_data <- cbind(CellTraject_data,tsne_umap_data[,10:dim(tsne_umap_data)[2]])
      write.csv(ceCellTraject_data,paste0(output_path,dataset[d],"_",temp_unicelltype,"_ceCellTraject_data.csv"))
      
    }
  }else{
    print(dataset[d])
    print(temp_celltypes)
    pd <- data.frame(celltype = seuractobject@meta.data$celltype,row.names = colnames(seuractobject))
    pd <- pd[pd$celltype==temp_celltypes,, drop = FALSE]
    matrix <- GetAssayData(seuractobject,layer="data",assay="RNA")
    geneid <-gsub('[_]','-',row.names(matrix))
    row.names(matrix)<-geneid
    data <- matrix[,rownames(pd)]
    fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
    #
    pd <- new("AnnotatedDataFrame", data = pd)
    fd <- new("AnnotatedDataFrame", data = fData)
    
    library(monocle)
    library("dplyr",warn.conflicts = FALSE)
    monocle_cds <- newCellDataSet(data,
                                  phenoData = pd,
                                  featureData = fd,
                                  lowerDetectionLimit = 0.5,
                                  expressionFamily = negbinomial.size())
    head(pData(monocle_cds))
    monocle_cds <- estimateSizeFactors(monocle_cds)
    monocle_cds <- estimateDispersions(monocle_cds)
    monocle_cds <- detectGenes(monocle_cds, min_expr = 3 )
    print(head(fData(monocle_cds)))
    
    expressed_genes <- row.names(subset(fData(monocle_cds), num_cells_expressed >= 20))
    diff_test_res <- dispersionTable(monocle_cds)
    ordering_genes <- subset(diff_test_res, mean_expression >= 0.1  & dispersion_empirical >= 1 * dispersion_fit)$gene_id
    #write.table(ordering_genes,paste0(output_path,dataset[i],"_ordering_genes.csv"))
    
    monocle_cds <-setOrderingFilter(monocle_cds,ordering_genes)
    monocle_cds <- reduceDimension(monocle_cds, max_components = 2, reduction_method = "DDRTree")
    monocle_cds <- orderCells(monocle_cds)
    colnames(pData(monocle_cds))
    # [1] "celltype"            "Size_Factor"         "num_genes_expressed"
    # [4] "Pseudotime"          "State"

    reduced_coords <- as.data.frame(t(reducedDimS(monocle_cds)))
    colnames(reduced_coords) <- c("x","y")
    pseudotime <- pData(monocle_cds)$Pseudotime
    CellTraject_data <- cbind(data.frame(cell_id = rownames(pData(monocle_cds)),pseudotime = pseudotime,row.names=rownames(pData(monocle_cds))),
                              reduced_coords)
    temp_unicelltype <- unique(Imm_cell_uni[select_temp_celltypes,"uni_celtype"])
    save(monocle_cds,file = paste0(output_path,dataset[d],"_",temp_unicelltype,"_monoTraject.Rdata"))
    tsne_umap_data <- read.csv(paste0("../output/output_impy/singlecell/",dataset[d],"_tsne_umap_data.csv"),header=T,row.names=1)
    uni_cellID <- intersect(rownames(tsne_umap_data),rownames(CellTraject_data))
    tsne_umap_data <- tsne_umap_data[uni_cellID,]
    CellTraject_data <- CellTraject_data[uni_cellID,]
    ceCellTraject_data <- cbind(CellTraject_data,tsne_umap_data[,10:dim(tsne_umap_data)[2]])
    write.csv(ceCellTraject_data,paste0(output_path,dataset[d],"_",temp_unicelltype,"_ceCellTraject_data.csv"))
  }
}
