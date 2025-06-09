library(TCGAbiolinks)
library(SummarizedExperiment)
library(scRNAseq)
library(data.table)
library(limma)
library(dplyr)
library(DT)
library(psych)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(ggalluvial)
library(tidyverse)
library(igraph)
library(ggraph)
library(tidygraph)
#install.packages("extrafont")
library(extrafont)
#install.packages("showtext")
library(showtext)
#install.packages("ggpubr")
library(ggpubr)
#install.packages('ggorrplot')
#library(ggorrplot)
library(reshape2)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
#library(ggnewscale)
#library(topGO)
#library(pathview)
#require(scales)
#library(ggtree)
#library(aplot)
library(survival)
library(ggplot2)
library(ggpubr)
library(survminer)
showtext_auto()

cancer <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","HNSC","KICH",
            "KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
            "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")

for(c in 1:length(cancer)){
  input_data_path <- "../TCGA_33cancer_data"
  output_imcepair_path <- "../output/output_imcepair"
  output_list <- "../output/output_list"
  #####################################################################################################################
  ############################################# ######################################################
  mylist <- list( ImCePair="null")
  setwd(input_data_path)
  
  expData <- load(paste(cancer[c],"_exp.rda",sep=""))
  tpm_data <-as.data.frame(assay(data,i = "tpm_unstrand")) 
  dim(tpm_data)
  rowNames <-rownames(tpm_data)
  pattern <- "^ENSG[0-9]+\\.[0-9]+$"
  pa <- grepl(pattern,rownames(tpm_data),ignore.case ="True")
  table(pa)
  rownames<-rowNames[pa]
  tpm_data<-tpm_data[rownames,]
  dim(tpm_data)
  

  gene_ano <- data.frame(gene_id=data@rowRanges@elementMetadata@listData[["gene_id"]],gene_type=data@rowRanges@elementMetadata@listData[["gene_type"]],gene_name=data@rowRanges@elementMetadata@listData[["gene_name"]])
  PCG_names <- gene_ano[gene_ano$gene_type=="protein_coding",]
  lnc_names <- gene_ano[gene_ano$gene_type=="lncRNA",]
  miR_names <- gene_ano[gene_ano$gene_type=="miRNA",]
  
  PCG_tpm_data <- tpm_data[rownames(tpm_data) %in% PCG_names$gene_id,]
  dim(PCG_tpm_data)
  #[1] 19944  1226
  PCG_symbol_names <- PCG_names[PCG_names$gene_id %in% rownames(PCG_tpm_data),]$gene_name
  PCG_tpm_data_s <- cbind(PCG_symbol_names,PCG_tpm_data)
  #length(table(PCG_symbol_names))  [1] 19938
  #length(PCG_symbol_names)  [1] 19944
  PCG_tpm_exprmean <- aggregate(.~PCG_symbol_names,mean,data=PCG_tpm_data_s)
  dim(PCG_tpm_exprmean)
  #[1] 19938  1227
  
  lnc_tpm_data <- tpm_data[rownames(tpm_data) %in% lnc_names$gene_id,]
  dim(lnc_tpm_data)
  #[1] 16889  1226
  lnc_symbol_names <- lnc_names[lnc_names$gene_id %in% rownames(lnc_tpm_data),]$gene_name
  lnc_tpm_data_s <- cbind(lnc_symbol_names,lnc_tpm_data)
  length(unique(lnc_symbol_names))
  #[1] 16882
  lnc_tpm_exprmean <- aggregate(.~lnc_symbol_names,mean,data=lnc_tpm_data_s)
  dim(lnc_tpm_exprmean)
  #[1] 16882  1227
  
  miR_tpm_data <- tpm_data[rownames(tpm_data) %in% miR_names$gene_id,]
  dim(miR_tpm_data)
  #[1] 1879 1226
  miR_symbol_names <- miR_names[miR_names$gene_id %in% rownames(miR_tpm_data),]$gene_name
  miR_tpm_data_s <- cbind(miR_symbol_names,miR_tpm_data)
  length(unique(miR_symbol_names))
  #[1] 1879
  miR_tpm_exprmean <- aggregate(.~miR_symbol_names,mean,data=miR_tpm_data_s)
  lnc_exp <- lnc_tpm_exprmean[,-1]
  rownames(lnc_exp) <- lnc_tpm_exprmean[,1]
  dim(lnc_exp)
  #[1] 16882  1226

  PCG_exp <-PCG_tpm_exprmean[,-1]
  rownames(PCG_exp) <- PCG_tpm_exprmean[,1]
  dim(PCG_exp)

  miR_exp <-miR_tpm_exprmean[,-1]
  rownames(miR_exp) <- miR_tpm_exprmean[,1]
  dim(miR_exp)
  head(rownames(lnc_exp))
  head(rownames(PCG_exp))
  head(rownames(miR_exp))

  setwd(input_data_path)
  lnc_Pathway_All <-read.table("ImmReg_lncRNA_Pathway_sig.txt",header = T,sep = "\t")
  IRG_All <-read.csv("Immport_IRG_pathway.csv",header = T,sep = ",")
  miR_Pathway_All <-read.table("ImmReg_miR_Pathway_sig.txt",header = T,sep = "\t")
  MIR_convert <- read.csv("miRIdCon_finally.csv",header = T,row.names = 1)
  head(rownames(lnc_exp))#"A1BG-AS1"    "A2M-AS1"     "A2ML1-AS1"   "A2ML1-AS2"   "AA06"        "AADACL2-AS1"
  head(rownames(PCG_exp))# "A1BG"    "A1CF"    "A2M"     "A2ML1"   "A3GALT2" "A4GALT" 
  head(rownames(miR_exp))#[1] "MIR494"    "MIRLET7E"  "MIR375"    "MIR30E"    "MIRLET7A2" "MIR429"  

  lnc_Pathway_All_uni <-  lnc_Pathway_All %>% dplyr :: distinct(Cancer,lncRNA.Symbol) 
  miR_Pathway_All_uni <-  miR_Pathway_All %>% dplyr :: distinct(Cancer,miRNA.Symbol) 
  head(miR_Pathway_All_uni$miRNA.Symbol)

  imlnc_canceri_names <- lnc_Pathway_All_uni[lnc_Pathway_All_uni$Cancer %in% cancer[c],]$lncRNA.Symbol
  length(imlnc_canceri_names) #3800
  imlnc_canceri_exp <- lnc_exp[rownames(lnc_exp) %in% imlnc_canceri_names,]
  dim(imlnc_canceri_exp)
  #[1] 3238 1226
  IRG_canceri_exp <-PCG_exp[rownames(PCG_exp) %in% IRG_All[,1],]
  dim(IRG_canceri_exp)
  #[1] 1361 1226
  immiR_canceri_names <- miR_Pathway_All_uni[miR_Pathway_All_uni$Cancer %in% cancer[c],]$miRNA.Symbol
  
  immiR_canceri_exp_names <- MIR_convert[MIR_convert$hsa %in% immiR_canceri_names,]$MIR
  length(unique(immiR_canceri_exp_names)) 
  immiR_canceri_exp <- miR_exp[rownames(miR_exp) %in% immiR_canceri_exp_names,]
  dim(immiR_canceri_exp)
  #[1] 416 1207
  #corr(lnc,gene)> 0 P<0.05
  #corr(lnc,miR)< 0 P<0.05
  #corr(miR,gene)< 0 P<0.05
  dim(imlnc_canceri_exp)
  #[1] 3238 1226
  dim(IRG_canceri_exp)
  #[1] 1361 1226
  dim(immiR_canceri_exp)
  #[1] 416 1207
  imlnc_exp_t <-as.data.frame(t(imlnc_canceri_exp))
  IRG_exp_t <- as.data.frame(t(IRG_canceri_exp))
  immiR_exp_t <- as.data.frame(t(immiR_canceri_exp))

  
  ############################
  #corr(lnc,gene) >0  P<0.05
  corr_imlnc_IRG <- corr.test(imlnc_exp_t,IRG_exp_t,use="pairwise",method="pearson",adjust = "none",ci=FALSE)
  #######################################
  corr_imlnc_IRG_r <- as.data.frame(corr_imlnc_IRG$r)
  corr_imlnc_IRG_p <- as.data.frame(corr_imlnc_IRG$p)
  dim(corr_imlnc_IRG_r)
  #[1] [1] 3238 1361
  length(unique(rownames(corr_imlnc_IRG_r)))
  length(unique(colnames(corr_imlnc_IRG_r)))
  # colnames(corr_imlnc_IRG_r) ——gene
  # rownames(corr_imlnc_IRG_r) ——lncRNA
  max(corr_imlnc_IRG_r) 
  corr_imlnc_IRG_r[is.na(corr_imlnc_IRG_r)] <-  0
  max(corr_imlnc_IRG_r) 
  min(corr_imlnc_IRG_r) 

  imlnc_IRG_data <- data.frame(lncRNA="",gene="",corr_r="")
  temp_data <-data.frame()
  temp_row_data <- data.frame()
  for(i in 1:dim(corr_imlnc_IRG_r)[1]){
    temp_row_data <- corr_imlnc_IRG_r[i,]#[1] "data.frame"
    if(!any(temp_row_data >0)) next
    temp_r <- which(temp_row_data > 0)
    temp_lnc_name <- rownames(temp_row_data)
    temp_gene_name <- colnames(temp_row_data)[temp_r]
    temp_r_data <- as.numeric(temp_row_data[,temp_row_data>0])
    temp_data <- data.frame(temp_lnc_name,temp_gene_name,temp_r_data)
    names(temp_data) <- c("lncRNA","gene","corr_r")
    imlnc_IRG_data <- rbind(imlnc_IRG_data,temp_data)
  }
  imlnc_IRG_data <- imlnc_IRG_data[-1,]
  dim(imlnc_IRG_data)
  #1] 1053244       3(ACC)
  
  
  #################################################################################################
  ##################################（lncRNA-miRNA-mRNA）
  
  
  ceRNA_data_lncmiR <- data.frame(lncRNA="",miRNA="",gene="")
  for(i in 1:dim(imlnc_immiR_data_p)[1]){
    temp_immiR_names <- imlnc_immiR_data_p[i,2]
    if(!any(temp_immiR_names %in% immiR_IRG_data_p[,1])) next
    temp_IRG_names <- immiR_IRG_data_p[immiR_IRG_data_p$miRNA == temp_immiR_names,]$gene#?是否用t() 先test一下
    temp_data <- data.frame(imlnc_immiR_data_p[i,1],imlnc_immiR_data_p[i,2],temp_IRG_names)
    names(temp_data) <-c("lncRNA","miRNA","gene")
    ceRNA_data_lncmiR <- rbind(ceRNA_data_lncmiR,temp_data)
  }
  
 
  ceRNA_lncgene <- function(x){
    temp_lnc_name <- x[1]
    temp_gene_name <- x[3]
    return(any(temp_gene_name %in% imlnc_IRG_data_p[imlnc_IRG_data_p$lncRNA == temp_lnc_name,]$gene))
  }
  ceRNA_app <- apply(ceRNA_data_lncmiR,1,ceRNA_lncgene) #TRUE or FALSE 值
  ceRNA_app_data <- ceRNA_data_lncmiR[ceRNA_app,]
  dim(ceRNA_app_data)
  
  write.csv(mylist,paste0(cancer[c],"_mylist.csv"),fileEncoding = "GB18030")
  print(cancer[c])
  print(mylist)
  save.image(paste(cancer[c],"_allRdata.RData",sep = ""))
  rm(list = ls())
  cancer <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","HNSC","KICH",
              "KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
              "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
  
}


#################################################################################################
#################################################################################################
cancer <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","HNSC","KICH",
            "KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
            "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
for(c in 1:length(cancer)){
  mylist <- list()
  input_data_path <- "../TCGA_33cancer_data"
  output_imcepair_path <- "../output/output_imcepair"
  output_list <- "../output/output_list"
  output_expce_path <- "../output/output_expce"
  setwd(output_imcepair_path)
  ceRNA_app_data <- read.csv(paste0(cancer[c],"_imceRNApair_data_p.csv"),header = T,row.names = 1)
  #############miRNA-gene  mirTarBase
  setwd(input_data_path)
#write.csv(expvali_MG_convertID,"expvali_mirTarBase_IDconvert.csv")
  expvali_MG_convertID <- read.csv("expvali_MG_IDconvert.csv",header=T)

  expvali_miRgene <- function(x){
    temp_miR_name <- x[2]
    temp_gene_name <- x[3]
    return(any(temp_gene_name %in% expvali_MG_convertID[expvali_MG_convertID$MIRNA == temp_miR_name,]$gene))
  }
  expvali_miRgene_app <- apply(ceRNA_app_data,1,expvali_miRgene) #TRUE or FALSE
  expvali_miRgene_data <- ceRNA_app_data[expvali_miRgene_app,]
  dim(expvali_miRgene_data)
  #[1][1] 4040    3
  
  mylist[["expvaliImCe"]]["miRgene"] <- paste(cancer[c],"miRgene number：",dim(expvali_miRgene_data)[1],";lncRNA_num:",length(unique(expvali_miRgene_data[,1])),";miRNA_num:",length(unique(expvali_miRgene_data[,2])),";gene_num:",length(unique(expvali_miRgene_data[,3])),sep="")
  
  expvali_LM_convertID <- read.csv("expvali_LM_IDconvert.csv",header=T)

  expvali_lncmiR <- function(x){
    temp_lnc_name <- x[1]
    temp_miR_name <- x[2]
    return(any(temp_lnc_name %in% expvali_LM_convertID[expvali_LM_convertID$MIRNA == temp_miR_name,]$lncRNA))
  }
  expvali_lncmiR_app <- apply(ceRNA_app_data,1,expvali_lncmiR) #TRUE or FALSE 值
  expvali_lncmiR_data <- ceRNA_app_data[expvali_lncmiR_app,]
  dim(expvali_lncmiR_data)
  #[1] [1] 701   3
  setwd(output_expce_path)
  write.csv(expvali_lncmiR_data,paste(cancer[c],"_expvali_lncmiR_data.csv",sep=""))

  
  mylist[["expvaliImCe"]]["lncmiR"] <- paste(cancer[c],"lncmiR number：",dim(expvali_lncmiR_data)[1],paste(";lncRNA_num:",length(unique(expvali_lncmiR_data[,1])),";miRNA_num:",length(unique(expvali_lncmiR_data[,2])),";gene_num:",length(unique(expvali_lncmiR_data[,3])),sep=""),sep="")
  
  paste_ceRNA <- function(x){
    return(paste(x[1],x[2],x[3],sep = ";"))
  }
  
  expvali_ceRNA <- intersect(apply(expvali_lncmiR_data,1,paste_ceRNA),apply(expvali_miRgene_data,1,paste_ceRNA))
  expvali_ceRNA_data <- as.data.frame(t(as.data.frame(strsplit(expvali_ceRNA,";"))))
  if(dim(expvali_ceRNA_data)[1]!=0) {
    colnames(expvali_ceRNA_data) <- c("lncRNA","miRNA","gene")
    rownames(expvali_ceRNA_data) <- 1:dim(expvali_ceRNA_data)[1]
    dim(expvali_ceRNA_data)
  }
  
  #[1] 26  3
  length(unique(expvali_ceRNA_data$lncRNA))#lncRNA 3
  length(unique(expvali_ceRNA_data$miRNA))#miRNA 1
  length(unique(expvali_ceRNA_data$gene))#gene 11
  write.csv(expvali_ceRNA_data,paste(cancer[c],"_expimceRNA_data.csv",sep=""))
  mylist[["expvaliImCe"]]["ceRNA"] <- paste(cancer[c],"number：",dim(expvali_ceRNA_data)[1],paste("  lncRNA_num:",length(unique(expvali_ceRNA_data[,1])),"  miRNA_num:",length(unique(expvali_ceRNA_data[,2])),"  gene_num:",length(unique(expvali_ceRNA_data[,3])),sep=""),sep="")
  print(cancer[c])
  print(mylist)
  write.csv(mylist,paste0(cancer[c],"_mylist_exp.csv"),fileEncoding = "GB18030")
  #save.image(paste(cancer[c],"_allRdata.RData",sep = ""))
  rm(list = ls())
  cancer <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","HNSC","KICH",
              "KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG",
              "PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
}



