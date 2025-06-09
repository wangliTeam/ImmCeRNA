#########immune cell infiltration

RNA_imcell_scc_fun <- function(imlnc_canceri_exp,IRG_canceri_exp,immiR_canceri_exp,im_infil_data,expvali_ceRNA_data){
  imlnc_canceri_exp_im <- imlnc_canceri_exp
  colnames(imlnc_canceri_exp_im) <- substr(colnames(imlnc_canceri_exp_im),1,15)
  inter_imfilexp_col <- intersect(colnames(imlnc_canceri_exp_im),rownames(im_infil_data))
  length(inter_imfilexp_col)#[1] 1212

  vali_ce_lnc_imdata <- imlnc_canceri_exp_im[unique(expvali_ceRNA_data$lncRNA),inter_imfilexp_col]

  cancer_immcell_data <- im_infil_data[inter_imfilexp_col,]
  cancer_immcell_data <- cancer_immcell_data[, colSums(is.na(cancer_immcell_data)) != nrow(cancer_immcell_data)]

  scc_lnc_imcell <- corr.test(as.data.frame(t(vali_ce_lnc_imdata)),cancer_immcell_data,method="spearman",adjust = "none",ci=FALSE)
  scc_lnc_imcell_R <- as.data.frame(scc_lnc_imcell$r)

  scc_lnc_imcell_R <- scc_lnc_imcell_R[, colSums(is.na(scc_lnc_imcell_R)) != nrow(scc_lnc_imcell_R)]


  lnc_imcell_data <- data.frame(lncRNA="",imcell="",scc_r="")
  temp_data <-data.frame()
  temp_row_data <- data.frame()
  for(i in 1:dim(scc_lnc_imcell_R)[1]){
    temp_row_data <- scc_lnc_imcell_R[i,]#[1] "data.frame"
    temp_lnc_name <- rownames(temp_row_data)
    temp_imcell_name <- colnames(temp_row_data)
    temp_r_data <- as.numeric(temp_row_data)
    temp_data <- data.frame(temp_lnc_name,temp_imcell_name,temp_r_data)
    names(temp_data) <- c("lncRNA","imcell","scc_r")
    lnc_imcell_data <- rbind(lnc_imcell_data,temp_data)
  }
  lnc_imcell_data <- lnc_imcell_data[-1,]
  dim(lnc_imcell_data)
  #[1] 159       3

  lnc_imcell_p <- c()
  for(i in 1:dim(scc_lnc_imcell_p)[1]){
    temp_row_data <- scc_lnc_imcell_p[i,]#[1] "data.frame"
    temp_p_data <- as.numeric(temp_row_data)
    lnc_imcell_p <- c(lnc_imcell_p,temp_p_data)
  }
  #lncRNA imcell scc_R scc_p
  lnc_imcell_data$scc_p <- lnc_imcell_p
  dim(lnc_imcell_data)
  #[1] 159   4
  lnc_imcell_data$scc_r <- as.numeric(lnc_imcell_data$scc_r)
 
  immiR_canceri_exp_im <- immiR_canceri_exp
  colnames(immiR_canceri_exp_im) <- substr(colnames(immiR_canceri_exp_im),1,15)
  inter_imfilexp_col <- intersect(colnames(immiR_canceri_exp_im),rownames(im_infil_data))
  length(inter_imfilexp_col)#[1] 1212
  
  vali_ce_miR_imdata <- immiR_canceri_exp_im[unique(expvali_ceRNA_data$miRNA),inter_imfilexp_col]

  cancer_immcell_data <- im_infil_data[inter_imfilexp_col,]
  cancer_immcell_data <- cancer_immcell_data[, colSums(is.na(cancer_immcell_data)) != nrow(cancer_immcell_data)]
  #

  scc_miR_imcell <- corr.test(as.data.frame(t(vali_ce_miR_imdata)),cancer_immcell_data,method="spearman",adjust = "none",ci=FALSE)
  scc_miR_imcell_R <- as.data.frame(scc_miR_imcell$r)
  scc_miR_imcell_p <- as.data.frame(scc_miR_imcell$p)

  scc_miR_imcell_R <- scc_miR_imcell_R[, colSums(is.na(scc_miR_imcell_R)) != nrow(scc_miR_imcell_R)]
  scc_miR_imcell_p <- scc_miR_imcell_p[, colSums(is.na(scc_miR_imcell_p)) != nrow(scc_miR_imcell_p)]
  #

  miR_imcell_data <- data.frame(miRNA="",imcell="",scc_r="")
  temp_data <-data.frame()
  temp_row_data <- data.frame()
  for(i in 1:dim(scc_miR_imcell_R)[1]){
    temp_row_data <- scc_miR_imcell_R[i,]#[1] "data.frame"
    temp_miR_name <- rownames(temp_row_data)
    temp_imcell_name <- colnames(temp_row_data)
    temp_r_data <- as.numeric(temp_row_data)
    temp_data <- data.frame(temp_miR_name,temp_imcell_name,temp_r_data)
    names(temp_data) <- c("miRNA","imcell","scc_r")
    miR_imcell_data <- rbind(miR_imcell_data,temp_data)
  }
  miR_imcell_data <- miR_imcell_data[-1,]
  dim(miR_imcell_data)
  #[1] 53       3
  

  miR_imcell_p <- c()
  for(i in 1:dim(scc_miR_imcell_p)[1]){
    temp_row_data <- scc_miR_imcell_p[i,]#[1] "data.frame"
    temp_p_data <- as.numeric(temp_row_data)
    miR_imcell_p <- c(miR_imcell_p,temp_p_data)
  }
  #lncRNA imcell scc_R scc_p
  miR_imcell_data$scc_p <- miR_imcell_p
  dim(miR_imcell_data)
  #[1] 53   4

  miR_imcell_data$scc_r <- as.numeric(miR_imcell_data$scc_r)
  

  miR_imcell_data <- miR_imcell_data[miR_imcell_data$scc_p < 0.05 ,]
  dim(miR_imcell_data)
  
  
  IRG_canceri_exp_im <- IRG_canceri_exp
  colnames(IRG_canceri_exp_im) <- substr(colnames(IRG_canceri_exp_im),1,15)
  inter_imfilexp_col <- intersect(colnames(IRG_canceri_exp_im),rownames(im_infil_data))
  length(inter_imfilexp_col)#[1] 1212

  vali_ce_IRG_imdata <- IRG_canceri_exp_im[unique(expvali_ceRNA_data$gene),inter_imfilexp_col]

  cancer_immcell_data <- im_infil_data[inter_imfilexp_col,]
  cancer_immcell_data <- cancer_immcell_data[, colSums(is.na(cancer_immcell_data)) != nrow(cancer_immcell_data)]

  scc_IRG_imcell <- corr.test(as.data.frame(t(vali_ce_IRG_imdata)),cancer_immcell_data,method="spearman",adjust = "none",ci=FALSE)
  scc_IRG_imcell_R <- as.data.frame(scc_IRG_imcell$r)
  scc_IRG_imcell_p <- as.data.frame(scc_IRG_imcell$p)

  scc_IRG_imcell_R <- scc_IRG_imcell_R[, colSums(is.na(scc_IRG_imcell_R)) != nrow(scc_IRG_imcell_R)]
  scc_IRG_imcell_p <- scc_IRG_imcell_p[, colSums(is.na(scc_IRG_imcell_p)) != nrow(scc_IRG_imcell_p)]

  max(scc_IRG_imcell_R) #0.7896078
  min(scc_IRG_imcell_R) #-0.41925
  IRG_imcell_data <- data.frame(gene="",imcell="",scc_r="")
  temp_data <-data.frame()
  temp_row_data <- data.frame()
  for(i in 1:dim(scc_IRG_imcell_R)[1]){
    temp_row_data <- scc_IRG_imcell_R[i,]#[1] "data.frame"
    temp_gene_name <- rownames(temp_row_data)
    temp_imcell_name <- colnames(temp_row_data)
    temp_r_data <- as.numeric(temp_row_data)
    temp_data <- data.frame(temp_gene_name,temp_imcell_name,temp_r_data)
    names(temp_data) <- c("gene","imcell","scc_r")
    IRG_imcell_data <- rbind(IRG_imcell_data,temp_data)
  }
  IRG_imcell_data <- IRG_imcell_data[-1,]
  dim(IRG_imcell_data)
  #[1] 583       3

  IRG_imcell_p <- c()
  for(i in 1:dim(scc_IRG_imcell_p)[1]){
    temp_row_data <- scc_IRG_imcell_p[i,]#[1] "data.frame"
    temp_p_data <- as.numeric(temp_row_data)
    IRG_imcell_p <- c(IRG_imcell_p,temp_p_data)
  }
  #lncRNA imcell scc_R scc_p
  IRG_imcell_data$scc_p <- IRG_imcell_p
  dim(IRG_imcell_data)
  #[1] 583   4

  IRG_imcell_data$scc_r <- as.numeric(IRG_imcell_data$scc_r)
  
  IRG_imcell_data <- IRG_imcell_data[IRG_imcell_data$scc_p < 0.05 ,]
  dim(IRG_imcell_data)
  out <- list(lnc_imcell_data=lnc_imcell_data,miR_imcell_data=miR_imcell_data,IRG_imcell_data=IRG_imcell_data)
  return(out)
}


