impycor <- function(cancer_ceRNA,exp_data,sample_anno,filename){
  ############1.#############
  lnc_exp <- exp_data[intersect(cancer_ceRNA$lncRNA,rownames(exp_data)),]
  lnc_exp <- lnc_exp[which(rowSums(lnc_exp) > 0),]
  gene_exp <- exp_data[intersect(cancer_ceRNA$gene,rownames(exp_data)),]
  gene_exp <- gene_exp[which(rowSums(gene_exp) > 0),]
  library(psych)
  library(tidyr)
  library(dplyr)
  corr_ce <- corr.test(as.data.frame(t(lnc_exp)),as.data.frame(t(gene_exp)),method="pearson",adjust = "fdr",ci=FALSE)
  corr_ce_r <- as.data.frame(cbind(lncRNA=rownames(corr_ce$r),corr_ce$r))
  corr_ce_r <- gather(corr_ce_r, key = "gene",value = "scc_r",-`lncRNA`)
  corr_ce_p <- as.data.frame(cbind(lncRNA=rownames(corr_ce$p),corr_ce$p))
  corr_ce_p <- gather(corr_ce_p, key = "gene",value = "pvalue",-`lncRNA`)
  corr_ce_fdr <- as.data.frame(cbind(lncRNA=rownames(corr_ce$p.adj),corr_ce$p.adj))
  corr_ce_fdr <- gather(corr_ce_fdr, key = "gene",value = "fdr",-`lncRNA`)
  corr_ce_r$pvalue <- corr_ce_p[,3]
  corr_ce_r$fdr <- corr_ce_fdr[,3]
  corr_ce_r$condi <- "all"
  corr_ce_r_all <- corr_ce_r
  
  #############2############
  if(colnames(sample_anno)[2]=="RNR"){RorT <- "R"}else{RorT <- "T"}
  lnc_exp <- exp_data[intersect(cancer_ceRNA$lncRNA,rownames(exp_data)),sample_anno[sample_anno[,2]==RorT,1]]
  lnc_exp <- lnc_exp[which(rowSums(lnc_exp) > 0),]
  filtered_data <- data.frame(num=apply(lnc_exp, 1, function(x){length(unique(x))}))
  filtered_data$name <- rownames(filtered_data)
  lnc_exp <- lnc_exp[filtered_data[filtered_data$num!=1,"name"],]
  gene_exp <- exp_data[intersect(cancer_ceRNA$gene,rownames(exp_data)),sample_anno[sample_anno[,2]==RorT,1]]
  gene_exp <- gene_exp[which(rowSums(gene_exp) > 0),]
  filtered_data <- data.frame(num=apply(gene_exp, 1, function(x){length(unique(x))}))
  filtered_data$name <- rownames(filtered_data)
  gene_exp <- gene_exp[filtered_data[filtered_data$num!=1,"name"],]
  library(psych)
  library(tidyr)
  corr_ce <- corr.test(as.data.frame(t(lnc_exp)),as.data.frame(t(gene_exp)),method="pearson",adjust = "fdr",ci=FALSE)
  corr_ce_r <- as.data.frame(cbind(lncRNA=rownames(corr_ce$r),corr_ce$r))
  corr_ce_r <- gather(corr_ce_r, key = "gene",value = "scc_r",-`lncRNA`)
  corr_ce_p <- as.data.frame(cbind(lncRNA=rownames(corr_ce$p),corr_ce$p))
  corr_ce_p <- gather(corr_ce_p, key = "gene",value = "pvalue",-`lncRNA`)
  corr_ce_fdr <- as.data.frame(cbind(lncRNA=rownames(corr_ce$p.adj),corr_ce$p.adj))
  corr_ce_fdr <- gather(corr_ce_fdr, key = "gene",value = "fdr",-`lncRNA`)
  corr_ce_r$pvalue <- corr_ce_p[,3]
  corr_ce_r$fdr <- corr_ce_fdr[,3]
  corr_ce_r$condi <- RorT
  corr_ce_r_R <- corr_ce_r
  #############3.############
  if(colnames(sample_anno)[2]=="RNR"){NRorT <- "NR"}else{NRorT <- "NT"}
  lnc_exp <- exp_data[intersect(cancer_ceRNA$lncRNA,rownames(exp_data)),sample_anno[sample_anno[,2]==NRorT,1]]
  lnc_exp <- lnc_exp[which(rowSums(lnc_exp) > 0),]
  filtered_data <- data.frame(num=apply(lnc_exp, 1, function(x){length(unique(x))}))
  filtered_data$name <- rownames(filtered_data)
  lnc_exp <- lnc_exp[filtered_data[filtered_data$num!=1,"name"],]
  gene_exp <- exp_data[intersect(cancer_ceRNA$gene,rownames(exp_data)),sample_anno[sample_anno[,2]==NRorT,1]]
  gene_exp <- gene_exp[which(rowSums(gene_exp) > 0),]
  filtered_data <- data.frame(num=apply(gene_exp, 1, function(x){length(unique(x))}))
  filtered_data$name <- rownames(filtered_data)
  gene_exp <- gene_exp[filtered_data[filtered_data$num!=1,"name"],]
  library(psych)
  library(tidyr)
  corr_ce <- corr.test(as.data.frame(t(lnc_exp)),as.data.frame(t(gene_exp)),method="pearson",adjust = "fdr",ci=FALSE)
  corr_ce_r <- as.data.frame(cbind(lncRNA=rownames(corr_ce$r),corr_ce$r))
  corr_ce_r <- gather(corr_ce_r, key = "gene",value = "scc_r",-`lncRNA`)
  corr_ce_p <- as.data.frame(cbind(lncRNA=rownames(corr_ce$p),corr_ce$p))
  corr_ce_p <- gather(corr_ce_p, key = "gene",value = "pvalue",-`lncRNA`)
  corr_ce_fdr <- as.data.frame(cbind(lncRNA=rownames(corr_ce$p.adj),corr_ce$p.adj))
  corr_ce_fdr <- gather(corr_ce_fdr, key = "gene",value = "fdr",-`lncRNA`)
  corr_ce_r$pvalue <- corr_ce_p[,3]
  corr_ce_r$fdr <- corr_ce_fdr[,3]
  corr_ce_r$condi <- NRorT
  corr_ce_r_NR <- corr_ce_r
  
  #summary
  corr_rbind <- rbind(corr_ce_r_all,rbind(corr_ce_r_R,corr_ce_r_NR))
  corr_sig <- subset(corr_rbind,scc_r > 0 & pvalue < 0.05)
  inte_ce <- merge(cancer_ceRNA,corr_rbind,by=c("lncRNA","gene"))
  inte_ce$imd_cancer <- strsplit(filename,"_")[[1]][2]
  inte_ce$Accession <- strsplit(filename,"_")[[1]][1]
  inte_ce <- inte_ce[,c(1:4,9:14)]
  # inte_ce_sig <- inte_ce_sig[,c(1:4,9:14)]
  out=list(corr_ce=corr_rbind,corr_ce_sig=corr_sig,sum_TCGA_impy=inte_ce)
  return(out)
}
####
