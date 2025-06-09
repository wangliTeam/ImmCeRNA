library(doParallel)
library(foreach)
library(dplyr)
library(psych)
library(tidyverse)
drug_DESPCC_fun <- function(cancer,sample_HighLow,sample_Riskscore,sample_drug_data,train_database){
  #
  ceRNA_names <- colnames(sample_HighLow)
  DE_drug_data <- data.frame()

  
  for(i in 1:length(ceRNA_names)){
    sample_drug_data_temp <- cbind(group=sample_HighLow[rownames(sample_drug_data),ceRNA_names[i]],
                                   sample_drug_data)
    #FC
    group_mean <- aggregate(sample_drug_data_temp[,-1], list(sample_drug_data_temp$group), mean)
    rownames(group_mean) <- as.character(group_mean[,1])
    group_mean <- group_mean[,-1]
    
    log2FC <- log2(as.numeric(group_mean["high-risk",])/as.numeric(group_mean["low-risk",]))
    
    pvalue <- foreach(j = 2 : ncol(sample_drug_data_temp), .combine = "c") %dopar% {wilcox.test( sample_drug_data_temp[,j] ~ group, data = sample_drug_data_temp)$p.value}
    pvalue <- cbind(drug_names, cbind(log2FC=as.numeric(t(log2FC)),pvalue))
    
    DE_drug_data <- rbind(DE_drug_data,
                          cbind(data.frame(cancer=cancer,ceRNA=ceRNA_names[i],train_database=train_database,method="Wilcoxon rank-sum test"),
                                pvalue))
    
  }
 
  corr_ce_drug <- corr.test(sample_Riskscore,sample_drug_data,use="pairwise",method="spearman",adjust = "fdr",ci=FALSE)
  corr_ce_drug_r <- as.data.frame(cbind(ceRNA=rownames(corr_ce_drug$r),corr_ce_drug$r))
  corr_ce_drug_r <- gather(corr_ce_drug_r, key = "drug",value = "scc_r",-`ceRNA`)
  corr_ce_drug_p <- as.data.frame(cbind(ceRNA=rownames(corr_ce_drug$p),corr_ce_drug$p))
  corr_ce_drug_p <- gather(corr_ce_drug_p, key = "drug",value = "pvalue",-`ceRNA`)
  corr_ce_drug_fdr <- as.data.frame(cbind(ceRNA=rownames(corr_ce_drug$p.adj),corr_ce_drug$p.adj))
  corr_ce_drug_fdr <- gather(corr_ce_drug_fdr, key = "drug",value = "fdr",-`ceRNA`)
  corr_ce_drug_r$pvalue <- corr_ce_drug_p[,3]
  corr_ce_drug_r$fdr <- corr_ce_drug_fdr[,3]
  #汇总
  SCC_drug_data <- cbind(data.frame(cancer=cancer,ceRNA=corr_ce_drug_r[,1],train_database=train_database,method="spearman"),
                         corr_ce_drug_r[,2:5])
  out <- list(DE_drug_data=DE_drug_data,SCC_drug_data=SCC_drug_data)

  return(out)
}

