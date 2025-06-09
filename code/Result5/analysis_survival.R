survival_fun <- function(lnc_exp,miR_exp,PCG_exp,clinical,expvali_ceRNA_data,cancer_type){

  group_list <- ifelse(as.numeric(substr(colnames(lnc_exp),14,15)) < 10,'tumor','normal')
  table(group_list)
  lnc_exp_tumor <- lnc_exp[,group_list=='tumor']
  dim(lnc_exp_tumor)#[1] 16882  1113

  colnames(lnc_exp_tumor) <- substr(colnames(lnc_exp_tumor),1,12) #1113
  length(unique(colnames(lnc_exp_tumor)))#1095
  colName <- as.data.frame(table(colnames(lnc_exp_tumor)))

  lnc_exp_tumor_uni <-  lnc_exp_tumor[,which(colName$Freq==1)]#16882  1082
  rep_colname <- c()
  for(i in 1:dim(colName)[1]){
    if(as.numeric(colName[i,2]) < 2) next
    rep_colname <- c(rep_colname,as.character(colName[i,1]))
    temp_sample <- apply(lnc_exp_tumor[,which(colnames(lnc_exp_tumor)==as.character(colName[i,1]))],1,mean)
    lnc_exp_tumor_uni <- cbind(lnc_exp_tumor_uni,temp_sample)
  }
  colnames(lnc_exp_tumor_uni) <- c(colnames(lnc_exp_tumor[,which(colName$Freq==1)]),rep_colname)

  
  meta = clinical
  
  meta <- meta[intersect(colnames(lnc_exp_tumor_uni),rownames(meta)),]
  lnc_exp_tumor_uni <- lnc_exp_tumor_uni[,intersect(colnames(lnc_exp_tumor_uni),rownames(meta))]
  dim(lnc_exp_tumor_uni)#[1] 16882  1064
  dim(meta)#[1] 1064    6
  
  miR_exp_tumor_uni <- miR_exp_tumor_uni[,intersect(colnames(miR_exp_tumor_uni),rownames(meta))]
  PCG_exp_tumor_uni <- PCG_exp_tumor_uni[,intersect(colnames(PCG_exp_tumor_uni),rownames(meta))]
 
  meta$days_to_death[is.na(meta$days_to_death)] <- 0   #
  meta$days_to_last_follow_up[is.na(meta$days_to_last_follow_up)] <- 0
  meta[meta$days_to_death=="'--","days_to_death"] <- 0   
  meta[meta$days_to_last_follow_up=="'--","days_to_last_follow_up"] <- 0

  
  ########################################################################################################

  dim(lnc_exp_tumor_uni)#[1] 16882  1064
  dim(meta)#[1] 1064    9
  imlnc_exp_tumor <- lnc_exp_tumor_uni[unique(expvali_ceRNA_data$lncRNA),]
  dim(imlnc_exp_tumor)#[1]3  1095
  
 
  lnc_coxtable<-as.data.frame(matrix(0,dim(imlnc_exp_tumor)[1],5))
  colnames(lnc_coxtable) <- c("lncRNA", "coef","HR", "95% CI for HR", "p.value")
  #
  for (i in 1:dim(imlnc_exp_tumor)[1]) {
    meta$gene<-t(imlnc_exp_tumor[i,rownames(meta)])
    meta$group<-factor(ifelse(meta$gene>median(meta$gene),"high","low"),levels = c("low","high"))
   
    cox_enh <- coxph(Surv(time,event)~ group,data=meta)
    tmp <- summary(cox_enh)
    lnc_coxtable[i,1] <- rownames(imlnc_exp_tumor)[i];
    lnc_coxtable[i,2]<-tmp[[7]][1][1]
    lnc_coxtable[i,5] <- signif(tmp$coefficients[ ,5],4);
    lnc_coxtable[i,3] <- signif(tmp$coefficients[ ,2],4);
    HR.confint.lower <- signif(tmp$conf.int[ ,"lower .95"], 4);
    HR.confint.upper <- signif(tmp$conf.int[ ,"upper .95"], 4);
    lnc_coxtable[i,4] <- paste0(HR.confint.lower, "-", HR.confint.upper);	
  }
  lnc_coxtable$coef<-as.numeric(lnc_coxtable$coef)
 
  coxtable <- rbind(lnc_coxtable,rbind(miR_coxtable,gene_coxtable))
  colnames(coxtable) <- c("RNA","coef","HR","95% CI for HR","p.value")
  #gene_coxtable
  #
  clin_riskScore <- meta[,c("time","event")]
  clin_HighLow <- meta[,c("time","event")]
  ceRNA_names <- c()
  for(i in 1:dim(expvali_ceRNA_data)[1]){
    temp_ce <- expvali_ceRNA_data[i,]
    temp_riskScore <- imlnc_exp_tumor[temp_ce$lncRNA,]*coxtable[which(coxtable$RNA==temp_ce$lncRNA),"coef"]+
      immiR_exp_tumor[temp_ce$miRNA,]*coxtable[which(coxtable$RNA==temp_ce$miRNA),"coef"]+
      iRG_exp_tumor[temp_ce$gene,]*coxtable[which(coxtable$RNA==temp_ce$gene),"coef"]
    temp_HighLow <- factor(ifelse(temp_riskScore>median(t(temp_riskScore)),"high-risk","low-risk"),levels = c("low-risk","high-risk"))
    clin_riskScore <- cbind(clin_riskScore,t(temp_riskScore))
    clin_HighLow <- cbind(clin_HighLow,temp_HighLow)
    ceRNA_names <- c(ceRNA_names,paste(temp_ce$lncRNA,temp_ce$miRNA,temp_ce$gene,sep = "-"))
  }
  colnames(clin_riskScore) <- c("time","event",ceRNA_names)
  colnames(clin_HighLow) <- c("time","event",ceRNA_names)
  dim(clin_HighLow)
  #[1] 1064   28
  out <- list(clin_riskScore=clin_riskScore,clin_HighLow=clin_HighLow)
  return(out)
}

