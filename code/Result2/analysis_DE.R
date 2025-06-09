#differentally expression

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- plyr::rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}


  DE_sig_fun <- function(lnc_exp){
    is_0_num <- function(m,n){
      m1_0_num <- apply(m[,1:n],1,function(x){return(sum(x==0))})
      m2_0_num <- apply(m[,(n+1):ncol(m)],1,function(x){return(sum(x==0))})
      return(cbind(m1_0_num,m2_0_num))
    }
    p_t_test<-function(m,n){
      m<-log2(m+0.01)
      tpvalue<-apply(m,1,function(x){t.test(x[1:n],x[(n+1):ncol(m)],paired = FALSE, alternative ="two.sided")$p.value})
    }
    FoldChange<-function(m,n){
      fc<-rowMeans(m[,(n+1):ncol(m)])/rowMeans(m[,1:n])
    }
    p_fisher<-function(m){
      table<-matrix(m,byrow = F,nrow = 2)
      fpvalue<-fisher.test(table,alternative = "two.sided")$p.value
    }
    
    colnames(lnc_exp) <- substr(colnames(lnc_exp),1,15)
    group_list <- factor(ifelse(as.numeric(substr(colnames(lnc_exp),14,15)) < 10,'tumor','normal'))
    table(group_list)
    colData_lnc <- data.frame(TCGA_ID = colnames(lnc_exp),group_list= group_list)
    exp_lnc_cancer <- lnc_exp[,colData_lnc[colData_lnc$group_list=="tumor",1]]
    exp_lnc_normal <- lnc_exp[,colData_lnc[colData_lnc$group_list=="normal",1]]
    dim(exp_lnc_cancer) #dim(exp_lnc_normal)
    #[1] 16882  1113    #[1] 16882   113
    exp_lnc_all<-as.matrix(cbind(exp_lnc_normal,exp_lnc_cancer))
    num_normal<-ncol(exp_lnc_normal)
    num_cancer<-ncol(exp_lnc_cancer)
    colnames(exp_lnc_all)<-c(rep(0,num_normal),rep(1,num_cancer))
    lnc_all_0_num<-is_0_num(exp_lnc_all,num_normal)
    judge<-((lnc_all_0_num[,1]/num_normal)<0.3) & ((lnc_all_0_num[,2]/num_cancer)<0.3)  
    lnc_t_test<-exp_lnc_all[judge,]
    lnc_t_test_p<-p_t_test(lnc_t_test,num_normal)
    lnc_fc<-FoldChange(lnc_t_test,num_normal)
    lnc_t_test_fdr<-p.adjust(lnc_t_test_p,method="BH")
    lnc_t_test_p_fdr_fc<-data.frame(lnc=names(lnc_t_test_p),p.value=unname(lnc_t_test_p),
                                    fdr=unname(lnc_t_test_fdr),ratio=unname(lnc_fc),
                                    test_type="t_test")
    lnc_fisher<-lnc_all_0_num[!judge,]
    lnc_fisher_table<- cbind(lnc_fisher,num_normal-lnc_fisher[,1],num_cancer-lnc_fisher[,2])
    colnames(lnc_fisher_table) <- c("normal_0_num","cancer_0_num","normal_not0_num","cancer_not0_num")
    lnc_fisher_p<-apply(lnc_fisher_table,1,p_fisher)
    lnc_fisher_fdr<-p.adjust(lnc_fisher_p)
    lnc_fisher_ratio<-data.frame((lnc_fisher[,1]/num_normal)/(lnc_fisher[,2]/num_cancer))
    lnc_fisher_p_fdr_ratio<-data.frame(lnc=names(lnc_fisher_p),p.value=unname(lnc_fisher_p),
                                       fdr=unname(lnc_fisher_fdr),ratio=unname(lnc_fisher_ratio),
                                       test_type="fisher")
    lnc_res_all <- rbind(lnc_t_test_p_fdr_fc,lnc_fisher_p_fdr_ratio)
    dim(lnc_res_all)#[1] 16882     5
    return(lnc_res_all)
  }
  
  
  DE_ce_fun <- function(lnc_res_sig,miR_res_sig,gene_res_sig,expvali_ceRNA_data){
    colnames(lnc_res_sig) <- c("lnc","p.value","fdr","ratio","test_type")
    colnames(miR_res_sig) <- c("miR","p.value","fdr","ratio","test_type")
    colnames(gene_res_sig) <- c("gene","p.value","fdr","ratio","test_type")
    DE_ceRNA <- data.frame()
    for(i in 1:dim(expvali_ceRNA_data)[1]){
      if(!((expvali_ceRNA_data[i,1] %in% lnc_res_sig$lnc)|(expvali_ceRNA_data[i,2] %in% miR_res_sig$miR)|(expvali_ceRNA_data[i,3] %in% gene_res_sig$gene))) next
      temp_data <- expvali_ceRNA_data[i,]
      DE_ceRNA <- rbind(DE_ceRNA,temp_data)
    }
    dim(DE_ceRNA)
    DE_ceRNA <- DE_ceRNA %>% mutate(DE_type=case_when(lncRNA %in% lnc_res_sig$lnc & miRNA %in% miR_res_sig$miR & gene %in% gene_res_sig$gene ~ "lnc_miR_gene_sig",
                                                      lncRNA %in% lnc_res_sig$lnc & miRNA %in% miR_res_sig$miR ~ "lnc_miR_sig",
                                                      lncRNA %in% lnc_res_sig$lnc & gene %in% gene_res_sig$gene ~ "lnc_gene_sig",
                                                      miRNA %in% miR_res_sig$miR & gene %in% gene_res_sig$gene ~ "miR_gene_sig",
                                                      lncRNA %in% lnc_res_sig$lnc ~ "lnc_sig",
                                                      miRNA %in% miR_res_sig$miR ~ "miR_sig",
                                                      gene %in% gene_res_sig$gene ~ "gene_sig"))  
    return(DE_ceRNA)
  }
  