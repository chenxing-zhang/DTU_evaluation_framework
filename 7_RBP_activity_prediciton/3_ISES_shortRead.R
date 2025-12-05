# ISES_ENCODE_NGS
library(tidyverse)
library(data.table)
library(fgsea)
library(parallel)
library(metap)

id2idTitle <- function(vec) {
  return(unlist(lapply(vec, function(x) unlist(str_split(x, '\\.'))[1])))
}

jaccard <- function(a, b) {
  if (length(a) == 0 | length(b) == 0) {
    return(0)
  } else {
    intersection = length(intersect(a, b))
    union = length(a) + length(b) - intersection
    return(intersection / union)
  }
}
dataType= 'ENCODE_NGS'
path_DTU_NGS = paste('~/result', paste('DTU',dataType,sep='_'),sep='/')
path_save = '~/result/RBP_activity_prediction/ISES_ENCODE_NGS'

path_support_NGS = '~/result/support_file_shortRead'
vec_dataSubTypeTermRBP = readRDS(file.path(path_support_NGS,'vec_bioTypeTerm1Term2',  "ENCODE_NGS.rds"))
list_vec_regulated_trans = readRDS(file.path(path_support_NGS,'list_vec_regulated_trans',  "ENCODE_NGS.rds"))
list_vec_common_trans = readRDS(file.path(path_support_NGS,'list_vec_common_trans',  "ENCODE_NGS.rds"))
list_vec_overlapped_trans = readRDS(file.path(path_support_NGS,'list_vec_overlapped_trans',  "ENCODE_NGS.rds"))
list_vec_trans_tpm = readRDS(file.path(path_support_NGS,'vec_trans_tpm',  paste(dataType,".rds",sep='')))

vec_dataSubTypeTermRBP1_dataSubTypeTermRBP2 = unlist(lapply(vec_dataSubTypeTermRBP, function(x) lapply(vec_dataSubTypeTermRBP, function(y) paste(x,y,sep='_'))))


isoformGSEA_NGS <-function(combin,vec_dataSubTypeTermRBP, vec_trans_tpm, vec_common_trans, vec_regulated_trans){
  vecSplit_dataSubTypeTermRBP_a <- unlist(str_split(combin, '_'))
  dataSubType_a <- vecSplit_dataSubTypeTermRBP_a[1]
  term_a <- vecSplit_dataSubTypeTermRBP_a[2]
  RBP_name_a <- vecSplit_dataSubTypeTermRBP_a[3]
  
  vec_m_DTU = c('satuRn')
  # 
  list_tab_DTU = list()
  for(combin_b in vec_dataSubTypeTermRBP){
    vecSplit_dataSubTypeTermRBP_b <- unlist(str_split(combin_b, '_'))
    dataSubType_b <- vecSplit_dataSubTypeTermRBP_b[1]
    term_b <- vecSplit_dataSubTypeTermRBP_b[2]
    RBP_name_b <- vecSplit_dataSubTypeTermRBP_b[3]
    
    for(m_DTU in vec_m_DTU){
      tab_DTU_ori = read.table(paste(path_DTU_NGS, dataSubType_b, paste('tab_DTU',term_b, paste(RBP_name_b,'human',sep='-'),m_DTU,sep='_'),sep='/'),sep='\t',header=T)
      if(m_DTU=='DRIMSeq'){
        tab_DTU = tab_DTU_ori %>% filter(grepl('^ENS',feature_id)&(!grepl('-[0-9]$',feature_id)))
        tab_DTU = tab_DTU %>% select(feature_id,prop_log2fc,pvalue,adj_pvalue)
      }
      if(m_DTU=='DEXSeq'){
        tab_DTU = tab_DTU_ori %>% filter(grepl('^ENS',featureID)&(!grepl('-[0-9]$',featureID)))
        tab_DTU = tab_DTU %>% select(featureID,grep('log2fold',colnames(tab_DTU),value=T),pvalue,padj)
      }
      if(m_DTU=='limma'){
        tab_DTU_ori['feature_id'] = rownames(tab_DTU_ori)
        tab_DTU = tab_DTU_ori %>% filter(grepl('^ENS',feature_id)&(!grepl('-[0-9]$',feature_id)))
        tab_DTU = tab_DTU %>% select(feature_id,logFC,P.Value,FDR)
      }
      if(m_DTU=='edgeR'){
        tab_DTU = tab_DTU_ori %>% filter(grepl('^ENS',feature_id)&(!grepl('-[0-9]$',feature_id)))
        tab_DTU = tab_DTU %>% select(feature_id,logFC,P.Value,FDR)
      }
      if(m_DTU=='satuRn'){
        tab_DTU_ori['feature_id'] = rownames(tab_DTU_ori)
        tab_DTU = tab_DTU_ori %>% filter(grepl('^ENS',feature_id)&(!grepl('-[0-9]$',feature_id)))
        tab_DTU = tab_DTU %>% select(feature_id,estimates,empirical_pval,empirical_FDR)
      }
      if(m_DTU=='DTUrtle'){
        tab_DTU = tab_DTU_ori %>% filter(grepl('^ENS',feature_id)&(!grepl('-[0-9]$',feature_id)))
        tab_DTU = tab_DTU %>% select(feature_id,lr,pvalue,adj_pvalue)
      }
      if(m_DTU=='NBSplice'){
        tab_DTU = tab_DTU_ori %>% filter(grepl('^ENS',iso)&(!grepl('-[0-9]$',iso)))
        tab_DTU = tab_DTU %>% select(iso,odd,pval,FDR)
      }
      if(m_DTU=='RATs'){
        tab_DTU = tab_DTU_ori %>% filter(grepl('^ENS',target_id)&(!grepl('-[0-9]$',target_id)))
        tab_DTU = tab_DTU %>% select(target_id,log2FC,pval,pval_corr)
      }
      if(m_DTU=='SPIT'){
        tab_DTU = tab_DTU_ori %>% filter(grepl('^ENS',transcript_id)&(!grepl('-[0-9]$',transcript_id)))
        tab_DTU = tab_DTU %>% select(transcript_id,likelihood,pvalue)
        tab_DTU['FDR'] = tab_DTU['pvalue']
      }
      if(m_DTU=='SUPPA2'){
        tab_DTU = tab_DTU_ori %>% filter(grepl('^ENS',transcript)&(!grepl('-[0-9]$',transcript)))
        tab_DTU = tab_DTU %>% select(transcript,dPSI,p_value)
        tab_DTU['FDR'] = tab_DTU['p_value']
      }
      colnames(tab_DTU) = c('transcript_id','logFC','pvalue','FDR')
      tab_DTU['transcript_id_title'] = id2idTitle(tab_DTU$transcript_id)
      
      # NA
      # tab_DTU = tab_DTU %>% filter(!is.na(pvalue))
      tab_DTU <- tab_DTU %>% mutate(pvalue = ifelse(is.na(pvalue), 1, pvalue),
                                    FDR = ifelse(is.na(FDR), 1, FDR))
      
      tab_DTU = tab_DTU[,c(grep('transcript',colnames(tab_DTU),value=TRUE),grep('gene',colnames(tab_DTU),value=TRUE),colnames(tab_DTU)[!grepl('transcript|gene',colnames(tab_DTU))])]
      tab_DTU['-log10pvalue'] = -log10(tab_DTU$pvalue)
      
      vec_pvalue = tab_DTU %>% pull(pvalue)
      vec_pvalue_unique = sort(unique(vec_pvalue))
      if(vec_pvalue_unique[1]==0){vec_pvalue[vec_pvalue==0] = vec_pvalue_unique[2]/2}
      if(vec_pvalue_unique[length(vec_pvalue_unique)]==1){vec_pvalue[vec_pvalue==1] = (vec_pvalue_unique[length(vec_pvalue_unique)-1] + 1)/2}
      tab_DTU['pvalue_noInf'] = vec_pvalue
      tab_DTU['-log10pvalue_noInf'] = -log10(tab_DTU$pvalue_noInf)
      tab_DTU['-lnpvalue_noInf'] = -log(tab_DTU$pvalue_noInf)
      
      list_tab_DTU[paste(m_DTU, combin_b,sep='_')] = list(tab_DTU)
    }
  }
  
  
  list_tab_jaccard_temp <- list()
  # create targets isoforms
  for (m_DTU in vec_m_DTU) {
      for(combin_b in vec_dataSubTypeTermRBP){
        vecSplit_dataSubTypeTermRBP_b <- unlist(str_split(combin_b, '_'))
        dataSubType_b <- vecSplit_dataSubTypeTermRBP_b[1]
        term_b <- vecSplit_dataSubTypeTermRBP_b[2]
        RBP_name_b <- vecSplit_dataSubTypeTermRBP_b[3]
        
        tab_DTU <- list_tab_DTU[[paste(m_DTU, combin_b,sep='_')]]
        vec_trans_sig_pvalue <- tab_DTU %>% filter(pvalue < 0.05) %>% pull(transcript_id)
        vec_trans_sig_FDR <- tab_DTU %>% filter(FDR < 0.05) %>% pull(transcript_id)
        vec_trans_DTU <- tab_DTU %>% pull(transcript_id)
        
        ## common -> AUROC, AUPRC
        if(length(vec_common_trans)!=0){
          labels_com <- as.numeric(tab_DTU$transcript_id %in% vec_common_trans)  
          scores_com <- 1 - tab_DTU$pvalue_noInf 
          roc_obj_com <- roc(response = labels_com, predictor = scores_com) 
          auroc_com <- as.numeric(auc(roc_obj_com)) 
          pr_obj_com <- pr.curve(scores.class0 = scores_com, weights.class0 = labels_com, curve = TRUE)  
          auprc_com <- pr_obj_com$auc.integral  
        }else{auroc_com=0;auprc_com=0}
        
        ## regulated_trans -> AUROC, AUPRC
        if(length(vec_regulated_trans)!=0){
          labels_reg <- as.numeric(tab_DTU$transcript_id %in% vec_regulated_trans) 
          scores_reg <- 1 - tab_DTU$pvalue_noInf  
          roc_obj_reg <- roc(response = labels_reg, predictor = scores_reg)  
          auroc_reg <- as.numeric(auc(roc_obj_reg))  
          pr_obj_reg <- pr.curve(scores.class0 = scores_reg, weights.class0 = labels_reg, curve = TRUE)  
          auprc_reg <- pr_obj_reg$auc.integral  
        }else{auroc_reg=0;auprc_reg=0}
        
        ## common & regulated_trans -> AUROC, AUPRC
        vec_comReg_trans = intersect(vec_common_trans, vec_regulated_trans)
        if(length(vec_comReg_trans)!=0){
          labels_comReg <- as.numeric(tab_DTU$transcript_id %in% vec_comReg_trans)  
          scores_comReg <- 1 - tab_DTU$pvalue_noInf 
          roc_obj_comReg <- roc(response = labels_comReg, predictor = scores_comReg) 
          auroc_comReg <- as.numeric(auc(roc_obj_comReg))  
          pr_obj_comReg <- pr.curve(scores.class0 = scores_comReg, weights.class0 = labels_comReg, curve = TRUE)  
          auprc_comReg <- pr_obj_comReg$auc.integral  
        }else{auroc_comReg=0;auprc_comReg=0}
        
        ## ES significant, DTU pvalue<0.05 isoforms
        vec_trans_sig = intersect(vec_trans_sig_pvalue, vec_trans_tpm)
        list_trans_sig = list(vec_trans_sig);names(list_trans_sig) = 'significant'
        
        ## ES regulated, RBP-regulated isoforms
        vec_trans_reg = intersect(vec_regulated_trans, vec_trans_tpm)
        list_trans_reg = list(vec_trans_reg);names(list_trans_reg) = 'regulated'
        ## DTU 1-pvalue list
        ranks_DTU = tab_DTU %>% arrange(desc(`-log10pvalue_noInf`)) %>% pull(`-log10pvalue_noInf`)
        names(ranks_DTU) = tab_DTU %>% arrange(desc(`-log10pvalue_noInf`)) %>% pull(transcript_id)
        
        ## ES common, common iosforms
        vec_trans_com = intersect(vec_common_trans, vec_trans_tpm)
        list_trans_com = list(vec_trans_com);names(list_trans_com) = 'common'
        
        ## ES comReg
        vec_trans_comReg = intersect(vec_comReg_trans, vec_trans_tpm)
        list_trans_comReg = list(vec_trans_comReg);names(list_trans_comReg) = 'comReg'
        
        
        gseaParam = 1
        
        ## ES regulated
        if(length(intersect(vec_trans_reg, names(ranks_DTU)))>1){
          GSEA_reg = fgsea(list_trans_reg, ranks_DTU,scoreType='pos',gseaParam=gseaParam,maxSize=300000, nperm=1000)
          ES_reg = GSEA_reg$ES; NES_reg = GSEA_reg$NES; pval_reg = GSEA_reg$pval
        }else{ES_reg = 0; NES_reg=0; pval_reg=1}
        
        ## ES common
        if(length(intersect(vec_trans_com, names(ranks_DTU)))>1){
          GSEA_com = fgsea(list_trans_com, ranks_DTU,scoreType='pos',gseaParam=gseaParam,maxSize=300000, nperm=1000)
          ES_com = GSEA_com$ES; NES_com = GSEA_com$NES; pval_com = GSEA_com$pval
        }else{ES_com = 0; NES_com=0; pval_com=1}
        
        ## ES comReg
        if(length(intersect(vec_trans_comReg, names(ranks_DTU)))>1){
          GSEA_comReg = fgsea(list_trans_comReg, ranks_DTU,scoreType='pos',gseaParam=gseaParam,maxSize=300000, nperm=1000)
          ES_comReg = GSEA_comReg$ES; NES_comReg = GSEA_comReg$NES; pval_comReg = GSEA_comReg$pval
        }else{ES_comReg = 0; NES_comReg=0; pval_comReg=1}
        
        tab_jaccard_temp_temp =
          data.frame(dataSubType_a = dataSubType_a,
                     term1_a = term_a,
                     term2_a = RBP_name_a,
                     
                     dataSubType_b = dataSubType_b,
                     term1_b = term_b,
                     term2_b = RBP_name_b,
                     
                     m_DTU = m_DTU,
                     gseaParam = gseaParam,
                     
                     num_trans_sig_pvalue = length(vec_trans_sig_pvalue),
                     num_trans_sig_FDR = length(vec_trans_sig_FDR),
                     
                     num_trans_com = length(vec_common_trans),
                     ES_com = ES_com, NES_com = NES_com, pval_com = pval_com,
                     jaccard_com = jaccard(vec_common_trans, vec_trans_sig_pvalue),
                     auroc_com = auroc_com, auprc_com = auprc_com,
                     
                     num_trans_reg = length(vec_regulated_trans),
                     ES_reg = ES_reg, NES_reg = NES_reg, pval_reg = pval_reg,
                     jaccard_reg = jaccard(vec_regulated_trans, vec_trans_sig_pvalue),
                     auroc_reg = auroc_reg, auprc_reg = auprc_reg,
                     
                     num_trans_comReg = length(vec_comReg_trans),
                     ES_comReg = ES_comReg, NES_comReg = NES_comReg, pval_comReg = pval_comReg,
                     jaccard_comReg = jaccard(vec_comReg_trans, vec_trans_sig_pvalue),
                     auroc_comReg = auroc_comReg, auprc_comReg = auprc_comReg)
        
        list_tab_jaccard_temp[[paste(m_DTU, combin_b,sep='_')]] = tab_jaccard_temp_temp
    }
  }
  
  tab_jaccard_temp = bind_rows(list_tab_jaccard_temp)
  return(tab_jaccard_temp)
}


list_tab_jaccard = mclapply(vec_dataSubTypeTermRBP, function(combin) {
  suppressPackageStartupMessages({
    library(tidyverse)
    library(fgsea)
    library(pROC)    
    library(PRROC)   
  })
  isoformGSEA_NGS(combin,vec_dataSubTypeTermRBP, list_vec_trans_tpm[[combin]], list_vec_common_trans[[combin]], list_vec_regulated_trans[[combin]])
}, 
mc.cores = 30, 
mc.preschedule = FALSE,  
mc.cleanup = TRUE        
)
tab_jaccard <- bind_rows(list_tab_jaccard)
gc()
tab_jaccard = tab_jaccard %>% rename(dataSubType_isoformSet = dataSubType_a,
                                     term1_isoformSet = term1_a,
                                     term2_isoformSet =term2_a, 
                                     dataSubType_DTU = dataSubType_b,
                                     term1_DTU = term1_b,
                                     term2_DTU =term2_b
                                     )
write.table(tab_jaccard, file.path(path_save, paste('tab_jaccard', 'ENCODE_NGS',sep='_')),
            row.names = FALSE, quote = FALSE, sep = '\t')

ranked_data <- tab_jaccard %>%
  group_by(dataSubType_isoformSet, term1_isoformSet, term2_isoformSet) %>%
  mutate(
    NES_comReg_rank = rank(-NES_comReg, ties.method = "min")
  ) %>%
  ungroup() 
