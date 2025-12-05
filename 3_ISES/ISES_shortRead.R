library(tidyverse)
library(data.table)
library(fgsea)
library(parallel)
library(metap)
library(pROC)    # AUROC
library(PRROC)   #AUPRC

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

# Paths
path_support = '~/result/support_file_shortRead'
path_DTU = '~/result'
path_save ='~/result'

dataType= 'ENCODE_NGS' # ENCODE_NGS, RNAWG_long-read

# load support 
vec_dataSubTypeTermRBP = readRDS(file.path(path_support,'vec_bioTypeTerm1Term2',  paste(dataType,".rds",sep='')))
list_ranks_FCFC = readRDS(file.path(path_support,'list_ranks_FCFC',  paste(dataType,".rds",sep='')))
list_vec_trans_tpm = readRDS(file.path(path_support,'vec_trans_tpm',  paste(dataType,".rds",sep='')))
list_vec_regulated_trans = readRDS(file.path(path_support,'list_vec_regulated_trans',  paste(dataType,".rds",sep='')))

# sig -> FCFC
# reg -> 1-pvalue
# com -> 1-pvalue
# comReg -> 1-pvalue
isoformGSEA_NGS <-function(combin,vec_trans_tpm, ranks_FCFC, vec_regulated_trans){
  vecSplit_dataSubTypeTermRBP <- unlist(str_split(combin, '_'))
  dataSubType <- vecSplit_dataSubTypeTermRBP[1]
  term <- vecSplit_dataSubTypeTermRBP[2]
  RBP_name <- vecSplit_dataSubTypeTermRBP[3]
  vec_m_DTU = c('DRIMSeq', 'DEXSeq','limma', 'edgeR','satuRn',
                'DTUrtle','NBSplice','RATs','SPIT','SUPPA2')
  # create common isoforms
  list_tab_DTU = list()
  for(m_DTU in vec_m_DTU){
    tab_DTU_ori = read.table(paste(path_DTU, paste('DTU',dataType,sep='_'), dataSubType, paste('tab_DTU',term, paste(RBP_name,'human',sep='-'),m_DTU,sep='_'),sep='/'),sep='\t',header=T)
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
    
    list_tab_DTU[m_DTU] = list(tab_DTU)
  }
  tab_num_DTU_sig <- as.data.frame(table(unlist(lapply(list_tab_DTU, function(x) { x %>% filter(pvalue < 0.05) %>% pull(transcript_id) }))))
  vec_trans_common_DTU <- tab_num_DTU_sig %>% filter(Freq >= (length(vec_m_DTU)/2)) %>% pull(1)
  tab_num_DTU <- as.data.frame(table(unlist(lapply(list_tab_DTU, function(x) { x %>% pull(transcript_id) }))))
  vec_trans_all_DTU <- tab_num_DTU %>% pull(1)
  
  list_tab_DTU_rand <- list()
  list_vec_trans_common_DTU <- list()
  list_tab_DTU_rand[[1]] <- list_tab_DTU
  list_vec_trans_common_DTU[[1]] <- vec_trans_common_DTU

  for (random_i in 2:2) {
    list_tab_DTU_rand[[random_i]] <- list()
    for (m_DTU in vec_m_DTU) {
      set.seed(random_i)
      list_tab_DTU_rand[[random_i]][[m_DTU]] <- list_tab_DTU[[m_DTU]] %>% mutate(transcript_id = sample(transcript_id))
    }
    list_vec_trans_common_DTU[[random_i]] <- as.data.frame(table(unlist(lapply(list_tab_DTU_rand[[random_i]],
                                                                               function(x) { x %>% filter(pvalue < 0.05) %>% pull(transcript_id) })))) %>%
      filter(Freq >= (length(vec_m_DTU)/2)) %>% pull(1)
  }


  list_tab_jaccard_temp <- list()
  # create targets isoforms
  for (random_i in 1:1) {
    ## common -> (roc,prc) -> AUC
    vec_trans_common_DTU = list_vec_trans_common_DTU[[random_i]]

    for (m_DTU in vec_m_DTU) {
      tab_DTU <- list_tab_DTU_rand[[random_i]][[m_DTU]]
      vec_trans_sig_pvalue <- tab_DTU %>% filter(pvalue < 0.05) %>% pull(transcript_id)
      vec_trans_sig_FDR <- tab_DTU %>% filter(FDR < 0.05) %>% pull(transcript_id)
      vec_trans_DTU <- tab_DTU %>% pull(transcript_id)

      ## common -> AUROC, AUPRC
      if(length(vec_trans_common_DTU)!=0){
        labels_com <- as.numeric(tab_DTU$transcript_id %in% vec_trans_common_DTU) 
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
      vec_comReg_trans = intersect(vec_trans_common_DTU, vec_regulated_trans)
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
      vec_trans_com = intersect(vec_trans_common_DTU, vec_trans_tpm)
      list_trans_com = list(vec_trans_com);names(list_trans_com) = 'common'

      ## ES comReg
      vec_trans_comReg = intersect(vec_comReg_trans, vec_trans_tpm)
      list_trans_comReg = list(vec_trans_comReg);names(list_trans_comReg) = 'comReg'

      for(gseaParam in c(1)){

        ##  ES target
        if(length(intersect(vec_trans_sig, names(ranks_FCFC)))>1){
          GSEA_sig = fgsea(list_trans_sig, ranks_FCFC,scoreType='pos',gseaParam=gseaParam,maxSize=300000, nperm=1000)
          ES_sig = GSEA_sig$ES; NES_sig = GSEA_sig$NES; pval_sig = GSEA_sig$pval
        }else{ES_sig = 0; NES_sig=0; pval_sig=1}

        ##  ES regulated
        if(length(intersect(vec_trans_reg, names(ranks_DTU)))>1){
          GSEA_reg = fgsea(list_trans_reg, ranks_DTU,scoreType='pos',gseaParam=gseaParam,maxSize=300000, nperm=1000)
          ES_reg = GSEA_reg$ES; NES_reg = GSEA_reg$NES; pval_reg = GSEA_reg$pval
        }else{ES_reg = 0; NES_reg=0; pval_reg=1}

        ##  ES common
        if(length(intersect(vec_trans_com, names(ranks_DTU)))>1){
          GSEA_com = fgsea(list_trans_com, ranks_DTU,scoreType='pos',gseaParam=gseaParam,maxSize=300000, nperm=1000)
          ES_com = GSEA_com$ES; NES_com = GSEA_com$NES; pval_com = GSEA_com$pval
        }else{ES_com = 0; NES_com=0; pval_com=1}

        ##  ES comReg
        if(length(intersect(vec_trans_comReg, names(ranks_DTU)))>1){
          GSEA_comReg = fgsea(list_trans_comReg, ranks_DTU,scoreType='pos',gseaParam=gseaParam,maxSize=300000, nperm=1000)
          ES_comReg = GSEA_comReg$ES; NES_comReg = GSEA_comReg$NES; pval_comReg = GSEA_comReg$pval
        }else{ES_comReg = 0; NES_comReg=0; pval_comReg=1}

        tab_jaccard_temp_temp =
          data.frame(dataSubType = dataSubType,
                     term1 = term,
                     term2 = RBP_name,
                     random_i = random_i,
                     m_DTU = m_DTU,
                     gseaParam = gseaParam,

                     num_trans_sig_pvalue = length(vec_trans_sig_pvalue),
                     num_trans_sig_FDR = length(vec_trans_sig_FDR),
                     ES_sig = ES_sig,NES_sig = NES_sig, pval_sig = pval_sig,

                     num_trans_com = length(vec_trans_common_DTU),
                     ES_com = ES_com, NES_com = NES_com, pval_com = pval_com,
                     jaccard_com = jaccard(vec_trans_common_DTU, vec_trans_sig_pvalue),
                     auroc_com = auroc_com, auprc_com = auprc_com,

                     num_trans_reg = length(vec_regulated_trans),
                     ES_reg = ES_reg, NES_reg = NES_reg, pval_reg = pval_reg,
                     jaccard_reg = jaccard(vec_regulated_trans, vec_trans_sig_pvalue),
                     auroc_reg = auroc_reg, auprc_reg = auprc_reg,

                     num_trans_comReg = length(vec_comReg_trans),
                     ES_comReg = ES_comReg, NES_comReg = NES_comReg, pval_comReg = pval_comReg,
                     jaccard_comReg = jaccard(vec_comReg_trans, vec_trans_sig_pvalue),
                     auroc_comReg = auroc_comReg, auprc_comReg = auprc_comReg)

        list_tab_jaccard_temp[[paste(as.character(random_i),as.character(gseaParam),m_DTU)]] = tab_jaccard_temp_temp
      }

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
  isoformGSEA_NGS(combin,list_vec_trans_tpm[[combin]],list_ranks_FCFC[[combin]],list_vec_regulated_trans[[combin]])
}, 
mc.cores = 30, 
mc.preschedule = FALSE, 
mc.cleanup = TRUE       
)
tab_jaccard <- bind_rows(list_tab_jaccard)
gc()

# Save results
path_save_result <- file.path(path_save, paste0('ISES_', dataType))
dir.create(path_save_result, showWarnings = FALSE)
write.table(tab_jaccard, file.path(path_save_result, paste('tab', dataType, 'FCFC_gseaPara_Simple',sep='_')),
            row.names = FALSE, quote = FALSE, sep = '\t')

