library(dplyr)
library(stringi)
library(data.table)
library(tidyverse)

id2idTitle <- function(vec){
  return(unlist(lapply(vec, function(x) unlist(str_split(x,'\\.'))[1])))
}

tab_DTU_ori_2_tab_DTU <- function(tab_DTU_ori, m_DTU){
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
  return(tab_DTU)
}

jaccard <- function (a, b) {
  if(length(a)==0 | length(b)==0){return(0)}
  intersection = length ( intersect (a,b))
  union = length (a) + length (b) - intersection
  return (intersection/union)
}

dataType= 'ENCODE_NGS' # ENCODE_NGS, RNAWG_long-read
# dataSubType = 'shRNA'  # CRISPR, shRNA, mouse, human, 


# load exp and DTU
path_exp = '~/data'
path_DTU = '~/result'
path_save = '~/result/result_overlap_betweenMethod'
path_support = '~/result/support_file_shortRead'

vec_m_DTU = c('DRIMSeq', 'DEXSeq','limma', 'edgeR','satuRn',
              'DTUrtle','NBSplice','RATs','SPIT','SUPPA2')

vec_dataSubTypeTermRBP = readRDS(file.path(path_support,'vec_bioTypeTerm1Term2',  paste(dataType,".rds",sep='')))

tab_overlap = data.frame()
for(combin in vec_dataSubTypeTermRBP){
  vecSplit_dataSubTypeTermRBP <- unlist(str_split(combin, '_'))
  dataSubType <- vecSplit_dataSubTypeTermRBP[1]
  term <- vecSplit_dataSubTypeTermRBP[2]
  RBP_name <- vecSplit_dataSubTypeTermRBP[3]
  
  for(m_DTU_i in 1:(length(vec_m_DTU)-1)){
    m_DTU_1 = vec_m_DTU[m_DTU_i]
    tab_DTU_ori_1 = read.table(paste(path_DTU, paste('DTU',dataType,sep='_'), dataSubType, paste('tab_DTU',term, paste(RBP_name,'human',sep='-'),m_DTU_1,sep='_'),sep='/'),sep='\t',header=T)
    tab_DTU_1 = tab_DTU_ori_2_tab_DTU(tab_DTU_ori_1, m_DTU_1)
    # vec_iosform_0.05_1 = tab_DTU_1 %>% arrange(pvalue) %>% filter(pvalue<0.05) %>% pull(transcript_id)
    vec_iosform_0.05_1 = tab_DTU_1 %>% arrange(pvalue) %>% pull(transcript_id)
    for(m_DTU_j in (m_DTU_i+1):length(vec_m_DTU)){
      m_DTU_2 = vec_m_DTU[m_DTU_j]
      print(paste(combin,m_DTU_1,m_DTU_2))
      tab_DTU_ori_2 = read.table(paste(path_DTU, paste('DTU',dataType,sep='_'), dataSubType, paste('tab_DTU',term, paste(RBP_name,'human',sep='-'),m_DTU_2,sep='_'),sep='/'),sep='\t',header=T)
      tab_DTU_2 = tab_DTU_ori_2_tab_DTU(tab_DTU_ori_2, m_DTU_2)
      # vec_iosform_0.05_2 = tab_DTU_2 %>% arrange(pvalue) %>% filter(pvalue<0.05) %>% pull(transcript_id)
      vec_iosform_0.05_2 = tab_DTU_2 %>% arrange(pvalue) %>% pull(transcript_id)
      
      tab_overlap_temp = data.frame(
        dataType = dataType,
        dataSubType = dataSubType,
        term1 = term,
        term2 = RBP_name,
        m_DTU_1 = m_DTU_1,
        m_DTU_2 = m_DTU_2,
        top = c(10,50,100,200,500),
        jaccard = c(jaccard(vec_iosform_0.05_1[1:10], vec_iosform_0.05_2[1:10]),
                    jaccard(vec_iosform_0.05_1[1:50], vec_iosform_0.05_2[1:50]),
                    jaccard(vec_iosform_0.05_1[1:100], vec_iosform_0.05_2[1:100]),
                    jaccard(vec_iosform_0.05_1[1:200], vec_iosform_0.05_2[1:200]),
                    jaccard(vec_iosform_0.05_1[1:500], vec_iosform_0.05_2[1:500])))
      
      if(dim(tab_overlap)[1]==0){tab_overlap = tab_overlap_temp}else{tab_overlap = rbind(tab_overlap,tab_overlap_temp)}
    }
  }
}
tab_overlap = tab_overlap %>% mutate(bioType='cell line',.before=term1)
write.table(tab_overlap,paste(path_save,paste('tab_overlapBetween',dataType,sep='_'),sep='/'), row.names = FALSE,quote=FALSE,sep='\t')


###############################################################################
# long-read, single-cell and sptial
path_support = '~/result/support_file_longRead'
path_DTU <- '~/result'
for(dataType in c('RNAWG_long-read', 'single_cell', 'spatial')){
  tab_overlap = data.frame()
  if(dataType=='RNAWG_long-read'){vec_dataSubType = c('mouse','human')}
  if(dataType=='single_cell'){vec_dataSubType = c('PromethION_5cl_rep1', 'PromethION_5cl_rep2', 'PromethION_MSC')}
  if(dataType=='spatial'){vec_dataSubType = c('GSE153859_CBS1', 'GSE153859_CBS2', 'GSE153859_MOB')}
  for(dataSubType in vec_dataSubType){
    vec_bioTypeTerm1Term2 = readRDS(file.path(path_support,'vec_bioTypeTerm1Term2', paste(paste(dataType, dataSubType,sep='_'),".rds",sep='')))
    for(bioTypeTerm1Term2 in vec_bioTypeTerm1Term2){
      sep_char <- ifelse(dataType == 'spatial', '-', '_')
      vecSplit_bioTypeTerm1Term2 <- unlist(str_split(bioTypeTerm1Term2, sep_char))
      bioType <- vecSplit_bioTypeTerm1Term2[1]
      term1 <- vecSplit_bioTypeTerm1Term2[2]
      term2 <- vecSplit_bioTypeTerm1Term2[3]
      
      for(m_DTU_i in 1:(length(vec_m_DTU)-1)){
        m_DTU_1 = vec_m_DTU[m_DTU_i]
        result_try = try({
          if (dataType == 'RNAWG_long-read') {
            tab_DTU_ori_1 <- read.table(paste(path_DTU, paste('DTU', dataType, sep = '_'), dataSubType, paste('tab_DTU', bioType, term1, term2, m_DTU_1, sep = '_'), sep = '/'), sep = '\t', header = TRUE)
          } else {
            tab_DTU_ori_1 <- read.table(paste(path_DTU, paste('DTU', dataType, sep = '_'), dataSubType, paste('tab_DTU', term1, term2, m_DTU_1, sep = '_'), sep = '/'), sep = '\t', header = TRUE)
          }
          tab_DTU_1 = tab_DTU_ori_2_tab_DTU(tab_DTU_ori_1, m_DTU_1)
          # vec_iosform_0.05_1 = tab_DTU_1 %>% arrange(pvalue) %>% filter(pvalue<0.05) %>% pull(transcript_id)
          vec_iosform_0.05_1 = tab_DTU_1 %>% arrange(pvalue) %>% pull(transcript_id)
        })
        if(inherits(result_try,'try-error')){
          vec_iosform_0.05_1 = c()
        }
        
        for(m_DTU_j in (m_DTU_i+1):length(vec_m_DTU)){
          m_DTU_2 = vec_m_DTU[m_DTU_j]
          print(paste(bioTypeTerm1Term2,m_DTU_1,m_DTU_2))
          result_try = try({
            if (dataType == 'RNAWG_long-read') {
              tab_DTU_ori_2 <- read.table(paste(path_DTU, paste('DTU', dataType, sep = '_'), dataSubType, paste('tab_DTU', bioType, term1, term2, m_DTU_2, sep = '_'), sep = '/'), sep = '\t', header = TRUE)
            } else {
              tab_DTU_ori_2 <- read.table(paste(path_DTU, paste('DTU', dataType, sep = '_'), dataSubType, paste('tab_DTU', term1, term2, m_DTU_2, sep = '_'), sep = '/'), sep = '\t', header = TRUE)
            }
            tab_DTU_2 = tab_DTU_ori_2_tab_DTU(tab_DTU_ori_2, m_DTU_2)
            # vec_iosform_0.05_2 = tab_DTU_2 %>% arrange(pvalue) %>% filter(pvalue<0.05) %>% pull(transcript_id)
            vec_iosform_0.05_2 = tab_DTU_2 %>% arrange(pvalue) %>% pull(transcript_id)
          })
          if(inherits(result_try,'try-error')){
            vec_iosform_0.05_2 = c()
          }
          
          tab_overlap_temp = data.frame(
            dataType = dataType,
            dataSubType = dataSubType,
            bioType = bioType, 
            term1 = term1,
            term2 = term2,
            m_DTU_1 = m_DTU_1,
            m_DTU_2 = m_DTU_2,
            top = c(10,50,100,200,500),
            jaccard = c(jaccard(vec_iosform_0.05_1[1:10], vec_iosform_0.05_2[1:10]),
                        jaccard(vec_iosform_0.05_1[1:50], vec_iosform_0.05_2[1:50]),
                        jaccard(vec_iosform_0.05_1[1:100], vec_iosform_0.05_2[1:100]),
                        jaccard(vec_iosform_0.05_1[1:200], vec_iosform_0.05_2[1:200]),
                        jaccard(vec_iosform_0.05_1[1:500], vec_iosform_0.05_2[1:500])))
          
          if(dim(tab_overlap)[1]==0){tab_overlap = tab_overlap_temp}else{tab_overlap = rbind(tab_overlap,tab_overlap_temp)}
        }
      }
      
    }
  }
  write.table(tab_overlap,paste(path_save,paste('tab_overlapBetween',dataType,sep='_'),sep='/'), row.names = FALSE,quote=FALSE,sep='\t')
}



