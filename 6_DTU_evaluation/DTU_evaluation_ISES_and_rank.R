library(tidyverse)

# glm, gseaPara=1, NES_comReg

path_main <- '~/result'
################################################################################
# merge all result
list_tab_jaccard_gseaPara1 = list()
# long-read
for(dataType in c('RNAWG_long-read','single_cell','spatial')){
  if(dataType=='RNAWG_long-read'){vec_dataSubType = c('mouse','human')}
  if(dataType=='single_cell'){vec_dataSubType = c('PromethION_5cl_rep1', 'PromethION_5cl_rep2', 'PromethION_MSC')}
  if(dataType=='spatial'){vec_dataSubType = c('GSE153859_CBS1', 'GSE153859_CBS2', 'GSE153859_MOB')}
  path_save_result <- file.path(path_main, paste0('ISES_', dataType))
  list_tab_jaccard_gseaPara1[[dataType]] = list()
  for(dataSubType in vec_dataSubType){
    # jaccard
    tab_jaccard_temp = read.table(file.path(path_save_result, 
                                            paste('tab_rand2_jaccard', dataType, dataSubType, 'classification2feature','gseaPara','Simple',sep='_')),header=T,sep='\t')
    list_tab_jaccard_gseaPara1[[dataType]][[dataSubType]] = 
      tab_jaccard_temp %>% 
      filter(gseaParam == 1, m_identify=='glm') %>% 
      mutate(dataType=dataType, .before = 1)
    
  }
}
# ENCODE_NGS
dataType= 'ENCODE_NGS' # ENCODE_NGS
list_tab_jaccard_gseaPara1[[dataType]] = list()
# load result
tab_jaccard = read.table(file.path(path_main,paste0('ISES_', dataType), paste('tab', dataType, 'FCFC_gseaPara','Simple',sep='_')),header=T,sep='\t')
tab_jaccard_gseaPara1 = tab_jaccard %>% filter(gseaParam == 1)
tab_jaccard_gseaPara1['dataSubType'] = paste(tab_jaccard_gseaPara1 %>% pull(dataSubType),tab_jaccard_gseaPara1 %>% pull(term1),sep = '_')
for(dataSubType in c("CRISPR_HepG2","CRISPR_K562","shRNA_HepG2","shRNA_K562")){
    list_tab_jaccard_gseaPara1[[dataType]][[dataSubType]] = 
      tab_jaccard_gseaPara1 %>% filter(dataSubType == !!dataSubType)
}

tab_jaccard = bind_rows(flatten(list_tab_jaccard_gseaPara1))
tab_jaccard_mean = tab_jaccard %>% group_by(dataType,dataSubType,m_DTU) %>% summarise(mean_NES_comReg=mean(NES_comReg))
tab_jaccard_mean = tab_jaccard_mean %>% group_by(dataType,dataSubType) %>% arrange(desc(mean_NES_comReg)) %>% mutate(rank = 1:10, rank_draw=10:1) %>% arrange(dataType,dataSubType,m_DTU)
tab_jaccard_mean_draw = as.data.frame(tab_jaccard_mean) %>% select(dataSubType,m_DTU,rank_draw) %>%
  pivot_wider(names_from = dataSubType, values_from = rank_draw)

