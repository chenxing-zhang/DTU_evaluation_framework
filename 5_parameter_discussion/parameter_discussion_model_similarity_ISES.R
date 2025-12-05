library(tidyverse)
library(ggpubr) # ggarrange
library(ggplot2)

path_save <- '~/result/parameter_discussion_gseaPara'
################################################################################
# merge all result
list_tab_jaccard_gseaPara1 = list()
for(dataType in c('RNAWG_long-read','single_cell','spatial')){
  if(dataType=='RNAWG_long-read'){vec_dataSubType = c('mouse','human')}
  if(dataType=='single_cell'){vec_dataSubType = c('PromethION_5cl_rep1', 'PromethION_5cl_rep2', 'PromethION_MSC')}
  if(dataType=='spatial'){vec_dataSubType = c('GSE153859_CBS1', 'GSE153859_CBS2', 'GSE153859_MOB')}
  for(dataSubType in vec_dataSubType){
    # jaccard
    # classification2feature, filterTopTop
    tab_jaccard_temp = read.table(file.path(path_save, 
                                            paste('tab_paraDiscuss', dataType, dataSubType,sep='_')),header=T,sep='\t')
    tab_jaccard_temp = tab_jaccard_temp %>% 
      filter(gseaParam == 1) %>% 
      mutate(dataType=dataType, .before = 1)
    
    # gseaPara to column name
    tab_jaccard_temp <- tab_jaccard_temp %>%
      pivot_wider(
        names_from = m_identify,
        values_from = c(num_trans_comReg, ES_comReg, NES_comReg, pval_comReg),
        names_sep = "_"
      )
    
    list_tab_jaccard_gseaPara1[[paste(dataType,dataSubType,sep='_')]] = tab_jaccard_temp
    
  }
}

# merge jaccard and rand
tab_jaccard_gseaPara1 = bind_rows(list_tab_jaccard_gseaPara1)
tab_jaccard_gseaPara1 = tab_jaccard_gseaPara1 %>% filter(NES_comReg_lda!=0 , NES_comReg_glm!=0)

## cor
list_cor <- tab_jaccard_gseaPara1 %>%
  group_by(dataType, dataSubType) %>%
  summarise({
    df <- cur_data()  
    cols <- c("NES_comReg_glm", "NES_comReg_lda")
    combs <- combn(cols, 2, simplify = FALSE)
    
    cor_list <- lapply(combs, function(vars) {
      ct <- cor.test(df[[vars[1]]], df[[vars[2]]],
                     method = "spearman",
                     use = "pairwise.complete.obs")
      tibble(
        var1 = vars[1],
        var2 = vars[2],
        spearman = round(ct$estimate, 2),
        p.value = format.pval(ct$p.value, digits = 2)
      )
    })
    bind_rows(cor_list)
  }, .groups = "drop")

