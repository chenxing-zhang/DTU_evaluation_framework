# function of DE

# add combin columns
add_combin_columns_DE <- function(tab_DE){
  tab_DE_logFC = tab_DE %>% select(grep('_logFC',colnames(tab_DE),value=TRUE))
  tab_DE_fdr = tab_DE %>% select(grep('_FDR',colnames(tab_DE),value=TRUE))
  tab_DE_pvalue = tab_DE %>% select(grep('_pvalue',colnames(tab_DE),value=TRUE))
  tab_DE['combinedLogFC'] = rowSums(tab_DE_logFC)
  tab_DE['combinedFDR'] = metapod::parallelFisher(as.list(tab_DE_fdr))$p.value
  tab_DE['combinedPvalue'] = metapod::parallelFisher(as.list(tab_DE_pvalue))$p.value
  tab_DE['combinedFDRFDR'] = p.adjust(tab_DE$combinedFDR,method='BH')
  tab_DE['combinedPvalueFDR'] = p.adjust(tab_DE$combinedPvalue,method='BH')
  tab_DE = tab_DE %>% arrange(combinedPvalue)
  return(tab_DE)
}

# filter significant gene
filter_DE_by_rate_significant <- function(tab_DE,upDown,threRate){
  tab_DE_pvalue = tab_DE[,grep('_pvalue',colnames(tab_DE),value=TRUE)]
  tab_DE_thre = tab_DE[rowSums(tab_DE_pvalue < 0.05) >= length(colnames(tab_DE_pvalue))*threRate,] # tab_DTU_12_filter_thre4_up
  if(upDown=='up'){tab_DE_filter = tab_DE_thre %>% filter(combinedLogFC>0)}
  if(upDown=='down'){tab_DE_filter = tab_DE_thre %>% filter(combinedLogFC<0)}
  if(upDown=='all') {tab_DE_filter = tab_DE_thre}
  tab_DE_filter = tab_DE_filter %>% arrange(combinedPvalue)
  return(tab_DE_filter)
}