# count to FPKM,RPKM and TPM
# https://zhuanlan.zhihu.com/p/513391213

countToTpm <- function(counts, effLen)
{
  rate <- log(counts) - log(effLen)
  denom <- log(sum(exp(rate)))
  exp(rate - denom + log(1e6))
}

countToFpkm <- function(counts, effLen)
{
  N <- sum(counts)
  exp( log(counts) + log(1e9) - log(effLen) - log(N) )
}

fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

countToEffCounts <- function(counts, len, effLen)
{
  counts * (len / effLen)
}

################################################################################An example
################################################################################
# exp_count = tab_count_temp
exp_count_2_exp_tpm <- function(exp_count, geneOrTrans){
  path_efflen = '~/data/Ref/GENECODE_r45' # '~/data/Ref/GENECODE_r45', '~/data/Ref/ENSEMBL_r111'
  if(geneOrTrans=='gene'){tab_efflen = read.csv(paste(path_efflen,'tab_efflen_gene',sep='/'),sep='\t',check.names = F)}
  if(geneOrTrans=='transcript'){tab_efflen = read.csv(paste(path_efflen,'tab_efflen_trans',sep='/'),sep='\t',check.names = F)}
  if(geneOrTrans=='all'){
    tab_efflen_gene = read.csv(paste(path_efflen,'tab_efflen_gene',sep='/'),sep='\t',check.names = F)
    tab_efflen_trans = read.csv(paste(path_efflen,'tab_efflen_trans',sep='/'),sep='\t',check.names = F)
    colnames(tab_efflen_gene) = c('all_id','efflen')
    colnames(tab_efflen_trans) = c('all_id','efflen')
    tab_efflen = rbind(tab_efflen_gene,tab_efflen_trans)
  }
  if(length(grep('gene|trans',colnames(exp_count),value = TRUE))==0){
    exp_count[paste(geneOrTrans,'id',sep='_')] = rownames(exp_count)
  }
  list_sample = colnames(exp_count)[!grepl('gene|trans|all',colnames(exp_count))]
  tab_count_efflen = merge(exp_count,tab_efflen,by=paste(geneOrTrans,'id',sep='_'),all.x=TRUE)
  exp_tpm = tab_count_efflen[paste(geneOrTrans,'id',sep='_')]
  for(sample in list_sample){
    exp_tpm[sample] = countToTpm(tab_count_efflen[sample],tab_count_efflen['efflen'])
  }
  return(exp_tpm)
}

exp_count_2_exp_tpm_path <- function(exp_count, geneOrTrans, path_efflen){
  # path_efflen = '~/data/Ref/GENECODE_r45' # '~/data/Ref/GENECODE_r45', '~/data/Ref/ENSEMBL_r111'
  if(geneOrTrans=='gene'){tab_efflen = read.csv(paste(path_efflen,'tab_efflen_gene',sep='/'),sep='\t',check.names = F)}
  if(geneOrTrans=='transcript'){tab_efflen = read.csv(paste(path_efflen,'tab_efflen_trans',sep='/'),sep='\t',check.names = F)}
  if(geneOrTrans=='all'){
    tab_efflen_gene = read.csv(paste(path_efflen,'tab_efflen_gene',sep='/'),sep='\t',check.names = F)
    tab_efflen_trans = read.csv(paste(path_efflen,'tab_efflen_trans',sep='/'),sep='\t',check.names = F)
    colnames(tab_efflen_gene) = c('all_id','efflen')
    colnames(tab_efflen_trans) = c('all_id','efflen')
    tab_efflen = rbind(tab_efflen_gene,tab_efflen_trans)
  }
  if(length(grep('gene|trans',colnames(exp_count),value = TRUE))==0){
    exp_count[paste(geneOrTrans,'id',sep='_')] = rownames(exp_count)
  }
  list_sample = colnames(exp_count)[!grepl('gene|trans|all',colnames(exp_count))]
  tab_count_efflen = merge(exp_count,tab_efflen,by=paste(geneOrTrans,'id',sep='_'),all.x=TRUE)
  exp_tpm = tab_count_efflen[paste(geneOrTrans,'id',sep='_')]
  for(sample in list_sample){
    exp_tpm[sample] = countToTpm(tab_count_efflen[sample],tab_count_efflen['efflen'])
  }
  return(exp_tpm)
}

id2idTitle <- function(vec){
  return(unlist(lapply(vec, function(x) unlist(str_split(x,'\\.'))[1])))
}

exp_count_2_exp_tpm_path_idTitle <- function(exp_count, geneOrTrans, path_efflen){
  # path_efflen = '~/data/Ref/GENECODE_r45' # '~/data/Ref/GENECODE_r45', '~/data/Ref/ENSEMBL_r111'
  if(geneOrTrans=='gene'){tab_efflen = read.csv(paste(path_efflen,'tab_efflen_gene',sep='/'),sep='\t',check.names = F)}
  if(geneOrTrans=='transcript'){tab_efflen = read.csv(paste(path_efflen,'tab_efflen_trans',sep='/'),sep='\t',check.names = F)}
  if(geneOrTrans=='all'){
    tab_efflen_gene = read.csv(paste(path_efflen,'tab_efflen_gene',sep='/'),sep='\t',check.names = F)
    tab_efflen_trans = read.csv(paste(path_efflen,'tab_efflen_trans',sep='/'),sep='\t',check.names = F)
    colnames(tab_efflen_gene) = c('all_id','efflen')
    colnames(tab_efflen_trans) = c('all_id','efflen')
    tab_efflen = rbind(tab_efflen_gene,tab_efflen_trans)
  }
  tab_efflen = tab_efflen[!grepl('PAR_Y',tab_efflen[,1]),]
  tab_efflen[paste(colnames(tab_efflen)[1],'title',sep='_')] = id2idTitle(tab_efflen[[colnames(tab_efflen)[1]]])
  if(length(grep('gene|trans',colnames(exp_count),value = TRUE))==0){
    exp_count[paste(geneOrTrans,'id','title',sep='_')] = rownames(exp_count)
  }
  list_sample = colnames(exp_count)[!grepl('gene|trans|all',colnames(exp_count))]
  tab_count_efflen = merge(exp_count,tab_efflen,by=paste(geneOrTrans,'id','title',sep='_'),all.x=TRUE)
  exp_tpm = tab_count_efflen[paste(geneOrTrans,'id','title',sep='_')]
  for(sample in list_sample){
    exp_tpm[sample] = countToTpm(tab_count_efflen[sample],tab_count_efflen['efflen'])
  }
  return(exp_tpm)
}
# cnts <- c(4250, 3300, 200, 1750, 50, 0)
# lens <- c(900, 1020, 2000, 770, 3000, 1777)
# countDf <- data.frame(count = cnts, length = lens)
# 
# #assume a mean(FLD) = 203.7
# countDf$effLength <- countDf$length - 203.7 + 1
# countDf$tpm <- with(countDf, countToTpm(count, effLength))
# countDf$fpkm <- with(countDf, countToFpkm(count, effLength))
# with(countDf, all.equal(tpm, fpkmToTpm(fpkm)))
# countDf$effCounts <- with(countDf, countToEffCounts(count, length, effLength))
