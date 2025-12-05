library(parallel)
library(edgeR)
library(limma)
library(DESeq2)
library(SummarizedExperiment) # using rowData
library(rtracklayer) # load gtf
library(tidyverse)

library(conflicted)
library(dplyr)
conflicts_prefer(dplyr::filter())
conflicts_prefer(base::intersect)

path_data = '~/data'
path_save = '~/result/RBP_activity_prediction'
dataType= 'ENCODE_NGS' 

# ENCODE_NGS 
list_tab_report_ep = list()
list_tab_meta_ep_gene = list()
list_tab_count_ep_gene = list()
list_tab_report_cl = list()
list_tab_meta_cl_gene = list()
list_tab_count_cl_gene = list()
for(dataSubType in c('CRISPR','shRNA')){
  # experiment
  tab_report_ep = read.csv(paste(path_data,dataType,dataSubType,'experiment_report.tsv',sep='/'),check.names = F,sep='\t',skip=1)
  tab_meta_ep = read.csv(paste(path_data,dataType,dataSubType,'metadata.tsv',sep='/'),check.names = F,sep='\t')
  tab_meta_ep_gene = tab_meta_ep %>% filter(`Genome annotation`=='V29') %>% 
    filter(`Output type` == 'gene quantifications') 
  tab_meta_ep_gene = tab_meta_ep_gene %>% arrange(`Biosample type`,`Biosample term name`,`Biosample term id`,`Experiment target`,`Experiment accession`,`Biological replicate(s)`,`Technical replicate(s)`)
  tab_count_ep_gene = read.csv(paste(path_data,dataType,dataSubType,'tab_count_merge_gene',sep='/'),sep='\t',header=T,check.names = F)
  tab_count_ep_gene[is.na(tab_count_ep_gene)]=0
  tab_count_ep_gene = tab_count_ep_gene %>% filter(!grepl('PAR_Y',gene_id))
  
  list_tab_report_ep[[dataSubType]] = tab_report_ep
  list_tab_meta_ep_gene[[dataSubType]] = tab_meta_ep_gene
  list_tab_count_ep_gene[[dataSubType]] = tab_count_ep_gene
  
  # control
  tab_report_cl = read.csv(paste(path_data,dataType,paste('cl',dataSubType,sep='_'),'experiment_report.tsv',sep='/'),check.names = F,sep='\t',skip=1)
  tab_meta_cl = read.csv(paste(path_data,dataType,paste('cl',dataSubType,sep='_'),'metadata.tsv',sep='/'),check.names = F,sep='\t')
  tab_meta_cl_gene = tab_meta_cl %>% filter(`Genome annotation`=='V29') %>% 
    filter(`Output type` == 'gene quantifications')
  tab_meta_cl_gene = tab_meta_cl_gene %>% arrange(`Biosample type`,`Biosample term name`,`Biosample term id`,`Experiment target`,`Experiment accession`,`Biological replicate(s)`,`Technical replicate(s)`)
  tab_count_cl_gene = read.csv(paste(path_data,dataType,paste('cl',dataSubType,sep='_'),'tab_count_merge_gene',sep='/'),sep='\t',header=T,check.names = F)
  tab_count_cl_gene[is.na(tab_count_cl_gene)]=0
  tab_count_cl_gene = tab_count_cl_gene %>% filter(!grepl('PAR_Y',gene_id))
  
  list_tab_report_cl[[dataSubType]] = tab_report_cl
  list_tab_meta_cl_gene[[dataSubType]] = tab_meta_cl_gene
  list_tab_count_cl_gene[[dataSubType]] = tab_count_cl_gene
}

# load gene annotation
path_anno = '~/data/Ref/GENECODE_r29_ENCODE'
table_anno_gene = read.table(paste(path_anno,'tab_anno_gene',sep='/'),sep='\t',header=T)

path_support_NGS = '~/result/support_file_shortRead'
vec_dataSubTypeTermRBP = readRDS(file.path(path_support_NGS,'vec_bioTypeTerm1Term2',  "ENCODE_NGS.rds"))
vec_RBP = unique(unlist(lapply(vec_dataSubTypeTermRBP, function(x) unlist(str_split(x, '_'))[3])))

# combin = combinList[1]
DE_multicore <-function(combin, table_anno_gene, vec_RBP){
  conflicts_prefer(dplyr::filter())
  conflicts_prefer(base::intersect)
  
  combins = unlist(strsplit(combin,split='_'))
  dataSubType = combins[1]
  term1 = combins[2]
  term2 = combins[3]
  
  tab_report_ep = list_tab_report_ep[[dataSubType]]
  tab_meta_ep_gene = list_tab_meta_ep_gene[[dataSubType]] 
  tab_count_ep_gene = list_tab_count_ep_gene[[dataSubType]]
  
  tab_report_cl = list_tab_report_cl[[dataSubType]]
  tab_meta_cl_gene = list_tab_meta_cl_gene[[dataSubType]] 
  tab_count_cl_gene = list_tab_count_cl_gene[[dataSubType]]
  gc()
  
  tab_meta_ep_gene_term_RBP = tab_meta_ep_gene %>% 
    filter(`Biosample term name`==term1 & `Experiment target`== term2)
  ## ep
  vec_s_ep_gene_term_RBP_ori = tab_meta_ep_gene_term_RBP %>% pull(`File accession`)
  ## cl
  access_ep_gene = tab_meta_ep_gene_term_RBP$`Experiment accession`[1]
  series_ep_gene = str_extract(tab_report_ep %>% filter(Accession == access_ep_gene) %>% pull(`Related series`),
                                "/gene-silencing-series/[^,]+")
  access_cl_gene = tab_report_cl %>% filter(grepl(series_ep_gene,`Related series`)) %>% pull(Accession)
  vec_s_cl_gene_term_RBP_ori = tab_meta_cl_gene %>% filter(`Experiment accession` == access_cl_gene) %>% pull(`File accession`)
  
  vec_sample1 = vec_s_cl_gene_term_RBP_ori
  vec_sample2 = vec_s_ep_gene_term_RBP_ori
  list_samples = c(vec_sample1,vec_sample2)
  
  tab_count = merge(tab_count_cl_gene[,c('gene_id',vec_sample1)],
                    tab_count_ep_gene[,c('gene_id',vec_sample2)],by=c('gene_id'),all=T)
  tab_count[is.na(tab_count)]=0
  vec_gene_id = tab_count$gene_id
  # to int
  tab_count <- as.data.frame(lapply(tab_count[,list_samples], as.integer))
  rownames(tab_count) = vec_gene_id
  
  # filter, important !! 
  tab_count = tab_count[rowSums(tab_count[,-c(1:2)]>5)>=2,]
  
  # create label
  term1_change = gsub('-| ','_',term1)
  term2_change = gsub('-| ','_',term2)
  tab_label = data.frame(sample_id = list_samples, group = relevel(factor(c(rep(term1_change,length(vec_sample1)),
                                                                            rep(term2_change,length(vec_sample2)))),ref=term1_change))
  # tab_label = data.frame(sample_id = list_samples, group = factor(c(rep(term1,length(vec_sample1)),rep(term2,length(vec_sample2)))))
  rownames(tab_label) = tab_label$sample_id
  group <- tab_label$group
  
  ################################################################################
  # different expression
  path_DE = paste(path_save,paste('DE',dataType,sep='_'),dataSubType,sep='/')
  
  ################################################################################
  # edgeR differential expression
  ################################################################################
  # logFC = log(MUT/WT)
  dgelist <- DGEList(counts = tab_count, group = group)
  # filter low expression
  # keep <- rowSums(cpm(dgelist) > 1 ) >= 2
  keep <- filterByExpr(dgelist)
  dgelist <- dgelist[keep, , keep.lib.sizes = FALSE]
  dgelist <- calcNormFactors(dgelist) # norm
  design <- model.matrix(~group) # dis
  dgelist <- estimateDisp(dgelist,design)
  fit <- glmFit(dgelist, design, robust = TRUE)# fit
  lrt <- glmLRT(fit,coef=2)
  tab_DE_edgeR = topTags(lrt, n = nrow(dgelist$counts))$table # output
  write.table(tab_DE_edgeR, file=paste(path_DE,paste('tab_DE',term1,term2,'edgeR',sep='_'),sep='/'), row.names = TRUE,sep="\t",quote=FALSE)
  
  ################################################################################
  # limma differential expression
  # https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf
  # Page71
  ################################################################################
  design <- model.matrix(~factor(group)+0)  
  colnames(design) <- c("WT", "KO")
  dgelist <- DGEList(counts = tab_count, group = group)
  keep <- filterByExpr(dgelist) # edgeR
  dgelist <- dgelist[keep, , keep.lib.sizes = FALSE]
  dgelist <- calcNormFactors(dgelist)
  mat_limma<- voom(dgelist, design) # voom norm
  df.fit <- lmFit(mat_limma, design)
  df.matrix <- makeContrasts(KO - WT , levels = design)
  fit <- contrasts.fit(df.fit, df.matrix)
  fit <- eBayes(fit)
  tab_DE_limma <- topTable(fit,n = Inf,adjust="fdr")
  write.table(tab_DE_limma, file=paste(path_DE,paste('tab_DE',term1,term2,'limma',sep='_'),sep='/'), row.names = TRUE,sep="\t",quote=FALSE)
  
  ################################################################################
  # DEseq2 differential expression
  ################################################################################
  smallestGroupSize <- 2
  keep <- rowSums(tab_count >= 10) >=smallestGroupSize
  cts <- tab_count[keep,]
  coldata = data.frame(condition=group)
  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = coldata,
                                design = ~ condition)
  dds <- DESeq(dds)
  res <- results(dds)
  tab_DE_DEseq2 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
  tab_DE_DEseq2 <- tab_DE_DEseq2[order(tab_DE_DEseq2$pvalue,decreasing = FALSE), ]
  write.table(tab_DE_DEseq2, file=paste(path_DE,paste('tab_DE',term1,term2,'DEseq2',sep='_'),sep='/'), row.names = TRUE,sep="\t",quote=FALSE)
  
  
  ################################################################################
  # merge result of edgeR and limma, and add expresion
  ################################################################################
  list_m_DE = c('edgeR','limma','DEseq2')
  source('~/utils/public/different_expression_func.R')
  
  tab_DE = data.frame()
  for(m_DE in list_m_DE){
    str_title = paste(m_DE,sep='_')
    tab_DE_ori = read.table(paste(path_DE,paste('tab_DE',term1, term2, m_DE,sep='_'),sep='/'),sep='\t',header=T)
    # tab_DE = tab_DE_ori
    if(m_DE=='edgeR'){tab_DE_temp = tab_DE_ori %>% select(logFC, PValue, FDR)}
    if(m_DE=='limma'){tab_DE_temp = tab_DE_ori %>% select(logFC, P.Value, adj.P.Val)}
    if(m_DE=='DEseq2'){tab_DE_temp = tab_DE_ori %>% select(log2FoldChange, pvalue, padj)}
    colnames(tab_DE_temp) = paste(str_title,c('logFC','pvalue','FDR'),sep='_')
    tab_DE_temp['name'] = rownames(tab_DE_temp)
    if(dim(tab_DE)[2]==0){tab_DE = tab_DE_temp}else{tab_DE = merge(tab_DE,tab_DE_temp,by='name',all=TRUE)}
  }
  tab_DE = column_to_rownames(tab_DE, var = "name")
  tab_DE[,grep('logFC', colnames(tab_DE),value = TRUE)] = tab_DE[,grep('logFC', colnames(tab_DE),value = TRUE)] %>% replace(is.na(.), 0)
  tab_DE[,grep('pvalue', colnames(tab_DE),value = TRUE)] = tab_DE[,grep('pvalue', colnames(tab_DE),value = TRUE)] %>% replace(is.na(.), 1)
  tab_DE[,grep('FDR', colnames(tab_DE),value = TRUE)] = tab_DE[,grep('FDR', colnames(tab_DE),value = TRUE)] %>% replace(is.na(.), 1)
  
  ################################################################################
  # gene and trans annotation
  ################################################################################
  
  # add annotation
  tab_DE['gene_id'] = rownames(tab_DE)
  tab_DE = merge(tab_DE,table_anno_gene, by='gene_id', all.x=T)
  tab_DE = tab_DE[,c(grep('gene',colnames(tab_DE),value=TRUE),colnames(tab_DE)[!grepl('gene',colnames(tab_DE))])]
  tab_DE = add_combin_columns_DE(tab_DE)
  tab_DE = tab_DE %>% arrange(combinedPvalue)
  
  write.table(tab_DE, file=paste(path_DE,paste('tab_DE',term1,term2,'all','anno',sep='_'),sep='/'), row.names = FALSE,sep="\t",quote=FALSE)
  
  #############################################################################
  # filter RBP
  tab_DE_RBP = tab_DE %>% filter(gene_name %in% vec_RBP)
  write.table(tab_DE_RBP, file=paste(path_DE,paste('tab_DE',term1,term2,'all','anno','RBP',sep='_'),sep='/'), row.names = FALSE,sep="\t",quote=FALSE)
  
  gc()
}


# Load RBP info
tab_RBPinfo <- read.csv('~/data/RBP_RNA/CLIPdb/POSTAR3_CLIPdb_module_browse_RBP_sub_info.csv', check.names = FALSE)
colnames(tab_RBPinfo) <- c('RBPname', 'geneID', 'Domain', 'Location', 'Function')
tab_RBPinfo_splice <- tab_RBPinfo[grepl('Splic', tab_RBPinfo$Function), ]
vec_RBP_RBPinfo <- tab_RBPinfo_splice[['RBPname']]

combinList = c()
for(dataSubType in c('CRISPR','shRNA')){
  for(bioCell in c('HepG2', 'K562')){
    vec_term = sort(tab_meta_ep %>% filter(`Biosample term name`==bioCell) %>% pull(`Experiment target`) %>% unique())
    term1 = bioCell
    
    # load meta trans
    tab_meta_ep_gene_term = list_tab_meta_ep_gene[[dataSubType]] %>% filter(`Biosample term name`==term1)
    vec_RBP_ep_gene = unlist(lapply(tab_meta_ep_gene_term[['Experiment target']], function(x) unlist(strsplit(x,'-'))[1]))
    
    # load ENCORE meta data
    tab_meta_ENCORE = readxl::read_excel(paste0('~/data/RBP_RNA/ENCORE/',term1,'/metadata_mergedbed.xlsx'))
    vec_RBP_ENCORE = unlist(lapply(tab_meta_ENCORE[['Experiment target']], function(x) unlist(strsplit(x,'-'))[1]))
    
    # intersect among vec_RBP_RBPinfo, vec_RBP_ENCORE and vec_RBP_ep
    vec_RBP = intersect(intersect(vec_RBP_RBPinfo,vec_RBP_ENCORE),vec_RBP_ep_gene)
    vec_term = paste(vec_RBP,'-human',sep='')
    
    for(term2 in vec_term){
      combinList = c(combinList, paste(dataSubType, term1,term2,sep='_')) 
    }
  }
}


mclapply(combinList, function(combin) {
  suppressPackageStartupMessages({
    library(parallel)
    library(edgeR)
    library(limma)
    library(DESeq2)
    library(SummarizedExperiment) # using rowData
    library(rtracklayer) # load gtf
    library(tidyverse)
    library(conflicted)
    library(dplyr)
  })
  DE_multicore(combin, table_anno_gene, vec_RBP)
},
mc.cores = 30, 
mc.preschedule = FALSE,  
mc.cleanup = TRUE        
)
gc()



# rank 
path_save = '~/result/RBP_activity_prediction'
dataType= 'ENCODE_NGS' 

path_support_NGS = '~/result/support_file_shortRead'
vec_dataSubTypeTermRBP = readRDS(file.path(path_support_NGS,'vec_bioTypeTerm1Term2',  "ENCODE_NGS.rds"))
vec_RBP = unique(unlist(lapply(vec_dataSubTypeTermRBP, function(x) unlist(str_split(x, '_'))[3])))

list_tab_DE_RBP = list()
for(combin in vec_dataSubTypeTermRBP){
  combins = unlist(strsplit(combin,split='_'))
  dataSubType = combins[1]
  term1 = combins[2]
  term2 = combins[3]
  path_DE = paste(path_save,paste('DE',dataType,sep='_'),dataSubType,sep='/')
  tab_DE_RBP = read.table(paste(path_DE,paste('tab_DE',term1,paste0(term2,'-human'),'all','anno','RBP',sep='_'),sep='/'),header=T)
  tab_DE_RBP = tab_DE_RBP %>% mutate(dataSubType=dataSubType, term1=term1, term2=term2, .before=1)
  list_tab_DE_RBP[[combin]] = tab_DE_RBP
}
tab_DE_RBP_all = bind_rows(list_tab_DE_RBP)

ranked_data <- tab_DE_RBP_all %>%
  group_by(dataSubType, term1, term2) %>%
  mutate(
    edgeR_rank = rank(edgeR_pvalue, ties.method = "min"),
    limma_rank = rank(limma_pvalue, ties.method = "min"),
    DEseq2_rank = rank(DEseq2_pvalue, ties.method = "min"),
    combined_rank = rank(combinedPvalue, ties.method = "min")
  ) %>%
  ungroup()

