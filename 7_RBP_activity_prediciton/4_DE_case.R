# DE 
library(edgeR)
library(limma)
library(DESeq2)
library(SummarizedExperiment) # using rowData
library(rtracklayer) # load gtf
library(tidyverse)

GEO_id = "GSE216825" # case_ISES:  GSE216825, GSE225633, GSE225637, GSE294883_NCI-H23
path_data = paste0('~/data/case/',GEO_id,'/count')
mat_exp_count = read.table(paste(path_data,'tab_gene_count.txt',sep='/'))
mat_exp_count = mat_exp_count[,c(2:dim(mat_exp_count)[2])]
vec_gene_id = row.names(mat_exp_count)
mat_exp_count <- as.data.frame(lapply(mat_exp_count, as.integer))
row.names(mat_exp_count) = vec_gene_id

# save path
path_DE = paste('~/result/RBP_activity_prediction/DE_case',GEO_id,sep='/')

if(GEO_id == 'GSE216825'){
  vec_sample1 = c('SRR22085277', 'SRR22085278', 'SRR22085279', 'SRR22085280', 'SRR22085281') # WT
  vec_sample2 = c('SRR22085273', 'SRR22085274', 'SRR22085275', 'SRR22085276') # MATR3 knockout
}
if(GEO_id == 'GSE225633'){
  vec_sample1 = c('SRR23557769', 'SRR23557770', 'SRR23557771') # WT
  vec_sample2 = c('SRR23557766', 'SRR23557767', 'SRR23557768') # KO
}
if(GEO_id == 'GSE225637'){
  vec_sample1 = c('SRR23558487', 'SRR23558489') # WT
  vec_sample2 = c('SRR23558483', 'SRR23558484', 'SRR23558485', 'SRR23558486') # KO
}
if(GEO_id == 'GSE294883_NCI-H23'){
  vec_sample1 = c('SRR33188891', 'SRR33188892', 'SRR33188893') # WT
  vec_sample2 = c('SRR33188882', 'SRR33188883', 'SRR33188890') # TRA2A Knockout
}
list_samples = c(vec_sample1,vec_sample2)


# create label and group
term1 = 'WT'
term2 = 'EXP'
tab_label = data.frame(sample_id = list_samples, group = relevel(factor(c(rep('WT',length(vec_sample1)),
                                                                          rep('EXP',length(vec_sample2)))),ref='WT'))
rownames(tab_label) = tab_label$sample_id
group <- tab_label$group

################################################################################
# edgeR differential expression
################################################################################
# logFC = log(MUT/WT)
dgelist <- DGEList(counts = mat_exp_count, group = group)
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
write.table(tab_DE_edgeR, file=paste(path_DE,paste('tab_DE', 'edgeR',sep='_'),sep='/'), row.names = TRUE,sep="\t",quote=FALSE)

################################################################################
# limma differential expression
# https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf
# Page71
################################################################################
design <- model.matrix(~factor(group)+0)  
colnames(design) <- c("WT", "KO")
dgelist <- DGEList(counts = mat_exp_count, group = group)
keep <- filterByExpr(dgelist) # edgeR
dgelist <- dgelist[keep, , keep.lib.sizes = FALSE]
dgelist <- calcNormFactors(dgelist)
mat_limma<- voom(dgelist, design) # voom norm
df.fit <- lmFit(mat_limma, design)
df.matrix <- makeContrasts(KO - WT , levels = design)
fit <- contrasts.fit(df.fit, df.matrix)
fit <- eBayes(fit)
tab_DE_limma <- topTable(fit,n = Inf,adjust="fdr")
write.table(tab_DE_limma, file=paste(path_DE,paste('tab_DE', 'limma',sep='_'),sep='/'), row.names = TRUE,sep="\t",quote=FALSE)

################################################################################
# DEseq2 differential expression
################################################################################
smallestGroupSize <- 2
keep <- rowSums(mat_exp_count >= 10) >=smallestGroupSize
cts <- mat_exp_count[keep,]
coldata = data.frame(condition=group)
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)
tab_DE_DEseq2 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
tab_DE_DEseq2 <- tab_DE_DEseq2[order(tab_DE_DEseq2$pvalue,decreasing = FALSE), ]
write.table(tab_DE_DEseq2, file=paste(path_DE,paste('tab_DE', 'DEseq2',sep='_'),sep='/'), row.names = TRUE,sep="\t",quote=FALSE)


################################################################################
# merge result of edgeR and limma, and add expresion
################################################################################
path_DE = paste('~/result/RBP_activity_prediction/DE_case',GEO_id,sep='/')
list_m_DE = c('edgeR','limma','DEseq2')

# load gene annotation
path_anno = '~/data/Ref/GENECODE_r45'
table_anno_gene = read.table(paste(path_anno,'tab_anno_gene',sep='/'),sep='\t',header=T)
source('~/utils/public/different_expression_func.R')

tab_DE = data.frame()
for(m_DE in list_m_DE){
  str_title = paste(m_DE,sep='_')
  tab_DE_ori = read.table(paste(path_DE,paste('tab_DE',m_DE,sep='_'),sep='/'),sep='\t',header=T)
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

write.table(tab_DE, file=paste(path_DE,paste('tab_DE','all','anno',sep='_'),sep='/'), row.names = FALSE,sep="\t",quote=FALSE)

#############################################################################
# filter RBP
path_support_NGS = '~/result/support_file_shortRead'
vec_dataSubTypeTermRBP = readRDS(file.path(path_support_NGS,'vec_bioTypeTerm1Term2',  "ENCODE_NGS.rds"))
vec_RBP = unique(unlist(lapply(vec_dataSubTypeTermRBP, function(x) unlist(str_split(x, '_'))[3])))
tab_DE_RBP = tab_DE %>% filter(gene_name %in% vec_RBP)
write.table(tab_DE_RBP, file=paste(path_DE,paste('tab_DE','all','anno','RBP',sep='_'),sep='/'), row.names = FALSE,sep="\t",quote=FALSE)
