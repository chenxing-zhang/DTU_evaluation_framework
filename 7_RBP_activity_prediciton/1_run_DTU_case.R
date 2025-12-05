library(tidyverse)
library(DRIMSeq)
library(DEXSeq)
library(limma)
library(edgeR)
library(satuRn)
library(parallel)
library(data.table) # for RATs
library(rats) # RATs
library(fishpond) # Swish
library(BiocParallel) # for DTUrtle
library(DTUrtle) # DTUrtle
library(NBSplice) # NBSplice
library(reticulate) 
use_condaenv("python", required = TRUE)
source_python("~/utils/DTU/SUPPA2/suppa_helpers.py") # SUPPA2

id2idTitle <- function(vec) {
  return(unlist(lapply(vec, function(x) unlist(str_split(x, '\\.'))[1])))
}

GEO_id = "GSE216825" # case_ISES: GSE216825, GSE225633, GSE225637, GSE294883_NCI-H23

# merge featureCounts count
library(data.table)
path_data = paste0('~/data/case/',GEO_id,'/count')

# gene, count
files = list.files(paste0(path_data), "*.genes.results$", full.names = T)
countData = data.frame(fread(files[1]),check.names = FALSE)[c(1,5)] # count
for(i in 2:length(files)) {countData = cbind(countData, data.frame(fread(files[i]))[5]) }
colnames(countData) = c("gene_id", gsub(paste0(path_data,'/'), "", files))
colnames(countData) = gsub(".genes.results", "", colnames(countData))
rownames(countData) = countData$gene_id
write.table(countData, file=paste0(path_data,"/tab_gene_count.txt"), quote=F,sep='\t',row.names=T, col.names=T)

# gene, tpm
files = list.files(paste0(path_data), "*.genes.results$", full.names = T)
countData = data.frame(fread(files[1]),check.names = FALSE)[c(1,6)]
for(i in 2:length(files)) {countData = cbind(countData, data.frame(fread(files[i]))[6]) }
colnames(countData) = c("gene_id", gsub(paste0(path_data,'/'), "", files))
colnames(countData) = gsub(".genes.results", "", colnames(countData))
rownames(countData) = countData$gene_id
write.table(countData, file=paste0(path_data,"/tab_gene_TPM.txt"), quote=F,sep='\t',row.names=T, col.names=T)

# isoform, count
files = list.files(paste0(path_data), "*.isoforms.results$", full.names = T)
countData = data.frame(fread(files[1]),check.names = FALSE)[c(1,2,5)] # count
for(i in 2:length(files)) {countData = cbind(countData, data.frame(fread(files[i]))[5]) }
colnames(countData) = c("transcript_id",'gene_id', gsub(paste0(path_data,'/'), "", files))
colnames(countData) = gsub(".isoforms.results", "", colnames(countData))
rownames(countData) = countData$transcript_id
write.table(countData, file=paste0(path_data,"/tab_isoform_count.txt"), quote=F,sep='\t',row.names=T, col.names=T)

# isoform, TPM
files = list.files(paste0(path_data), "*.isoforms.results$", full.names = T)
countData = data.frame(fread(files[1]),check.names = FALSE)[c(1,2,6)] # TPM
for(i in 2:length(files)) {countData = cbind(countData, data.frame(fread(files[i]))[6]) }
colnames(countData) = c("transcript_id",'gene_id', gsub(paste0(path_data,'/'), "", files))
colnames(countData) = gsub(".isoforms.results", "", colnames(countData))
rownames(countData) = countData$transcript_id
write.table(countData, file=paste0(path_data,"/tab_isoform_TPM.txt"), quote=F,sep='\t',row.names=T, col.names=T)

################################################################################
# # load count
path_save = paste0('~/result/RBP_activity_prediction/DTU_case/',GEO_id)
if(!dir.exists(path_save)){dir.create(path_save)}

path_data = paste0('~/data/case/',GEO_id,'/count')

tab_count = read.table(paste(path_data,'tab_isoform_count.txt',sep='/'))
tab_tpm_trans  = read.table(paste(path_data,'tab_isoform_TPM.txt',sep='/'))

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
tab_count[is.na(tab_count)]=0
colnames(tab_count) = c('feature_id','gene_id',list_samples)
colnames(tab_tpm_trans) = c('feature_id','gene_id',list_samples)
# filter, important !! 
tab_count = tab_count[rowSums(tab_count[,-c(1:2)]>5)>=2,]
tab_tpm_trans = tab_tpm_trans[rownames(tab_count),]
tab_count[is.na(tab_count)]=0
tab_tpm_trans[is.na(tab_tpm_trans)]=0

# create label
term1 = 'WT'
term2 = 'EXP'
tab_label = data.frame(sample_id = list_samples, group = relevel(factor(c(rep('WT',length(vec_sample1)),
                                                                          rep('EXP',length(vec_sample2)))),ref='WT'))
rownames(tab_label) = tab_label$sample_id


################################################################################
# different transcriptome
DTU <- function(){
  use_condaenv("python", required = TRUE)
  dataType = GEO_id
  
  # DRIMSeq
  if(!file.exists(paste(path_save,paste('tab_DTU',term1,term2,'DRIMSeq',sep='_'),sep='/'))){
    start_time <- Sys.time()
    tab_count_temp = tab_count
    tab_count_temp <- tab_count_temp %>%
      mutate(across(3:ncol(.), round)) # 20000
    d <- DRIMSeq::dmDSdata(counts = tab_count_temp, samples = tab_label)
    # d <- DRIMSeq::dmFilter(d, min_samps_gene_expr = 2, min_samps_feature_expr = 2,
    #                        min_gene_expr = 5, min_feature_expr = 5)
    group <- relevel(factor(tab_label$group),ref=term1)
    design <- model.matrix(~group)
    d <- DRIMSeq::dmPrecision(d, design = design)
    d <- DRIMSeq::dmFit(d, design = design, verbose = 1)
    d <- DRIMSeq::dmTest(d, coef = paste0("group",term2), verbose = 1)
    tab_prop_DRIMSeq = DRIMSeq::proportions(d)
    list_prop_term1 = rowMeans(tab_prop_DRIMSeq[,3:(3+length(vec_sample1)-1)])
    list_prop_term2 = rowMeans(tab_prop_DRIMSeq[,(3+length(vec_sample1)):(3+length(vec_sample1)+length(vec_sample2)-1)])
    list_prop_log2fc = log2(list_prop_term2/list_prop_term1)
    tab_result_DRIMSeq = DRIMSeq::results(d, level = "feature") # 'feature','gene'
    tab_result_DRIMSeq = tab_result_DRIMSeq %>% mutate(prop_term1=list_prop_term1,
                                                       prop_term2=list_prop_term2,
                                                       prop_log2fc = list_prop_log2fc,.before=lr)
    
    end_time <- Sys.time()
    all_time = as.numeric(end_time) - as.numeric(start_time)
    tab_time_DRIMSeq = data.frame(dataType=dataType,
                                  term1=term1,term2=term2,
                                  m_DTU = 'DRIMSeq',time = all_time)
    write.table(tab_time_DRIMSeq, paste(path_save,paste('tab_time',term1,term2,'DRIMSeq',sep='_'),sep='/'), row.names = TRUE,sep="\t",quote=FALSE)
    
    write.table(tab_result_DRIMSeq, paste(path_save,paste('tab_DTU',term1,term2,'DRIMSeq',sep='_'),sep='/'), row.names = TRUE,sep="\t",quote=FALSE)
  }
  gc()
  
  # SPIT
  if(!file.exists(paste(path_save,paste('tab_DTU',term1,term2,'SPIT',sep='_'),sep='/'))){
    source_python("~/utils/DTU/SPIT/dtu_detection_zcx.py")
    start_time <- Sys.time()
    list_result_SPIT = compute_double_stats_zcx(tab_count[,c('gene_id',vec_sample1,vec_sample2)],
                                                tab_count[,c('feature_id','gene_id')],
                                                vec_sample1, vec_sample2)
    tab_result_SPIT = as.data.frame(list_result_SPIT[1:4])
    colnames(tab_result_SPIT) = c('transcript_id','gene_id','likelihood','pvalue')
    end_time <- Sys.time()
    all_time = as.numeric(end_time) - as.numeric(start_time)
    tab_time_SPIT = data.frame(dataType=dataType,
                               term1=term1,term2=term2,
                               m_DTU = 'SPIT',time = all_time)
    write.table(tab_time_SPIT, paste(path_save,paste('tab_time',term1,term2,'SPIT',sep='_'),sep='/'), row.names = TRUE,sep="\t",quote=FALSE)
    write.table(tab_result_SPIT, paste(path_save,paste('tab_DTU',term1,term2,'SPIT',sep='_'),sep='/'), row.names = TRUE,sep="\t",quote=FALSE)
  }
  
  # SUPPA2
  if(!file.exists(paste(path_save,paste('tab_DTU',term1,term2,'SUPPA2',sep='_'),sep='/'))){
    source_python("~/utils/DTU/SUPPA2/suppa_helpers.py") # SUPPA2
    start_time <- Sys.time()
    ioi_data = tab_tpm_trans[,c('feature_id','gene_id')]
    colnames(ioi_data) = c('transcript', 'gene')
    psi_output_A <- psiPerIsoform_direct(tab_tpm_trans[,vec_sample1], ioi_data)
    psi_output_B <- psiPerIsoform_direct(tab_tpm_trans[,vec_sample2], ioi_data)
    tab_result_SUPPA2 <- diffSplice_transcript_direct(psi_output_A, psi_output_B)
    tab_result_SUPPA2['transcript'] = rownames(tab_result_SUPPA2)
    tab_result_SUPPA2 = merge(tab_result_SUPPA2 ,ioi_data, by='transcript', all.x=T)
    tab_result_SUPPA2 = tab_result_SUPPA2[,c('transcript','gene','dPSI','p_value')]
    end_time <- Sys.time()
    all_time = as.numeric(end_time) - as.numeric(start_time)
    tab_time_SUPPA2 = data.frame(dataType=dataType,
                                 term1=term1,term2=term2,
                                 m_DTU = 'SUPPA2',time = all_time)
    write.table(tab_time_SUPPA2, paste(path_save,paste('tab_time',term1,term2,'SUPPA2',sep='_'),sep='/'), row.names = TRUE,sep="\t",quote=FALSE)
    write.table(tab_result_SUPPA2, paste(path_save,paste('tab_DTU',term1,term2,'SUPPA2',sep='_'),sep='/'), row.names = TRUE,sep="\t",quote=FALSE)
  }
  
  # RATs
  if(!file.exists(paste(path_save,paste('tab_DTU',term1,term2,'RATs',sep='_'),sep='/'))){
    start_time <- Sys.time()
    mydtu <- rats::call_DTU(annot= tab_tpm_trans[,c('feature_id','gene_id')],
                            TARGET_COL = "feature_id", PARENT_COL = "gene_id",
                            count_data_A= setDT(tab_tpm_trans[,c('feature_id',vec_sample1)]),
                            count_data_B= setDT(tab_tpm_trans[,c('feature_id',vec_sample2)]),
                            abund_thresh = 0, dprop_thresh=0, verbose= FALSE)
    tab_result_RATs = mydtu[["Transcripts"]]
    end_time <- Sys.time()
    all_time = as.numeric(end_time) - as.numeric(start_time)
    tab_time_RATs = data.frame(dataType=dataType,
                               term1=term1,term2=term2,
                               m_DTU = 'RATs',time = all_time)
    write.table(tab_time_RATs, paste(path_save,paste('tab_time',term1,term2,'RATs',sep='_'),sep='/'), row.names = TRUE,sep="\t",quote=FALSE)
    write.table(tab_result_RATs, paste(path_save,paste('tab_DTU',term1,term2,'RATs',sep='_'),sep='/'), row.names = TRUE,sep="\t",quote=FALSE)
  }
  
  # DTUrtle
  if(!file.exists(paste(path_save,paste('tab_DTU',term1,term2,'DTUrtle',sep='_'),sep='/'))){
    start_time <- Sys.time()
    biocpar <- BiocParallel::SerialParam()
    dturtle <- DTUrtle::run_drimseq(counts = as.matrix(tab_count[,c(vec_sample1,vec_sample2)]),
                                    tx2gene = tab_count[,c('feature_id','gene_id')],
                                    pd=tab_label, id_col = "sample_id",cond_col = "group",
                                    filtering_strategy = "own",
                                    BPPARAM = biocpar,
                                    min_samps_gene_expr=0, min_samps_feature_expr=0,
                                    min_samps_feature_prop=0,min_gene_expr=0,
                                    min_feature_expr=0, min_feature_prop=0,run_gene_twice=F)
    tab_result_DTUrtle = dturtle[["drim"]]@results_feature
    end_time <- Sys.time()
    all_time = as.numeric(end_time) - as.numeric(start_time)
    tab_time_DTUrtle = data.frame(dataType=dataType,
                                  term1=term1,term2=term2,
                                  m_DTU = 'DTUrtle',time = all_time)
    write.table(tab_time_DTUrtle, paste(path_save,paste('tab_time',term1,term2,'DTUrtle',sep='_'),sep='/'), row.names = TRUE,sep="\t",quote=FALSE)
    write.table(tab_result_DTUrtle, paste(path_save,paste('tab_DTU',term1,term2,'DTUrtle',sep='_'),sep='/'), row.names = TRUE,sep="\t",quote=FALSE)
  }
  
  # NBSplice
  if(!file.exists(paste(path_save,paste('tab_DTU',term1,term2,'NBSplice',sep='_'),sep='/'))){
    start_time <- Sys.time()
    geneIso = tab_count[,c('gene_id','feature_id')]
    colnames(geneIso) = c('gene_id','isoform_id')
    colName = 'group'
    myIsoDataSet<-NBSplice::IsoDataSet(tab_count[,c(vec_sample1,vec_sample2)],
                                       tab_label, colName, geneIso)
    myIsoDataSet<-buildLowExpIdx(myIsoDataSet, colName, ratioThres = 0.01,countThres = 1)
    myDSResults<-NBTest(myIsoDataSet, colName, test="F")
    tab_result_NBSplice = myDSResults@results
    end_time <- Sys.time()
    all_time = as.numeric(end_time) - as.numeric(start_time)
    tab_time_NBSplice = data.frame(dataType=dataType,
                                   term1=term1,term2=term2,
                                   m_DTU = 'NBSplice',time = all_time)
    write.table(tab_time_NBSplice, paste(path_save,paste('tab_time',term1,term2,'NBSplice',sep='_'),sep='/'), row.names = TRUE,sep="\t",quote=FALSE)
    write.table(tab_result_NBSplice, paste(path_save,paste('tab_DTU',term1,term2,'NBSplice',sep='_'),sep='/'), row.names = TRUE,sep="\t",quote=FALSE)
    
  }
  
  
  # DEXSeq
  if(!file.exists(paste(path_save,paste('tab_DTU',term1,term2,'DEXSeq',sep='_'),sep='/'))){
    start_time <- Sys.time()
    group <- relevel(factor(tab_label$group),ref=term1)
    dxd <- DEXSeq::DEXSeqDataSet(countData = round(as.matrix(tab_count[,-c(1:2)])),
                                 sampleData=tab_label,
                                 design=~sample+exon+group:exon,
                                 featureID=tab_count$feature_id,
                                 groupID=tab_count$gene_id)
    dxd <- DEXSeq::estimateSizeFactors(dxd)
    dxd <- DEXSeq::estimateDispersions(dxd)
    dxd <- DEXSeq::testForDEU(dxd, reducedModel = ~sample+exon)
    dxd = DEXSeq::estimateExonFoldChanges(dxd,fitExpToVar = "group")
    dxr = DEXSeq::DEXSeqResults(dxd, independentFiltering = FALSE)
    tab_result_DEXSeq = as.data.frame(dxr[,c("featureID", "groupID", term2, term1 ,paste("log2fold",term2,term1,sep='_'),"pvalue", "padj")])
    
    end_time <- Sys.time()
    all_time = as.numeric(end_time) - as.numeric(start_time)
    tab_time_DEXSeq = data.frame(dataType=dataType,
                                 term1=term1,term2=term2,
                                 m_DTU = 'DEXSeq',time = all_time)
    write.table(tab_time_DEXSeq, paste(path_save,paste('tab_time',term1,term2,'DEXSeq',sep='_'),sep='/'), row.names = TRUE,sep="\t",quote=FALSE)
    
    write.table(tab_result_DEXSeq, paste(path_save,paste('tab_DTU',term1,term2,'DEXSeq',sep='_'),sep='/'), row.names = TRUE,sep="\t",quote=FALSE)
  }
  gc()
  
  # limma-DTU
  if(!file.exists(paste(path_save,paste('tab_DTU',term1,term2,'limma',sep='_'),sep='/'))){
    start_time <- Sys.time()
    y_limma <- DGEList(counts = as.matrix(tab_count[,-c(1:2)]),
                       samples = tab_label,
                       genes = tab_count[,c(1,2)])
    rownames(y_limma) <- tab_count$feature_id
    y_limma <- calcNormFactors(y_limma)
    group <- relevel(factor(tab_label$group),ref=term1)
    design <- model.matrix(~group)
    y_limma <- voom(y_limma, design, plot=F)
    fit_limma <- lmFit(y_limma,design)
    ex <- diffSplice(fit_limma, geneid = tab_count$gene_id, exonid=tab_count$feature_id)
    tab_result_limma = topSplice(ex, number = Inf,test='t') # test='t':trans, test='F'or'simes':gene
    
    end_time <- Sys.time()
    all_time = as.numeric(end_time) - as.numeric(start_time)
    tab_time_limma = data.frame(dataType=dataType,
                                term1=term1,term2=term2,
                                m_DTU = 'limma',time = all_time)
    write.table(tab_time_limma, paste(path_save,paste('tab_time',term1,term2,'limma',sep='_'),sep='/'), row.names = TRUE,sep="\t",quote=FALSE)
    
    write.table(tab_result_limma, paste(path_save,paste('tab_DTU',term1,term2,'limma',sep='_'),sep='/'), row.names = TRUE,sep="\t",quote=FALSE)
  }
  gc()
  
  # edgeR-DTU
  if(!file.exists(paste(path_save,paste('tab_DTU',term1,term2,'edgeR',sep='_'),sep='/'))){
    start_time <- Sys.time()
    y_edgeR <- DGEList(counts = as.matrix(tab_count[,-c(1:2)]),
                       samples = tab_label,
                       genes = tab_count[,c(1,2)])
    rownames(y_edgeR) <- tab_count$feature_id
    group <- relevel(factor(tab_label$group),ref=term1)
    design <- model.matrix(~group)
    y_edgeR <- estimateDisp(y_edgeR, design)
    fit_edgeR <- glmQLFit(y_edgeR, design)
    ds <- diffSpliceDGE(fit_edgeR, geneid = tab_count$gene_id, exonid=tab_count$feature_id)
    tab_result_edgeR <- topSpliceDGE(ds, test = "exon", number = Inf)
    
    end_time <- Sys.time()
    all_time = as.numeric(end_time) - as.numeric(start_time)
    tab_time_edgeR = data.frame(dataType=dataType,
                                term1=term1,term2=term2,
                                m_DTU = 'edgeR',time = all_time)
    write.table(tab_time_edgeR, paste(path_save,paste('tab_time',term1,term2,'edgeR',sep='_'),sep='/'), row.names = TRUE,sep="\t",quote=FALSE)
    
    write.table(tab_result_edgeR, paste(path_save,paste('tab_DTU',term1,term2,'edgeR',sep='_'),sep='/'), row.names = TRUE,sep="\t",quote=FALSE)
  }
  gc()
  
  # # satuRn
  if(!file.exists(paste(path_save,paste('tab_DTU',term1,term2,'satuRn',sep='_'),sep='/'))){
    start_time <- Sys.time()
    counts <- tab_count[,-c(1,2)]
    rownames(counts) <- tab_count[,1]
    tab_rowData = data.frame(gene_id = tab_count$gene_id, isoform_id = tab_count$feature_id)
    rownames(tab_rowData) <- tab_count$feature_id
    tab_rowData <- subset(tab_rowData, duplicated(gene_id) | duplicated(gene_id, fromLast = TRUE))
    counts = counts[which(rownames(counts)%in% tab_rowData$isoform_id), ]
    tab_colData = tab_label
    colnames(tab_colData) = c('sample_name','group')
    sumExp <- SummarizedExperiment(
      assays = list(counts = counts),colData = tab_colData,rowData = tab_rowData)
    metadata(sumExp)$formula <- ~0 + as.factor(colData(sumExp)$group)
    group <- relevel(factor(tab_colData$group),ref=term1)
    design <- model.matrix(~0+group)
    colnames(design) <- levels(group)
    sumExp <- fitDTU(object = sumExp, formula=~0+group,verbose = FALSE)
    contr <- makeContrasts(contrasts=paste(term2, "-", term1, sep = ""), levels=design)
    sumExp <- testDTU(object = sumExp, contrasts = contr,diagplot1 = FALSE,diagplot2 = FALSE)
    tab_result_satuRn = rowData(sumExp)[[paste("fitDTUResult",paste(term2, "-", term1, sep = ""),sep='_')]]
    
    end_time <- Sys.time()
    all_time = as.numeric(end_time) - as.numeric(start_time)
    tab_time_satuRn = data.frame(dataType=dataType,
                                 term1=term1,term2=term2,
                                 m_DTU = 'satuRn',time = all_time)
    write.table(tab_time_satuRn, paste(path_save,paste('tab_time',term1,term2,'satuRn',sep='_'),sep='/'), row.names = TRUE,sep="\t",quote=FALSE)
    
    write.table(tab_result_satuRn, paste(path_save,paste('tab_DTU',term1,term2,'satuRn',sep='_'),sep='/'), row.names = TRUE,sep="\t",quote=FALSE)
  }
  gc()
}
DTU()