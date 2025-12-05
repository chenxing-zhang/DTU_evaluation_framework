# run_DTU_shortRead_bulk
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

dataType= 'ENCODE_NGS' 
path_data = '~/data'
path_save = paste('~/result',paste('DTU',dataType,sep='_'),sep='/')


list_tab_report_ep = list()
list_tab_meta_ep_trans = list()
list_tab_count_ep_trans = list()
list_tab_report_cl = list()
list_tab_meta_cl_trans = list()
list_tab_count_cl_trans = list()
for(dataSubType in c('CRISPR','shRNA')){
  # experiment
  tab_report_ep = read.csv(paste(path_data,dataSubType,'experiment_report.tsv',sep='/'),check.names = F,sep='\t',skip=1)
  tab_meta_ep = read.csv(paste(path_data,dataSubType,'metadata.tsv',sep='/'),check.names = F,sep='\t')
  tab_meta_ep_trans = tab_meta_ep %>% filter(`Genome annotation`=='V29') %>% 
    filter(`Output type` == 'transcript quantifications') %>% filter(Size>30000000) # RSEM quantification
  tab_meta_ep_trans = tab_meta_ep_trans %>% arrange(`Biosample type`,`Biosample term name`,`Biosample term id`,`Experiment target`,`Experiment accession`,`Biological replicate(s)`,`Technical replicate(s)`)
  tab_count_ep_trans = read.csv(paste(path_data,dataSubType,'tab_count_merge_trans',sep='/'),sep='\t',header=T,check.names = F)
  tab_count_ep_trans[is.na(tab_count_ep_trans)]=0
  tab_count_ep_trans = tab_count_ep_trans %>% filter(!grepl('PAR_Y',transcript_id))
  
  list_tab_report_ep[[dataSubType]] = tab_report_ep
  list_tab_meta_ep_trans[[dataSubType]] = tab_meta_ep_trans
  list_tab_count_ep_trans[[dataSubType]] = tab_count_ep_trans
  
  # control
  tab_report_cl = read.csv(paste(path_data,paste('cl',dataSubType,sep='_'),'experiment_report.tsv',sep='/'),check.names = F,sep='\t',skip=1)
  tab_meta_cl = read.csv(paste(path_data,paste('cl',dataSubType,sep='_'),'metadata.tsv',sep='/'),check.names = F,sep='\t')
  tab_meta_cl_trans = tab_meta_cl %>% filter(`Genome annotation`=='V29') %>% 
    filter(`Output type` == 'transcript quantifications')%>% filter(Size>30000000) # RSEM quantification
  tab_meta_cl_trans = tab_meta_cl_trans %>% arrange(`Biosample type`,`Biosample term name`,`Biosample term id`,`Experiment target`,`Experiment accession`,`Biological replicate(s)`,`Technical replicate(s)`)
  tab_count_cl_trans = read.csv(paste(path_data,paste('cl',dataSubType,sep='_'),'tab_count_merge_trans',sep='/'),sep='\t',header=T,check.names = F)
  tab_count_cl_trans[is.na(tab_count_cl_trans)]=0
  tab_count_cl_trans = tab_count_cl_trans %>% filter(!grepl('PAR_Y',transcript_id))
  
  list_tab_report_cl[[dataSubType]] = tab_report_cl
  list_tab_meta_cl_trans[[dataSubType]] = tab_meta_cl_trans
  list_tab_count_cl_trans[[dataSubType]] = tab_count_cl_trans
}

# combin = combinList[1]
DTU_multicore <-function(combin){
  combins = unlist(strsplit(combin,split='_'))
  dataSubType = combins[1]
  term1 = combins[2]
  term2 = combins[3]
  
  tab_report_ep = list_tab_report_ep[[dataSubType]]
  tab_meta_ep_trans = list_tab_meta_ep_trans[[dataSubType]] 
  tab_count_ep_trans = list_tab_count_ep_trans[[dataSubType]]
  
  tab_report_cl = list_tab_report_cl[[dataSubType]]
  tab_meta_cl_trans = list_tab_meta_cl_trans[[dataSubType]] 
  tab_count_cl_trans = list_tab_count_cl_trans[[dataSubType]]
  gc()
  
  # filter samples
  tab_meta_ep_trans_term_RBP = tab_meta_ep_trans %>% 
    filter(`Biosample term name`==term1 & `Experiment target`== term2)
  ## ep
  vec_s_ep_trans_term_RBP_ori = tab_meta_ep_trans_term_RBP %>% pull(`File accession`)
  ## cl
  access_ep_trans = tab_meta_ep_trans_term_RBP$`Experiment accession`[1]
  series_ep_trans = str_extract(tab_report_ep %>% filter(Accession == access_ep_trans) %>% pull(`Related series`),
                                "/gene-silencing-series/[^,]+")
  access_cl_trans = tab_report_cl %>% filter(grepl(series_ep_trans,`Related series`)) %>% pull(Accession)
  vec_s_cl_trans_term_RBP_ori = tab_meta_cl_trans %>% filter(`Experiment accession` == access_cl_trans) %>% pull(`File accession`)
  
  vec_sample1 = vec_s_cl_trans_term_RBP_ori
  vec_sample2 = vec_s_ep_trans_term_RBP_ori
  list_samples = c(vec_sample1,vec_sample2)
  
  # filter expression
  tab_count_ori = merge(tab_count_cl_trans[,c('transcript_id','gene_id',vec_sample1)],
                    tab_count_ep_trans[,c('transcript_id','gene_id',vec_sample2)],by=c('transcript_id','gene_id'),all=T)
  tab_count_ori[is.na(tab_count_ori)]=0
  colnames(tab_count_ori) = c('feature_id','gene_id',list_samples)
  rownames(tab_count_ori) = tab_count_ori$feature_id
  tab_count = tab_count[rowSums(tab_count[,-c(1:2)]>5)>=2,] # filter, important !! 
  
  # normalization
  source('~/utils/public/count_2_FPKM_RPKM_TPM.R') 
  path_efflen = '~/data/Ref/GENECODE_r29_ENCODE'
  tab_tpm_trans = exp_count_2_exp_tpm_path(tab_count[,list_samples],'transcript',path_efflen)
  tab_tpm_trans = cbind(tab_count[,c('feature_id','gene_id')],tab_tpm_trans[,list_samples])
  
  # create label
  term1_change = gsub('-| ','_',term1)
  term2_change = gsub('-| ','_',term2)
  tab_label = data.frame(sample_id = list_samples, group = relevel(factor(c(rep(term1_change,length(vec_sample1)),
                                                                            rep(term2_change,length(vec_sample2)))),ref=term1_change))
  rownames(tab_label) = tab_label$sample_id

  ################################################################################
  # different transcriptome
  suppressPackageStartupMessages({
    library(reticulate)
  })
  if (!reticulate::py_available(initialize = FALSE)) {
    reticulate::use_condaenv("python", required = TRUE)
  }
  
  # SPIT
  if(!file.exists(paste(path_save,dataSubType,paste('tab_DTU',term1,term2,'SPIT',sep='_'),sep='/'))){
    source_python("~/utils/DTU/SPIT/dtu_detection_zcx.py")
    start_time <- Sys.time()
    list_result_SPIT = compute_double_stats_zcx(tab_count[,c('gene_id',vec_sample1,vec_sample2)],
                                                tab_count[,c('feature_id','gene_id')],
                                                vec_sample1, vec_sample2)
    tab_result_SPIT = as.data.frame(list_result_SPIT[1:4])
    colnames(tab_result_SPIT) = c('transcript_id','gene_id','likelihood','pvalue')
    end_time <- Sys.time()
    all_time = as.numeric(end_time) - as.numeric(start_time)
    tab_time_SPIT = data.frame(dataType=dataType,dataSubType=dataSubType,term1=term1,term2=term2,m_DTU = 'SPIT',time = all_time)
    write.table(tab_time_SPIT, paste(path_save,dataSubType,paste('tab_time',term1,term2,'SPIT',sep='_'),sep='/'), row.names = TRUE,sep="\t",quote=FALSE)
    write.table(tab_result_SPIT, paste(path_save,dataSubType,paste('tab_DTU',term1,term2,'SPIT',sep='_'),sep='/'), row.names = TRUE,sep="\t",quote=FALSE)
  }
  gc()
  
  # SUPPA2
  if(!file.exists(paste(path_save,dataSubType,paste('tab_DTU',term1,term2,'SUPPA2',sep='_'),sep='/'))){
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
    tab_time_SUPPA2 = data.frame(dataType=dataType,dataSubType=dataSubType,term1=term1,term2=term2,m_DTU = 'SUPPA2',time = all_time)
    write.table(tab_time_SUPPA2, paste(path_save,dataSubType,paste('tab_time',term1,term2,'SUPPA2',sep='_'),sep='/'), row.names = TRUE,sep="\t",quote=FALSE)
    write.table(tab_result_SUPPA2, paste(path_save,dataSubType,paste('tab_DTU',term1,term2,'SUPPA2',sep='_'),sep='/'), row.names = TRUE,sep="\t",quote=FALSE)
  }
  gc()
  
  # RATs
  if(!file.exists(paste(path_save,dataSubType,paste('tab_DTU',term1,term2,'RATs',sep='_'),sep='/'))){
    start_time <- Sys.time()
    mydtu <- rats::call_DTU(annot= tab_tpm_trans[,c('feature_id','gene_id')],
                            TARGET_COL = "feature_id", PARENT_COL = "gene_id",
                            count_data_A= setDT(tab_tpm_trans[,c('feature_id',vec_sample1)]),
                            count_data_B= setDT(tab_tpm_trans[,c('feature_id',vec_sample2)]),
                            abund_thresh = 0, dprop_thresh=0, verbose= FALSE)
    tab_result_RATs = mydtu[["Transcripts"]]
    end_time <- Sys.time()
    all_time = as.numeric(end_time) - as.numeric(start_time)
    tab_time_RATs = data.frame(dataType=dataType,dataSubType=dataSubType,term1=term1,term2=term2,m_DTU = 'RATs',time = all_time)
    write.table(tab_time_RATs, paste(path_save,dataSubType,paste('tab_time',term1,term2,'RATs',sep='_'),sep='/'), row.names = TRUE,sep="\t",quote=FALSE)
    write.table(tab_result_RATs, paste(path_save,dataSubType,paste('tab_DTU',term1,term2,'RATs',sep='_'),sep='/'), row.names = TRUE,sep="\t",quote=FALSE)
  }
  gc()
  
  ################################################################################
  # DTUrtle
  if(!file.exists(paste(path_save,dataSubType,paste('tab_DTU',term1,term2,'DTUrtle',sep='_'),sep='/'))){
    start_time <- Sys.time()
    dturtle <- DTUrtle::run_drimseq(counts = as.matrix(tab_count[,c(vec_sample1,vec_sample2)]),
                                    tx2gene = tab_count[,c('feature_id','gene_id')],
                                    pd=tab_label, id_col = "sample_id",cond_col = "group",
                                    filtering_strategy = "own",
                                    min_samps_gene_expr=0, min_samps_feature_expr=0,
                                    min_samps_feature_prop=0,min_gene_expr=0,
                                    min_feature_expr=0, min_feature_prop=0,run_gene_twice=F)
    tab_result_DTUrtle = dturtle[["drim"]]@results_feature
    end_time <- Sys.time()
    all_time = as.numeric(end_time) - as.numeric(start_time)
    tab_time_DTUrtle = data.frame(dataType=dataType,dataSubType=dataSubType,term1=term1,term2=term2,m_DTU = 'DTUrtle',time = all_time)
    write.table(tab_time_DTUrtle, paste(path_save,dataSubType,paste('tab_time',term1,term2,'DTUrtle',sep='_'),sep='/'), row.names = TRUE,sep="\t",quote=FALSE)
    write.table(tab_result_DTUrtle, paste(path_save,dataSubType,paste('tab_DTU',term1,term2,'DTUrtle',sep='_'),sep='/'), row.names = TRUE,sep="\t",quote=FALSE)
  }
  gc()
  
  # NBSplice
  if(!file.exists(paste(path_save,dataSubType,paste('tab_DTU',term1,term2,'NBSplice',sep='_'),sep='/'))){
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
    tab_time_NBSplice = data.frame(dataType=dataType,dataSubType=dataSubType,term1=term1,term2=term2,m_DTU = 'NBSplice',time = all_time)
    write.table(tab_time_NBSplice, paste(path_save,dataSubType,paste('tab_time',term1,term2,'NBSplice',sep='_'),sep='/'), row.names = TRUE,sep="\t",quote=FALSE)
    write.table(tab_result_NBSplice, paste(path_save,dataSubType,paste('tab_DTU',term1,term2,'NBSplice',sep='_'),sep='/'), row.names = TRUE,sep="\t",quote=FALSE)
  }
  gc()
  
  # limma-DTU
  if(!file.exists(paste(path_save,dataSubType,paste('tab_DTU',term1,term2,'limma',sep='_'),sep='/'))){
    start_time <- Sys.time()
    y_limma <- DGEList(counts = as.matrix(tab_count[,-c(1:2)]),
                       samples = tab_label,
                       genes = tab_count[,c(1,2)])
    rownames(y_limma) <- tab_count$feature_id
    y_limma <- calcNormFactors(y_limma)
    group <- relevel(factor(tab_label$group),ref=term1_change)
    design <- model.matrix(~group)
    y_limma <- voom(y_limma, design, plot=F)
    fit_limma <- lmFit(y_limma,design)
    ex <- diffSplice(fit_limma, geneid = tab_count$gene_id, exonid=tab_count$feature_id)
    tab_result_limma = topSplice(ex, number = Inf,test='t') # test='t':trans, test='F'or'simes':gene
    end_time <- Sys.time()
    all_time = as.numeric(end_time) - as.numeric(start_time)
    tab_time_limma = data.frame(dataType=dataType,dataSubType=dataSubType,term1=term1,term2=term2,m_DTU = 'limma',time = all_time)
    write.table(tab_time_limma, paste(path_save,dataSubType,paste('tab_time',term1,term2,'limma',sep='_'),sep='/'), row.names = TRUE,sep="\t",quote=FALSE)
    write.table(tab_result_limma, paste(path_save,dataSubType,paste('tab_DTU',term1,term2,'limma',sep='_'),sep='/'), row.names = TRUE,sep="\t",quote=FALSE)
  }
  gc()
  
  # edgeR-DTU
  if(!file.exists(paste(path_save,dataSubType,paste('tab_DTU',term1,term2,'edgeR',sep='_'),sep='/'))){
    start_time <- Sys.time()
    y_edgeR <- DGEList(counts = as.matrix(tab_count[,-c(1:2)]),
                       samples = tab_label,
                       genes = tab_count[,c(1,2)])
    rownames(y_edgeR) <- tab_count$feature_id
    group <- relevel(factor(tab_label$group),ref=term1_change)
    design <- model.matrix(~group)
    y_edgeR <- estimateDisp(y_edgeR, design)
    fit_edgeR <- glmQLFit(y_edgeR, design)
    ds <- diffSpliceDGE(fit_edgeR, geneid = tab_count$gene_id, exonid=tab_count$feature_id)
    tab_result_edgeR <- topSpliceDGE(ds, test = "exon", number = Inf)
    end_time <- Sys.time()
    all_time = as.numeric(end_time) - as.numeric(start_time)
    tab_time_edgeR = data.frame(dataType=dataType,dataSubType=dataSubType,term1=term1,term2=term2,m_DTU = 'edgeR',time = all_time)
    write.table(tab_time_edgeR, paste(path_save,dataSubType,paste('tab_time',term1,term2,'edgeR',sep='_'),sep='/'), row.names = TRUE,sep="\t",quote=FALSE)
    write.table(tab_result_edgeR, paste(path_save,dataSubType,paste('tab_DTU',term1,term2,'edgeR',sep='_'),sep='/'), row.names = TRUE,sep="\t",quote=FALSE)
  }
  gc()
  
  # # satuRn
  if(!file.exists(paste(path_save,dataSubType,paste('tab_DTU',term1,term2,'satuRn',sep='_'),sep='/'))){
    options(device = function(...) pdf(file = NULL))
    invisible({
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
      group <- relevel(factor(tab_colData$group),ref=term1_change)
      design <- model.matrix(~0+group)
      colnames(design) <- levels(group)
      sumExp <- fitDTU(object = sumExp, formula=~0+group,verbose = FALSE)
      contr <- makeContrasts(contrasts=paste(term2_change, "-", term1_change, sep = ""), levels=design)
      sumExp <- testDTU(object = sumExp, contrasts = contr,diagplot1 = FALSE,diagplot2 = FALSE)
      tab_result_satuRn = rowData(sumExp)[[paste("fitDTUResult",paste(term2_change, "-", term1_change, sep = ""),sep='_')]]
      end_time <- Sys.time()
      all_time = as.numeric(end_time) - as.numeric(start_time)
      tab_time_satuRn = data.frame(dataType=dataType,dataSubType=dataSubType,term1=term1,term2=term2,m_DTU = 'satuRn',time = all_time)
      write.table(tab_time_satuRn, paste(path_save,dataSubType,paste('tab_time',term1,term2,'satuRn',sep='_'),sep='/'), row.names = TRUE,sep="\t",quote=FALSE)
      write.table(tab_result_satuRn, paste(path_save,dataSubType,paste('tab_DTU',term1,term2,'satuRn',sep='_'),sep='/'), row.names = TRUE,sep="\t",quote=FALSE)
    })
  }
  gc()
  
  # ################################################################################
  # DRIMSeq
  if(!file.exists(paste(path_save,dataSubType,paste('tab_DTU',term1,term2,'DRIMSeq',sep='_'),sep='/'))){
    # try({
    start_time <- Sys.time()
    d <- DRIMSeq::dmDSdata(counts = tab_count, samples = tab_label)
    group <- relevel(factor(tab_label$group),ref=term1_change)
    design <- model.matrix(~group)
    d <- DRIMSeq::dmPrecision(d, design = design,verbose=0)
    d <- DRIMSeq::dmFit(d, design = design, verbose = 0)
    d <- DRIMSeq::dmTest(d, coef = paste0("group",term2_change), verbose = 0)
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
    tab_time_DRIMSeq = data.frame(dataType=dataType,dataSubType=dataSubType,term1=term1,term2=term2,m_DTU = 'DRIMSeq',time = all_time)
    write.table(tab_time_DRIMSeq, paste(path_save,dataSubType,paste('tab_time',term1,term2,'DRIMSeq',sep='_'),sep='/'), row.names = TRUE,sep="\t",quote=FALSE)
    write.table(tab_result_DRIMSeq, paste(path_save,dataSubType,paste('tab_DTU',term1,term2,'DRIMSeq',sep='_'),sep='/'), row.names = TRUE,sep="\t",quote=FALSE)
  }
  gc()
  
  # DEXSeq
  if(!file.exists(paste(path_save,dataSubType,paste('tab_DTU',term1,term2,'DEXSeq',sep='_'),sep='/'))){
    # try({
    start_time <- Sys.time()
    group <- relevel(factor(tab_label$group),ref=term1_change)
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
    tab_result_DEXSeq = as.data.frame(dxr[,c("featureID", "groupID", term2_change, term1_change ,paste("log2fold",term2_change,term1_change,sep='_'),"pvalue", "padj")])
    end_time <- Sys.time()
    all_time = as.numeric(end_time) - as.numeric(start_time)
    tab_time_DEXSeq = data.frame(dataType=dataType,term1=term1,term2=term2,dataSubType=dataSubType,m_DTU = 'DEXSeq',time = all_time)
    write.table(tab_time_DEXSeq, paste(path_save,dataSubType,paste('tab_time',term1,term2,'DEXSeq',sep='_'),sep='/'), row.names = TRUE,sep="\t",quote=FALSE)
    write.table(tab_result_DEXSeq, paste(path_save,dataSubType,paste('tab_DTU',term1,term2,'DEXSeq',sep='_'),sep='/'), row.names = TRUE,sep="\t",quote=FALSE)
  }
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
    tab_meta_ep_trans_term = list_tab_meta_ep_trans[[dataSubType]] %>% filter(`Biosample term name`==term1)
    vec_RBP_ep_trans = unlist(lapply(tab_meta_ep_trans_term[['Experiment target']], function(x) unlist(strsplit(x,'-'))[1]))
    
    # load ENCORE meta data
    tab_meta_ENCORE = readxl::read_excel(paste0('~/data/RBP_RNA/ENCORE/',term1,'/metadata_mergedbed.xlsx'))
    vec_RBP_ENCORE = unlist(lapply(tab_meta_ENCORE[['Experiment target']], function(x) unlist(strsplit(x,'-'))[1]))
    
    # intersect among vec_RBP_RBPinfo, vec_RBP_ENCORE and vec_RBP_ep
    vec_RBP = intersect(intersect(vec_RBP_RBPinfo,vec_RBP_ENCORE),vec_RBP_ep_trans)
    vec_term = paste(vec_RBP,'-human',sep='')
    
    for(term2 in vec_term){
      combinList = c(combinList, paste(dataSubType, term1,term2,sep='_')) 
    }
  }
}


mclapply(combinList, function(combin) {
  suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
    library(DRIMSeq)
    library(DEXSeq)
    library(limma)
    library(edgeR)
    library(satuRn)
    library(rats)
    library(BiocParallel) # for DTUrtle
    library(DTUrtle) # DTUrtle
    library(NBSplice) # NBSplice
    library(reticulate) 
  })
  DTU_multicore(combin)
}, 
mc.cores = 30, # 32-2
mc.preschedule = FALSE,  
mc.cleanup = TRUE       
)
gc()

