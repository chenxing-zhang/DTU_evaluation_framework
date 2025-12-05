# run_DTU_longRead_singleCell
library(tidyverse)
library(DRIMSeq)
library(DEXSeq)
library(limma)
library(edgeR)
library(satuRn)
library(parallel)
library(rats) # RATs
library(fishpond) # Swish
library(BiocParallel) # for DTUrtle
library(DTUrtle) # DTUrtle
library(NBSplice) # NBSplice
library(data.table)
library(reticulate)

dataType = 'single_cell'
path_data = '~/data'
path_save = paste('~/result',paste('DTU',dataType,sep='_'),sep='/')

list_vec_term = list()
list_tab_count_trans = list()
list_tab_anno_cell = list()
for(dataSubType in c('PromethION_5cl_rep1', 'PromethION_5cl_rep2', 'PromethION_MSC')){
  
  tab_count_trans_ori <- read.csv(file.path(path_data, dataSubType,"FLTSA_output","transcript_count.csv.gz"), stringsAsFactors=FALSE)
  isoform_FSM_annotation <- read.csv(file.path(path_data, dataSubType,"FLTSA_output","isoform_FSM_annotation.csv") , stringsAsFactors=FALSE)
  tab_anno_cell <- read.csv(file.path(path_data, dataSubType,"cluster_annotation.csv"), stringsAsFactors=FALSE)
  
  # load expression and merge trans information
  rownames(tab_count_trans_ori) = tab_count_trans_ori$transcript_id
  tab_count_trans = tab_count_trans_ori[sort(intersect(tab_count_trans_ori$transcript_id, isoform_FSM_annotation$transcript_id)),
                                        c("transcript_id",sort(intersect(colnames(tab_count_trans_ori),tab_anno_cell$barcode_seq)))]
  rownames(isoform_FSM_annotation) = isoform_FSM_annotation$transcript_id
  tab_count_trans = merge(tab_count_trans,isoform_FSM_annotation[,c('transcript_id','gene_id','FSM_match')],by='transcript_id',all.x = T)
  tab_anno_trans = tab_count_trans[,c("transcript_id","gene_id","FSM_match")]
  vec_cells = colnames(tab_count_trans)[!(colnames(tab_count_trans) %in% c("transcript_id","gene_id","FSM_match"))]
  tab_count_trans = tab_count_trans %>% group_by(gene_id,FSM_match) %>% summarise_at(vec_cells,sum)
  tab_count_trans = tab_count_trans[,c(c("FSM_match","gene_id"),vec_cells)]
  colnames(tab_count_trans) = c('feature_id','gene_id',vec_cells)
  tab_count_trans = as.data.frame(tab_count_trans)
  
  # load sample information
  tab_anno_cell = tab_anno_cell %>% filter(barcode_seq %in% sort(intersect(colnames(tab_count_trans_ori),tab_anno_cell$barcode_seq))) %>% arrange(barcode_seq)
  
  # load terms
  vec_term = sort(unique(tab_anno_cell$groups))
  
  list_vec_term[[dataSubType]] = vec_term
  list_tab_count_trans[[dataSubType]] = tab_count_trans
  list_tab_anno_cell[[dataSubType]] = tab_anno_cell
}
 

DTU_multicore <-function(combin){
  combins = unlist(strsplit(combin,split='-'))
  dataSubType = combins[1]
  i_term = as.integer(combins[2])
  j_term = as.integer(combins[3])
  
  # load terms
  vec_term = list_vec_term[[dataSubType]]
  term1 = vec_term[i_term]
  term2 = vec_term[j_term]
  
  # load expression and sample information
  tab_count_trans = list_tab_count_trans[[dataSubType]]
  tab_anno_cell = list_tab_anno_cell[[dataSubType]] 
  
  # filter samples
  vec_sample1 = sort(tab_anno_cell %>% filter(groups== term1)%>% pull(barcode_seq))
  vec_sample2 = sort(tab_anno_cell %>% filter(groups== term2)%>% pull(barcode_seq))
  list_samples = c(vec_sample1,vec_sample2)
  
  # filter expression
  tab_count = tab_count_trans[,c('feature_id','gene_id',list_samples)]
  tab_count = tab_count[rowSums(tab_count[,-c(1:2)]>=3)>=3,]# filter, important !!
  rownames(tab_count) = tab_count$feature_id
  
  # normalization
  source('~/utils/public/count_2_FPKM_RPKM_TPM.R') 
  if(dataSubType=='PromethION_5cl_rep1' | dataSubType=='PromethION_5cl_rep2'){path_efflen = '~/data/Ref/GENECODE_r29_ENCODE'}
  if(dataSubType=='PromethION_MSC'){path_efflen = '~/data/Ref_mouse/GENECODE_M21_ENCODE'}
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
  
  # SUPPA2, use U-test
  if(!file.exists(paste(path_save,dataSubType,paste('tab_DTU',term1,term2,'SUPPA2',sep='_'),sep='/'))){
    use_condaenv("python", required = TRUE)
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
  
  # RATs
  if(!file.exists(paste(path_save,dataSubType,paste('tab_DTU', term1,term2,'RATs',sep='_'),sep='/'))){
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
  
  ################################################################################
  # DTUrtle
  if(!file.exists(paste(path_save,dataSubType,paste('tab_DTU',term1,term2,'DTUrtle',sep='_'),sep='/'))){
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
    tab_time_DTUrtle = data.frame(dataType=dataType,dataSubType=dataSubType,term1=term1,term2=term2,m_DTU = 'DTUrtle',time = all_time)
    write.table(tab_time_DTUrtle, paste(path_save,dataSubType,paste('tab_time',term1,term2,'DTUrtle',sep='_'),sep='/'), row.names = TRUE,sep="\t",quote=FALSE)
    write.table(tab_result_DTUrtle, paste(path_save,dataSubType,paste('tab_DTU',term1,term2,'DTUrtle',sep='_'),sep='/'), row.names = TRUE,sep="\t",quote=FALSE)
  }
  
  # NBSplice
  if(!file.exists(paste(path_save,dataSubType,paste('tab_DTU',term1,term2,'NBSplice',sep='_'),sep='/'))){
    start_time <- Sys.time()
    geneIso = tab_count[,c('gene_id','feature_id')]
    colnames(geneIso) = c('gene_id','isoform_id')
    colName = 'group'
    myIsoDataSet<-NBSplice::IsoDataSet(tab_count[,c(vec_sample1,vec_sample2)],tab_label, colName, geneIso)
    myIsoDataSet<-buildLowExpIdx(myIsoDataSet, colName, ratioThres = 0.01,countThres = 1)
    myDSResults<-NBTest(myIsoDataSet, colName, test="F")
    tab_result_NBSplice = myDSResults@results
    end_time <- Sys.time()
    all_time = as.numeric(end_time) - as.numeric(start_time)
    tab_time_NBSplice = data.frame(dataType=dataType,dataSubType=dataSubType,term1=term1,term2=term2,m_DTU = 'NBSplice',time = all_time)
    write.table(tab_time_NBSplice, paste(path_save,dataSubType,paste('tab_time',term1,term2,'NBSplice',sep='_'),sep='/'), row.names = TRUE,sep="\t",quote=FALSE)
    write.table(tab_result_NBSplice, paste(path_save,dataSubType,paste('tab_DTU',term1,term2,'NBSplice',sep='_'),sep='/'), row.names = TRUE,sep="\t",quote=FALSE)
  }
  
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
    tab_time_limma = data.frame(dataType=dataType,dataSubType=dataSubType,term1=term1, term2=term2,m_DTU = 'limma',time = all_time)
    write.table(tab_time_limma, paste(path_save,dataSubType,paste('tab_time',term1,term2,'limma',sep='_'),sep='/'), row.names = TRUE,sep="\t",quote=FALSE)
    write.table(tab_result_limma, paste(path_save,dataSubType,paste('tab_DTU',term1,term2,'limma',sep='_'),sep='/'), row.names = TRUE,sep="\t",quote=FALSE)
  }
  
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
  
  # satuRn
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
      sumExp <- testDTU(object = sumExp, contrasts = contr,diagplot1 = 0,diagplot2 = 0)
      tab_result_satuRn = rowData(sumExp)[[paste("fitDTUResult",paste(term2_change, "-", term1_change, sep = ""),sep='_')]]
      end_time <- Sys.time()
      all_time = as.numeric(end_time) - as.numeric(start_time)
      tab_time_satuRn = data.frame(dataType=dataType,dataSubType=dataSubType,term1=term1, term2=term2,m_DTU = 'satuRn',time = all_time)
      write.table(tab_time_satuRn, paste(path_save,dataSubType,paste('tab_time',term1,term2,'satuRn',sep='_'),sep='/'), row.names = TRUE,sep="\t",quote=FALSE)
      write.table(tab_result_satuRn, paste(path_save,dataSubType,paste('tab_DTU',term1,term2,'satuRn',sep='_'),sep='/'), row.names = TRUE,sep="\t",quote=FALSE)
    })
  }
  
  ################################################################################
  # DRIMSeq 
  if(!file.exists(paste(path_save,dataSubType,paste('tab_DTU',term1,term2,'DRIMSeq',sep='_'),sep='/'))){
    start_time <- Sys.time()
    d <- DRIMSeq::dmDSdata(counts = tab_count, samples = tab_label)
    group <- relevel(factor(tab_label$group),ref=term1_change)
    design <- model.matrix(~group)
    d <- DRIMSeq::dmPrecision(d, design = design, verbose = 0)
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
  
  # DEXSeq
  if(!file.exists(paste(path_save,dataSubType,paste('tab_DTU',term1,term2,'DEXSeq',sep='_'),sep='/'))){
    start_time <- Sys.time()
    group <- relevel(factor(tab_label$group),ref=term1_change)
    dxd <- DEXSeq::DEXSeqDataSet(countData = as.matrix(tab_count[,-c(1:2)]),
                                 sampleData=tab_label,
                                 design=~sample+exon+group:exon,
                                 featureID=tab_count$feature_id,
                                 groupID=tab_count$gene_id)
    try_idx = 1
    try_idx = try({dxd <- DEXSeq::estimateSizeFactors(dxd)})
    if('try-error' %in% class(try_idx)){
      dxd <- DEXSeq::DEXSeqDataSet(countData = as.matrix(tab_count[,-c(1:2)])+1,
                                   sampleData=tab_label,
                                   design=~sample+exon+group:exon,
                                   featureID=tab_count$feature_id,
                                   groupID=tab_count$gene_id)
      dxd <- DEXSeq::estimateSizeFactors(dxd)
    }
    dxd <- DEXSeq::estimateDispersions(dxd)
    dxd <- DEXSeq::testForDEU(dxd, reducedModel = ~sample+exon)
    dxd = DEXSeq::estimateExonFoldChanges(dxd,fitExpToVar = "group")
    dxr = DEXSeq::DEXSeqResults(dxd, independentFiltering = FALSE)
    tab_result_DEXSeq = as.data.frame(dxr[,c("featureID", "groupID", term2_change, term1_change ,paste("log2fold",term2_change,term1_change,sep='_'),"pvalue", "padj")])
    end_time <- Sys.time()
    all_time = as.numeric(end_time) - as.numeric(start_time)
    tab_time_DEXSeq = data.frame(dataType=dataType,dataSubType=dataSubType,term1=term1, term2=term2,m_DTU = 'DEXSeq',time = all_time)
    write.table(tab_time_DEXSeq, paste(path_save,dataSubType,paste('tab_time',term1,term2,'DEXSeq',sep='_'),sep='/'), row.names = TRUE,sep="\t",quote=FALSE)
    write.table(tab_result_DEXSeq, paste(path_save,dataSubType,paste('tab_DTU',term1,term2,'DEXSeq',sep='_'),sep='/'), row.names = TRUE,sep="\t",quote=FALSE)
  }
  
}

combinList = c()
for(dataSubType in c('PromethION_5cl_rep1', 'PromethION_5cl_rep2', 'PromethION_MSC')){
  vec_term = list_vec_term[[dataSubType]]
  for(i_term in 1:(length(vec_term)-1)){
    for(j_term in (i_term+1):length(vec_term)){
      combinList = c(combinList, paste(dataSubType, as.character(i_term),as.character(j_term),sep='-'))    
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
mc.cores = 30, 
mc.preschedule = FALSE,  
mc.cleanup = TRUE        
)
gc()



