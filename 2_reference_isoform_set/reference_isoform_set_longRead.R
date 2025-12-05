library(tidyverse)
library(bedtoolsr)
options(bedtools.path = "~/anaconda3/envs/py_fastq/bin")
library(data.table)
library(parallel)

id2idTitle <- function(vec) {
  return(unlist(lapply(vec, function(x) unlist(str_split(x, '\\.'))[1])))
}

path_data <- '~/data'
path_save <- '~/result'
path_support = '~/result/support_file_longRead'


# load regulated prediction model
path_model = '~/result/predict_regulated_isoform'
lda_model  = readRDS(file.path(path_model, "lda_model_allData_2feature.rds"))
glm_model  = readRDS(file.path(path_model, "glm_model_allData_2feature.rds"))

# get support file: vec_trans_tpm, vec_bioTypeTerm1Term2 and regulated isoform set
# For RNAWG_long-read: 'mouse', 'human'
# For single_cell: 'PromethION_5cl_rep1', 'PromethION_5cl_rep2', 'PromethION_MSC'
# For spatial: 'GSE153859_CBS1', 'GSE153859_CBS2', 'GSE153859_MOB'
for(dataType in c('RNAWG_long-read', 'single_cell', 'spatial')){
  if(dataType=='RNAWG_long-read'){vec_dataSubType = c('mouse','human')}
  if(dataType=='single_cell'){vec_dataSubType = c('PromethION_5cl_rep1', 'PromethION_5cl_rep2', 'PromethION_MSC')}
  if(dataType=='spatial'){vec_dataSubType = c('GSE153859_CBS1', 'GSE153859_CBS2', 'GSE153859_MOB')}
  for(dataSubType in vec_dataSubType){
    # Determine species based on dataSubType
    if (dataSubType %in% c('mouse', 'human')) {
      species <- dataSubType
    } else if (dataSubType %in% c("PromethION_5cl_rep1", "PromethION_5cl_rep2")) {
      species <- 'human'
    } else if (dataSubType %in% c("GSE153859_CBS1", "GSE153859_CBS2", "GSE153859_MOB", "PromethION_MSC")) {
      species <- 'mouse'
    } else {
      stop("Unknown dataSubType")
    }
    
    # Load RBP info
    tab_RBPinfo <- read.csv('~/data/RBP_RNA/CLIPdb/POSTAR3_CLIPdb_module_browse_RBP_sub_info.csv', check.names = FALSE)
    colnames(tab_RBPinfo) <- c('RBPname', 'geneID', 'Domain', 'Location', 'Function')
    tab_RBPinfo_splice <- tab_RBPinfo[grepl('Splic', tab_RBPinfo$Function), ]
    vec_RBP_RBPinfo <- tab_RBPinfo_splice[['RBPname']]
    
    # Load RBP-RNA data
    if (species == 'human') {
      tab_CLIPdb <- data.table::fread('~/data/RBP_RNA/CLIPdb/human.txt', sep = '\t')
    } else if (species == 'mouse') {
      tab_CLIPdb <- data.table::fread('~/data/RBP_RNA/CLIPdb/mouse.txt', sep = '\t')
    }
    colnames(tab_CLIPdb) <- c('seqnames', 'start', 'end', 'peak_id', 'strand', 'RBP_name',
                              'experiment_method', 'sample', 'dataID', 'confidence_score')
    
    # Filter RBP data
    if (species == 'human') {
      A.bed <- tab_CLIPdb %>% filter(RBP_name %in% vec_RBP_RBPinfo) %>%
        select(seqnames, start, end, RBP_name, experiment_method, sample) %>%
        filter(experiment_method == 'eCLIP' & sample != 'adrenal gland')
    } else if (species == 'mouse') {
      A.bed <- tab_CLIPdb %>% filter(RBP_name %in% vec_RBP_RBPinfo) %>%
        select(seqnames, start, end, RBP_name, experiment_method, sample)
      A.bed$RBP_name <- str_to_title(A.bed$RBP_name)
    }
    
    # Load annotation data
    if (species == 'human') {
      tab_anno_trans <- read.table('~/data/Ref/GENECODE_r29_ENCODE/tab_anno_trans', sep = '\t', header = TRUE)
      tab_anno_gene <- read.table('~/data/Ref/GENECODE_r29_ENCODE/tab_anno_gene', sep = '\t', header = TRUE)
    } else if (species == 'mouse') {
      if (dataType == 'spatial') {
        tab_anno_trans <- read.table('~/data/Ref_mouse/GENECODE_M24/tab_anno_trans', sep = '\t', header = TRUE)
        tab_anno_gene <- read.table('~/data/Ref_mouse/GENECODE_M24/tab_anno_gene', sep = '\t', header = TRUE)
      } else {
        tab_anno_trans <- read.table('~/data/Ref_mouse/GENECODE_M21_ENCODE/tab_anno_trans', sep = '\t', header = TRUE)
        tab_anno_gene <- read.table('~/data/Ref_mouse/GENECODE_M21_ENCODE/tab_anno_gene', sep = '\t', header = TRUE)
      }
    }
    tab_anno_trans <- tab_anno_trans %>% mutate(transcript_id_title = id2idTitle(transcript_id), .before = transcript_id) %>%
      filter(!grepl('PAR', transcript_id) & !grepl('PAR', gene_id))
    tab_anno_gene <- tab_anno_gene %>% mutate(gene_id_title = id2idTitle(gene_id), .before = gene_id) %>%
      filter(!grepl('PAR', gene_id))
    
    B.bed <- tab_anno_trans %>%
      select(seqnames, start, end, transcript_id, transcript_name, gene_id, gene_name) %>%
      arrange(seqnames, start, end)
    
    # Intersect A and B
    inter <- function(A, B, spe) {
      vec_merge_col <- c('RNA_seqnames', 'RNA_start', 'RNA_end', 'RBP_name',
                         'experiment_method', 'sample','target_seqnames', 'target_start', 'target_end',
                         'target_trans_id', 'target_trans_name','target_gene_id', 'target_gene_name','num_overlap_seq')
      tab_inter_temp <- as.data.table(bedtoolsr::bt.intersect(a = A, b = B, wo = TRUE, f = 1))
      if (dim(tab_inter_temp)[1] == 0) {return(data.frame())}
      colnames(tab_inter_temp) <- vec_merge_col
      tab_inter_temp <- tab_inter_temp %>% select(RBP_name, experiment_method, sample,
                                                  target_trans_id, target_trans_name,
                                                  target_gene_id, target_gene_name) %>% unique()
      tab_inter_temp_num <- tab_inter_temp %>%
        group_by(RBP_name, target_trans_id, target_trans_name, target_gene_id, target_gene_name) %>%
        unique() %>%
        mutate(num_experiment_method = length(unique(experiment_method)),
               num_sample = length(unique(sample)))
      if (spe == 'mouse') {
        tab_inter_temp_num_filer <- tab_inter_temp_num %>%
          filter(num_sample >= 1, num_experiment_method >= 1) %>%
          arrange(RBP_name, target_trans_id)
      } else if (spe == 'human') {
        tab_inter_temp_num_filer <- tab_inter_temp_num %>%
          filter(num_sample >= 2) %>%
          arrange(RBP_name, target_trans_id)
      }
      tab_inter_temp_num_filer_uni <- tab_inter_temp_num_filer %>%
        select(RBP_name, target_trans_id, target_trans_name, target_gene_id, target_gene_name) %>% unique() %>%
        arrange(RBP_name, target_trans_id, target_trans_name, target_gene_id, target_gene_name)
      return(tab_inter_temp_num_filer_uni)
    }
    
    vec_seqnames <- c("chr1", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
                      "chr16", "chr17", "chr18", "chr19", "chr2", "chr20", "chr21",
                      "chr22", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
                      "chr9", "chrM", "chrX", "chrY")
    
    cl <- makeCluster(5, type = "FORK")
    list_tab_inter <- parLapplyLB(cl, vec_seqnames,
                                  function(seqName) inter(A.bed %>% filter(seqnames == seqName),
                                                          B.bed %>% filter(seqnames == seqName),
                                                          species), chunk.size = 1)
    stopCluster(cl)
    tab_inter <- bind_rows(list_tab_inter)
    rm(list_tab_inter);gc()
    vec_RBP_inter <- unique(tab_inter$RBP_name)
    
    
    # Load and process count data based on dataType
    if (dataType == 'RNAWG_long-read') {
      # Bulk data processing
      tab_meta <- read.csv(paste(path_data, dataType, dataSubType, 'metadata.tsv', sep = '/'), check.names = FALSE, sep = '\t')
      tab_meta <- tab_meta %>% filter(`Output type` == 'transcript quantifications')
      tab_meta <- tab_meta %>% arrange(`Biosample type`, `Biosample term name`, `Biosample term id`)
      
      tab_meta_count <- tab_meta %>% group_by(`Biosample type`, `Biosample term name`, `Biosample term id`) %>%
        summarise(count = length(`Biosample term name`))
      tab_meta_filter <- tab_meta %>% filter(`Biosample type` %in% c('cell line', 'tissue') &
                                               `Biosample term name` %in% (tab_meta_count %>% filter(count >= 3) %>% pull(`Biosample term name`)))
      
      tab_count_trans <- read.csv(paste(path_data, dataType, dataSubType, 'tab_count_merge', sep = '/'), sep = '\t', header = TRUE, check.names = FALSE)
      tab_count_trans[is.na(tab_count_trans)] <- 0
      vec_samples <- grep('ENCFF', colnames(tab_count_trans), value = TRUE)
      tab_count_trans <- tab_count_trans[rowSums(tab_count_trans[, vec_samples] > 5) >= 3, ]
      rownames(tab_count_trans) <- tab_count_trans$transcript_id
      vec_columns_trans <- c('transcript_id', 'transcript_name', 'gene_id', 'gene_name')
      
      tab_count_gene <- read.csv(paste(path_data, dataType, dataSubType, 'tab_count_merge_gene', sep = '/'), sep = '\t', header = TRUE, check.names = FALSE)
      tab_count_gene[is.na(tab_count_gene)] <- 0
      tab_count_gene <- tab_count_gene[rowSums(tab_count_gene[, vec_samples] > 5) >= 3, ]
      rownames(tab_count_gene) <- tab_count_gene$gene_id
      vec_columns_gene <- c('gene_id', 'gene_name')
      
    } else if (dataType == 'single_cell') {
      # Single-cell data processing
      tab_meta <- read.csv(file.path(path_data, dataType, dataSubType, "cluster_annotation.csv"), stringsAsFactors = FALSE)
      tab_count_trans <- read.csv(paste(path_data, dataType, dataSubType, 'tab_count_trans', sep = '/'), sep = '\t', header = TRUE, check.names = FALSE)
      tab_count_trans <- tab_count_trans %>%
        mutate(transcript_id_title = id2idTitle(feature_id), .before = 1) %>%
        rename(transcript_id = feature_id)
      tab_count_trans[is.na(tab_count_trans)] <- 0
      vec_samples <- colnames(tab_count_trans)[!grepl('transcript|gene', colnames(tab_count_trans))]
      tab_count_trans <- tab_count_trans[rowSums(tab_count_trans[, vec_samples] > 3) >= 2, ]
      tab_count_trans <- merge(tab_count_trans[, c('transcript_id_title', vec_samples)],
                               tab_anno_trans[, c('transcript_id_title', 'transcript_id', 'transcript_name', 'gene_id', 'gene_name')],
                               by = 'transcript_id_title', all = FALSE) %>% unique()
      rownames(tab_count_trans) <- tab_count_trans$transcript_id
      vec_columns_trans <- c('transcript_id', 'transcript_name', 'gene_id', 'gene_name')
      tab_count_trans <- tab_count_trans[, c(vec_columns_trans, vec_samples)]
      
      tab_count_gene <- read.csv(paste(path_data, dataType, dataSubType, 'tab_count_gene', sep = '/'), sep = '\t', header = TRUE, check.names = FALSE)
      tab_count_gene <- tab_count_gene %>% mutate(gene_id_title = id2idTitle(gene_id), .before = gene_id)
      tab_count_gene[is.na(tab_count_gene)] <- 0
      tab_count_gene <- tab_count_gene[rowSums(tab_count_gene[, vec_samples] > 3) >= 2, ]
      tab_count_gene <- merge(tab_count_gene[, c('gene_id_title', vec_samples)],
                              tab_anno_gene[, c('gene_id_title', 'gene_id', 'gene_name')],
                              by = 'gene_id_title', all = FALSE) %>% unique()
      rownames(tab_count_gene) <- tab_count_gene$gene_id
      vec_columns_gene <- c('gene_id', 'gene_name')
      tab_count_gene <- tab_count_gene[, c(vec_columns_gene, vec_samples)]
      
      tab_meta_filter <- tab_meta %>% filter(barcode_seq %in% vec_samples)
      
    } else if (dataType == 'spatial') {
      # Spatial data processing
      seu <- read_rds(paste(path_data, dataType, paste0(dataSubType, '.rds'), sep = '/'))
      colnames(seu) <- gsub('-| |/|\\(|\\)|\\+', '_', colnames(seu))
      tab_meta <- seu@meta.data
      tab_meta_filter <- tab_meta
      
      tab_count_trans_ori <- as.data.frame(seu@assays[["ISO"]]@counts)
      vec_trans_ori <- rownames(seu@assays[["ISO"]]@counts)
      vec_samples <- colnames(seu@assays[["ISO"]]@counts)
      vec_trans <- unlist(lapply(vec_trans_ori, function(x) unlist(strsplit(x, '\\.\\.'))[2]))
      tab_count_trans_ori['transcript_id'] <- vec_trans
      tab_count_trans <- merge(tab_count_trans_ori,
                               tab_anno_trans[, c('transcript_id', 'transcript_name', 'gene_id', 'gene_name')],
                               by = 'transcript_id', all = FALSE)
      vec_columns_trans <- c('transcript_id', 'transcript_name', 'gene_id', 'gene_name')
      tab_count_trans <- tab_count_trans[, c(vec_columns_trans, vec_samples)]
      rownames(tab_count_trans) <- tab_count_trans$transcript_id
      tab_count_trans <- tab_count_trans[rowSums(tab_count_trans[, vec_samples] > 1) >= 2, ]
      
      tab_count_gene <- tab_count_trans %>% group_by(gene_id, gene_name) %>% summarise_at(vec_samples, sum)
      tab_count_gene <- tab_count_gene[rowSums(tab_count_gene[, vec_samples] > 1) >= 2, ]
      vec_columns_gene <- c('gene_id', 'gene_name')
      tab_count_gene <- as.data.frame(tab_count_gene[, c(vec_columns_gene, vec_samples)])
      rownames(tab_count_gene) <- tab_count_gene$gene_id
    }
    
    # Normalize count data to TPM
    source('~/utils/public/count_2_FPKM_RPKM_TPM.R')
    if (species == 'human') {
      path_efflen <- '~/data/Ref/GENECODE_r29_ENCODE'
    } else if (species == 'mouse') {
      if (dataType == 'spatial') {
        path_efflen <- '~/data/Ref_mouse/GENECODE_M24'
      } else {
        path_efflen <- '~/data/Ref_mouse/GENECODE_M21_ENCODE'
      }
    }
    tab_tpm_trans <- exp_count_2_exp_tpm_path(tab_count_trans[, vec_samples], 'transcript', path_efflen)
    tab_tpm_trans <- merge(tab_count_trans[, vec_columns_trans], tab_tpm_trans, by = 'transcript_id')
    rownames(tab_tpm_trans) <- tab_tpm_trans[, 'transcript_id']
    vec_trans_tpm <- rownames(tab_tpm_trans)
    
    tab_tpm_gene <- exp_count_2_exp_tpm_path(tab_count_gene[, vec_samples], 'gene', path_efflen)
    tab_tpm_gene <- merge(tab_count_gene[, vec_columns_gene], tab_tpm_gene, by = 'gene_id')
    rownames(tab_tpm_gene) <- tab_tpm_gene[, 'gene_id']
    
    # support file: vec_trans_tpm
    saveRDS(vec_trans_tpm, file.path(paste(path_support,'vec_trans_tpm',sep='/'),
                                     paste(paste(dataType, dataSubType,sep='_'),".rds",sep='')))
    
    # Group samples and prepare data 
    vec_bioTypeTerm1Term2 <- c()
    # list_tab_inter_term <- list()
    list_list_vec_regulated_trans_id  = list()
    if (dataType == 'RNAWG_long-read') {
      for (bioType in c('cell line', 'tissue')) {
        vec_term <- sort(tab_meta_filter %>% filter(`Biosample type` == bioType) %>% pull(`Biosample term name`) %>% unique())
        for (i_term in 1:(length(vec_term) - 1)) {
          for (j_term in (i_term + 1):length(vec_term)) {
            term1 <- vec_term[i_term]
            term2 <- vec_term[j_term]
            
            bioTypeTerm1Term2 <- paste(bioType, term1, term2, sep = '_')
            vec_bioTypeTerm1Term2 <- c(vec_bioTypeTerm1Term2, bioTypeTerm1Term2)
            
            vec_sample1 <- tab_meta_filter %>% filter(`Biosample term name` == term1) %>% pull(`File accession`)
            vec_sample2 <- tab_meta_filter %>% filter(`Biosample term name` == term2) %>% pull(`File accession`)
            list_samples <- c(vec_sample1, vec_sample2)
            
            tab_tpm_gene_term <- tab_tpm_gene[rowSums(tab_tpm_gene[, list_samples] > 1) >= 2, c(vec_columns_gene, list_samples)]
            tab_tpm_gene_term['abs_logFC'] <- abs(log2((rowMeans(tab_tpm_gene_term[, vec_sample1]) + 1) / (rowMeans(tab_tpm_gene_term[, vec_sample2]) + 1)))
            tab_tpm_gene_term <- tab_tpm_gene_term %>%
              filter(gene_name %in% vec_RBP_inter) %>% select(c(vec_columns_gene, 'abs_logFC')) %>%
              arrange(desc(abs_logFC))
            colnames(tab_tpm_gene_term) <- c(vec_columns_gene, 'RBP_abs_logFC')
            
            vec_columns_saved <- c('transcript_id', 'transcript_name', 'gene_id', 'gene_name')
            tab_count_trans_term_RBP <- tab_tpm_trans[, c(vec_columns_saved, list_samples)]
            tab_count_trans_term_RBP <- tab_count_trans_term_RBP[rowSums(tab_count_trans_term_RBP[, list_samples] > 1) >= 2, ]
            tab_tpm_trans_mean <- cbind(tab_count_trans_term_RBP[, vec_columns_saved],
                                        data.frame(cl = rowMeans(tab_count_trans_term_RBP[, vec_sample1]),
                                                   ep = rowMeans(tab_count_trans_term_RBP[, vec_sample2])))
            tab_tpm_trans_mean_group <- tab_tpm_trans_mean %>% group_by(gene_id) %>% summarise(sum_cl = sum(cl), sum_ep = sum(ep))
            tab_tpm_trans_mean_merge <- merge(tab_tpm_trans_mean, tab_tpm_trans_mean_group, by = 'gene_id', all.x = TRUE)
            tab_tpm_trans_mean_merge['rate_cl'] <- tab_tpm_trans_mean_merge['cl'] / tab_tpm_trans_mean_merge['sum_cl']
            tab_tpm_trans_mean_merge['rate_ep'] <- tab_tpm_trans_mean_merge['ep'] / tab_tpm_trans_mean_merge['sum_ep']
            tab_tpm_trans_mean_merge[is.na(tab_tpm_trans_mean_merge)] <- 0
            tab_tpm_trans_mean_merge['logFC'] <- log2((tab_tpm_trans_mean_merge['rate_ep'] + 0.01) / (tab_tpm_trans_mean_merge['rate_cl'] + 0.01))
            tab_tpm_trans_mean_merge['abs_logFC'] <- abs(tab_tpm_trans_mean_merge['logFC'])
            tab_tpm_trans_mean_merge['diff_rate'] <- tab_tpm_trans_mean_merge['rate_ep'] - tab_tpm_trans_mean_merge['rate_cl']
            tab_tpm_trans_mean_merge['abs_diff_rate'] <- abs(tab_tpm_trans_mean_merge['diff_rate'])
            tab_tpm_trans_mean_merge <- tab_tpm_trans_mean_merge %>% arrange(desc(abs_diff_rate))
            tab_trans_diffRate <- tab_tpm_trans_mean_merge %>% select(transcript_id, abs_logFC, abs_diff_rate)
            colnames(tab_trans_diffRate) <- c('transcript_id', 'trans_abs_logFC', 'trans_abs_diff_rate')
            
            tab_inter_term <- merge(tab_inter, tab_tpm_gene_term, by.x = 'RBP_name', by.y = 'gene_name', all = FALSE) %>%
              filter(target_trans_id %in% tab_tpm_trans$transcript_id)
            tab_inter_term <- merge(tab_inter_term, tab_trans_diffRate, by.x = 'target_trans_id', by.y = 'transcript_id')
            
            # calculate rate
            num_RBP = length(unique(tab_inter_term$RBP_name))
            tab_inter_term  = tab_inter_term %>% 
              group_by(target_trans_id) %>% 
              mutate(num_RBP_perTrans = n(),rate_RBP_perTrans = num_RBP_perTrans/num_RBP) %>% ungroup() 
            num_trans = length(unique(tab_inter_term$target_trans_id))
            tab_inter_term  = tab_inter_term %>% 
              group_by(RBP_name) %>% 
              mutate(num_trans_perRBP = n(), rate_trans_perRBP = num_trans_perRBP/num_trans) %>% ungroup()
            
            # calculate max1_abs_logFC
            tab_inter_term = tab_inter_term %>% 
              mutate(RBP_max1_abs_logFC = RBP_abs_logFC/max(RBP_abs_logFC),
                     trans_max1_abs_logFC = trans_abs_logFC/max(trans_abs_logFC))
            
            # pred
            list_list_vec_regulated_trans_id[[bioTypeTerm1Term2]] = list()
            for(m_identify in c('lda','glm')){
              if(m_identify=='lda'){pred = predict(lda_model, newdata = tab_inter_term, type = "prob")[, 2]}
              if(m_identify=='glm'){pred = predict(glm_model, newdata = tab_inter_term, type = "prob")[, 2]}
              
              tab_inter_term[paste('predClass',m_identify,sep='_')] = ifelse(pred > 0.5, 'untarget', 'target')
              list_list_vec_regulated_trans_id[[bioTypeTerm1Term2]][[m_identify]] = tab_inter_term %>%
                filter(!!sym(paste('predClass',m_identify,sep='_'))=='target') %>%
                pull(target_trans_id) %>% unique()
            }
          }
        }
      }
    } else if (dataType == 'single_cell') {
      for (bioType in c('single cell')) {
        vec_term <- sort(tab_meta_filter %>% pull(groups) %>% unique())
        for (i_term in 1:(length(vec_term) - 1)) {
          for (j_term in (i_term + 1):length(vec_term)) {
            term1 <- vec_term[i_term]
            term2 <- vec_term[j_term]
            
            bioTypeTerm1Term2 <- paste(bioType, term1, term2, sep = '_')
            vec_bioTypeTerm1Term2 <- c(vec_bioTypeTerm1Term2, bioTypeTerm1Term2)
            
            vec_sample1 <- tab_meta_filter %>% filter(groups == term1) %>% pull(barcode_seq)
            vec_sample2 <- tab_meta_filter %>% filter(groups == term2) %>% pull(barcode_seq)
            list_samples <- c(vec_sample1, vec_sample2)
            
            tab_tpm_gene_term <- tab_tpm_gene[rowSums(tab_tpm_gene[, list_samples] > 1) >= 2, c(vec_columns_gene, list_samples)]
            tab_tpm_gene_term['abs_logFC'] <- abs(log2((rowMeans(tab_tpm_gene_term[, vec_sample1]) + 1) / (rowMeans(tab_tpm_gene_term[, vec_sample2]) + 1)))
            tab_tpm_gene_term <- tab_tpm_gene_term %>%
              filter(gene_name %in% vec_RBP_inter) %>% select(c(vec_columns_gene, 'abs_logFC')) %>%
              arrange(desc(abs_logFC))
            colnames(tab_tpm_gene_term) <- c(vec_columns_gene, 'RBP_abs_logFC')
            
            vec_columns_saved <- c('transcript_id', 'transcript_name', 'gene_id', 'gene_name')
            tab_count_trans_term_RBP <- tab_tpm_trans[, c(vec_columns_saved, list_samples)]
            tab_count_trans_term_RBP <- tab_count_trans_term_RBP[rowSums(tab_count_trans_term_RBP[, list_samples] > 1) >= 2, ]
            tab_tpm_trans_mean <- cbind(tab_count_trans_term_RBP[, vec_columns_saved],
                                        data.frame(cl = rowMeans(tab_count_trans_term_RBP[, vec_sample1]),
                                                   ep = rowMeans(tab_count_trans_term_RBP[, vec_sample2])))
            tab_tpm_trans_mean_group <- tab_tpm_trans_mean %>% group_by(gene_id) %>% summarise(sum_cl = sum(cl), sum_ep = sum(ep))
            tab_tpm_trans_mean_merge <- merge(tab_tpm_trans_mean, tab_tpm_trans_mean_group, by = 'gene_id', all.x = TRUE)
            tab_tpm_trans_mean_merge['rate_cl'] <- tab_tpm_trans_mean_merge['cl'] / tab_tpm_trans_mean_merge['sum_cl']
            tab_tpm_trans_mean_merge['rate_ep'] <- tab_tpm_trans_mean_merge['ep'] / tab_tpm_trans_mean_merge['sum_ep']
            tab_tpm_trans_mean_merge[is.na(tab_tpm_trans_mean_merge)] <- 0
            tab_tpm_trans_mean_merge['logFC'] <- log2((tab_tpm_trans_mean_merge['rate_ep'] + 0.01) / (tab_tpm_trans_mean_merge['rate_cl'] + 0.01))
            tab_tpm_trans_mean_merge['abs_logFC'] <- abs(tab_tpm_trans_mean_merge['logFC'])
            tab_tpm_trans_mean_merge['diff_rate'] <- tab_tpm_trans_mean_merge['rate_ep'] - tab_tpm_trans_mean_merge['rate_cl']
            tab_tpm_trans_mean_merge['abs_diff_rate'] <- abs(tab_tpm_trans_mean_merge['diff_rate'])
            tab_tpm_trans_mean_merge <- tab_tpm_trans_mean_merge %>% arrange(desc(abs_diff_rate))
            tab_trans_diffRate <- tab_tpm_trans_mean_merge %>% select(transcript_id, abs_logFC, abs_diff_rate)
            colnames(tab_trans_diffRate) <- c('transcript_id', 'trans_abs_logFC', 'trans_abs_diff_rate')
            
            tab_inter_term <- merge(tab_inter, tab_tpm_gene_term, by.x = 'RBP_name', by.y = 'gene_name', all = FALSE) %>%
              filter(target_trans_id %in% tab_tpm_trans$transcript_id)
            tab_inter_term <- merge(tab_inter_term, tab_trans_diffRate, by.x = 'target_trans_id', by.y = 'transcript_id')
            
            # calculate rate
            num_RBP = length(unique(tab_inter_term$RBP_name))
            tab_inter_term  = tab_inter_term %>% 
              group_by(target_trans_id) %>% 
              mutate(num_RBP_perTrans = n(),rate_RBP_perTrans = num_RBP_perTrans/num_RBP) %>% ungroup() 
            num_trans = length(unique(tab_inter_term$target_trans_id))
            tab_inter_term  = tab_inter_term %>% 
              group_by(RBP_name) %>% 
              mutate(num_trans_perRBP = n(), rate_trans_perRBP = num_trans_perRBP/num_trans) %>% ungroup()
            
            # calculate max1_abs_logFC
            tab_inter_term = tab_inter_term %>% 
              mutate(RBP_max1_abs_logFC = RBP_abs_logFC/max(RBP_abs_logFC),
                     trans_max1_abs_logFC = trans_abs_logFC/max(trans_abs_logFC))
            
            # pred
            list_list_vec_regulated_trans_id[[bioTypeTerm1Term2]] = list()
            for(m_identify in c('lda','glm')){
              if(m_identify=='lda'){pred = predict(lda_model, newdata = tab_inter_term, type = "prob")[, 2]}
              if(m_identify=='glm'){pred = predict(glm_model, newdata = tab_inter_term, type = "prob")[, 2]}
              
              tab_inter_term[paste('predClass',m_identify,sep='_')] = ifelse(pred > 0.5, 'untarget', 'target')
              list_list_vec_regulated_trans_id[[bioTypeTerm1Term2]][[m_identify]] = tab_inter_term %>%
                filter(!!sym(paste('predClass',m_identify,sep='_'))=='target') %>%
                pull(target_trans_id) %>% unique()
            }
            
          }
        }
      }
    } else if (dataType == 'spatial') {
      for (bioType in c('spatial')) {
        vec_term <- sort(as.character(unique(tab_meta_filter$ClusterName)))
        for (i_term in 1:(length(vec_term) - 1)) {
          for (j_term in (i_term + 1):length(vec_term)) {
            term1 <- vec_term[i_term]
            term2 <- vec_term[j_term]
            if (dataSubType == 'GSE153859_CBS2' && term1 == 'Hippocampus area' && term2 == 'Isocortex-1') next
            if (dataSubType == 'GSE153859_CBS2' && term1 == 'Isocortex-1' && term2 == 'Midbrain') next
            
            term1_change <- gsub('-| |/|\\(|\\)|\\+', '_', term1)
            term2_change <- gsub('-| |/|\\(|\\)|\\+', '_', term2)
            bioTypeTerm1Term2 <- paste(bioType, term1_change, term2_change, sep = '-')
            vec_bioTypeTerm1Term2 <- c(vec_bioTypeTerm1Term2, bioTypeTerm1Term2)
            
            vec_sample1 <- sort(rownames(tab_meta_filter %>% filter(ClusterName == term1)))
            vec_sample2 <- sort(rownames(tab_meta_filter %>% filter(ClusterName == term2)))
            list_samples <- c(vec_sample1, vec_sample2)
            
            tab_tpm_gene_term <- tab_tpm_gene[rowSums(tab_tpm_gene[, list_samples] > 1) >= 2, c(vec_columns_gene, list_samples)]
            tab_tpm_gene_term['abs_logFC'] <- abs(log2((rowMeans(tab_tpm_gene_term[, vec_sample1]) + 1) / (rowMeans(tab_tpm_gene_term[, vec_sample2]) + 1)))
            tab_tpm_gene_term <- tab_tpm_gene_term %>%
              filter(gene_name %in% vec_RBP_inter) %>% select(c(vec_columns_gene, 'abs_logFC')) %>%
              arrange(desc(abs_logFC))
            colnames(tab_tpm_gene_term) <- c(vec_columns_gene, 'RBP_abs_logFC')
            
            vec_columns_saved <- c('transcript_id', 'transcript_name', 'gene_id', 'gene_name')
            tab_count_trans_term_RBP <- tab_tpm_trans[, c(vec_columns_saved, list_samples)]
            tab_count_trans_term_RBP <- tab_count_trans_term_RBP[rowSums(tab_count_trans_term_RBP[, list_samples] > 1) >= 2, ]
            tab_tpm_trans_mean <- cbind(tab_count_trans_term_RBP[, vec_columns_saved],
                                        data.frame(cl = rowMeans(tab_count_trans_term_RBP[, vec_sample1]),
                                                   ep = rowMeans(tab_count_trans_term_RBP[, vec_sample2])))
            tab_tpm_trans_mean_group <- tab_tpm_trans_mean %>% group_by(gene_id) %>% summarise(sum_cl = sum(cl), sum_ep = sum(ep))
            tab_tpm_trans_mean_merge <- merge(tab_tpm_trans_mean, tab_tpm_trans_mean_group, by = 'gene_id', all.x = TRUE)
            tab_tpm_trans_mean_merge['rate_cl'] <- tab_tpm_trans_mean_merge['cl'] / tab_tpm_trans_mean_merge['sum_cl']
            tab_tpm_trans_mean_merge['rate_ep'] <- tab_tpm_trans_mean_merge['ep'] / tab_tpm_trans_mean_merge['sum_ep']
            tab_tpm_trans_mean_merge[is.na(tab_tpm_trans_mean_merge)] <- 0
            tab_tpm_trans_mean_merge['logFC'] <- log2((tab_tpm_trans_mean_merge['rate_ep'] + 0.01) / (tab_tpm_trans_mean_merge['rate_cl'] + 0.01))
            tab_tpm_trans_mean_merge['abs_logFC'] <- abs(tab_tpm_trans_mean_merge['logFC'])
            tab_tpm_trans_mean_merge['diff_rate'] <- tab_tpm_trans_mean_merge['rate_ep'] - tab_tpm_trans_mean_merge['rate_cl']
            tab_tpm_trans_mean_merge['abs_diff_rate'] <- abs(tab_tpm_trans_mean_merge['diff_rate'])
            tab_tpm_trans_mean_merge <- tab_tpm_trans_mean_merge %>% arrange(desc(abs_diff_rate))
            tab_trans_diffRate <- tab_tpm_trans_mean_merge %>% select(transcript_id, abs_logFC, abs_diff_rate)
            colnames(tab_trans_diffRate) <- c('transcript_id', 'trans_abs_logFC', 'trans_abs_diff_rate')
            
            tab_inter_term <- merge(tab_inter, tab_tpm_gene_term, by.x = 'RBP_name', by.y = 'gene_name', all = FALSE) %>%
              filter(target_trans_id %in% tab_tpm_trans$transcript_id)
            tab_inter_term <- merge(tab_inter_term, tab_trans_diffRate, by.x = 'target_trans_id', by.y = 'transcript_id')
            
            # calculate rate
            num_RBP = length(unique(tab_inter_term$RBP_name))
            tab_inter_term  = tab_inter_term %>% 
              group_by(target_trans_id) %>% 
              mutate(num_RBP_perTrans = n(),rate_RBP_perTrans = num_RBP_perTrans/num_RBP) %>% ungroup() 
            num_trans = length(unique(tab_inter_term$target_trans_id))
            tab_inter_term  = tab_inter_term %>% 
              group_by(RBP_name) %>% 
              mutate(num_trans_perRBP = n(), rate_trans_perRBP = num_trans_perRBP/num_trans) %>% ungroup()
            
            # calculate max1_abs_logFC
            tab_inter_term = tab_inter_term %>% 
              mutate(RBP_max1_abs_logFC = RBP_abs_logFC/max(RBP_abs_logFC),
                     trans_max1_abs_logFC = trans_abs_logFC/max(trans_abs_logFC))
            
            # pred
            list_list_vec_regulated_trans_id[[bioTypeTerm1Term2]] = list()
            for(m_identify in c('lda','glm')){
              if(m_identify=='lda'){pred = predict(lda_model, newdata = tab_inter_term, type = "prob")[, 2]}
              if(m_identify=='glm'){pred = predict(glm_model, newdata = tab_inter_term, type = "prob")[, 2]}
              
              tab_inter_term[paste('predClass',m_identify,sep='_')] = ifelse(pred > 0.5, 'untarget', 'target')
              list_list_vec_regulated_trans_id[[bioTypeTerm1Term2]][[m_identify]] = tab_inter_term %>%
                filter(!!sym(paste('predClass',m_identify,sep='_'))=='target') %>%
                pull(target_trans_id) %>% unique()
            }
            
          }
        }
      }
    }
    
    # support file: vec_bioTypeTerm1Term2
    saveRDS(vec_bioTypeTerm1Term2, file.path(paste(path_support,'vec_bioTypeTerm1Term2',sep='/'),
                                             paste(paste(dataType, dataSubType,sep='_'),".rds",sep='')))
    # support file: regulated isoform set
    saveRDS(list_list_vec_regulated_trans_id, file.path(paste(path_support,'list_list_vec_regulated_trans_id',sep='/'),
                                                     paste(paste(dataType, dataSubType,sep='_'),".rds",sep='')))
  
    # support file: common isoform set
    list_vec_trans_common_DTU = mclapply(vec_bioTypeTerm1Term2, function(bioTypeTerm1Term2) {
      suppressPackageStartupMessages({library(tidyverse)})
      get_common_trans(combin)
    }, 
    mc.cores = 15, 
    mc.preschedule = FALSE,  
    mc.cleanup = TRUE       
    )
    gc()
    names(list_vec_trans_common_DTU) = vec_dataSubTypeTermRBP
    saveRDS(list_vec_trans_common_DTU, file.path(paste(path_support,'list_vec_common_trans',sep='/'), 
                                                 paste(dataType,dataSubType,".rds",sep='')))
    # support file: overlapped isoform set
    list_list_vec_overlapped_trans = list()
    for(combin in vec_dataSubTypeTermRBP){
      list_list_vec_overlapped_trans[[combin]] = list()
      list_list_vec_overlapped_trans[[combin]][['lda']] = intersect(list_vec_trans_common_DTU[[combin]],list_list_vec_regulated_trans_id[[combin]][['lda']]) 
      list_list_vec_overlapped_trans[[combin]][['glm']] = intersect(list_vec_trans_common_DTU[[combin]],list_list_vec_regulated_trans_id[[combin]][['glm']]) 
    }
    saveRDS(list_list_vec_overlapped_trans, file.path(paste(path_support,'list_list_vec_overlapped_trans',sep='/'), 
                                                 paste(dataType,dataSubType,".rds",sep='')))
    
  }
}

