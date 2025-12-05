library(tidyverse)
library(bedtoolsr)
options(bedtools.path = "~/anaconda3/envs/py_fastq/bin")
library(data.table)
library(fgsea)
library(caret)
library(parallel)
library(pROC)    # AUROC
library(PRROC)   # AUPRC

# Utility functions
id2idTitle <- function(vec) {
  return(unlist(lapply(vec, function(x) unlist(str_split(x, '\\.'))[1])))
}

jaccard <- function(a, b) {
  if (length(a) == 0 | length(b) == 0) {
    return(0)
  } else {
    intersection = length(intersect(a, b))
    union = length(a) + length(b) - intersection
    return(intersection / union)
  }
}

# Paths
path <- '~/data'
path_save <- '~/result/parameter_discussion_gseaPara'
path_support = '~/result/support_file_longRead'

# Define isoform GSEA function
isoformGSEA_func <- function(bioTypeTerm1Term2, dataType, dataSubType, list_vec_target_trans_id, vec_trans_tpm) {
  sep_char <- ifelse(dataType == 'spatial', '-', '_')
  vecSplit_bioTypeTerm1Term2 <- unlist(str_split(bioTypeTerm1Term2, sep_char))
  bioType <- vecSplit_bioTypeTerm1Term2[1]
  term1 <- vecSplit_bioTypeTerm1Term2[2]
  term2 <- vecSplit_bioTypeTerm1Term2[3]
  vec_m_DTU = c('DRIMSeq', 'DEXSeq','limma', 'edgeR','satuRn',
                'DTUrtle','NBSplice','RATs','SPIT','SUPPA2')
  
  path_DTU <- '~/result'
  
  list_tab_DTU <- list()
  for (m_DTU in vec_m_DTU) {
    result_try = try({
      if (dataType == 'RNAWG_long-read') {
        tab_DTU_ori <- read.table(paste(path_DTU, paste('DTU', dataType, sep = '_'), dataSubType, paste('tab_DTU', bioType, term1, term2, m_DTU, sep = '_'), sep = '/'), sep = '\t', header = TRUE)
      } else {
        tab_DTU_ori <- read.table(paste(path_DTU, paste('DTU', dataType, sep = '_'), dataSubType, paste('tab_DTU', term1, term2, m_DTU, sep = '_'), sep = '/'), sep = '\t', header = TRUE)
      }
      if (m_DTU == 'DRIMSeq') {
        tab_DTU <- tab_DTU_ori %>% filter(grepl('^ENS', feature_id) & (!grepl('-[0-9]$', feature_id)))
        tab_DTU <- tab_DTU %>% select(feature_id, prop_log2fc, pvalue, adj_pvalue)
      } else if (m_DTU == 'DEXSeq') {
        tab_DTU <- tab_DTU_ori %>% filter(grepl('^ENS', featureID) & (!grepl('-[0-9]$', featureID)))
        tab_DTU <- tab_DTU %>% select(featureID, grep('log2fold', colnames(tab_DTU), value = TRUE), pvalue, padj)
      } else if (m_DTU == 'limma') {
        tab_DTU_ori['feature_id'] <- rownames(tab_DTU_ori)
        tab_DTU <- tab_DTU_ori %>% filter(grepl('^ENS', feature_id) & (!grepl('-[0-9]$', feature_id)))
        tab_DTU <- tab_DTU %>% select(feature_id, logFC, P.Value, FDR)
      } else if (m_DTU == 'edgeR') {
        tab_DTU <- tab_DTU_ori %>% filter(grepl('^ENS', feature_id) & (!grepl('-[0-9]$', feature_id)))
        tab_DTU <- tab_DTU %>% select(feature_id, logFC, P.Value, FDR)
      } else if (m_DTU == 'satuRn') {
        tab_DTU_ori['feature_id'] <- rownames(tab_DTU_ori)
        tab_DTU <- tab_DTU_ori %>% filter(grepl('^ENS', feature_id) & (!grepl('-[0-9]$', feature_id)))
        if (sum(is.na(tab_DTU$empirical_pval)) == dim(tab_DTU)[1]) {
          tab_DTU <- tab_DTU %>% select(feature_id, estimates, pval, regular_FDR)
        } else {
          tab_DTU <- tab_DTU %>% select(feature_id, estimates, empirical_pval, empirical_FDR)
        }
      }else if(m_DTU=='DTUrtle'){
        tab_DTU = tab_DTU_ori %>% filter(grepl('^ENS',feature_id)&(!grepl('-[0-9]$',feature_id)))
        tab_DTU = tab_DTU %>% select(feature_id,lr,pvalue,adj_pvalue)
      }else if(m_DTU=='NBSplice'){
        tab_DTU = tab_DTU_ori %>% filter(grepl('^ENS',iso)&(!grepl('-[0-9]$',iso)))
        tab_DTU = tab_DTU %>% select(iso,odd,pval,FDR)
      }else if(m_DTU=='RATs'){
        tab_DTU = tab_DTU_ori %>% filter(grepl('^ENS',target_id)&(!grepl('-[0-9]$',target_id)))
        tab_DTU = tab_DTU %>% select(target_id,log2FC,pval,pval_corr)
      }else if(m_DTU=='SPIT'){
        tab_DTU = tab_DTU_ori %>% filter(grepl('^ENS',transcript_id)&(!grepl('-[0-9]$',transcript_id)))
        tab_DTU = tab_DTU %>% select(transcript_id,likelihood,pvalue)
        tab_DTU['FDR'] = tab_DTU['pvalue']
      }else if(m_DTU=='SUPPA2'){
        tab_DTU = tab_DTU_ori %>% filter(grepl('^ENS',transcript)&(!grepl('-[0-9]$',transcript)))
        tab_DTU = tab_DTU %>% select(transcript,dPSI,p_value)
        tab_DTU['FDR'] = tab_DTU['p_value']
      }
      colnames(tab_DTU) <- c('transcript_id', 'logFC', 'pvalue', 'FDR')
      tab_DTU['transcript_id_title'] <- id2idTitle(tab_DTU$transcript_id)
      tab_DTU <- tab_DTU %>% filter(!is.na(pvalue))
      tab_DTU <- tab_DTU[, c(grep('transcript', colnames(tab_DTU), value = TRUE), grep('gene', colnames(tab_DTU), value = TRUE), colnames(tab_DTU)[!grepl('transcript|gene', colnames(tab_DTU))])]
      tab_DTU['-log10pvalue'] <- -log10(tab_DTU$pvalue)
      
      vec_pvalue <- tab_DTU %>% pull(pvalue)
      vec_pvalue_unique <- sort(unique(vec_pvalue))
      if(vec_pvalue_unique[1]==0){vec_pvalue[vec_pvalue==0] = ifelse(vec_pvalue_unique[2]/2==0,vec_pvalue_unique[2],vec_pvalue_unique[2]/2)} # 如果/2过小显示为0,则取不/2的值
      if (vec_pvalue_unique[length(vec_pvalue_unique)] == 1) {vec_pvalue[vec_pvalue == 1] <- (vec_pvalue_unique[length(vec_pvalue_unique) - 1] + 1) / 2 }
      tab_DTU['pvalue_noInf'] <- vec_pvalue
      tab_DTU['-log10pvalue_noInf'] <- -log10(tab_DTU$pvalue_noInf)
      tab_DTU['-lnpvalue_noInf'] <- -log(tab_DTU$pvalue_noInf)
    })
    if(inherits(result_try,'try-error')){
      list_tab_DTU[[m_DTU]] <- data.frame(
        transcript_id = character(),
        transcript_id_title = character(),
        logFC = numeric(),
        pvalue = numeric(),
        FDR = numeric(),
        `-log10pvalue` = numeric(),
        pvalue_noInf = numeric(),
        `-log10pvalue_noInf` = numeric(),
        `-lnpvalue_noInf` = numeric(),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }else{
      list_tab_DTU[[m_DTU]] <- tab_DTU
    }
  }
  tab_num_DTU_sig <- as.data.frame(table(unlist(lapply(list_tab_DTU, function(x) { x %>% filter(pvalue < 0.05) %>% pull(transcript_id) }))))
  vec_trans_common_DTU <- tab_num_DTU_sig %>% filter(Freq >= (length(vec_m_DTU)/2)) %>% pull(1)
  tab_num_DTU <- as.data.frame(table(unlist(lapply(list_tab_DTU, function(x) { x %>% pull(transcript_id) }))))
  vec_trans_all_DTU <- tab_num_DTU %>% pull(1)
  
  list_tab_DTU_rand <- list()
  list_vec_trans_common_DTU <- list()
  list_tab_DTU_rand[[1]] <- list_tab_DTU
  list_vec_trans_common_DTU[[1]] <- vec_trans_common_DTU
  
  for (random_i in 2:2) {
    list_tab_DTU_rand[[random_i]] <- list()
    for (m_DTU in vec_m_DTU) {
      set.seed(random_i)
      list_tab_DTU_rand[[random_i]][[m_DTU]] <- list_tab_DTU[[m_DTU]] %>% mutate(transcript_id = sample(transcript_id))
    }
    list_vec_trans_common_DTU[[random_i]] <- as.data.frame(table(unlist(lapply(list_tab_DTU_rand[[random_i]],
                                                                               function(x) { x %>% filter(pvalue < 0.05) %>% pull(transcript_id) })))) %>%
      filter(Freq >= (length(vec_m_DTU)/2)) %>% pull(1)
  }
  
  # tab_jaccard_temp <- data.frame()
  list_tab_jaccard_temp <- list()
  for(m_identify in names(list_vec_target_trans_id)){
      vec_regulated_trans = list_vec_target_trans_id[[m_identify]]
    
    for (random_i in 1:1) {
      vec_trans_common_DTU = list_vec_trans_common_DTU[[random_i]]
      for (m_DTU in vec_m_DTU) {
        # print(paste(dataType, dataSubType, bioType, term1, term2,m_identify,random_i,gseaParam, m_DTU, sep = '_'))
        tab_DTU <- list_tab_DTU_rand[[random_i]][[m_DTU]]
        vec_trans_sig_pvalue <- tab_DTU %>% filter(pvalue < 0.05) %>% pull(transcript_id)
        vec_trans_sig_FDR <- tab_DTU %>% filter(FDR < 0.05) %>% pull(transcript_id)
        vec_trans_DTU <- tab_DTU %>% pull(transcript_id)
        
        ## DTU 1-pvalue list
        ranks_DTU = tab_DTU %>% arrange(desc(`-log10pvalue_noInf`)) %>% pull(`-log10pvalue_noInf`)
        names(ranks_DTU) = tab_DTU %>% arrange(desc(`-log10pvalue_noInf`)) %>% pull(transcript_id)
        
        ## common & regulated_trans
        vec_comReg_trans = intersect(vec_trans_common_DTU, vec_regulated_trans)
        
        ## ES comReg
        vec_trans_comReg = intersect(vec_comReg_trans, vec_trans_tpm)
        list_trans_comReg = list(vec_trans_comReg);names(list_trans_comReg) = 'comReg'
        
        for(gseaParam in c(10,9,8,7,6,5,4,3,2,1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0)){
          
          ## ES comReg
          if(length(intersect(vec_trans_comReg, names(ranks_DTU)))>1){
            GSEA_comReg = fgsea(list_trans_comReg, ranks_DTU,scoreType='pos',gseaParam=gseaParam,maxSize=300000,nperm=1000)
            if(dim(GSEA_comReg)[1]==0){ES_comReg = 0; NES_comReg=0; pval_comReg=1}else{
              ES_comReg = GSEA_comReg$ES; NES_comReg = GSEA_comReg$NES; pval_comReg = GSEA_comReg$pval
            }
          }else{ES_comReg = 0; NES_comReg=0; pval_comReg=1}
          
          tab_jaccard_temp_temp = data.frame(dataSubType = dataSubType,
                                             bioType = bioType,
                                             term1 = term1,
                                             term2 = term2,
                                             random_i = random_i,
                                             m_DTU = m_DTU,
                                             m_identify = m_identify,
                                             gseaParam = gseaParam,
                                             
                                             num_trans_sig_pvalue = length(vec_trans_sig_pvalue),
                                             num_trans_sig_FDR = length(vec_trans_sig_FDR),
                                             
                                             num_trans_comReg = length(vec_comReg_trans),
                                             ES_comReg = ES_comReg, NES_comReg = NES_comReg, pval_comReg = pval_comReg)
          # print(tab_jaccard_temp_temp)
          list_tab_jaccard_temp[[paste(m_identify,as.character(random_i),as.character(gseaParam),m_DTU,sep='_')]] = tab_jaccard_temp_temp
        }
        
      }
    }
  }
  
  tab_jaccard_temp = bind_rows(list_tab_jaccard_temp)
  return(tab_jaccard_temp)
}

# Run GSEA in parallel
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
    
    # load support 
    vec_bioTypeTerm1Term2 = readRDS(file.path(path_support,'vec_bioTypeTerm1Term2', paste(paste(dataType, dataSubType,sep='_'),".rds",sep='')))
    list_list_vec_target_trans_id =
      readRDS(file.path(path_support,'list_list_vec_target_trans_id',
                        paste(paste(dataType, dataSubType,sep='_'),".rds",sep='')))
    vec_trans_tpm = readRDS(file.path(path_support,'vec_trans_tpm', paste(paste(dataType, dataSubType,sep='_'),".rds",sep='')))
    
    # run isoformGSEA_func
    list_tab_jaccard = mclapply(vec_bioTypeTerm1Term2, function(bioTypeTerm1Term2) {
      suppressPackageStartupMessages({library(tidyverse);library(fgsea);library(pROC);library(PRROC)})
      isoformGSEA_func(bioTypeTerm1Term2, dataType, dataSubType, list_list_vec_target_trans_id[[bioTypeTerm1Term2]], vec_trans_tpm)
    },
    mc.cores = 30, # 32-2
    mc.preschedule = FALSE,  
    mc.cleanup = TRUE        
    )
    tab_jaccard <- bind_rows(list_tab_jaccard)
    gc()
    
    # Save results
    write.table(tab_jaccard,
                file.path(path_save,
                          paste0('tab_paraDiscuss_', dataType, '_', dataSubType)),
                row.names = FALSE, quote = FALSE, sep = '\t')
  }
}


################################################################################
# Effectiveness, Unbiasedness, and Stability
# Effective, unbiased, and stable
library(tidyverse)
path_save <- '~/result/parameter_discussion_gseaPara'
################################################################################
# merge all result
list_tab_jaccard_gseaPara = list()
for(dataType in c('RNAWG_long-read','single_cell','spatial')){
  if(dataType=='RNAWG_long-read'){vec_dataSubType = c('mouse','human')}
  if(dataType=='single_cell'){vec_dataSubType = c('PromethION_5cl_rep1', 'PromethION_5cl_rep2', 'PromethION_MSC')}
  if(dataType=='spatial'){vec_dataSubType = c('GSE153859_CBS1', 'GSE153859_CBS2', 'GSE153859_MOB')}
  for(dataSubType in vec_dataSubType){
    # jaccard
    tab_jaccard_temp = read.table(file.path(path_save, 
                                            paste('tab_paraDiscuss', dataType, dataSubType,sep='_')),header=T,sep='\t')
    tab_jaccard_temp = tab_jaccard_temp %>% 
      filter(m_identify == 'glm') %>% 
      mutate(dataType=dataType, .before = 1)
    
    # gseaPara to column name
    tab_jaccard_temp <- tab_jaccard_temp %>%
      mutate(gseaParam = as.character(gseaParam)) %>%
      pivot_wider(
        names_from = gseaParam,
        values_from = c(ES_comReg, NES_comReg, pval_comReg),
        names_sep = "_"
      )
    
    list_tab_jaccard_gseaPara[[paste(dataType,dataSubType,sep='_')]] = tab_jaccard_temp

  }
}
################################################################################
# merge jaccard and rand
tab_jaccard_gseaPara = bind_rows(list_tab_jaccard_gseaPara)
tab_jaccard_and_rand_mIdentify <- tab_jaccard_gseaPara 

################################################################################
# analysis 1. Effectiveness 
vec_pval <- c(
  "pval_comReg_10","pval_comReg_9","pval_comReg_8","pval_comReg_7","pval_comReg_6",
  "pval_comReg_5","pval_comReg_4","pval_comReg_3","pval_comReg_2",
  "pval_comReg_1","pval_comReg_0.9","pval_comReg_0.8","pval_comReg_0.7","pval_comReg_0.6",
  "pval_comReg_0.5","pval_comReg_0.4","pval_comReg_0.3","pval_comReg_0.2","pval_comReg_0.1",
  "pval_comReg_0"
)
vec_scoreName = c(
  "NES_comReg_10","NES_comReg_9","NES_comReg_8","NES_comReg_7","NES_comReg_6",
  "NES_comReg_5","NES_comReg_4","NES_comReg_3","NES_comReg_2",
  "NES_comReg_1","NES_comReg_0.9","NES_comReg_0.8","NES_comReg_0.7","NES_comReg_0.6",
  "NES_comReg_0.5","NES_comReg_0.4","NES_comReg_0.3","NES_comReg_0.2","NES_comReg_0.1",
  "NES_comReg_0"
)

# compute medians - result is a single-row dataframe
median_df <- data.frame(t(sapply(tab_jaccard_and_rand_mIdentify[vec_pval], median, na.rm = TRUE)))
colnames(median_df) <- vec_pval
rownames(median_df) <- "median"

# compute means - result is a single-row dataframe
mean_df <- data.frame(t(sapply(tab_jaccard_and_rand_mIdentify[vec_pval], mean, na.rm = TRUE)))
colnames(mean_df) <- vec_pval
rownames(mean_df) <- "mean"

# transpose results so variable names become row names
median_result <- as.data.frame(t(median_df))
mean_result <- as.data.frame(t(mean_df))

# rename columns
colnames(median_result) <- "median"
colnames(mean_result) <- "mean"

tab_effective <- cbind(median_result, mean_result)
rownames(tab_effective) = vec_scoreName

###############################################################################
# analysis2. Unbiasedness
library(tidyr)    
library(purrr)    
library(utils)    

# define columns of interest
columns_of_interest <- c(
  "num_trans_sig_pvalue", 
  "NES_comReg_10","NES_comReg_9","NES_comReg_8","NES_comReg_7","NES_comReg_6",
  "NES_comReg_5","NES_comReg_4","NES_comReg_3","NES_comReg_2",
  "NES_comReg_1","NES_comReg_0.9","NES_comReg_0.8","NES_comReg_0.7","NES_comReg_0.6",
  "NES_comReg_0.5","NES_comReg_0.4","NES_comReg_0.3","NES_comReg_0.2","NES_comReg_0.1",
  "NES_comReg_0"
)

# generate all ordered pairs
ordered_pairs_df <- expand.grid(columns_of_interest, columns_of_interest)
combinations <- apply(ordered_pairs_df, 1, as.character) 
combinations <- lapply(1:ncol(combinations), function(i) combinations[,i]) 
unique_pairs <- sapply(combinations, function(x) paste(x, collapse = "-"))  

# group by dataSubType, term1, term2 and compute correlation for each group
grouped_data <- tab_jaccard_and_rand_mIdentify %>%
  group_by(dataSubType, bioType,term1, term2) %>%
  nest()  

# compute pairwise correlations for each group
grouped_data <- grouped_data %>%
  mutate(
    correlations = map(data, ~ {
      df_sub <- .x[, columns_of_interest, drop = FALSE]  
      cor_list <- lapply(combinations, function(pair) {
        col1 <- pair[1]
        col2 <- pair[2]
        if (all(c(col1, col2) %in% names(df_sub))) {  
          cor(df_sub[[col1]], df_sub[[col2]], use = "complete.obs", method = "spearman")  
        } else {NA  
        }
      })
      names(cor_list) <- unique_pairs  
      return(cor_list)
    })
  )

# unnest the correlation results (convert list to wide dataframe)
cor_df <- grouped_data %>%
  unnest_wider(correlations)

# build final results dataframe with pair1 and pair2 labels
final_results <- data.frame(
  pair1 = sapply(strsplit(unique_pairs, "-"), `[`, 1),  
  pair2 = sapply(strsplit(unique_pairs, "-"), `[`, 2),  
  median_correlation = sapply(unique_pairs, function(p) median(abs(cor_df[[p]]), na.rm = TRUE)),  
  mean_correlation = sapply(unique_pairs, function(p) mean(abs(cor_df[[p]]), na.rm = TRUE))  
)
tab_unbiased = final_results %>% filter(pair2 == 'num_trans_sig_pvalue')


################################################################################
# analysis3. stability (stability within the same large dataset (long, sc, spatial),
# and stability across smaller datasets (like MOB))

library(wCorr)
score_cols <- c(
  "NES_comReg_10","NES_comReg_9","NES_comReg_8","NES_comReg_7","NES_comReg_6",
  "NES_comReg_5","NES_comReg_4","NES_comReg_3","NES_comReg_2",
  "NES_comReg_1","NES_comReg_0.9","NES_comReg_0.8","NES_comReg_0.7","NES_comReg_0.6",
  "NES_comReg_0.5","NES_comReg_0.4","NES_comReg_0.3","NES_comReg_0.2","NES_comReg_0.1",
  "NES_comReg_0"
)
median_scores_all <- tab_jaccard_and_rand_mIdentify %>%
  group_by(m_DTU) %>%
  summarise(
    across(all_of(score_cols), ~ median(.x, na.rm = TRUE)),
    .groups = 'drop' 
  )

# compute ranks for each column; larger value => larger rank
median_scores_all_rank <- as.data.frame(lapply(median_scores_all[,score_cols], function(x) rank(x, ties.method = "average")))
median_scores_all_rank = median_scores_all_rank %>% mutate(m_DTU=median_scores_all[['m_DTU']], .before=1)

# compute median scores within each (dataSubType, term1, m_DTU) across term2
median_scores <- tab_jaccard_and_rand_mIdentify %>%
  group_by(dataType, dataSubType,m_DTU) %>%
  summarise(
    across(all_of(score_cols), ~ median(.x, na.rm = TRUE)),
    .groups = 'drop' 
  )

# create unique group identifier
median_scores <- median_scores %>%
  mutate(group_id = paste(dataType, dataSubType, sep = "-"))

# compute Spearman correlation grouped by score columns
correlation_results_by_score <- list()
for (score_col in score_cols) {
  cat("Calculating correlation for score:", score_col, "\n")
  wide_data_for_score <- median_scores %>%
    select(group_id, m_DTU, !!sym(score_col)) %>% 
    filter(complete.cases(.)) %>%
    pivot_wider(
      names_from = m_DTU,
      values_from = !!sym(score_col)
    )
  numeric_data_for_corr <- wide_data_for_score %>%
    select(-group_id) %>%
    as.data.frame() # 
  rownames(numeric_data_for_corr) <- wide_data_for_score$group_id
  if (nrow(numeric_data_for_corr) < 2 || ncol(numeric_data_for_corr) < 2) {
    cat("Skipping correlation for score '", score_col, "' - Need at least 2 (dataSubType, term1) groups and 2 m_DTU methods with data.\n\n", sep = "")
    correlation_results_by_score[[score_col]] <- NA 
    next 
  }
  
  # compute Spearman correlation matrix
  correlation_matrix <- cor(
    t(numeric_data_for_corr),
    method = "spearman",
    use = "pairwise.complete.obs" 
  )
  
  # store results
  correlation_results_by_score[[score_col]] <- correlation_matrix
}


library(reshape2)
library(purrr) 
all_correlations_long_by_score <- purrr::map_dfr(correlation_results_by_score, ~{
  if(is.matrix(.x)){
    reshape2::melt(.x, varnames = c("Group1", "Group2"), value.name = "Spearman_Correlation")
  } else {
    data.frame(Group1=NA, Group2=NA, Spearman_Correlation=NA)
  }
}, .id = "Score_Column")

all_correlations_long_by_score <- all_correlations_long_by_score %>%
  mutate(
    prefix1 = sub("_.*", "", Group1),
    prefix2 = sub("_.*", "", Group2)
  ) %>%
  filter(prefix1 == prefix2) %>%
  select(-prefix1, -prefix2) 

all_correlations_long_by_score_mean  = 
  all_correlations_long_by_score %>% 
  filter(Group1!=Group2) %>%
  group_by(Score_Column) %>% 
  summarise(mean_Spearman_Correlation = mean(abs(Spearman_Correlation)),
            median_Spearman_Correlation = median(abs(Spearman_Correlation))) %>% 
  arrange(desc(median_Spearman_Correlation))

tab_stable = all_correlations_long_by_score_mean
################################################################################
tab_effective_draw = tab_effective %>% mutate(Effectiveness=1-median) %>% select(Effectiveness)  %>% mutate(ScoreType = rownames(.))
tab_unbiased_draw = tab_unbiased %>% mutate(Unbiasedness=1-abs(median_correlation)) %>% select(Unbiasedness,pair1) %>% rename(ScoreType=pair1)
tab_stable_draw = tab_stable %>% mutate(Stability = median_Spearman_Correlation)%>% select(Stability,Score_Column)%>% rename(ScoreType=Score_Column)
tab_merge_effUnbSta = merge(merge(tab_effective_draw,tab_unbiased_draw,by='ScoreType'),tab_stable_draw, by='ScoreType')
tab_merge_effUnbSta['mean_score'] = (tab_merge_effUnbSta$Effectiveness + tab_merge_effUnbSta$Unbiasedness + tab_merge_effUnbSta$Stability)/3
tab_merge_effUnbSta['geomean_score'] = (tab_merge_effUnbSta$Effectiveness * tab_merge_effUnbSta$Unbiasedness * tab_merge_effUnbSta$Stability)^(1/3)

