# Effectiveness, Unbiasedness, and Stability
# Effective, unbiased, and stable
library(tidyverse)
path_save <- '~/result'
################################################################################
# merge all result
list_tab_jaccard_gseaPara1 = list()
list_tab_rand = list()
for(dataType in c('RNAWG_long-read','single_cell','spatial')){
  if(dataType=='RNAWG_long-read'){vec_dataSubType = c('mouse','human')}
  if(dataType=='single_cell'){vec_dataSubType = c('PromethION_5cl_rep1', 'PromethION_5cl_rep2', 'PromethION_MSC')}
  if(dataType=='spatial'){vec_dataSubType = c('GSE153859_CBS1', 'GSE153859_CBS2', 'GSE153859_MOB')}
  path_save_result <- file.path(path_save, paste0('ISES_', dataType))
  for(dataSubType in vec_dataSubType){
    # jaccard
    # classification2feature, filterTopTop
    tab_jaccard_temp = read.table(file.path(path_save_result, 
                                            paste('tab_rand2_jaccard', dataType, dataSubType, 'classification2feature','gseaPara','Simple',sep='_')),header=T,sep='\t')
    list_tab_jaccard_gseaPara1[[paste(dataType,dataSubType,sep='_')]] = 
      tab_jaccard_temp %>% 
      filter(gseaParam == 1) %>% 
      mutate(dataType=dataType, .before = 1)
    
    # rand
    tab_rand_temp = read.table(file.path(path_save_result, 
                                         paste('tab_jaccard', dataType, dataSubType, 'classification2feature','noNES_rand',sep='_')),header=T,sep='\t')
    list_tab_rand[[paste(dataType,dataSubType,sep='_')]] = 
      tab_rand_temp %>% 
      mutate(dataType=dataType, .before = 1)
  }
}

tab_jaccard_gseaPara1 = bind_rows(list_tab_jaccard_gseaPara1)
tab_rand = bind_rows(list_tab_rand)

################################################################################
# Effectiveness
# Define grouping variables and the metric columns for p-value calculation
group_vars <- c("dataType", "dataSubType", "bioType", "term1", "term2", "m_DTU", "m_identify")
metric_vars <- c(
  "jaccard_com", "auroc_com", "auprc_com",
  "jaccard_reg", "auroc_reg", "auprc_reg",
  "jaccard_comReg", "auroc_comReg", "auprc_comReg"
)

# Perform grouped calculations
result_df <- tab_rand %>%
  group_by(across(all_of(group_vars))) %>% 
  summarise(
    .groups = "drop",
    p_values = {
      current_group_data <- cur_data()
      baseline_row <- current_group_data %>% filter(random_i == 1)
      random_comparison_rows <- current_group_data %>% filter(random_i != 1)
      p_val_list <- list()
      if (nrow(baseline_row) == 1 && nrow(random_comparison_rows) > 0) {
        n <- nrow(random_comparison_rows) 
        for (metric in metric_vars) {
          baseline_value <- baseline_row[[metric]][1] 
          # compute m: count of random_i != 1 rows where metric >= baseline
          m <- sum(random_comparison_rows[[metric]] >= baseline_value, na.rm = TRUE)
          # compute p-value using (m+1)/(n+1)
          p_val_list[[paste0("pval_", metric)]] <- (m + 1) / (n + 1)
        }
      } else if (nrow(baseline_row) == 1 && nrow(random_comparison_rows) == 0) {
        # if only the baseline row exists and there are no random rows to compare
        # treat as (m+1)/(n+1) = (0+1)/(0+1) = 1
        n <- 0
        m <- 0
        for (metric in metric_vars) {
          p_val_list[[paste0("p_", metric)]] <- (m + 1) / (n + 1)
        }
      } else {
        warning(paste("Group with", paste(current_group_data[1, group_vars], collapse=", "), 
                      "does not have a unique baseline (random_i == 1) or has issues."))
        for (metric in metric_vars) {
          p_val_list[[paste0("p_", metric)]] <- NA_real_
        }
      }
      as_tibble(p_val_list)
    }
  ) %>%
  tidyr::unnest(p_values)

# merge jaccard and rand 
tab_jaccard_and_rand = merge(tab_jaccard_gseaPara1, result_df, by=group_vars, all.x=T)

# filter
m_identify = 'glm'
gseaParam = 1
random_i = 1
tab_jaccard_and_rand_mIdentify <- tab_jaccard_and_rand %>%
  filter(m_identify == !!m_identify,gseaParam ==!!gseaParam, random_i ==!!random_i) 


vec_pval <- c(
  "pval_com", "pval_jaccard_com", "pval_auroc_com", "pval_auprc_com",
  "pval_reg", "pval_jaccard_reg", "pval_auroc_reg", "pval_auprc_reg",
  "pval_comReg", "pval_jaccard_comReg", "pval_auroc_comReg", "pval_auprc_comReg"
)
vec_scoreName = c(
  "NES_com", "jaccard_com", "auroc_com", "auprc_com",
  "NES_reg", "jaccard_reg", "auroc_reg", "auprc_reg",
  "NES_comReg", "jaccard_comReg", "auroc_comReg", "auprc_comReg"
)

# compute medians - result is a single-row data frame
median_df <- data.frame(t(sapply(tab_jaccard_and_rand_mIdentify[vec_pval], median, na.rm = TRUE)))
colnames(median_df) <- vec_pval
rownames(median_df) <- "median"

# compute means - result is a single-row data frame
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
  "num_trans_sig_pvalue", "NES_com", "jaccard_com",
  "auroc_com", "auprc_com", "NES_reg", "jaccard_reg",
  "auroc_reg", "auprc_reg", "NES_comReg", "jaccard_comReg",
  "auroc_comReg", "auprc_comReg"
)

# generate all ordered pairs (combinations)
ordered_pairs_df <- expand.grid(columns_of_interest, columns_of_interest)
combinations <- apply(ordered_pairs_df, 1, as.character) 
combinations <- lapply(1:ncol(combinations), function(i) combinations[,i]) 

unique_pairs <- sapply(combinations, function(x) paste(x, collapse = "-"))  
# group by dataSubType, term1, term2 and compute correlations for each group
grouped_data <- tab_jaccard_and_rand_mIdentify %>%
  group_by(dataSubType, bioType,term1, term2) %>%
  nest()  

# compute all pairwise correlations for each nested group
grouped_data <- grouped_data %>%
  mutate(
    correlations = map(data, ~ {
      df_sub <- .x[, columns_of_interest, drop = FALSE]  
      cor_list <- lapply(combinations, function(pair) {
        col1 <- pair[1]
        col2 <- pair[2]
        if (all(c(col1, col2) %in% names(df_sub))) {  
          cor(df_sub[[col1]], df_sub[[col2]], use = "complete.obs", method = "spearman")  # compute Spearman correlation
        } else {NA  
        }
      })
      names(cor_list) <- unique_pairs  # name the list elements
      return(cor_list)
    })
  )

# unnest correlation results (convert list to wide data frame)
cor_df <- grouped_data %>%
  unnest_wider(correlations)  

# create final results data frame, ensuring pair1 and pair2 are shown
final_results <- data.frame(
  pair1 = sapply(strsplit(unique_pairs, "-"), `[`, 1),  # first column name, e.g. "num_trans_sig_pvalue"
  pair2 = sapply(strsplit(unique_pairs, "-"), `[`, 2),  # second column name, e.g. "NES_com"
  median_correlation = sapply(unique_pairs, function(p) median(abs(cor_df[[p]]), na.rm = TRUE)),  # median of absolute correlations
  mean_correlation = sapply(unique_pairs, function(p) mean(abs(cor_df[[p]]), na.rm = TRUE))  # mean of absolute correlations
)
tab_unbiased = final_results %>% filter(pair2 == 'num_trans_sig_pvalue')


################################################################################
# analysis3. stability (stability among small datasets within the same large dataset: e.g., MOB among spatial)
library(wCorr)
score_cols <- c(
  "NES_com", "jaccard_com",
  "auroc_com", "auprc_com", "NES_reg", "jaccard_reg",
  "auroc_reg", "auprc_reg", "NES_comReg", "jaccard_comReg",
  "auroc_comReg", "auprc_comReg"
)
median_scores_all <- tab_jaccard_and_rand_mIdentify %>%
  group_by(m_DTU) %>%
  summarise(
    across(all_of(score_cols), ~ median(.x, na.rm = TRUE)),
    .groups = 'drop' 
  )
# compute ranks for each column; larger values get larger ranks
median_scores_all_rank <- as.data.frame(lapply(median_scores_all[,score_cols], function(x) rank(x, ties.method = "average")))
median_scores_all_rank = median_scores_all_rank %>% mutate(m_DTU=median_scores_all[['m_DTU']], .before=1)


# compute median scores per dataSubType, term1, and m_DTU
median_scores <- tab_jaccard_and_rand_mIdentify %>%
  group_by(dataType, dataSubType,m_DTU) %>%
  summarise(
    across(all_of(score_cols), ~ median(.x, na.rm = TRUE)),
    .groups = 'drop' # ungroup after calculation
  )

# create a unique group identifier
median_scores <- median_scores %>%
  mutate(group_id = paste(dataType, dataSubType, sep = "-"))

# compute Spearman correlations per score column
# create a list to store correlation matrices per score column
correlation_results_by_score <- list()

# loop through each score column
for (score_col in score_cols) {
  cat("Calculating correlation for score:", score_col, "\n")
  # reshape to wide format: rows = group_id, columns = m_DTU, values = score_col median
  wide_data_for_score <- median_scores %>%
    select(group_id, m_DTU, !!sym(score_col)) %>% 
    filter(complete.cases(.)) %>%
    pivot_wider(
      names_from = m_DTU,
      values_from = !!sym(score_col)
    )
  
  numeric_data_for_corr <- wide_data_for_score %>%
    select(-group_id) %>%
    as.data.frame() # ensure a data frame
  rownames(numeric_data_for_corr) <- wide_data_for_score$group_id
  
  
  # check if there are enough groups and methods to compute correlation
  # need at least 2 group_id (rows) and 2 m_DTU (columns)
  if (nrow(numeric_data_for_corr) < 2 || ncol(numeric_data_for_corr) < 2) {
    cat("Skipping correlation for score '", score_col, "' - Need at least 2 (dataSubType, term1) groups and 2 m_DTU methods with data.\n\n", sep = "")
    correlation_results_by_score[[score_col]] <- NA # or NULL or a message
    next # go to next score column
  }
  
  # compute Spearman correlation matrix
  correlation_matrix <- cor(
    t(numeric_data_for_corr),
    method = "spearman",
    use = "pairwise.complete.obs" # handle possible NAs in rows
  )
  
  # store results in list
  correlation_results_by_score[[score_col]] <- correlation_matrix
  
}

# if needed, further process correlation_results_by_score
# e.g., combine all results into a long-format data frame
library(reshape2)
library(purrr) 
all_correlations_long_by_score <- purrr::map_dfr(correlation_results_by_score, ~{
  if(is.matrix(.x)){
    reshape2::melt(.x, varnames = c("Group1", "Group2"), value.name = "Spearman_Correlation")
  } else {
    data.frame(Group1=NA, Group2=NA, Spearman_Correlation=NA)
  }
}, .id = "Score_Column")

# extract the prefix before the first "_" from Group1 and Group2
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
# tab_merge_effUnbSta['mean_score'] = (tab_merge_effUnbSta$Effectiveness * tab_merge_effUnbSta$Unbiasedness * tab_merge_effUnbSta$Stability)^(1/3)
