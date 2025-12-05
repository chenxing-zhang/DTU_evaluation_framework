# Effectiveness, Unbiasedness, and Stability
# Effective, unbiased, and stable
library(tidyverse)
dataType= 'ENCODE_NGS' # ENCODE_NGS
path_save <- '~/result'
path_save_result <- file.path(path_save, paste0('ISES_', dataType))
################################################################################
# load result
tab_jaccard = read.table(file.path(path_save_result, paste('tab', dataType, 'FCFC_gseaPara','Simple',sep='_')),header=T,sep='\t')
tab_jaccard_gseaPara1 = tab_jaccard %>% filter(gseaParam == 1)

tab_rand = read.table(file.path(path_save_result, paste('tab', dataType, 'noNES_rand',sep='_')),header=T,sep='\t')

################################################################################
# analysis 1. Effectiveness
# Define grouping variables and the metric columns to compute p-values
group_vars <- c("dataSubType", "term1", "term2", "m_DTU")
metric_vars <- c(
  "jaccard_com", "auroc_com", "auprc_com",
  "jaccard_reg", "auroc_reg", "auprc_reg",
  "jaccard_comReg", "auroc_comReg", "auprc_comReg"
)

# calculate p-value
result_df <- tab_rand %>%
  group_by(across(all_of(group_vars))) %>% 
  summarise(
    .groups = "drop", 
    p_values = {
      current_group_data <- cur_data() 
      baseline_row <- current_group_data %>% filter(random_i == 1)
      random_comparison_rows <- current_group_data %>% filter(random_i != 1)
      p_val_list <- list()
      # Case: exactly 1 baseline row and valid random rows
      if (nrow(baseline_row) == 1 && nrow(random_comparison_rows) > 0) {
        n <- nrow(random_comparison_rows) 
        for (metric in metric_vars) {
          baseline_value <- baseline_row[[metric]][1] 
          m <- sum(random_comparison_rows[[metric]] >= baseline_value, na.rm = TRUE)
          p_val_list[[paste0("pval_", metric)]] <- (m + 1) / (n + 1)
        }
      
      # Case: only baseline row, no random rows
      } else if (nrow(baseline_row) == 1 && nrow(random_comparison_rows) == 0) {
        n <- 0
        m <- 0
        for (metric in metric_vars) {
          p_val_list[[paste0("p_", metric)]] <- (m + 1) / (n + 1)
        }
      # Unexpected group: missing or duplicated baseline rows
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
  tidyr::unnest(p_values) # expand list-column of p-values

# merge jaccard and rand 
tab_jaccard_and_rand = merge(tab_jaccard_gseaPara1, result_df, by=group_vars, all.x=T)

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

# compute medians (one-row dataframe)
median_df <- data.frame(t(sapply(tab_jaccard_and_rand[vec_pval], median, na.rm = TRUE)))
colnames(median_df) <- vec_pval
rownames(median_df) <- "median"

# compute means (one-row dataframe)
mean_df <- data.frame(t(sapply(tab_jaccard_and_rand[vec_pval], mean, na.rm = TRUE)))
colnames(mean_df) <- vec_pval
rownames(mean_df) <- "mean"

# transpose: variable names become rows
median_result <- as.data.frame(t(median_df))
mean_result <- as.data.frame(t(mean_df))

# 重命名列
colnames(median_result) <- "median"
colnames(mean_result) <- "mean"

tab_effective <- cbind(median_result, mean_result)
rownames(tab_effective) = vec_scoreName


################################################################################
# analysis 2. Unbiasedness
library(tidyr)  
library(purrr)    
library(utils) 

# columns to evaluate
columns_of_interest <- c(
  "num_trans_sig_pvalue", "NES_sig", "NES_com", "jaccard_com",
  "auroc_com", "auprc_com", "NES_reg", "jaccard_reg",
  "auroc_reg", "auprc_reg", "NES_comReg", "jaccard_comReg",
  "auroc_comReg", "auprc_comReg"
)

# generate all column pairs
ordered_pairs_df <- expand.grid(columns_of_interest, columns_of_interest)
combinations <- apply(ordered_pairs_df, 1, as.character) 
combinations <- lapply(1:ncol(combinations), function(i) combinations[,i])
unique_pairs <- sapply(combinations, function(x) paste(x, collapse = "-")) 

# group by dataset and compute correlations
grouped_data <- tab_jaccard_gseaPara1 %>%
  group_by(dataSubType, term1, term2) %>%
  nest()  

# compute pairwise correlations per group
grouped_data <- grouped_data %>%
  mutate(
    correlations = map(data, ~ {
      df_sub <- .x[, columns_of_interest, drop = FALSE]  
      cor_list <- lapply(combinations, function(pair) {
        col1 <- pair[1]
        col2 <- pair[2]
        if (all(c(col1, col2) %in% names(df_sub))) { 
          cor(df_sub[[col1]], df_sub[[col2]], use = "complete.obs", method = "spearman")  # 计算Pearson相关性
        } else {
          NA  
        }
      })
      names(cor_list) <- unique_pairs  
      return(cor_list)
    })
  )

# expand results
cor_df <- grouped_data %>%
  unnest_wider(correlations) 

# summarize unbiasedness metrics
final_results <- data.frame(
  pair1 = sapply(strsplit(unique_pairs, "-"), `[`, 1),  
  pair2 = sapply(strsplit(unique_pairs, "-"), `[`, 2), 
  median_correlation = sapply(unique_pairs, function(p) median(abs(cor_df[[p]]), na.rm = TRUE)), 
  mean_correlation = sapply(unique_pairs, function(p) mean(abs(cor_df[[p]]), na.rm = TRUE)) 
)
tab_unbiased = final_results %>% filter(pair2 == 'num_trans_sig_pvalue')

################################################################################
# analysis 3. stability (stability across small datasets within a dataset: long / sc / spatial)

library(wCorr)
score_cols <- c(
  "NES_sig", "NES_com", "jaccard_com",
  "auroc_com", "auprc_com", "NES_reg", "jaccard_reg",
  "auroc_reg", "auprc_reg", "NES_comReg", "jaccard_comReg",
  "auroc_comReg", "auprc_comReg"
)

# median per m_DTU
median_scores_all <- tab_jaccard_gseaPara1 %>%
  group_by(m_DTU) %>%
  summarise(
    across(all_of(score_cols), ~ median(.x, na.rm = TRUE)),
    .groups = 'drop' 
  )

# compute rank for each score column
median_scores_all_rank <- as.data.frame(lapply(median_scores_all[,score_cols], function(x) rank(x, ties.method = "average")))
median_scores_all_rank = median_scores_all_rank %>% mutate(m_DTU=median_scores_all[['m_DTU']], .before=1)


# compute median score per (dataSubType, term1, m_DTU)
median_scores <- tab_jaccard_gseaPara1 %>%
  group_by(dataSubType, term1, m_DTU) %>%
  summarise(
    across(all_of(score_cols), ~ median(.x, na.rm = TRUE)),
    .groups = 'drop' # 计算完后解除分组
  )

# create unique group id
median_scores <- median_scores %>%
  mutate(group_id = paste(dataSubType, term1, sep = "_"))

# compute Spearman correlation for each score column
correlation_results_by_score <- list()
for (score_col in score_cols) {
  wide_data_for_score <- median_scores %>%
    select(group_id, m_DTU, !!sym(score_col)) %>% 
    filter(complete.cases(.)) %>%
    pivot_wider(
      names_from = m_DTU,
      values_from = !!sym(score_col)
    )
  numeric_data_for_corr <- wide_data_for_score %>%
    select(-group_id) %>%
    as.data.frame() 
  rownames(numeric_data_for_corr) <- wide_data_for_score$group_id
  if (nrow(numeric_data_for_corr) < 2 || ncol(numeric_data_for_corr) < 2) {
    cat("Skipping correlation for score '", score_col, "' - Need at least 2 (dataSubType, term1) groups and 2 m_DTU methods with data.\n\n", sep = "")
    correlation_results_by_score[[score_col]] <- NA 
    next 
  }
  correlation_matrix <- cor(
    t(numeric_data_for_corr),
    method = "spearman",
    use = "pairwise.complete.obs" 
  )
  weighted_correlation_matrix = correlation_matrix
  for(group_id_1 in rownames(weighted_correlation_matrix)){
    for(group_id_2 in rownames(weighted_correlation_matrix)){
      weighted_correlation_matrix[group_id_1,group_id_2] = 
        weightedCorr(as.numeric(numeric_data_for_corr[group_id_1,]), 
                     as.numeric(numeric_data_for_corr[group_id_2,]), 
                     method='Spearman',weights=median_scores_all_rank[[score_col]])
    }
  }
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
print(all_correlations_long_by_score)

all_correlations_long_by_score_mean  = 
  all_correlations_long_by_score %>% 
  filter(Group1!=Group2) %>%
  group_by(Score_Column) %>% 
  summarise(mean_Spearman_Correlation = mean(abs(Spearman_Correlation)),
            median_Spearman_Correlation = median(abs(Spearman_Correlation))) %>% 
  arrange(desc(median_Spearman_Correlation))
tab_stable = all_correlations_long_by_score_mean

################################################################################
# drawing tables for final plotting
tab_effective_draw = tab_effective %>% mutate(Effectiveness=1-median) %>% select(Effectiveness)  %>% mutate(ScoreType = rownames(.))
tab_unbiased_draw = tab_unbiased %>% mutate(Unbiasedness=1-abs(median_correlation)) %>% select(Unbiasedness,pair1) %>% rename(ScoreType=pair1)
tab_stable_draw = tab_stable %>% mutate(Stability = median_Spearman_Correlation)%>% select(Stability,Score_Column)%>% rename(ScoreType=Score_Column)
tab_merge_effUnbSta = merge(merge(tab_effective_draw,tab_unbiased_draw,by='ScoreType',all.x=T),tab_stable_draw, by='ScoreType',all.x=T)
tab_merge_effUnbSta['mean_score'] = (tab_merge_effUnbSta$Effectiveness + tab_merge_effUnbSta$Unbiasedness + tab_merge_effUnbSta$Stability)/3
