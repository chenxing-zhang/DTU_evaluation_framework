library(tidyverse)
library(data.table)
library(bedtoolsr)
options(bedtools.path = "~/anaconda3/envs/py_fastq/bin")

id2idTitle <- function(vec) {
  return(unlist(lapply(vec, function(x) unlist(str_split(x, '\\.'))[1])))
}

## predict regulated isoform
# load model
path_model = '~/result/predict_regulated_isoform'
glm_model  = readRDS(file.path(path_model, "glm_model_allData_2feature.rds"))

# Load RBP info
tab_RBPinfo <- read.csv('~/data/RBP_RNA/CLIPdb/POSTAR3_CLIPdb_module_browse_RBP_sub_info.csv', check.names = FALSE)
colnames(tab_RBPinfo) <- c('RBPname', 'geneID', 'Domain', 'Location', 'Function')
tab_RBPinfo_splice <- tab_RBPinfo[grepl('Splic', tab_RBPinfo$Function), ]
vec_RBP_RBPinfo <- tab_RBPinfo_splice[['RBPname']]

species = 'human'

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
    dplyr::select(seqnames, start, end, RBP_name, experiment_method, sample) %>%
    filter(experiment_method == 'eCLIP' & sample != 'adrenal gland')
} else if (species == 'mouse') {
  A.bed <- tab_CLIPdb %>% filter(RBP_name %in% vec_RBP_RBPinfo) %>%
    dplyr::select(seqnames, start, end, RBP_name, experiment_method, sample)
  A.bed$RBP_name <- str_to_title(A.bed$RBP_name)
}

# Load annotation data
if (species == 'human') {
  tab_anno_trans <- read.table('~/data/Ref/GENECODE_r45/tab_anno_trans', sep = '\t', header = TRUE)
  tab_anno_gene <- read.table('~/data/Ref/GENECODE_r45/tab_anno_gene', sep = '\t', header = TRUE)
} else {
  tab_anno_trans <- read.table('~/data/Ref_mouse/GENECODE_M34/tab_anno_trans', sep = '\t', header = TRUE)
  tab_anno_gene <- read.table('~/data/Ref_mouse/GENECODE_M34/tab_anno_gene', sep = '\t', header = TRUE)
}

tab_anno_trans <- tab_anno_trans %>% mutate(transcript_id_title = id2idTitle(transcript_id), .before = transcript_id) %>%
  filter(!grepl('PAR', transcript_id) & !grepl('PAR', gene_id))
tab_anno_gene <- tab_anno_gene %>% mutate(gene_id_title = id2idTitle(gene_id), .before = gene_id) %>%
  filter(!grepl('PAR', gene_id))

B.bed <- tab_anno_trans %>% 
  dplyr::select(seqnames, start, end, transcript_id, transcript_name, gene_id, gene_name) %>%
  arrange(seqnames, start, end)

# Intersect A and B
inter <- function(A, B, spe) {
  vec_merge_col <- c('RNA_seqnames', 'RNA_start', 'RNA_end', 'RBP_name',
                     'experiment_method', 'sample',
                     'target_seqnames', 'target_start', 'target_end',
                     'target_trans_id', 'target_trans_name',
                     'target_gene_id', 'target_gene_name',
                     'num_overlap_seq')
  tab_inter_temp <- as.data.table(bedtoolsr::bt.intersect(a = A, b = B, wo = TRUE, f = 1))
  if (dim(tab_inter_temp)[1] == 0) {
    return(data.frame())
  }
  colnames(tab_inter_temp) <- vec_merge_col
  tab_inter_temp <- tab_inter_temp %>% dplyr::select(RBP_name, experiment_method, sample,
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
    dplyr::select(RBP_name, target_trans_id, target_trans_name, target_gene_id, target_gene_name) %>% unique() %>%
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


# reload tpm tab_gene and trans
GEO_id = "GSE294883_NCI-H23" # case_ISES: GSE216825, GSE225633, GSE225637, GSE294883_NCI-H23

path_data = paste0('~/data/case/',GEO_id,'/count')
path_support = paste0('~/result/RBP_activity_prediction/support_file_case/',GEO_id)

if(!dir.exists(path_support)){dir.create(path_support)}

# # load count
tab_tpm_gene = read.table(paste(path_data,'tab_gene_TPM.txt',sep='/'))
tab_tpm_trans  = read.table(paste(path_data,'tab_isoform_TPM.txt',sep='/'))

# add gene name and trans name
tab_tpm_gene = merge(tab_tpm_gene, tab_anno_gene[,c('gene_id','gene_name')],by='gene_id',all.x=T)
tab_tpm_trans = merge(tab_tpm_trans, tab_anno_trans[,c('transcript_id','transcript_name','gene_name')],by='transcript_id',all.x=T)

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
list_samples = c(vec_sample1, vec_sample2)

vec_columns_gene <- c('gene_id', 'gene_name')
tab_tpm_gene_term <- tab_tpm_gene[rowSums(tab_tpm_gene[, list_samples] > 1) >= 2, c(vec_columns_gene, list_samples)]
tab_tpm_gene_term['abs_logFC'] <- abs(log2((rowMeans(tab_tpm_gene_term[, vec_sample1]) + 1) / (rowMeans(tab_tpm_gene_term[, vec_sample2]) + 1)))
tab_tpm_gene_term <- tab_tpm_gene_term %>%
  filter(gene_name %in% vec_RBP_inter) %>% dplyr::select(c(vec_columns_gene, 'abs_logFC')) %>%
  arrange(desc(abs_logFC))
colnames(tab_tpm_gene_term) <- c(vec_columns_gene, 'RBP_abs_logFC')

vec_columns_trans <- c('transcript_id', 'transcript_name', 'gene_id', 'gene_name')
tab_count_trans_term_RBP <- tab_tpm_trans[, c(vec_columns_trans, list_samples)]
tab_count_trans_term_RBP <- tab_count_trans_term_RBP[rowSums(tab_count_trans_term_RBP[, list_samples] > 1) >= 2, ]
tab_tpm_trans_mean <- cbind(tab_count_trans_term_RBP[, vec_columns_trans],
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
tab_trans_diffRate <- tab_tpm_trans_mean_merge %>% dplyr::select(transcript_id, abs_logFC, abs_diff_rate)
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
pred = predict(glm_model, newdata = tab_inter_term, type = "prob")[, 2]
tab_inter_term['predClass'] = ifelse(pred > 0.5, 'untarget', 'target')
vec_regulated_trans_id = tab_inter_term %>%
  filter(predClass=='target') %>%
  pull(target_trans_id) %>% unique()
saveRDS(vec_regulated_trans_id, file.path(paste(path_support,'vec_regulated_trans_id.rds',sep='/')))


##############################################################################
## common isoforms
GEO_id = "GSE294883_NCI-H23" # case_ISES: GSE216825, GSE225633, GSE225637, GSE294883_NCI-H23

path_DTU = paste0('~/result/RBP_activity_prediction/DTU_case/',GEO_id)
path_support = paste('~/result/case_ISES',GEO_id,'support_file',sep='/')

if(!dir.exists(path_support)){dir.create(path_support)}

term1 = 'WT'; term2 = 'EXP'
vec_m_DTU = c('DRIMSeq', 'DEXSeq','limma', 'edgeR','satuRn',
              'DTUrtle','NBSplice','RATs','SPIT','SUPPA2')
list_tab_DTU <- list()
for (m_DTU in vec_m_DTU) {
  result_try = try({
    tab_DTU_ori <- read.table(paste(path_DTU, paste('tab_DTU', term1, term2, m_DTU, sep = '_'), sep = '/'), sep = '\t', header = TRUE)
    
    if (m_DTU == 'DRIMSeq') {
      tab_DTU <- tab_DTU_ori %>% filter(grepl('^ENS', feature_id) & (!grepl('-[0-9]$', feature_id)))
      tab_DTU <- tab_DTU %>% dplyr::select(feature_id, prop_log2fc, pvalue, adj_pvalue)
    } else if (m_DTU == 'DEXSeq') {
      tab_DTU <- tab_DTU_ori %>% filter(grepl('^ENS', featureID) & (!grepl('-[0-9]$', featureID)))
      tab_DTU <- tab_DTU %>% dplyr::select(featureID, grep('log2fold', colnames(tab_DTU), value = TRUE), pvalue, padj)
    } else if (m_DTU == 'limma') {
      tab_DTU_ori['feature_id'] <- rownames(tab_DTU_ori)
      tab_DTU <- tab_DTU_ori %>% filter(grepl('^ENS', feature_id) & (!grepl('-[0-9]$', feature_id)))
      tab_DTU <- tab_DTU %>% dplyr::select(feature_id, logFC, P.Value, FDR)
    } else if (m_DTU == 'edgeR') {
      tab_DTU <- tab_DTU_ori %>% filter(grepl('^ENS', feature_id) & (!grepl('-[0-9]$', feature_id)))
      tab_DTU <- tab_DTU %>% dplyr::select(feature_id, logFC, P.Value, FDR)
    } else if (m_DTU == 'satuRn') {
      tab_DTU_ori['feature_id'] <- rownames(tab_DTU_ori)
      tab_DTU <- tab_DTU_ori %>% filter(grepl('^ENS', feature_id) & (!grepl('-[0-9]$', feature_id)))
      if (sum(is.na(tab_DTU$empirical_pval)) == dim(tab_DTU)[1]) {
        tab_DTU <- tab_DTU %>% dplyr::select(feature_id, estimates, pval, regular_FDR)
      } else {
        tab_DTU <- tab_DTU %>% dplyr::select(feature_id, estimates, empirical_pval, empirical_FDR)
      }
    }else if(m_DTU=='DTUrtle'){
      tab_DTU = tab_DTU_ori %>% filter(grepl('^ENS',feature_id)&(!grepl('-[0-9]$',feature_id)))
      tab_DTU = tab_DTU %>% dplyr::select(feature_id,lr,pvalue,adj_pvalue)
    }else if(m_DTU=='NBSplice'){
      tab_DTU = tab_DTU_ori %>% filter(grepl('^ENS',iso)&(!grepl('-[0-9]$',iso)))
      tab_DTU = tab_DTU %>% dplyr::select(iso,odd,pval,FDR)
    }else if(m_DTU=='RATs'){
      tab_DTU = tab_DTU_ori %>% filter(grepl('^ENS',target_id)&(!grepl('-[0-9]$',target_id)))
      tab_DTU = tab_DTU %>% dplyr::select(target_id,log2FC,pval,pval_corr)
    }else if(m_DTU=='SPIT'){
      tab_DTU = tab_DTU_ori %>% filter(grepl('^ENS',transcript_id)&(!grepl('-[0-9]$',transcript_id)))
      tab_DTU = tab_DTU %>% dplyr::select(transcript_id,likelihood,pvalue)
      tab_DTU['FDR'] = tab_DTU['pvalue']
    }else if(m_DTU=='SUPPA2'){
      tab_DTU = tab_DTU_ori %>% filter(grepl('^ENS',transcript)&(!grepl('-[0-9]$',transcript)))
      tab_DTU = tab_DTU %>% dplyr::select(transcript,dPSI,p_value)
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
saveRDS(vec_trans_common_DTU, file.path(paste(path_support,'vec_common_trans_id.rds',sep='/')))

##############################################################################
# overlapped trans
vec_trans_common_DTU  = readRDS(file.path(path_support, "vec_common_trans_id.rds"))
vec_regulated_trans_id  = readRDS(file.path(path_support, "vec_regulated_trans_id.rds"))
vec_overlapped_trans_id = intersect(vec_trans_common_DTU,vec_regulated_trans_id)
saveRDS(vec_overlapped_trans_id, file.path(paste(path_support,'vec_overlapped_trans_id.rds',sep='/')))
