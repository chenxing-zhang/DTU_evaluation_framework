library(tidyverse)
library(data.table)
library(parallel)
library(bedtoolsr)
options(bedtools.path = "~/anaconda3/envs/py_fastq/bin")

id2idTitle <- function(vec) {
  return(unlist(lapply(vec, function(x) unlist(str_split(x, '\\.'))[1])))
}

path_support = '~/result/support_file_shortRead'
path_DTU = '~/result'
path_data = '~/data'

dataType= 'ENCODE_NGS' 
species = 'human'

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

B.bed <- tab_anno_trans %>% select(seqnames, start, end, transcript_id, transcript_name, gene_id, gene_name) %>%
  arrange(seqnames, start, end)

# Intersect A and B
inter <- function(A, B, spe, dataType) {
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
  } else if (spe == 'human' & dataType=='ENCODE_NGS'){
    tab_inter_temp_num_filer = tab_inter_temp_num %>% arrange(sample, RBP_name, target_trans_id)
  } else{
    tab_inter_temp_num_filer <- tab_inter_temp_num %>%
      filter(num_sample >= 2) %>%
      arrange(RBP_name, target_trans_id)
  }
  
  if(dataType=='ENCODE_NGS'){
    tab_inter_temp_num_filer_uni = tab_inter_temp_num_filer %>% select(sample, RBP_name, target_trans_id, target_trans_name, target_gene_id, target_gene_name)
  } else{
    tab_inter_temp_num_filer_uni <- tab_inter_temp_num_filer %>%
      select(RBP_name, target_trans_id, target_trans_name, target_gene_id, target_gene_name) %>% unique() %>%
      arrange(RBP_name, target_trans_id, target_trans_name, target_gene_id, target_gene_name)
  }
  return(tab_inter_temp_num_filer_uni)
}

vec_seqnames <- c("chr1", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
                  "chr16", "chr17", "chr18", "chr19", "chr2", "chr20", "chr21",
                  "chr22", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
                  "chr9", "chrM", "chrX", "chrY")

cl <- makeCluster(20, type = "FORK")
list_tab_inter <- parLapplyLB(cl, vec_seqnames,
                              function(seqName) inter(A.bed %>% filter(seqnames == seqName),
                                                      B.bed %>% filter(seqnames == seqName),
                                                      species,dataType), chunk.size = 1)
stopCluster(cl)
tab_inter <- bind_rows(list_tab_inter)
rm(list_tab_inter)
gc()

# create support data
list_vec_trans_tpm = list()
vec_dataSubTypeTermRBP = c()
list_ranks_FCFC = list() 
list_vec_regulated_trans = list()
# dataSubType = 'CRISPR'
for(dataSubType in c('CRISPR','shRNA')){
  # load exp, experiment, trans
  tab_meta_ep_trans = read.csv(paste(path_data,dataType,dataSubType,'metadata.tsv',sep='/'),check.names = F,sep='\t')
  tab_meta_ep_trans = tab_meta_ep_trans %>% filter(`Genome annotation`=='V29')%>% 
    filter(`Output type` == 'transcript quantifications')%>% filter(Size>30000000) # RSEM quantification
  tab_meta_ep_trans = tab_meta_ep_trans %>% arrange(`Biosample type`,`Biosample term name`,`Biosample term id`,`Experiment target`,`Experiment accession`,`Biological replicate(s)`,`Technical replicate(s)`)
  tab_count_ep_trans = read.csv(paste(path_data,dataType,dataSubType,'tab_count_merge_trans',sep='/'),sep='\t',header=T,check.names = F)
  tab_count_ep_trans[is.na(tab_count_ep_trans)]=0
  tab_count_ep_trans = tab_count_ep_trans %>% filter(!grepl('PAR_Y',transcript_id))
  tab_count_ep_trans['transcript_id_title'] = id2idTitle(tab_count_ep_trans$transcript_id)
  vec_s_ep_trans = grep('ENCFF',colnames(tab_count_ep_trans),value=T)
  
  # load exp, experiment, gene
  tab_meta_ep_gene = read.csv(paste(path_data,dataType,dataSubType,'metadata.tsv',sep='/'),check.names = F,sep='\t')
  tab_meta_ep_gene = tab_meta_ep_gene %>% filter(`Genome annotation`=='V29') %>% 
    filter(`Output type` == 'gene quantifications') # RSEM quantification
  tab_meta_ep_gene = tab_meta_ep_gene %>% arrange(`Biosample type`,`Biosample term name`,`Biosample term id`,`Experiment target`)
  tab_count_ep_gene = read.csv(paste(path_data,dataType,dataSubType,'tab_count_merge_gene',sep='/'),sep='\t',header=T,check.names = F)
  tab_count_ep_gene[is.na(tab_count_ep_gene)]=0
  tab_count_ep_gene = tab_count_ep_gene %>% filter(!grepl('PAR_Y',gene_id))
  vec_s_ep_gene = grep('ENCFF',colnames(tab_count_ep_gene),value=T)
  
  # load exp, control, trans
  tab_meta_cl_trans = read.csv(paste(path_data,dataType,'control','metadata.tsv',sep='/'),check.names = F,sep='\t')
  tab_meta_cl_trans = tab_meta_cl_trans %>% filter(`Genome annotation`=='V29') %>% 
    filter(`Output type` == 'transcript quantifications')%>% filter(Size>30000000) # RSEM quantification
  tab_meta_cl_trans = tab_meta_cl_trans %>% arrange(`Biosample type`,`Biosample term name`,`Biosample term id`,`Experiment target`,`Experiment accession`,`Biological replicate(s)`,`Technical replicate(s)`)
  tab_count_cl_trans = read.csv(paste(path_data,dataType,'control','tab_count_merge_trans',sep='/'),sep='\t',header=T,check.names = F)
  tab_count_cl_trans[is.na(tab_count_cl_trans)]=0
  tab_count_cl_trans = tab_count_cl_trans %>% filter(!grepl('PAR_Y',transcript_id))
  tab_count_cl_trans['transcript_id_title'] = id2idTitle(tab_count_cl_trans$transcript_id)
  vec_s_cl_trans = grep('ENCFF',colnames(tab_count_cl_trans),value=T)
  
  # load exp, control, gene
  tab_meta_cl_gene = read.csv(paste(path_data,dataType,'control','metadata.tsv',sep='/'),check.names = F,sep='\t')
  tab_meta_cl_gene = tab_meta_cl_gene %>% filter(`Genome annotation`=='V29') %>% 
    filter(`Output type` == 'gene quantifications') # RSEM quantification
  tab_meta_cl_gene = tab_meta_cl_gene %>% arrange(`Biosample type`,`Biosample term name`,`Biosample term id`,`Experiment target`,`Experiment accession`,`Biological replicate(s)`,`Technical replicate(s)`)
  tab_count_cl_gene = read.csv(paste(path_data,dataType,'control','tab_count_merge_gene',sep='/'),sep='\t',header=T,check.names = F)
  tab_count_cl_gene[is.na(tab_count_cl_gene)]=0
  tab_count_cl_gene = tab_count_cl_gene %>% filter(!grepl('PAR_Y',gene_id))
  vec_s_cl_gene = grep('ENCFF',colnames(tab_count_cl_gene),value=T)
  
  # merge expression, trans 
  vec_columns_saved_trans = c('transcript_id','transcript_name','gene_id','gene_name')
  tab_count_trans = merge(tab_count_cl_trans[,c(vec_columns_saved_trans, vec_s_cl_trans)],
                          tab_count_ep_trans[,c(vec_columns_saved_trans, vec_s_ep_trans)],
                          by=vec_columns_saved_trans, all=T)
  tab_count_trans[is.na(tab_count_trans)] = 0
  tab_count_trans = tab_count_trans[rowSums(tab_count_trans[,c(vec_s_cl_trans,vec_s_ep_trans)]>5)>=2,]
  rownames(tab_count_trans) = tab_count_trans$transcript_id
  
  # merge expression, gene 
  vec_columns_saved_gene = c('gene_id','gene_name')
  tab_count_gene = merge(tab_count_cl_gene[,c(vec_columns_saved_gene, vec_s_cl_gene)],
                         tab_count_ep_gene[,c(vec_columns_saved_gene, vec_s_ep_gene)],
                         by=vec_columns_saved_gene, all=T)
  tab_count_gene[is.na(tab_count_gene)] = 0
  tab_count_gene = tab_count_gene[rowSums(tab_count_gene[,c(vec_s_cl_gene,vec_s_ep_gene)]>5)>=2,]
  rownames(tab_count_gene) = tab_count_gene$gene_id
  
  # norm, trans
  source('~/utils/public/count_2_FPKM_RPKM_TPM.R') 
  path_efflen = '~/data/Ref/GENECODE_r29_ENCODE'
  tab_tpm_trans = exp_count_2_exp_tpm_path(tab_count_trans[,c(vec_s_cl_trans,vec_s_ep_trans)],'transcript',path_efflen)
  tab_tpm_trans = merge(tab_count_trans[,vec_columns_saved_trans],tab_tpm_trans,by='transcript_id')
  rownames(tab_tpm_trans) = tab_tpm_trans[,c('transcript_id')]
  tab_log_tpm_trans = log2(tab_tpm_trans[,c(vec_s_cl_trans,vec_s_ep_trans)]+1)
  tab_log_tpm_trans['transcript_id'] = rownames(tab_log_tpm_trans)
  tab_log_tpm_trans = merge(tab_count_trans[,vec_columns_saved_trans],tab_log_tpm_trans,by='transcript_id')
  rownames(tab_log_tpm_trans) = tab_log_tpm_trans[,c('transcript_id')]
  vec_transcript_id = rownames(tab_log_tpm_trans)
  
  # norm gene
  tab_tpm_gene = exp_count_2_exp_tpm_path(tab_count_gene[,c(vec_s_cl_gene,vec_s_ep_gene)],'gene',path_efflen)
  tab_tpm_gene = merge(tab_count_gene[,vec_columns_saved_gene],tab_tpm_gene,by='gene_id')
  rownames(tab_tpm_gene) = tab_tpm_gene[,c('gene_id')]
  tab_log_tpm_gene = log2(tab_tpm_gene[,c(vec_s_cl_gene,vec_s_ep_gene)]+1)
  tab_log_tpm_gene['gene_id'] = rownames(tab_log_tpm_gene)
  tab_log_tpm_gene = merge(tab_count_gene[,vec_columns_saved_gene],tab_log_tpm_gene,by='gene_id')
  rownames(tab_log_tpm_gene) = tab_log_tpm_gene[,c('gene_id')]
  
  for(term in c('HepG2', 'K562')){
    tab_meta_ep_trans_term = tab_meta_ep_trans %>% filter(`Biosample term name`==term)
    vec_RBP_ep_trans = unlist(lapply(tab_meta_ep_trans_term[['Experiment target']], function(x) unlist(strsplit(x,'-'))[1]))
    
    tab_meta_ep_gene_term = tab_meta_ep_gene %>% filter(`Biosample term name`==term)
    vec_RBP_ep_gene = unlist(lapply(tab_meta_ep_gene_term[['Experiment target']], function(x) unlist(strsplit(x,'-'))[1]))
    
    # load ENCORE meta data
    tab_meta_ENCORE = readxl::read_excel(paste0('~/data/RBP_RNA/ENCORE/',term,'/metadata_mergedbed.xlsx'))
    vec_RBP_ENCORE = unlist(lapply(tab_meta_ENCORE[['Experiment target']], function(x) unlist(strsplit(x,'-'))[1]))
    
    # intersect among vec_RBP_RBPinfo, vec_RBP_ENCORE and vec_RBP_ep
    vec_RBP = intersect(intersect(intersect(vec_RBP_RBPinfo,vec_RBP_ENCORE),vec_RBP_ep_trans),vec_RBP_ep_gene)
    
    # RBP_name = vec_RBP[2]
    for(RBP_name in vec_RBP){
      dataSubTypeTermRBP = paste(dataSubType, term, RBP_name, sep='_')
      vec_dataSubTypeTermRBP = c(vec_dataSubTypeTermRBP,dataSubTypeTermRBP)
      
      vec_s_ep_trans_term_RBP_ori = tab_meta_ep_trans %>% 
        filter(`Biosample term name`==term & `Experiment target`==paste0(RBP_name,'-human')) %>% 
        pull(`File accession`)
      vec_s_cl_trans_term_RBP_ori = tab_meta_cl_trans %>% 
        filter(`Biosample term name`==term) %>% 
        pull(`File accession`)
      vec_s_ep_gene_term_RBP_ori = tab_meta_ep_gene %>% 
        filter(`Biosample term name`==term & `Experiment target`==paste0(RBP_name,'-human')) %>% 
        pull(`File accession`)
      vec_s_cl_gene_term_RBP_ori = tab_meta_cl_gene %>% 
        filter(`Biosample term name`==term) %>% 
        pull(`File accession`)
      vec_s_trans_term_RBP_ori = c(vec_s_cl_trans_term_RBP_ori, vec_s_ep_trans_term_RBP_ori)
      vec_s_gene_term_RBP_ori = c(vec_s_cl_gene_term_RBP_ori, vec_s_ep_gene_term_RBP_ori)
      
      
      vec_s_ep_trans_term_RBP = vec_s_ep_trans_term_RBP_ori
      vec_s_cl_trans_term_RBP = vec_s_cl_trans_term_RBP_ori
      vec_s_ep_gene_term_RBP = vec_s_ep_gene_term_RBP_ori
      vec_s_cl_gene_term_RBP = vec_s_cl_gene_term_RBP_ori
      
      
      tab_log_tpm_trans_term_RBP = tab_log_tpm_trans[,c(vec_columns_saved_trans,vec_s_cl_trans_term_RBP,vec_s_ep_trans_term_RBP)]
      tab_log_tpm_trans_term_RBP = tab_log_tpm_trans_term_RBP[rowSums(tab_log_tpm_trans_term_RBP>2)>=2,]
      
      tab_log_tpm_gene_term_RBP = tab_log_tpm_gene[,c(vec_columns_saved_gene,vec_s_cl_gene_term_RBP,vec_s_ep_gene_term_RBP)]
      tab_log_tpm_gene_term_RBP = tab_log_tpm_gene_term_RBP[rowSums(tab_log_tpm_gene_term_RBP>2)>=2,]
      
      list_vec_trans_tpm[[dataSubTypeTermRBP]] = rownames(tab_log_tpm_trans_term_RBP)
      
      # logFC
      ## trans
      tab_logFC_trans_term_RBP <- tab_log_tpm_trans_term_RBP %>% 
        select(any_of(vec_columns_saved_trans)) %>%
        mutate(exp_trans_cl = rowMeans(tab_log_tpm_trans_term_RBP[, vec_s_cl_trans_term_RBP]),
               exp_trans_ep = rowMeans(tab_log_tpm_trans_term_RBP[, vec_s_ep_trans_term_RBP])) %>%
        group_by(gene_id) %>%
        mutate(
          rate_trans_cl = exp_trans_cl / sum(exp_trans_cl),
          rate_trans_ep = exp_trans_ep / sum(exp_trans_ep)
        ) %>%
        replace(is.na(.), 0) %>%
        mutate(trans_abs_logFC = abs(log2((rate_trans_ep + 0.01) / (rate_trans_cl + 0.01))),
               trans_abs_diffRate = abs(rate_trans_ep - rate_trans_cl)) %>% 
        as.data.frame() %>% 
        select(transcript_id, trans_abs_logFC, trans_abs_diffRate) %>% 
        arrange(desc(trans_abs_logFC))
      ## gene
      tab_logFC_gene_term_RBP = tab_log_tpm_gene_term_RBP %>% 
        select(any_of(vec_columns_saved_gene)) %>% 
        mutate(exp_gene_cl = rowMeans(tab_log_tpm_gene_term_RBP[, vec_s_cl_gene_term_RBP]),
               exp_gene_ep = rowMeans(tab_log_tpm_gene_term_RBP[, vec_s_ep_gene_term_RBP]),
               RBP_abs_logFC = abs(log2(((exp_gene_ep + 1) / (exp_gene_cl + 1))))) %>%
        arrange(desc(RBP_abs_logFC))
      
      # Merge interactive data
      # the sample column in tab_inter represents different cell lines (terms) and needs to be distinguished!!
      tab_inter_term <- merge(tab_inter %>% filter(sample==!!term), tab_logFC_gene_term_RBP, by.x = 'RBP_name', by.y = 'gene_name') %>%
        filter(target_trans_id %in% vec_transcript_id) %>%
        group_by(RBP_name) %>%
        mutate(num_trans_perRBP = n()) %>%
        merge(tab_logFC_trans_term_RBP, by.x = 'target_trans_id', by.y = 'transcript_id')
      
      # list_regulated_trans_id
      list_vec_regulated_trans[[dataSubTypeTermRBP]] = 
        tab_inter_term %>% filter(RBP_name == !!RBP_name) %>% pull(target_trans_id) %>% unique()
      
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
               trans_max1_abs_logFC = trans_abs_logFC/max(trans_abs_logFC),
               geoMean_max1_abs_logFC = sqrt(RBP_max1_abs_logFC * trans_max1_abs_logFC)) %>% 
        arrange(desc(geoMean_max1_abs_logFC)) %>% filter(!is.na(geoMean_max1_abs_logFC))
      
      
      tab_inter_term_sim <- tab_inter_term %>% 
        select(target_trans_id,geoMean_max1_abs_logFC) %>% 
        group_by(target_trans_id) %>%
        slice_max(geoMean_max1_abs_logFC, n = 1,with_ties = F) %>%
        ungroup() %>% 
        arrange(desc(geoMean_max1_abs_logFC))
      
      # list_ranks_FCFC
      ranks_FCFC = tab_inter_term_sim  %>% pull(geoMean_max1_abs_logFC)
      names(ranks_FCFC) = tab_inter_term_sim  %>% pull(target_trans_id)
      list_ranks_FCFC[[dataSubTypeTermRBP]] = ranks_FCFC
      
    }
  }
}


# save support file: vec_trans_tpm
saveRDS(list_vec_trans_tpm, file.path(paste(path_support,'vec_trans_tpm',sep='/'), 
                                      paste(dataType,".rds",sep='')))
# save support file: vec_bioTypeTerm1Term2, list_ranks_FCFC
saveRDS(vec_dataSubTypeTermRBP, file.path(paste(path_support,'vec_bioTypeTerm1Term2',sep='/'), 
                                          paste(dataType,".rds",sep='')))
saveRDS(list_ranks_FCFC, file.path(paste(path_support,'list_ranks_FCFC',sep='/'), 
                                   paste(dataType,".rds",sep='')))
# save support file: regulated isoform set
saveRDS(list_vec_regulated_trans, file.path(paste(path_support,'list_vec_regulated_trans',sep='/'), 
                                            paste(dataType,".rds",sep='')))


# get support file: common, overlapped
get_common_trans <-function(combin){
  vecSplit_dataSubTypeTermRBP <- unlist(str_split(combin, '_'))
  dataSubType <- vecSplit_dataSubTypeTermRBP[1]
  term <- vecSplit_dataSubTypeTermRBP[2]
  RBP_name <- vecSplit_dataSubTypeTermRBP[3]
  vec_m_DTU = c('DRIMSeq', 'DEXSeq','limma', 'edgeR','satuRn',
                'DTUrtle','NBSplice','RATs','SPIT','SUPPA2')
  # create common isoforms
  list_tab_DTU = list()
  for(m_DTU in vec_m_DTU){
    tab_DTU_ori = read.table(paste(path_DTU, paste('DTU',dataType,sep='_'), dataSubType, paste('tab_DTU',term, paste(RBP_name,'human',sep='-'),m_DTU,sep='_'),sep='/'),sep='\t',header=T)
    if(m_DTU=='DRIMSeq'){
      tab_DTU = tab_DTU_ori %>% filter(grepl('^ENS',feature_id)&(!grepl('-[0-9]$',feature_id)))
      tab_DTU = tab_DTU %>% select(feature_id,prop_log2fc,pvalue,adj_pvalue)
    }
    if(m_DTU=='DEXSeq'){
      tab_DTU = tab_DTU_ori %>% filter(grepl('^ENS',featureID)&(!grepl('-[0-9]$',featureID)))
      tab_DTU = tab_DTU %>% select(featureID,grep('log2fold',colnames(tab_DTU),value=T),pvalue,padj)
    }
    if(m_DTU=='limma'){
      tab_DTU_ori['feature_id'] = rownames(tab_DTU_ori)
      tab_DTU = tab_DTU_ori %>% filter(grepl('^ENS',feature_id)&(!grepl('-[0-9]$',feature_id)))
      tab_DTU = tab_DTU %>% select(feature_id,logFC,P.Value,FDR)
    }
    if(m_DTU=='edgeR'){
      tab_DTU = tab_DTU_ori %>% filter(grepl('^ENS',feature_id)&(!grepl('-[0-9]$',feature_id)))
      tab_DTU = tab_DTU %>% select(feature_id,logFC,P.Value,FDR)
    }
    if(m_DTU=='satuRn'){
      tab_DTU_ori['feature_id'] = rownames(tab_DTU_ori)
      tab_DTU = tab_DTU_ori %>% filter(grepl('^ENS',feature_id)&(!grepl('-[0-9]$',feature_id)))
      tab_DTU = tab_DTU %>% select(feature_id,estimates,empirical_pval,empirical_FDR)
    }
    if(m_DTU=='DTUrtle'){
      tab_DTU = tab_DTU_ori %>% filter(grepl('^ENS',feature_id)&(!grepl('-[0-9]$',feature_id)))
      tab_DTU = tab_DTU %>% select(feature_id,lr,pvalue,adj_pvalue)
    }
    if(m_DTU=='NBSplice'){
      tab_DTU = tab_DTU_ori %>% filter(grepl('^ENS',iso)&(!grepl('-[0-9]$',iso)))
      tab_DTU = tab_DTU %>% select(iso,odd,pval,FDR)
    }
    if(m_DTU=='RATs'){
      tab_DTU = tab_DTU_ori %>% filter(grepl('^ENS',target_id)&(!grepl('-[0-9]$',target_id)))
      tab_DTU = tab_DTU %>% select(target_id,log2FC,pval,pval_corr)
    }
    if(m_DTU=='SPIT'){
      tab_DTU = tab_DTU_ori %>% filter(grepl('^ENS',transcript_id)&(!grepl('-[0-9]$',transcript_id)))
      tab_DTU = tab_DTU %>% select(transcript_id,likelihood,pvalue)
      tab_DTU['FDR'] = tab_DTU['pvalue']
    }
    if(m_DTU=='SUPPA2'){
      tab_DTU = tab_DTU_ori %>% filter(grepl('^ENS',transcript)&(!grepl('-[0-9]$',transcript)))
      tab_DTU = tab_DTU %>% select(transcript,dPSI,p_value)
      tab_DTU['FDR'] = tab_DTU['p_value']
    }
    colnames(tab_DTU) = c('transcript_id','logFC','pvalue','FDR')
    tab_DTU['transcript_id_title'] = id2idTitle(tab_DTU$transcript_id)
    
    # NA
    # tab_DTU = tab_DTU %>% filter(!is.na(pvalue))
    tab_DTU <- tab_DTU %>% mutate(pvalue = ifelse(is.na(pvalue), 1, pvalue),
                                  FDR = ifelse(is.na(FDR), 1, FDR))
    
    tab_DTU = tab_DTU[,c(grep('transcript',colnames(tab_DTU),value=TRUE),grep('gene',colnames(tab_DTU),value=TRUE),colnames(tab_DTU)[!grepl('transcript|gene',colnames(tab_DTU))])]
    tab_DTU['-log10pvalue'] = -log10(tab_DTU$pvalue)
    
    vec_pvalue = tab_DTU %>% pull(pvalue)
    vec_pvalue_unique = sort(unique(vec_pvalue))
    if(vec_pvalue_unique[1]==0){vec_pvalue[vec_pvalue==0] = vec_pvalue_unique[2]/2}
    if(vec_pvalue_unique[length(vec_pvalue_unique)]==1){vec_pvalue[vec_pvalue==1] = (vec_pvalue_unique[length(vec_pvalue_unique)-1] + 1)/2}
    tab_DTU['pvalue_noInf'] = vec_pvalue
    tab_DTU['-log10pvalue_noInf'] = -log10(tab_DTU$pvalue_noInf)
    tab_DTU['-lnpvalue_noInf'] = -log(tab_DTU$pvalue_noInf)
    
    list_tab_DTU[m_DTU] = list(tab_DTU)
  }
  tab_num_DTU_sig <- as.data.frame(table(unlist(lapply(list_tab_DTU, function(x) { x %>% filter(pvalue < 0.05) %>% pull(transcript_id) }))))
  vec_trans_common_DTU <- tab_num_DTU_sig %>% filter(Freq >= (length(vec_m_DTU)/2)) %>% pull(1)
  tab_num_DTU <- as.data.frame(table(unlist(lapply(list_tab_DTU, function(x) { x %>% pull(transcript_id) }))))
  vec_trans_all_DTU <- tab_num_DTU %>% pull(1)
  return(vec_trans_common_DTU)
}

# get support file: common isoform set
list_vec_trans_common_DTU = mclapply(vec_dataSubTypeTermRBP, function(combin) {
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
                                             paste(dataType,".rds",sep='')))

# get support file: overlapped isoform set
list_vec_overlapped_trans = list()
for(combin in vec_dataSubTypeTermRBP){
  list_vec_overlapped_trans[[combin]] = intersect(list_vec_trans_common_DTU[[combin]],list_vec_regulated_trans[[combin]]) 
}
saveRDS(list_vec_overlapped_trans, file.path(paste(path_support,'list_vec_overlapped_trans',sep='/'), 
                                             paste(dataType,".rds",sep='')))
