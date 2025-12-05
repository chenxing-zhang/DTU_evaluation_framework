library(tidyverse)
library(bedtoolsr)
options(bedtools.path = "~/anaconda3/envs/py_fastq/bin")
library(data.table)
library(parallel)
library(caret) # confusionMatrix , glm

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

num_cores <- max(detectCores() - 2, 1)

path_save = '~/result/predict_regulated_isoform'
dataType= 'ENCODE_NGS' # ENCODE_NGS, RNAWG_long-read
path_data = '~/data'
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

cl <- makeCluster(30, type = "FORK")
list_tab_inter <- parLapplyLB(cl, vec_seqnames,
                              function(seqName) inter(A.bed %>% filter(seqnames == seqName),
                                                      B.bed %>% filter(seqnames == seqName),
                                                      species,dataType), chunk.size = 1)
stopCluster(cl)
tab_inter <- bind_rows(list_tab_inter)
rm(list_tab_inter)
gc()

vec_RBP_inter <- unique(tab_inter$RBP_name)
tab_inter <- tab_inter %>% group_by(RBP_name) %>% mutate(num_trans_perRBP = n())

list_tab_inter_term_temp = list()
vec_columns_saved_trans = c('transcript_id','transcript_name','gene_id','gene_name')
vec_columns_saved_gene = c('gene_id','gene_name')
# dataSubType = 'CRISPR'
for(dataSubType in c('CRISPR','shRNA')){
  # load control experiment_report.tsv
  tab_report_ep = read.csv(paste(path_data,dataType,dataSubType,'experiment_report.tsv',sep='/'),check.names = F,sep='\t',skip=1)
  
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
  
  # load control experiment_report.tsv
  tab_report_cl = read.csv(paste(path_data,dataType,paste('cl',dataSubType,sep='_'),'experiment_report.tsv',sep='/'),check.names = F,sep='\t',skip=1)
  
  # load exp, control, trans
  tab_meta_cl_trans = read.csv(paste(path_data,dataType,paste('cl',dataSubType,sep='_'),'metadata.tsv',sep='/'),check.names = F,sep='\t')
  tab_meta_cl_trans = tab_meta_cl_trans %>% filter(`Genome annotation`=='V29') %>% 
    filter(`Output type` == 'transcript quantifications')%>% filter(Size>30000000) # RSEM quantification
  tab_meta_cl_trans = tab_meta_cl_trans %>% arrange(`Biosample type`,`Biosample term name`,`Biosample term id`,`Experiment target`,`Experiment accession`,`Biological replicate(s)`,`Technical replicate(s)`)
  tab_count_cl_trans = read.csv(paste(path_data,dataType,paste('cl',dataSubType,sep='_'),'tab_count_merge_trans',sep='/'),sep='\t',header=T,check.names = F)
  tab_count_cl_trans[is.na(tab_count_cl_trans)]=0
  tab_count_cl_trans = tab_count_cl_trans %>% filter(!grepl('PAR_Y',transcript_id))
  tab_count_cl_trans['transcript_id_title'] = id2idTitle(tab_count_cl_trans$transcript_id)
  vec_s_cl_trans = grep('ENCFF',colnames(tab_count_cl_trans),value=T)
  
  # load exp, control, gene
  tab_meta_cl_gene = read.csv(paste(path_data,dataType,paste('cl',dataSubType,sep='_'),'metadata.tsv',sep='/'),check.names = F,sep='\t')
  tab_meta_cl_gene = tab_meta_cl_gene %>% filter(`Genome annotation`=='V29') %>% 
    filter(`Output type` == 'gene quantifications') # RSEM quantification
  tab_meta_cl_gene = tab_meta_cl_gene %>% arrange(`Biosample type`,`Biosample term name`,`Biosample term id`,`Experiment target`,`Experiment accession`,`Biological replicate(s)`,`Technical replicate(s)`)
  tab_count_cl_gene = read.csv(paste(path_data,dataType,paste('cl',dataSubType,sep='_'),'tab_count_merge_gene',sep='/'),sep='\t',header=T,check.names = F)
  tab_count_cl_gene[is.na(tab_count_cl_gene)]=0
  tab_count_cl_gene = tab_count_cl_gene %>% filter(!grepl('PAR_Y',gene_id))
  vec_s_cl_gene = grep('ENCFF',colnames(tab_count_cl_gene),value=T)
  
  # merge expression, trans 
  tab_count_trans = merge(tab_count_cl_trans[,c(vec_columns_saved_trans, vec_s_cl_trans)],
                          tab_count_ep_trans[,c(vec_columns_saved_trans, vec_s_ep_trans)],
                          by=vec_columns_saved_trans, all=T)
  tab_count_trans[is.na(tab_count_trans)] = 0
  tab_count_trans = tab_count_trans[rowSums(tab_count_trans[,c(vec_s_cl_trans,vec_s_ep_trans)]>5)>=2,]
  rownames(tab_count_trans) = tab_count_trans$transcript_id
  
  # merge expression, gene 
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
  # tab_log_tpm_trans = log2(tab_tpm_trans[,c(vec_s_cl_trans,vec_s_ep_trans)]+1)
  tab_log_tpm_trans = tab_tpm_trans[,c(vec_s_cl_trans,vec_s_ep_trans)]
  tab_log_tpm_trans['transcript_id'] = rownames(tab_log_tpm_trans)
  tab_log_tpm_trans = merge(tab_count_trans[,vec_columns_saved_trans],tab_log_tpm_trans,by='transcript_id')
  rownames(tab_log_tpm_trans) = tab_log_tpm_trans[,c('transcript_id')]
  vec_transcript_id = rownames(tab_log_tpm_trans)
  
  # norm gene
  tab_tpm_gene = exp_count_2_exp_tpm_path(tab_count_gene[,c(vec_s_cl_gene,vec_s_ep_gene)],'gene',path_efflen)
  tab_tpm_gene = merge(tab_count_gene[,vec_columns_saved_gene],tab_tpm_gene,by='gene_id')
  rownames(tab_tpm_gene) = tab_tpm_gene[,c('gene_id')]
  # tab_log_tpm_gene = log2(tab_tpm_gene[,c(vec_s_cl_gene,vec_s_ep_gene)]+1)
  tab_log_tpm_gene = tab_tpm_gene[,c(vec_s_cl_gene,vec_s_ep_gene)]
  tab_log_tpm_gene['gene_id'] = rownames(tab_log_tpm_gene)
  tab_log_tpm_gene = merge(tab_count_gene[,vec_columns_saved_gene],tab_log_tpm_gene,by='gene_id')
  rownames(tab_log_tpm_gene) = tab_log_tpm_gene[,c('gene_id')]
  
  # term = 'HepG2'
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
    
    for(RBP_name in vec_RBP){
      # sample of trans 
      tab_meta_ep_trans_term_RBP = tab_meta_ep_trans %>% 
        filter(`Biosample term name`==term & `Experiment target`==paste0(RBP_name,'-human'))
      ## ep
      vec_s_ep_trans_term_RBP_ori = tab_meta_ep_trans_term_RBP %>% pull(`File accession`)
      ## cl
      access_ep_trans = tab_meta_ep_trans_term_RBP$`Experiment accession`[1]
      series_ep_trans = str_extract(tab_report_ep %>% filter(Accession == access_ep_trans) %>% pull(`Related series`),
                                    "/gene-silencing-series/[^,]+")
      access_cl_trans = tab_report_cl %>% filter(grepl(series_ep_trans,`Related series`)) %>% pull(Accession)
      vec_s_cl_trans_term_RBP_ori = tab_meta_cl_trans %>% filter(`Experiment accession` == access_cl_trans) %>% pull(`File accession`)
      
      # sample of gene 
      tab_meta_ep_gene_term_RBP = tab_meta_ep_gene %>% 
        filter(`Biosample term name`==term & `Experiment target`==paste0(RBP_name,'-human'))
      ## ep
      vec_s_ep_gene_term_RBP_ori = tab_meta_ep_gene_term_RBP %>% pull(`File accession`)
      ## cl
      access_ep_gene = tab_meta_ep_gene_term_RBP$`Experiment accession`[1]
      series_ep_gene = str_extract(tab_report_ep %>% filter(Accession == access_ep_gene) %>% pull(`Related series`),
                                   "/gene-silencing-series/[^,]+")
      access_cl_gene = tab_report_cl %>% filter(grepl(series_ep_gene,`Related series`)) %>% pull(Accession)
      vec_s_cl_gene_term_RBP_ori = tab_meta_cl_gene %>% filter(`Experiment accession` == access_cl_gene) %>% pull(`File accession`)
      
      
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
        merge(tab_logFC_trans_term_RBP, by.x = 'target_trans_id', by.y = 'transcript_id') %>%
        mutate(geoMean_abs_logFC = sqrt(RBP_abs_logFC * trans_abs_logFC))%>%
        as.data.frame() %>% 
        mutate(RBP_max1_abs_logFC = RBP_abs_logFC/max(RBP_abs_logFC),
               trans_max1_abs_logFC = trans_abs_logFC/max(trans_abs_logFC),
               geoMean_max1_abs_logFC = sqrt(RBP_max1_abs_logFC * trans_max1_abs_logFC))%>% 
        arrange(desc(geoMean_max1_abs_logFC))
      
      # calculate rate_RBP_perTrans and rate_trans_perRBP
      num_RBP = length(unique(tab_inter_term$RBP_name))
      tab_inter_term  = tab_inter_term %>% 
        group_by(target_trans_name) %>% 
        mutate(num_RBP_perTrans = n(),rate_RBP_perTrans = num_RBP_perTrans/num_RBP) %>% ungroup() 
      num_trans = length(unique(tab_inter_term$target_trans_name))
      tab_inter_term  = tab_inter_term %>% 
        group_by(RBP_name) %>% 
        mutate(num_trans_perRBP = n(), rate_trans_perRBP = num_trans_perRBP/num_trans) %>% ungroup()
      
      # mark target and untarget
      tab_inter_term = tab_inter_term %>% mutate(class=ifelse(RBP_name==!!RBP_name,'target','untarget'))
      
      # select useful columns
      # tab_inter_term_temp = tab_inter_term %>% 
      #   mutate(sampleRBPTrans = paste(sample, RBP_name, target_trans_name,sep='_')) %>% 
      #   select(sampleRBPTrans, RBP_max1_abs_logFC, trans_max1_abs_logFC, rate_RBP_perTrans, rate_trans_perRBP, class)
      tab_inter_term_temp = tab_inter_term %>% 
        mutate(sampleRBPTrans = paste(sample, RBP_name, target_trans_name,sep='_')) %>% 
        select(sampleRBPTrans, sample, RBP_name, target_trans_name, RBP_max1_abs_logFC, trans_max1_abs_logFC, rate_RBP_perTrans, rate_trans_perRBP, class) %>% 
        mutate(class = as.factor(class))
      
      write.table(tab_inter_term_temp,paste(path_save,paste('tab_inter_term',dataSubType, term, RBP_name,sep='_'),sep='/'), row.names = FALSE,quote=FALSE,sep='\t')
      list_tab_inter_term_temp[[paste(dataSubType, term, RBP_name,sep='_')]] = tab_inter_term_temp
      
    }
  }
}

tab_inter_term_all = bind_rows(list_tab_inter_term_temp)
write.table(tab_inter_term_all,paste(path_save,'tab_inter_term_all',sep='/'), row.names = FALSE,quote=FALSE,sep='\t')

################################################################################

# prediction
library(tidyverse)
library(caret)
library(pROC)
library(doParallel)
library(data.table) 
path_save = '~/result/predict_regulated_isoform'

# Data preprocessing 
tab_inter_term_all = fread(paste(path_save,'tab_inter_term_all',sep='/'),sep='\t',header=T)

# Ensure the class variable is a factor
tab_inter_term_all[, class := as.factor(class)]


# Split training/test sets (efficient sampling)
set.seed(123)
tab_inter_term_all[, test := {
  sampled <- sample(.N, size = floor(.N * 0.2))
  idx <- rep(FALSE, .N)
  idx[sampled] <- TRUE
  idx
}, by = c('sample','RBP_name')]


# Split into training and test data
train_data <- tab_inter_term_all[test == FALSE]
test_data  <- tab_inter_term_all[test == TRUE]
 
# Remove temporary column
train_data[, test := NULL]
test_data[, test := NULL]

# Configure parallel computing
cl <- makePSOCKcluster(5)
registerDoParallel(cl)

# Define training control parameters
ctrl <- trainControl(
  method = "cv", 
  number = 5,
  classProbs = TRUE,   # Needed for probability output
  summaryFunction = twoClassSummary, # For binary classification
  allowParallel = TRUE
)

# Feature selection
features <- c("RBP_max1_abs_logFC", "trans_max1_abs_logFC") # FCFC
tab_inter_term_all <- tab_inter_term_all[,.(RBP_max1_abs_logFC, trans_max1_abs_logFC, class)]

# Free up memory (useful for large datasets)
# rm(tab_inter_term_all); gc()

# Model training
model_formula <- as.formula(paste("class ~", paste(features, collapse = "+")))

# GLM
start_time <- Sys.time()
set.seed(123)
glm_model <- train(
  model_formula,
  data = train_data,
  method = "glm",
  family = "binomial",
  trControl = ctrl
)
end_time <- Sys.time()
print(paste('glm',as.numeric(end_time) - as.numeric(start_time)))
saveRDS(glm_model, file.path(path_save, "glm_model_allData_2feature.rds"))

# LDA
set.seed(123)
lda_model <- train(
  model_formula,
  data = train_data,
  method = "lda",
  trControl = ctrl,
  preProcess = c("center", "scale") 
)
saveRDS(lda_model, file.path(path_save, "lda_model_allData_2feature.rds"))


# Stop parallel workers
stopCluster(cl)
