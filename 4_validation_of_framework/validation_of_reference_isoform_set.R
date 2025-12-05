library(tidyverse)
library(bedtoolsr)
options(bedtools.path = "~/anaconda3/envs/py_fastq/bin")

id2idTitle <- function(vec){
  return(unlist(lapply(vec, function(x) unlist(str_split(x,'\\.'))[1])))
}

# load anno data
B.bed = read.table('~/data/Ref/GENECODE_r29_ENCODE/tab_anno_trans',sep='\t',header=T)
B.bed = B.bed %>% select(seqnames,start,end,
                         transcript_id,transcript_type,transcript_name,
                         gene_id,gene_type,gene_name) %>% arrange(seqnames,start,end)

# load RBPinfo
tab_RBPinfo = read.csv('~/data/RBP_RNA/CLIPdb/POSTAR3_CLIPdb_module_browse_RBP_sub_info.csv',check.names = F)
colnames(tab_RBPinfo) = c('RBPname','geneID','Domain','Location','Function')
tab_RBPinfo_splice = tab_RBPinfo[grepl('Splic',tab_RBPinfo$Function),]
vec_RBP_RBPinfo = tab_RBPinfo_splice[['RBPname']]


# load exp 
dataType= 'ENCODE_NGS' # ENCODE_NGS
path_exp = '~/data'

# from ~/utils/20240910_DTU_benchmark_based_RNA-seq/isoformGSEA_NGS_rightCL_identifyTargetTrans_gseaPara_multi.R
path_support = '~/result/support_file_shortRead'
list_vec_trans_common_DTU = readRDS(paste(path_support,'list_vec_trans_common_DTU',sep='/', paste(dataType,".rds",sep='')))
list_vec_trans_regulated = readRDS(paste(path_support,'list_vec_regulated_trans',sep='/', paste(dataType,".rds",sep='')))

# dataSubType = 'CRISPR'  # CRISPR, shRNA
tab_compare = data.frame()
tab_logFC = data.frame()
for(dataSubType in c('CRISPR','shRNA')){
  # load exp, experiment
  tab_report_ep = read.csv(paste(path_exp,dataType,dataSubType,'experiment_report.tsv',sep='/'),check.names = F,sep='\t',skip=1)
  tab_meta_ep_trans = read.csv(paste(path_exp,dataType,dataSubType,'metadata.tsv',sep='/'),check.names = F,sep='\t')
  tab_meta_ep_trans = tab_meta_ep_trans %>% filter(`Genome annotation`=='V29')%>% 
    filter(`Output type` == 'transcript quantifications')%>% filter(Size>30000000) # RSEM quantification
  tab_meta_ep_trans = tab_meta_ep_trans %>% arrange(`Biosample type`,`Biosample term name`,`Biosample term id`,`Experiment target`,`Experiment accession`,`Biological replicate(s)`,`Technical replicate(s)`)
  tab_count_ep_trans = read.csv(paste(path_exp,dataType,dataSubType,'tab_count_merge_trans',sep='/'),sep='\t',header=T,check.names = F)
  tab_count_ep_trans[is.na(tab_count_ep_trans)]=0
  tab_count_ep_trans = tab_count_ep_trans %>% filter(!grepl('PAR_Y',transcript_id))
  tab_count_ep_trans['transcript_id_title'] = id2idTitle(tab_count_ep_trans$transcript_id)
  vec_s_ep_trans = grep('ENCFF',colnames(tab_count_ep_trans),value=T)
  
  # load exp, control
  tab_report_cl = read.csv(paste(path_exp,dataType,paste('cl',dataSubType,sep='_'),'experiment_report.tsv',sep='/'),check.names = F,sep='\t',skip=1)
  tab_meta_cl_trans = read.csv(paste(path_exp,dataType,paste('cl',dataSubType,sep='_'),'metadata.tsv',sep='/'),check.names = F,sep='\t')
  tab_meta_cl_trans = tab_meta_cl_trans %>% filter(`Genome annotation`=='V29') %>% 
    filter(`Output type` == 'transcript quantifications')%>% filter(Size>30000000) # RSEM quantification
  tab_meta_cl_trans = tab_meta_cl_trans %>% arrange(`Biosample type`,`Biosample term name`,`Biosample term id`,`Experiment target`,`Experiment accession`,`Biological replicate(s)`,`Technical replicate(s)`)
  tab_count_cl_trans = read.csv(paste(path_exp,dataType,paste('cl',dataSubType,sep='_'),'tab_count_merge_trans',sep='/'),sep='\t',header=T,check.names = F)
  tab_count_cl_trans[is.na(tab_count_cl_trans)]=0
  tab_count_cl_trans = tab_count_cl_trans %>% filter(!grepl('PAR_Y',transcript_id))
  tab_count_cl_trans['transcript_id_title'] = id2idTitle(tab_count_cl_trans$transcript_id)
  vec_s_cl_trans = grep('ENCFF',colnames(tab_count_cl_trans),value=T)

  
  # merge expression
  vec_columns_saved = c('transcript_id','transcript_name','gene_id','gene_name')
  tab_count_trans = merge(tab_count_cl_trans[,c(vec_columns_saved, vec_s_cl_trans)],
                          tab_count_ep_trans[,c(vec_columns_saved, vec_s_ep_trans)],
                          by=vec_columns_saved, all=T)
  tab_count_trans[is.na(tab_count_trans)] = 0
  tab_count_trans = tab_count_trans[rowSums(tab_count_trans[,c(vec_s_cl_trans,vec_s_ep_trans)]>5)>=2,]
  rownames(tab_count_trans) = tab_count_trans$transcript_id
  
  # norm
  source('~/utils/public/count_2_FPKM_RPKM_TPM.R') 
  path_efflen = '~/data/Ref/GENECODE_r29_ENCODE'
  tab_tpm_trans = exp_count_2_exp_tpm_path(tab_count_trans[,c(vec_s_cl_trans,vec_s_ep_trans)],'transcript',path_efflen)
  tab_tpm_trans = merge(tab_count_trans[,vec_columns_saved],tab_tpm_trans,by='transcript_id')
  rownames(tab_tpm_trans) = tab_tpm_trans[,c('transcript_id')]
  # tab_log_tpm_trans = log2(tab_tpm_trans[,-c(1)]+1)
  # rownames(tab_log_tpm_trans) = tab_tpm_trans[,c(1)]
  
  # tab_compare = data.frame()
  # tab_logFC = data.frame()
  # term = 'HepG2'
  for(term in c('HepG2', 'K562')){
    tab_meta_ep_trans_term = tab_meta_ep_trans %>% filter(`Biosample term name`==term)
    vec_RBP_ep = unlist(lapply(tab_meta_ep_trans_term[['Experiment target']], function(x) unlist(strsplit(x,'-'))[1]))
    
    # load ENCORE meta data
    tab_meta_ENCORE = readxl::read_excel(paste0('~/data/RBP_RNA/ENCORE/',term,'/metadata_mergedbed.xlsx'))
    vec_RBP_ENCORE = unlist(lapply(tab_meta_ENCORE[['Experiment target']], function(x) unlist(strsplit(x,'-'))[1]))
    
    # intersect among vec_RBP_RBPinfo, vec_RBP_ENCORE and vec_RBP_ep
    vec_RBP = intersect(intersect(vec_RBP_RBPinfo,vec_RBP_ENCORE),vec_RBP_ep)
    
    # RBP_name = vec_RBP[1]
    for(RBP_name in vec_RBP){
      
      tab_meta_ep_trans_term_RBP = tab_meta_ep_trans %>% 
        filter(`Biosample term name`==term & `Experiment target`== paste0(RBP_name,'-human'))
      ## ep
      vec_s_ep_trans_term_RBP = tab_meta_ep_trans_term_RBP %>% pull(`File accession`)
      ## cl
      access_ep_trans = tab_meta_ep_trans_term_RBP$`Experiment accession`[1]
      series_ep_trans = str_extract(tab_report_ep %>% filter(Accession == access_ep_trans) %>% pull(`Related series`),
                                    "/gene-silencing-series/[^,]+")
      access_cl_trans = tab_report_cl %>% filter(grepl(series_ep_trans,`Related series`)) %>% pull(Accession)
      vec_s_cl_trans_term_RBP = tab_meta_cl_trans %>% filter(`Experiment accession` == access_cl_trans) %>% pull(`File accession`)
      
      
      tab_count_trans_term_RBP = tab_tpm_trans[,c(vec_columns_saved,vec_s_cl_trans_term_RBP,vec_s_ep_trans_term_RBP)]
      tab_count_trans_term_RBP = tab_count_trans_term_RBP[rowSums(tab_count_trans_term_RBP>5)>=2,]
      tab_tpm_trans_mean = cbind(tab_count_trans_term_RBP[,vec_columns_saved],
                                 data.frame(cl=rowSums(tab_count_trans_term_RBP[,vec_s_cl_trans_term_RBP]),
                                            ep=rowSums(tab_count_trans_term_RBP[,vec_s_ep_trans_term_RBP])))
      tab_tpm_trans_mean_group = tab_tpm_trans_mean %>% group_by(gene_id) %>% summarise(sum_cl=sum(cl),sum_ep=sum(ep))
      tab_tpm_trans_mean_merge = merge(tab_tpm_trans_mean,tab_tpm_trans_mean_group,by='gene_id',all.x=T)
      tab_tpm_trans_mean_merge['rate_cl'] = tab_tpm_trans_mean_merge['cl']/tab_tpm_trans_mean_merge['sum_cl']
      tab_tpm_trans_mean_merge['rate_ep'] = tab_tpm_trans_mean_merge['ep']/tab_tpm_trans_mean_merge['sum_ep']
      tab_tpm_trans_mean_merge[is.na(tab_tpm_trans_mean_merge)] = 0
      tab_tpm_trans_mean_merge['logFC'] = log2(tab_tpm_trans_mean_merge['rate_ep']/tab_tpm_trans_mean_merge['rate_cl'])
      tab_tpm_trans_mean_merge[is.na(tab_tpm_trans_mean_merge)] = 0
      tab_tpm_trans_mean_merge['abs_logFC'] = abs(tab_tpm_trans_mean_merge['logFC'])
      tab_tpm_trans_mean_merge['diff_rate'] = tab_tpm_trans_mean_merge['rate_ep'] - tab_tpm_trans_mean_merge['rate_cl']
      tab_tpm_trans_mean_merge['abs_diff_rate'] = abs(tab_tpm_trans_mean_merge['diff_rate'])
      
      
      ## create vec_target_trans_id
    
      # regulated_isoforms
      
      # load ENCORE for each RBP
      # name_file_ENCORE = tab_meta_ENCORE %>% filter(`Experiment target`==paste0(RBP_name,'-human')) %>% pull(`File accession`)
      # A.bed = read.table(paste0('~/data/RBP_RNA/ENCORE/',term,'/',name_file_ENCORE,'.bed.gz'),sep='\t',header=F)
      # A.bed = A.bed[,c('V1','V2','V3','V7','V8')] %>% arrange(V1,V2,V3) # V7: geomean of the log2 fold changes; V8: minimum of the -log10 p-value between two replicates
      # colnames(A.bed) = c('seqnames','start','end','geomean','pvalue')
      # # intersect between A and B
      # tab_inter = bedtoolsr::bt.intersect(a = A.bed, b = B.bed, wo=T,f=1) # f=1: RBP fully covers the sequence
      # colnames(tab_inter) = c('RNA_seqnames','RNA_start','RNA_end','RNA_geomean','RNA_pvalue',
      #                         'target_seqnames','target_start','target_end',
      #                         'target_trans_id','target_trans_type','target_trans_name',
      #                         'target_gene_id','target_gene_type','target_gene_name',
      #                         'num_overlap_seq')
      # vec_target_trans_id  = unique(tab_inter[['target_trans_id']])
      
      # # regulated_isoforms
      # vec_target_trans_id = list_vec_trans_regulated[[paste(dataSubType,term,RBP_name, sep='_')]]
      # # common_isoforms
      # vec_target_trans_id = list_vec_trans_common_DTU[[paste(dataSubType,term,RBP_name, sep='_')]]
      # common and regulated
      vec_target_trans_id = intersect(list_vec_trans_regulated[[paste(dataSubType,term,RBP_name, sep='_')]],
                                      list_vec_trans_common_DTU[[paste(dataSubType,term,RBP_name, sep='_')]])
      
      # compare between "tab_tpm_trans_mean_merge in vec_target_trans_id" and "not in"
      vec_logFC_within = tab_tpm_trans_mean_merge %>% filter(abs_logFC!=Inf & transcript_id %in% vec_target_trans_id) %>% pull(abs_logFC)
      vec_logFC_without = tab_tpm_trans_mean_merge %>% filter(abs_logFC!=Inf & !transcript_id %in% vec_target_trans_id) %>% pull(abs_logFC)
      
      
      
      tab_logFC_temp = data.frame(dataType = dataType,
                                  dataSubType = dataSubType,
                                  term=term,
                                  RBP_name=RBP_name,
                                  withinOrOut = c(rep('within',length(vec_logFC_within)),rep('without',length(vec_logFC_without))),
                                  logFC = c(vec_logFC_within,vec_logFC_without))
      if(dim(tab_logFC)[1]==0){tab_logFC = tab_logFC_temp}else{tab_logFC = rbind(tab_logFC,tab_logFC_temp)}
      
      result_ttest = t.test(vec_logFC_within, vec_logFC_without,alternative = "greater")
      result_ranksum = wilcox.test(vec_logFC_within, vec_logFC_without, alternative = "greater")
      tab_compare_temp = data.frame(dataType = dataType,
                                    dataSubType = dataSubType,
                                    term=term,
                                    RBP_name=RBP_name,
                                    mean_logFC_within = mean(vec_logFC_within),
                                    mean_logFC_without = mean(vec_logFC_without),
                                    ttest_stat = result_ttest[['statistic']],
                                    ttest_pvalue = result_ttest[['p.value']],
                                    ransum_stat = result_ranksum[['statistic']],
                                    ransum_pvalue = result_ranksum[['p.value']])
      if(dim(tab_compare)[1]==0){tab_compare = tab_compare_temp}else{tab_compare = rbind(tab_compare,tab_compare_temp)}
    }
  }
  # write.table(tab_compare,paste('~/result/result_compare_logFC_rightCL/tab_compare',dataSubType,sep='_'),sep='\t',row.names = FALSE,quote=FALSE)
}
tab_compare = merge(tab_compare,tab_RBPinfo_splice[,c('RBPname','Function')],by.x='RBP_name',by.y='RBPname',all.x=T)
tab_compare = tab_compare %>% arrange(dataSubType,term,RBP_name)

write.table(tab_compare,'~/result/result_compare_logFC_rightCL/tab_compare',sep='\t',row.names = FALSE,quote=FALSE)
write.table(tab_logFC,'~/result/result_compare_logFC_rightCL/tab_logFC',sep='\t',row.names = FALSE,quote=FALSE)


