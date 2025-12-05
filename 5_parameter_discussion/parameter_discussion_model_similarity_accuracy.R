library(tidyverse)
library(caret)
library(pROC)    # AUROC
library(PRROC)   # AUPRC
library(ranger)
library(reticulate)
use_condaenv("r-tensorflow", required = TRUE)

path_model = '~/result/predict_regulated_isoform'
path_save = '~/result/predict_regulated_isoform/tab_inter_term'
vec_fileName = list.files(path_save)
path_save_pred = '~/result/predict_regulated_isoform/tab_inter_term_GLM_LDA'

# load model
lda_model  = readRDS(file.path(path_model, "lda_model_allData_2feature.rds"))
glm_model  = readRDS(file.path(path_model, "glm_model_allData_2feature.rds"))


for(fileName in vec_fileName){
  tab_inter_term = read.table(paste(path_save,fileName,sep='/'),header = T,sep='\t')
  tab_inter_term['class'] = factor(tab_inter_term$class)
  test_data = tab_inter_term
  # predicte
  pred_lda = predict(lda_model, newdata = test_data, type = "prob")[, 2]
  pred_glm = predict(glm_model, newdata = test_data, type = "prob")[, 2] 
  
  tab_inter_term['predScore_lda'] = pred_lda
  tab_inter_term['predClass_lda'] = factor(ifelse(pred_lda > 0.5, levels(test_data$class)[2], levels(test_data$class)[1]), levels = levels(test_data$class))
  tab_inter_term['predScore_glm'] = pred_glm
  tab_inter_term['predClass_glm'] = factor(ifelse(pred_glm > 0.5, levels(test_data$class)[2], levels(test_data$class)[1]), levels = levels(test_data$class))
  
  write.table(tab_inter_term,paste(path_save_pred,fileName,sep='/'), row.names = FALSE,quote=FALSE,sep='\t')
}

# accuracy
list_accuracy = list()
for(fileName in vec_fileName){
  test_data = read.table(paste(path_save_pred,fileName,sep='/'), header=T,sep='\t')
  test_data['class'] = factor(test_data$class)
  list_accuracy[[fileName]] = list()
  for(m in c('lda','glm')){
    pred = test_data[[paste('predScore',m,sep='_')]]

    # AUROC
    roc_obj <- roc(response = test_data$class, predictor = pred)  
    auroc <- as.numeric(auc(roc_obj))  
    # AUPRC
    pr_obj <- pr.curve(scores.class0 = pred, weights.class0 = as.numeric(test_data$class)-1, curve = TRUE)  
    auprc <- pr_obj$auc.integral  
    
    list_accuracy[[fileName]][[m]] = list(roc_obj=roc_obj, auroc=auroc, 
                                     pr_obj=pr_obj, auprc=auprc)
  }
}

# list_accuracy to tab_accuracy
list_tab_accuracy = list()
for(fileName in vec_fileName){
  vec_splited_fileName = unlist(str_split(fileName,'_'))
  for(m in c('lda','glm')){
    list_tab_accuracy[[paste(fileName,m,sep='_')]] = 
      data.frame(dataSubType = paste(vec_splited_fileName[4],vec_splited_fileName[5],sep='_'),
                 term1 = vec_splited_fileName[5],
                 term2 = vec_splited_fileName[6],
                 m_identify = m,
                 m_accuracy = c('auroc','auprc'),
                 accuracy = c(list_accuracy[[fileName]][[m]][['auroc']], list_accuracy[[fileName]][[m]][['auprc']]))
  }
}

tab_accuracy = bind_rows(list_tab_accuracy)

