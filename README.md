# DTU Evaluation Framework

A unbiased and real-data-driven framework for benchmarking differential transcript usage methods.

![Framework Overview](framework.png)

---

## Table of Contents

- [Overview](#overview)
- [Key Features](#key-features)
- [Framework Architecture](#framework-architecture)
- [Project Structure](#project-structure)
- [Environment](#environment)
- [Datasets and Results](#datasets-and-results)
  - [Datasets](#datasets)
  - [Results](#results)
- [Code](#code)
- [Contact](#contact)
- [Reference](#reference)

---

## Overview

Differential Transcript Usage (DTU) describes changes in the relative abundance of transcript isoforms of the same gene between conditions. Although many DTU analysis tools exist, a unbiased and real-data-based evaluation framework is still lacking. This repository provides such a framework, enabling unbiased comparison of DTU methods using real RNA-seq datasets, modular pipelines, and DTU-specific evaluation metrics.

---

## Key Features

- **Real RNA-seq–based benchmarking**
- **Unbiased DTU-specific evaluation metric**
- **Modular pipeline** covering DTU method execution, reference isoform sets construction, isoform set enrichment score calculation, framework validation, parameter discussion and DTU method evaluation.
- Support for **multiple DTU methods** (e.g., DRIMSeq, DEXSeq, satuRn, RATs, etc.)
- **RBP activity prediction method**

---

## Framework Architecture

The overall workflow (illustrated in `framework.png`) consists of:

1. **Datasets Collection**  
   Collect short-read bulk RNA-seq datasets from RBP KO/KD experiments, as well as long-read bulk, long-read single-cell, and long-read spatial RNA-seq datasets without KO/KD.

2. **Running DTU Methods**  
   Execute 10 representative differential transcript usage methods.

3. **Reference Isoform Set Construction**  
   Construct the common isoform set, the regulated isoform set and the overlapped isoform set.

4. **Isoform Set Enrichment Score Calculation**  
   Design of a DTU-specific evaluation metric.

---

## Project Structure
```
├── 0_overlap_between_DTU/           # Overlap analysis across DTU methods (after run DTU)
├── 1_run_DTU/                       # Scripts to run each DTU method
├── 2_reference_isoform_set/         # Generate reference isoform sets
├── 3_ISES/                          # Caculate isoform set enrichment score (ISES)
├── 4_validation_of_framework/       # Validation of reference isoform set and ISES
├── 5_parameter_discussion/          # Parameter discussion of framework
├── 6_DTU_evaluation/                # Evaluation result
├── 7_RBP_activity_prediction/       # RBP activity prediction 
├── utils/                           # Utility functions
├── framework.png                    # Framework diagram
└── README.md
```

---

## Environment

- **R ≥ 4.4**
```r
# 1. Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c(
  "DRIMSeq", "DEXSeq", "edgeR", "satuRn", "NBSplice", "limma", 
  "SummarizedExperiment", "GenomicRanges", "rtracklayer", "BiocParallel","DESeq2"
))

# 2. CRAN packages 
install.packages(c(
  "bedtoolsr", "Seurat", "reticulate",
  "tidyverse", "dplyr", "tidyr", "readr", "stringr",
  "stringi", "purrr", "reshape2", "data.table", "parallel",
  "doParallel", "caret", "ranger", "wCorr", "conflicted", "utils"
))

# 3. github pacakges
if(!requireNamespace("remotes", quietly = TRUE)){
    install.packages("remotes")
}
remotes::install_github("TobiTekath/DTUrtle")                 # https://github.com/TobiTekath/DTUrtle

if(!requireNamespace("devtools", quietly = TRUE)){
    install.packages("devtools")
}
devtools::install_github("bartongroup/rats", ref="master")    # https://github.com/bartongroup/Rats

# 4. Evaluation and statistics packages
install.packages(c(
  "pROC",     # ROC, AUROC
  "PRROC",    # PR curves
  "metap"     # p-value combination
))

```
- **Python ≥ 3.7**
```bash
conda create -n python python=3.11 -y
conda activate python
pip install numpy pandas scipy matplotlib seaborn tqdm collections
pip install SUPPA==2.3     # https://github.com/comprna/SUPPA
pip install spit           # https://github.com/berilerdogdu/SPIT
```

---

## Datasets and Results

All datasets used in this evaluation framework have been uploaded to Zenodo (**DOI link: 10.5281/zenodo.17698058**).

- ### Datasets

Please download the following archive files in order:
| File Name | Description |
|:-----|:-----|
| `data_Ref.zip`    | Reference sequence of human  |
| `data_Ref_mouse.zip`    | Reference sequence of mouse   |
| `data_RBP_RNA.zip`    | RBP-RNA interaction   |
| `data_ENCODE_NGS.zip`    | short-read RNA-seq data (RBP KO/KD)  |
| `data_RNAWG_long-read.zip`    | long-read bulk RNA-seq data   |
| `data_single_cell.zip`    | long-read single-cell RNA-seq data   |
| `data_spatial.zip`    | long-read spatial RNA-seq data   |
| `data_case.zip`    | RNA-seq data, case  |

After downloading, **extract each file into a folder named** `data` in the project directory.

Once extracted, the directory structure should become:
```
~/data/Ref/
~/data/Ref_mouse/
~/data/RBP_RNA/
~/data/ENCODE_NGS/
~/data/RNAWG_long-read/
~/data/single_cell/
~/data/spatial/
~/data/case/
```
These directories contain all real datasets required for reproducing our benchmarking workflow.



- ### Results

To ensure that all input directories match the expected structure and to bypass some of the computationally intensive intermediate steps (e.g., predicting the regulated isoform set), we also provide **all results**.
Please download the following archive files in order:
| File Name | Description |
|:-----|:-----|
|`result_DTU_ENCODE_NGS.zip` | DTU results in short-read RNA-seq data (RBP KO/KD) |
|`result_DTU_RNAWG_long-read.zip` | DTU results in long-read bulk RNA-seq data |
|`result_DTU_single_cell.zip` | DTU results in long-read single-cell RNA-seq data |
|`result_DTU_spatial.zip` | DTU results in long-read spatial RNA-seq data |
|`result_ISES_ENCODE_NGS.zip` | isoform enrichment score in short-read RNA-seq data (RBP KO/KD) |
|`result_ISES_RNAWG_long-read.zip` | isoform enrichment score in long-read bulk RNA-seq data |
|`result_ISES_single_cell.zip` | isoform enrichment score in long-read single-cell RNA-seq data |
|`result_ISES_spatial.zip` | isoform enrichment score in long-read spatial RNA-seq data |
|`result_parameter_discussion_gseaPara.zip` | parameter discussion of ISES parameter |
|**`result_predict_regulated_isoform.zip`** | regulated isoform prediction model and results|
|`result_support_file_longRead.zip` | support files in long-read data|
|`result_support_file_shortRead.zip` | support files in short-read data|
|`result_RBP_activity_prediction.zip` | RBP activity prediction|

Because `result_predict_regulated_isoform.zip` is very large, it is split into **13 part files**:

`result_predict_regulated_isoform_part_aa` ~ `result_predict_regulated_isoform_part_am`

After downloading **all 13 parts**, merge them using the following command:
```bash
cat result_predict_regulated_isoform_part_* > result_predict_regulated_isoform.zip      # Linux / macOS
copy /b result_predict_regulated_isoform_part_* result_predict_regulated_isoform.zip    # Windows
```

After downloading and merging, **extract all archives into a folder named** `result`.

the directory structure should become:

```
~/result/DTU_ENCODE_NGS/
~/result/DTU_RNAWG_long-read/
~/result/DTU_single_cell/
~/result/DTU_spatial/
~/result/ISES_ENCODE_NGS/
~/result/ISES_RNAWG_long-read/
~/result/ISES_single_cell/
~/result/ISES_spatial/
~/result/parameter_discussion_gseaPara/
~/result/predict_regulated_isoform/
~/result/support_file_longRead/
~/result/support_file_shortRead/
~/result/RBP_activity_prediction/
```
These directories contain all processed outputs used in our evaluation pipeline.

---

## Code

**1. Run differential transcript usage methods**
```
Rscript 1_run_DTU/run_DTU_longRead_bulk.R
Rscript 1_run_DTU/run_DTU_longRead_singleCell.R
Rscript 1_run_DTU/run_DTU_longRead_spatial.R
Rscript 1_run_DTU/run_DTU_shortRead_bulk.R
```

**2. Identify reference isoform set**
```
Rscript 2_reference_isoform_set/predict_regulated_isoforms.R
Rscript 2_reference_isoform_set/reference_isoform_set_longRead.R
Rscript 2_reference_isoform_set/reference_isoform_set_shortRead.R
```

**3. Caculate isoform enrichment score**
```
Rscript 3_ISES/ISES_longRead.R
Rscript 3_ISES/ISES_shortRead.R
```

**4. Validate framework (optional)**
```
Rscript 4_validation_of_framework/validation_of_reference_isoform_set.R
Rscript 4_validation_of_framework/validation_of_ISES_longRead.R
Rscript 4_validation_of_framework/validation_of_ISES_shortRead.R
```

**5. Discuss parameters (optional)**
```
Rscript 5_parameter_discussion/parameter_discussion_model_similarity_accuracy.R
Rscript 5_parameter_discussion/parameter_discussion_model_similarity_ISES.R
Rscript 5_parameter_discussion/parameter_discussion_gseaPara_longRead.R
Rscript 5_parameter_discussion/parameter_discussion_gseaPara_shortRead.R
```

**6. Evaluate differential transcript usage method**
```
Rscript 6_DTU_evaluation/DTU_evaluation_ISES_and_rank.R
```

**7. Predict RBP acticity (optional)**
```
Rscript 7_RBP_activity_prediciton/1_run_DTU_case.R
Rscript 7_RBP_activity_prediciton/2_reference_isoform_set_case.R
Rscript 7_RBP_activity_prediciton/3_ISES_case.R
Rscript 7_RBP_activity_prediciton/3_ISES_shortRead.R
Rscript 7_RBP_activity_prediciton/4_DE_case.R
Rscript 7_RBP_activity_prediciton/4_DE_shortRead.R
```

---


## Contact
- Author: Chenxing Zhang
- GitHub: https://github.com/chenxing-zhang/DTU_evaluation_framework
- Email: chenxing_zhang813@163.com

---

## Reference
Chenxing Zhang, Jun Liu, Qi Zhao, Huilong Yin, Angang Yang, Minhua Zheng*, Rui Zhang*. An Unbiased and Real-Data-Driven Framework for Benchmarking Differential Transcript Usage Methods. (under review)
