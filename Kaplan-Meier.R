#importing the required packages 
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TCGAbiolinks")

install.packages("survminer")
```
```{r}
install.packages("survival")
BiocManager::install("DESeq2")

#Loading the required packages 
library(tidyverse)
library(dplyr)
library(survival)
library(survminer)
library(TCGAbiolinks)
library(DESeq2)

#get the list of projects from TCGA
GDC_project <-getGDCprojects()
clinical_data <- GDCquery_clinic("TARGET-NBL")
any(colnames(clinical_data) %in% c("vital_status","days_to_last_follow_up","days_to_death"))
which(colnames(clinical_data)%in% c("vital_status","days_to_last_follow_up","days_to_death"))
selected_columns <- clinical_data[,c(25,35,41)]
head(selected_columns)

table(selected_columns$vital_status)

clinical_data$censored <- ifelse(clinical_data$vital_status== 'Dead', TRUE, FALSE)

table(clinical_data$censored)
clinical_data$overall_survival <- ifelse(clinical_data$vital_status=="Dead",clinical_data$days_to_death,clinical_data$days_to_last_follow_up)

query_NB <- GDCquery(project = "TARGET-NBL",
                     data.category = "Transcriptome Profiling", experimental.strategy = "RNA-Seq", workflow.type = "STAR - Counts", data.type = "Gene Expression Quantification",
                     sample.type = "Primary Tumor", access = "open")

output_NB <- getResults(query_NB)
head(output_NB)

Neuroblastoma <- output_NB$cases[1:10]

query_NB <- GDCquery(project = "TARGET-NBL",
                     data.category = "Transcriptome Profiling", experimental.strategy = "RNA-Seq", workflow.type = "STAR - Counts", data.type = "Gene Expression Quantification",
                     sample.type = "Primary Tumor", access = "open",
                     barcode = Neuroblastoma)

GDCdownload(query_NB, directory = "C:/Users/notay003/Downloads/")

tcga_NB <- GDCprepare(query_NB, summarizedExperiment = T, directory = "C:/Users/notay003/Downloads/")
tcga_mat <- assay(tcga_NB, 'unstranded')

gene_metadata <- as.data.frame(rowData(tcga_NB))
gene_metadata[1:5]
coldata <- as.data.frame(colData(tcga_NB))
coldata

dds <-DESeqDataSetFromMatrix(countData = tcga_mat,
                             colData = coldata,
                             design = ~ 1)

keep <- rowSums(counts(dds)) >=10
dds <- dds[keep,]

vsd <- vst(dds,blind = F)
NB_mat_vst <- assay(vsd)
NB_mat_vst[1:2,1:2]

NB_MYCN <- NB_mat_vst %>% as.data.frame()%>%
  rownames_to_column("gene_id")%>%
  gather(key = 'case', value='counts', -gene_id) %>% 
  left_join(.,gene_metadata, by='gene_id') %>% filter(gene_name == "MYCN") 

median_val <- median(NB_MYCN$counts)
NB_MYCN$exprs <- ifelse(NB_MYCN$counts >= median_val, "High","Low")

NB_MYCN$case_id <- gsub("-01.*","",NB_MYCN$case)
NB_MYCN <- merge(NB_MYCN, clinical_data_u, by.x="case_id", by.y= "submitter_id")

NB_MYCN

fit <- survfit(Surv(overall_survival,censored) ~ exprs,data=NB_MYCN)
fit
```
```{r,echo=T}
ggsurvplot(fit,NB_MYCN,pval = T,risk.table = T)









