# Survival_Analysis
This repo contians basic code for survival analysis using TCGA dataset

## R Packages 
TCGAbiolinks 
survival
survminer
DESeq2 
tidyverse
dplyr

### Basic functions 
GetGDCprojects() - to obtain a list of projects from TCGA
GDCquery_clinic() - to obtain clinical data for project of interest
Columns to check in the clinical data:
(i) vital_status   (ii) days_to_last_follow_up   (iii) days_to_death 
GDCquery() ,GDCdownload(), GDCprepare() to obtain data for related project
DESeqDataSetFromMatrix()
survfit()
ggsurvplot()




