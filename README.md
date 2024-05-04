# Survival_Analysis
This repo contians basic code for survival analysis using TCGA dataset

## R Packages 
1.TCGAbiolinks 
2.survival
3.survminer
4.DESeq2 
5.tidyverse
6.dplyr

### Basic functions 
1.GetGDCprojects() - to obtain a list of projects from TCGA
2.GDCquery_clinic() - to obtain clinical data for project of interest
3.Columns to check in the clinical data:
(i) vital_status   (ii) days_to_last_follow_up   (iii) days_to_death 
4.GDCquery() ,GDCdownload(), GDCprepare() to obtain data for related project
5.DESeqDataSetFromMatrix()
6.survfit()
7.ggsurvplot()




