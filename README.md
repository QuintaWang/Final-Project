# Final Project
# Title
## Differential Gene Expression in Different age Using StarCOUNT
# Author
## Jiaxin Wang
# Overview of Project
## I will identify differentially expressed genes between primary tumor for female age under 55 vs age over 55. This analysis will utilize the package DeSEQ2 (http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) and complete the entire vignette. For this analysis, I’ll utilize the TCGA-OV cohort, and have identified STAR files for tumors that fit within my cohort with 28 female under or at age of 55 and 36 female above age of 55. Within the analysis, I will control for race and ethnicity. 

## Vignette: http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

# Data
## I will use the data from https://portal.gdc.cancer.gov/repository. Examining clinical data, there are 238 primary tumor samples, I eliminated the outliar and 28 are defined by me as young patient and 36 are identified as senior patient. The specific files are available are here [link].
# Milestone 1
## I will pick the data from 238 clinical cases and make a new data sheet. And I will download the repository files individually and label them with numeric names. Based on the name I will merge total 63 cases into one text to have first column and unstranded data listed. 
# Milestone 2
## After all data loaded on R-Studio, I will make graph to show the result based on the method showen on the Deseq by first convert the deseq to starcount.
# Deliverable
## R MarkDown
