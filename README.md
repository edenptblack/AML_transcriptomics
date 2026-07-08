# Differential expression in AML: t-test vs limma

Differential expression analysis of acute myeloid leukaemia, calling DEGs using 2 independent methods: a per-gene t-test and a limma linear model. Both methods use identical filtering and Benjamini-Hochburg FDR correction. The two results are then compared. Adapted from Teesside University MSc Bioinformatics coursework.

## Description

For each of ~22,000 probes, expression is compared between AML and healthy bone marrow/peripheral blood samples using 2 methods: 1/ a per-gene two-sample Welch's t-test, and 2/ limma's moderated linear model. Both DEG lists are produced through identical filter steps, checking for reliable expression, non-saturation, a log fold change > 1, and a BH-adjusted p-value < 0.05. Therefore, the only separation arrises from how each method estimates variance. The two sets are compared by their overlap.

## Results and limitations

The t-test calls 507 DEGs and limma calls 512, with 490 shared. There are 17 t-test only DEGs and 22 limma-only DEGs. Age was determined to have a significant relationship with disease status; however, as 14/64 patients in the dataset had missing age data, it was not included in the limma model as a confounding variable.

## Stack & methods

R - tidyverse - limma (Bioconductor) - VennDiagram

## Reproducability

tidyverse, VennDiagram, and BioConductor can be installed via CRAN; limma must be installed using BioConductor (BiocManager::install("limma)). 

## Data

Stirewalt DL et al. Identification of genes with abnormal expression changes in acute myeloid leukemia. Genes Chromosomes Cancer 2008; 47(1):8-20. GEO accession GSE9476.