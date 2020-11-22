# Developing-an-infection-state-predictive-model

This repository contains code supporting the study from https://www.medrxiv.org/content/10.1101/2020.07.28.20163329v2. Scripts and functionality included are (1) multiple-step batch correction, (2) gene feature selection with backwards elimination implemented in VarSelRF with Ranger, (3) gene feature selection with a genetically inspired search algorithm, implemented with GALGO. 

## Abstract

A fundamental problem for disease treatment is that while antibiotics are a powerful counter to bacteria, they are ineffective against viruses. Often, bacterial and viral infections are confused due to their similar symptoms and lack of rapid diagnostics. With many clinicians relying primarily on symptoms for diagnosis, overuse and misuse of modern antibiotics are rife, contributing to the growing pool of antibiotic resistance. To ensure a given individual receives optimal treatment given their disease state and to reduce over-prescription of antibiotics, the host response can be measured quickly to distinguish between the two states. To establish a predictive biomarker panel of disease state (viral/bacterial/no-infection) we conducted a meta-analysis of human blood infection studies using Machine Learning (ML). We focused on publicly available gene expression data from two widely used platforms, Affymetrix and Illumina microarrays as they represented a significant proportion of the available data. We were able to develop multi-class models with high accuracies with our best model predicting 93% of bacterial and 89% viral samples correctly. To compare the selected features in each of the different technologies, we reverse engineered the underlying molecular regulatory network and explored the neighbourhood of the selected features. This highlighted that although on the gene-level the models differed, they did contain genes from the same areas of the network. Specifically, this convergence was to pathways including the Type I interferon Signalling Pathway, Chemotaxis, Apoptotic Processes, and Inflammatory / Innate Response. 

If used please cite: Myall, A.C., Perkins, S., Rushton, D., Jonathan, D., Spencer, P., Jones, A.R. and Antczak, P., 2020. Identifying robust biomarkers of infection through an omics-based meta-analysis. medRxiv.

## Batch correction

We developed a two-phase batch correction pipeline, which included the combination of platforms and studies. Results of batch correction were investigated for impact on the underlying biology by comparing dataset differentially expressed genes before and after batch correction using a Fisher's exact t-test and a Hyper Geometric Test. For a significant result, the overlap of differentially expressed genes was significant, and we inferred batch correction had not removed biological variation. Conversely, for a non-significant differentially expressed gene overlap, batch correction had likely removed biological variation and damaged the data quality. We also validated that the batch correction removed the previously observed clustering by study, by visualising before and after principle component plots:

![alt text](https://raw.githubusercontent.com/ashm97/Developing-an-infection-state-predictive-model/main/images/preview_batch_pca.png) 

## Backwards elimination

A 60/20/20 training/test/evaluation data split is used in our Backwards elimination, with 60 used for model training, 20 used to select trained models, then a final 20 as a held-out subset for final evaluation and reporting, a standard technique in machine learning. For each dataset we ran multiple search procedures, using Out-of-bag (OOB) error as the minimisation criterion. Each run generated a single optimal model which minimised OOB. For each dataset, a gene selection frequency is compared between runs, as well as optimal model size to accuracy, and OOB by model size within each run.

![alt text](https://raw.githubusercontent.com/ashm97/Developing-an-infection-state-predictive-model/main/images/example_backward_elim.png) 

## Genetically inspired search algorithm (GALGO)

** package currently unsuported by R 3.6

Victor Trevino, Francesco Falciani, GALGO: an R package for multivariate variable selection using genetic algorithms, Bioinformatics, Volume 22, Issue 9, 1 May 2006, Pages 1154–1156, https://doi.org/10.1093/bioinformatics/btl074

