# Longitudinal single-cell RNA-seq analysis of immune responses

This repository provides a structured computational workflow for the integrative analysis of longitudinal single-cell transcriptomic data to investigate dynamic immune responses and regulatory programs over time.

The pipeline combines differential expression analysis, gene module identification, functional enrichment and pathway activity quantification to uncover mechanisms underlying immune heterogeneity.

---

## Overview

Understanding immune dynamics requires integrating multiple analytical layers across time.  
This project implements a reproducible framework to:

- identify cell-type-specific transcriptional changes across time points  
- detect gene modules representing coordinated regulatory programs  
- interpret these modules using functional enrichment analysis  
- quantify pathway activity dynamics using single-cell AUC-based scoring  

---

## Analysis Workflow

The workflow is organised into four main steps:

### 1. Longitudinal differential expression analysis
- Paired differential expression using **MAST**
- Performed across multiple cell types
- Identifies dynamic transcriptional responses over time  

---

### 2. Monocyte-focused module analysis
- Construction of effect-size matrices across time points  
- Gene module identification using hierarchical clustering  
- Silhouette-based selection of optimal cluster number  
- Heatmap visualization of temporal expression patterns  

---

### 3. Functional enrichment analysis
Gene modules are interpreted using:

- Gene Ontology Biological Process (**GO_BP**)  
- Reactome pathway analysis  

Outputs include:
- full enrichment results  
- significance-filtered results  
- module-level pathway summaries  

---

### 4. Pathway activity dynamics (UCell-based)
- Selection of representative pathways based on adjusted p-values  
- Gene sets derived from enrichment results  
- Single-cell pathway scoring using **UCell**  
- Aggregation across time points to reveal temporal dynamics  

---

## Key Features

- Longitudinal design with paired statistical modeling  
- Module-based interpretation of dynamic gene programs  
- Integration of enrichment results into pathway activity scoring  
- Flexible and extensible workflow for multi-step analysis  

---

## Data and Usage Notes

- Data are not included due to privacy and data protection restrictions  
- Scripts are provided as representative examples of the analysis workflow  
- Gene sets used for pathway scoring are derived from enrichment results  
