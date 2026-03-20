# Longitudinal single-cell RNA-seq analysis of immune responses

This repository contains code for the integrative analysis of longitudinal single-cell transcriptomic data to investigate dynamic immune responses over time.

The workflow combines differential expression analysis, gene module identification, functional enrichment and pathway activity quantification to uncover regulatory programs underlying immune heterogeneity.

---

## Overview

Understanding how immune responses evolve over time requires integrating information across multiple analytical layers.  
In this project, I implement a structured computational pipeline to:

- identify cell-type-specific transcriptional changes across time points  
- detect gene modules capturing coordinated regulatory programs  
- interpret these modules using functional enrichment analysis  
- quantify pathway activity dynamics using AUC-based scoring  

---

## Analysis Workflow

The analysis is organised into four main steps:

### 1. Longitudinal differential expression screening
Paired differential expression analysis across time points using **MAST**, performed across multiple cell types to identify the most responsive populations.

### 2. Monocyte module analysis
Focused downstream analysis of monocytes, including:

- construction of effect-size matrices across time points  
- gene module identification using hierarchical clustering  
- heatmap visualization of temporal expression patterns  

### 3. Functional enrichment analysis
Gene modules are functionally interpreted using:

- Gene Ontology Biological Process (**GO_BP**)  
- Reactome pathway analysis  

Both full and significant enrichment results are retained for downstream interpretation.

### 4. Pathway activity dynamics (UCell-based scoring)
Representative enriched pathways are selected based on adjusted p-values.  
Pathway activity is quantified at the single-cell level using **UCell**, and summarized across time points to reveal temporal dynamics.

---

## Notes

- Data are not included due to privacy and data protection restrictions  
- Scripts are provided as representative examples of the computational workflow  
- Gene sets used for pathway scoring are derived from enrichment results
