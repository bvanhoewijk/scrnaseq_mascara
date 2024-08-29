# Minor project MSc: Elucidating the suberin pathway of Tomato using single-cell RNA-seq and MASCARA

## Abstract
Bulk RNA-seq allows for investigating the gene expression of samples between wild-type and test conditions and allows for the exploration of co-expressed genes. However, with this method, it is not possible to compensate for cell type and cell-cycle heterogeneity. However, with single-cell RNA-seq (scRNAseq) this is possible.
In a previous study, the MASCARA approach showed a novel approach for compensating within-sample variation, however, this method was not yet applied on sparse scRNA-seq data. MASCARA is a combination of ASCA (ANOVA Simultaneous Component Analysis) and Partial Least Squares (PLS). Here the ASCA part is used to to compensate for the gene expression differences between cell types and the PLS part is used to find genes associated with a known pathway. In this work, we show the MASCARA approach applied to scRNA-seq data in which we investigate the suberin biosynthesis pathway that is involved in drought stress response. Furthermore, we show that the MASCARA method works for scRNA-seq data and offer suggestions to further improve on this approach. We expect that this study will allow for a valuable and relevant alternative to other established methods for studying co-expression in scRNA-seq data.

## Repo contents
This repo contains all the custom code for my master project at the UvA (University of Amsterdam). 

Folder contents:
- 01_playground
  - Test playground for various scripts and data

- 02_scRNA_analysis
  - Primary analysis for scRNA data analysis
  - Reference gff files
  - gene signature lists

- 03_MASCARA
  - Markdown files for MASCARA analysis using results from 03

- 04_normal_WGCNA_tutorial
   - Tutorial I found about WGCNA (unfinished)

- 05_hdWGCNA
   - Playground for hdWGCNA and meta-cells approaches


This work relies heavily on the code of Fred White et al.:

- https://github.com/BiosystemsDataAnalysis/MASCARA
