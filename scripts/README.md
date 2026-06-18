This folder contains scripts used for analyzing transcriptomic, epigenomic, and cross-species conservation data from the in vitro wood-smoke exposure study. The number at the beginning of each script name indicates the order in which the analyses were run.

1. `01_WGCNA/1_Wood_Smoke_WGCNA_analysis.R` performs weighted gene co-expression network analysis of RNA-seq count data.
2. `02_Methylation_Expression_Correlation/2_Methylation_Expression_correlation.R` calculates gene-level correlations between DNA methylation and gene expression.
3. `03_Transcriptome_Similarity/3_Transcriptome_similarity_analysis.R` compares differential-expression signatures using overlap counts and Jaccard similarity.
4. `04_Rhesus_Human_Conservation/4_Rhesus_Human_conservation_analysis.R` evaluates rhesus-to-human orthology and genomic interval liftOver.

Each analysis subfolder contains a README describing its required input files.
