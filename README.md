# InVitroWoodSmoke

This repository contains code used for analyses accompanying the in vitro wood-smoke exposure manuscript from the Ji Lab at the University of California, Davis.

The `scripts/` folder contains numbered analysis scripts. The number at the beginning of each script name indicates the order in which the analyses were run. Each analysis folder contains a README listing the required input data.

## Repository structure

```text
InVitroWoodSmoke/
├── README.md
├── scripts/
│   ├── README.md
│   ├── 01_WGCNA/
│   │   ├── README.md
│   │   └── 1_Wood_Smoke_WGCNA_analysis.R
│   ├── 02_Methylation_Expression_Correlation/
│   │   ├── README.md
│   │   └── 2_Methylation_Expression_correlation.R
│   ├── 03_Transcriptome_Similarity/
│   │   ├── README.md
│   │   └── 3_Transcriptome_similarity_analysis.R
│   └── 04_Rhesus_Human_Conservation/
│       ├── README.md
│       └── 4_Rhesus_Human_conservation_analysis.R
├── data/
└── results/
```

Input data should be placed under `data/` using the paths specified near the beginning of each script. Generated tables, figures, and R objects are written under `results/`.

Large or controlled-access data should not be committed directly to GitHub. Provide accession numbers, download instructions, or an external archival link instead.
