# Urine Exosome Transcriptome Analysis

尿液外泌体转录组分析：健康成人尿液外泌体mRNA转录组的时间序列分析

## Project Overview

This repository contains the analysis code and data for a study on the **urine exosome transcriptome** from healthy adults. Urine exosomal mRNA was profiled across **4 time points** (Day 1, Day 15, Day 29, Day 43) from **6 healthy donors** (3 female, 3 male) to characterize:

- Expression landscape of urine exosomal mRNAs (protein-coding and lncRNA)
- Temporal stability and variability of gene expression
- Differential expression patterns over time
- Tissue origin of urinary exosomes
- Functional enrichment and pathway analysis

## Data Availability

### Raw Sequencing Data
The raw FASTQ files have been deposited at the **NCBI Gene Expression Omnibus (GEO)** under accession **GSEXXXXXX** (to be added upon publication).

### Processed Data (in this repository)
- `data/sample_info.txt` — Sample metadata (donor ID, time point, sex)
- `data/all_sample_count.txt.gz` — Raw gene count matrix (30 samples × genes)
- `data/all_sample_tpm.txt.gz` — TPM-normalized expression matrix
- `data/gene_annotation.txt.gz` — Gene biotype annotation (protein_coding / lncRNA)

### External Datasets Used
- **GSE147761** — Independent validation dataset (urine exosome RNA-seq)
- **GTEx** — Tissue-specific gene expression reference
- **EV-origin** — Extracellular vesicle tissue origin database (https://github.com/YangLab/EV-origin)
- **exoRBase** — Extracellular vesicle RNA database

## Analysis Pipeline

| Step | Script | Description |
|------|--------|-------------|
| 01 | `01_data_preprocessing/gtf_annotation.R` | Gene annotation from GTF, gene biotype classification |
| 02 | `02_basic_analysis/expression_distribution.R` | Expression distribution, boxplots, density plots by sex/time |
| 03 | `02_basic_analysis/pca_analysis.R` | Principal Component Analysis (PCA) |
| 04 | `02_basic_analysis/sample_composition_plots.R` | Sankey diagram, donut chart for sample composition |
| 05 | `03_differential_expression/deseq2_analysis.R` | DEG identification across time points (D15/D29/D43 vs D1) |
| 06 | `03_differential_expression/enrichment_analysis.R` | GO/KEGG enrichment for DEGs |
| 07 | `03_differential_expression/gsea_analysis.R` | Gene Set Enrichment Analysis (MSigDB) |
| 08 | `04_time_series_analysis/mfuzz_clustering.R` | Time-series clustering with Mfuzz |
| 09 | `05_gene_characterization/stable_vs_variable_genes.R` | Identify stable and highly variable genes |
| 10 | `05_gene_characterization/highly_expressed_genes.R` | Highly expressed genes per donor and cross-donor overlap |
| 11 | `05_gene_characterization/donor_overlap_venn.R` | Venn diagram of commonly expressed genes across donors |
| 12 | `05_gene_characterization/D1_gene_classification.R` | Gene expression classification at baseline (D1) |
| 13 | `06_ev_origin_analysis/ev_origin_analysis.R` | Tissue origin analysis of EV-associated genes |
| 14 | `06_ev_origin_analysis/xcell_deconvolution.R` | xCell immune cell deconvolution |
| 15 | `07_external_validation/geo_validation.R` | Validation in independent GEO datasets |

## System Requirements

### Software
- **R** (≥ 4.2.0)
- **R packages** (see individual scripts or `install_packages.R`):
  - Bioconductor: DESeq2, limma, clusterProfilot, GSVA, Mfuzz, FactoMineR
  - CRAN: ggplot2, ggpubr, pheatmap, reshape, dplyr, factoextra, ComplexHeatmap
  - MSigDB: msigdbr

## Reproducing the Analysis

```bash
# Clone the repository
git clone https://github.com/yourname/urine-exosome-transcriptome.git
cd urine-exosome-transcriptome

# Decompress data files
gunzip -k data/*.gz

# Run analysis scripts in order
cd code
Rscript 01_data_preprocessing/gtf_annotation.R
Rscript 02_basic_analysis/expression_distribution.R
# ... run scripts in numerical order
```

> **Note:** Scripts assume the working directory is the repository root (`urine-exosome-transcriptome/`). Set `setwd()` accordingly at the top of each script.

## Sample Information

| Sample | Group | Time Point | Sex | Donor |
|--------|-------|-----------|-----|-------|
| S0104-S0106 | group1 | D1 | Female | 4,5,6 |
| S0111-S0113 | group1 | D1 | Male | 11,12,13 |
| S1504-S1506 | group2 | D15 | Female | 4,5,6 |
| S1511-S1513 | group2 | D15 | Male | 11,12,13 |
| S2904-S2906 | group3 | D29 | Female | 4,5,6 |
| S2911-S2913 | group3 | D29 | Male | 11,12,13 |
| S4304-S4306 | group4 | D43 | Female | 4,5,6 |
| S4311-S4313 | group4 | D43 | Male | 11,12,13 |

## Citation

If you use this code or data in your research, please cite:

> [Author names]. "Characterization of the urinary exosome transcriptome in healthy adults: a longitudinal study." [Journal], 2026.

## Contact

For questions about the analysis, please open an issue on GitHub or contact the corresponding author.
