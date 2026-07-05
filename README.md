# mrps16-expression-analysis-
# Exploration of the Dichotomous Role of MRPS16 in Cancer

## Comparative Analysis of Gene Expression in Glioblastoma (GBM) and Hepatocellular Carcinoma (LIHC)

### Overview

This repository contains the source code, datasets and analysis results for the Master's mini-project in Bioinformatics at the University of Kinshasa (UNIKIN).

The study investigates the expression of the **MRPS16** gene in tumor and normal tissues from two cancer types:

- Glioblastoma (GBM)
- Liver Hepatocellular Carcinoma (LIHC)

using transcriptomic data from the **TCGA TARGET GTEx** cohort available through the **UCSC Xena Browser**.

---

## Objectives

The objectives of this study are to:

- Compare MRPS16 expression between tumor and normal tissues.
- Evaluate the statistical significance of the observed differences.
- Compare the expression profiles of GBM and LIHC.
- Discuss the biological significance of the findings based on the scientific literature.

---

## Data Source

The data were downloaded from:

**UCSC Xena Browser**

Dataset:

- **TCGA TARGET GTEx (Toil RNA-seq Recompute)**

This dataset combines:

- TCGA tumor samples
- GTEx normal tissue samples

using the same RNA-seq processing pipeline (Toil), allowing direct comparison.

---

## Study Design

### Glioblastoma (GBM)

| Group | Samples |
|--------|---------|
| Tumor | 171 |
| Normal | 1141 |

### Liver Hepatocellular Carcinoma (LIHC)

| Group | Samples |
|--------|---------|
| Tumor | 421 |
| Normal | 110 |

---

## Repository Structure

```
mrps16-expression-analysis/
├── gbm_mrps16.tsv
└── lihc_mrps16.tsv
└── mrps16_analysis.py
├── mrps16_boxplots.png
└── mrps16_statistiques.csv
├── README.md
```

---

## Software

The analysis was performed in Python using the following libraries:

- pandas
- scipy
- matplotlib
- seaborn

---

## Statistical Analysis

The comparison between tumor and normal tissues was performed using the **Mann–Whitney U test** (Wilcoxon rank-sum test).

A p-value lower than 0.05 was considered statistically significant.

---

## Main Results

| Cancer | log₂ Fold Change | p-value |
|---------|-----------------|----------|
| GBM | +0.39 | 4.86 × 10⁻⁴³ |
| LIHC | +0.69 | 1.79 × 10⁻³⁷ |

The analysis shows that **MRPS16 is significantly overexpressed in tumor tissues compared with normal tissues in both GBM and LIHC**.

The observed "dichotomous role" reported in the literature appears to be related to the biological function of MRPS16 rather than to opposite expression patterns.

---

## Reproducibility

Install the required Python packages:

```bash
pip install -r requirements.txt
```

Run the analysis:

```bash
python mrps16_analysis.py
```

The script generates:

- Statistical summary
- Boxplots
- CSV results

---

## Citation

If you use this repository, please cite:

- Goldman MJ et al. *The UCSC Xena Platform for Cancer Genomics Data Visualization and Interpretation*. Nature Biotechnology, 2020.
- GTEx Consortium. *The GTEx Atlas of Genetic Regulatory Effects Across Human Tissues*. Science, 2020.
- Wang Z et al. *MRPS16 facilitates tumor progression via PI3K/AKT/Snail signaling*. Journal of Cancer, 2020.

---

## Authors
TWELA MWAY Elton, KASONGO NJIMINY Landers, NGOLO MALONDA Kethia, MUELA BANANGA Miradel, OYEMA MBOLADINGA Dawen

Master's AI and Data science

Faculty of Sciences

University of Kinshasa (UNIKIN)

Academic Year: 2025–2026

---
