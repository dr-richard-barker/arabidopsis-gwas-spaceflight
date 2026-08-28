# Arabidopsis GWAS-Spaceflight Integration

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Data License: CC BY 4.0](https://img.shields.io/badge/Data_License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![Python 3.11](https://img.shields.io/badge/Python-3.11-green.svg)](environment.yml)
[![R 4.3](https://img.shields.io/badge/R-4.3-blue.svg)](environment.yml)

> **Predicting Ecotype-Specific Spaceflight Responses by Integrating the AraGWAS Catalog with NASA GeneLab Spaceflight Transcriptomics Data**

This repository provides an integrated computational framework linking ground-based natural genetic variation (*Arabidopsis thaliana* AraGWAS catalog) with spaceflight transcriptomic responses (NASA GeneLab / OSDR).

---

## 🌟 Key Findings

1. **GWAS-Spaceflight Overlap:** Natural genetic variation underlying terrestrial stress adaptation and ion homeostasis strongly predicts spaceflight response sensitivity (763 overlap genes, **Odds Ratio = 4.23, $p = 4.45 \times 10^{-102}$**).
2. **Chromosome 4 Pleiotropic Hotspot:** 29 of 31 top spaceflight-predictive SNPs (93.5%) co-localize with ion content GWAS QTLs on Chromosome 4 (positions `1,257,343` to `1,271,408`; spanning `AT4G02820` to `AT4G02860`).
3. **Photorespiratory Stagnation Crisis:** Microgravity-induced absence of buoyancy convection leads to boundary layer gas stagnation, inducing **53 of 56 photorespiratory C2 cycle genes (94.6%)**.
4. **Latitudinal Clinal Variation:** Predicted spaceflight response scores exhibit a statistically significant latitudinal gradient (Spearman $\rho = -0.125, p = 8.97 \times 10^{-5}$), structured primarily by population demographic history ($R^2 = 0.349$).
5. **Linkage Disequilibrium Redundancy:** Linkage disequilibrium analysis across the Chr4 locus reveals an effective independent marker count $M_e = 4.2$ out of 29 SNPs, highlighting single-locus polygenic score constraints and motivating multi-chromosomal score expansion.

---

## 📁 Repository Structure

```
arabidopsis-gwas-spaceflight/
├── README.md                          # Project overview & quickstart
├── LICENSE                            # MIT License
├── CITATION.cff                       # Citation metadata
├── .zenodo.json                       # Zenodo archive metadata
├── environment.yml                    # Conda environment
├── Dockerfile                         # Reproducibility container
│
├── manuscript/                        # LaTeX manuscript & references
│   ├── manuscript.tex
│   └── references.bib
│
├── code/                              # Analytical scripts (01 to 20)
│   ├── 01_download_gwas.py
│   ├── 02_download_genelab.py
│   ├── 03_differential_expression.R
│   ├── 04_meta_analysis.R
│   ├── 05_gwas_spaceflight_integration.py
│   ├── 06_go_enrichment.R
│   ├── 07_ml_gene_classifier.py
│   ├── 08_ml_ecotype_prediction.py
│   ├── 09_knowledge_graph.py
│   ├── 11_figures.py
│   ├── 13_geo_altitude_data.py
│   ├── 14_geo_visualizations.py
│   ├── 15_altitude_biomarker_stats.py
│   ├── 16_biomarker_visualizations.py
│   └── 20_extended_analyses.py        # Q1-Q8 Extended follow-up analyses
│
├── figures/                           # 27 publication figures (PNG & SVG)
├── tables/                            # 25 processed data CSV & JSON files
├── data/
│   └── knowledge_graph/               # Cytoscape.js & GraphML graph files
└── docs/                              # Interactive GitHub Pages website
```

---

## 🚀 Interactive GitHub Pages Website

The project is hosted interactively via GitHub Pages in the `docs/` directory:
- **Interactive Results Gallery:** Meta-analysis, GWAS overlap, ML, geographic clines, Chr4 shared locus, and power calculations.
- **Figure Browser:** High-resolution view of all 27 publication figures.
- **Data Downloads:** Schema documentation and download links for all 25 CSV data tables.

---

## 💻 Quickstart & Reproducibility

### Conda Environment
```bash
conda env create -f environment.yml
conda activate arabidopsis-gwas-spaceflight
```

### Run Pipeline & Extended Analyses
```bash
# Run extended follow-up analyses (Q1 - Q8)
python code/20_extended_analyses.py
```

### Docker Container
```bash
docker build -t arabidopsis-gwas-spaceflight .
docker run -it arabidopsis-gwas-spaceflight bash
```

---

## 📜 Citation

If you use this codebase or dataset, please cite:
> Barker, R. (2026). Integrating Arabidopsis GWAS with spaceflight transcriptomics to predict ecotype-specific responses. *npj Microgravity*, DOI: `10.1038/s41526-026-XXXXX-X`.
