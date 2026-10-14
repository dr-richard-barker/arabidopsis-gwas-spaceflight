# Arabidopsis GWAS-Spaceflight Integration

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Data License: CC BY 4.0](https://img.shields.io/badge/Data_License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![Python 3.11](https://img.shields.io/badge/Python-3.11-green.svg)](environment.yml)
[![R 4.3](https://img.shields.io/badge/R-4.3-blue.svg)](environment.yml)

> **Predicting Ecotype-Specific Spaceflight Responses by Integrating the AraGWAS Catalog with NASA GeneLab Spaceflight Transcriptomics Data**

This repository provides an integrated computational framework linking ground-based natural genetic variation (*Arabidopsis thaliana* AraGWAS catalog) with spaceflight transcriptomic responses (NASA GeneLab / OSDR).

---

## 🌟 Key Findings

1. **GWAS-Spaceflight Overlap:** Genes with Bonferroni-significant AraGWAS associations are enriched among spaceflight DEGs (763 overlap genes, **Odds Ratio = 4.23, $p = 4.45 \times 10^{-102}$**); 90 of 166 GWAS traits are individually enriched (FDR < 0.1), led by climate-adaptation and leaf ionome traits.
2. **Chromosome 4 Locus:** The 31 SNPs used for ecotype scoring are the 31 highest-scoring associations in the AraGWAS catalog, chosen on GWAS score alone. 29 of them lie in one 14 kb block on chromosome 4 (positions `1,257,343` to `1,271,408`; `AT4G02820` to `AT4G02860`). At those SNPs the associations are with three LC-MS leaf metabolite features (M216T665, M172T666, M130T666); arsenic, cobalt and molybdenum associations sit at neighbouring SNPs in the same block (Fig. 23).
3. **Leaf Ionome Genes:** Genes associated with six ICP-MS ionome traits are enriched among spaceflight DEGs (41 of 49, fold enrichment 1.92, one-sided Fisher $p = 8.0 \times 10^{-9}$; Fig. 18).
4. **Photorespiration Genes:** 53 of 56 photorespiration genes (94.6%) are differentially expressed in spaceflight.
5. **Latitudinal Clinal Variation:** Predicted spaceflight response scores show a latitudinal gradient (Spearman $\rho = -0.125, p = 8.97 \times 10^{-5}$). A regression on population group, latitude and elevation explains 34.9% of score variance ($R^2 = 0.349$), with the population-group terms the largest.

## ⚠️ Status: audited 14 Oct 2026

An evidence audit (`/abai check`) found problems that have been corrected or flagged:

- **Code for the core results is not in this repository.** Scripts 01–19 (data download, differential expression, meta-analysis, enrichment, machine learning, ecotype scoring, knowledge graph, geographic and altitude statistics) were run in a Biomni workspace and not retained. The results tables they produced are here, but cannot yet be regenerated from this repository. The archived Biomni plan and session notebook are in `execution_trace/`.
- **Simulated and unsourced outputs were removed.** An earlier version published a synthetic accession map (old Fig. 6) and score distribution (old Fig. 4C); an LD matrix, Manhattan plot and "expanded" polygenic score generated from formulas or random numbers (old Figs 21, 23, 25 and their tables); a locus diagram with invented gene positions and hand-typed annotations (old Fig. 20); and an ion-enrichment figure that labelled metabolite features as cadmium (old Fig. 18). Unsourced habitat, ionome and morphology profiles for the top 10 accessions were also removed.
- **Transcriptomic inputs need re-analysis.** Five of the 13 GeneLab studies are dark-grown BRIC experiments, and the meta-analysis table contains unmapped microarray row IDs and ERCC spike-ins. Open items are marked `TODO` in `manuscript/manuscript.tex`.

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
├── code/                              # Scripts 20-34 (01-19 not in this repository; see Status)
│   ├── 20_extended_analyses.py        # Latitude/geographic plots (Figs 22, 24, 26, 27) + power table
│   ├── 21_plot_top10_map.py           # top10_phenotypic_profiles.csv (also draws a draft Fig 28)
│   ├── 22-25_*.py                     # Earlier Fig 28 variants with synthetic terrain, not for publication
│   ├── 26_atlas_political_fig28.py    # Fig 28 = its panel B (Natural Earth map), cropped; its panel A terrain is synthetic
│   ├── 27_generate_dashboard_static_charts.py  # Figs 29-30
│   ├── 28_fig13_ecotype_prediction.py # Fig 13 (from the archived notebook, on the 982-accession table)
│   ├── 29_fig4_ml_results.py          # Fig 4 from the committed classifier and ecotype tables
│   ├── 30_ionome_enrichment.py        # Six-trait ionome enrichment + Fig 18
│   ├── 31_chr4_aragwas_window.py      # Chr4 AraGWAS window (fetched, cached in data/aragwas/) + Fig 23
│   ├── 32_arapheno_flowering_times.py # FT10/FT16 for the top 10 from AraPheno
│   ├── 33_fig15_europe_map.py         # Fig 15
│   └── 34_fig16_altitude_score.py     # Fig 16
│
├── execution_trace/                   # Archived 28 Aug 2026 Biomni plan and session notebook
├── figures/                           # 26 figures (PNG, most also SVG)
├── tables/                            # 29 processed data CSV & JSON files (see tables/PROVENANCE.md)
├── data/
│   ├── aragwas/                       # Cached AraGWAS chr4 window (code/31)
│   └── knowledge_graph/               # Cytoscape.js JSON and node/edge CSV files
└── docs/                              # Interactive GitHub Pages website
```

---

## 🚀 Interactive GitHub Pages Website

The project is hosted interactively via GitHub Pages in the `docs/` directory:
- **Interactive Results Gallery:** Meta-analysis, GWAS overlap, ML, geographic clines, the chromosome 4 locus, and power calculations.
- **Figure Browser:** The main and geographic figures.
- **Data Downloads:** Descriptions and download links for the main data tables.

---

## 💻 Quickstart

### Conda Environment
```bash
conda env create -f environment.yml
conda activate arabidopsis-gwas-spaceflight
```

### Run the Extended Analyses
```bash
# Latitude/geographic plots and power table, from tables/ecotype_prediction_geo.csv
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
