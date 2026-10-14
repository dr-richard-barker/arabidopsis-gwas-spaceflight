> **Historical record, not current results.** This is the plan the Biomni agent wrote on 28 Aug 2026
> for the geographic, altitude-biomarker and manuscript-revision step. It is kept unchanged
> below for provenance. Its counts are as stated at planning time. Note that it records an earlier
> version of the analysis that used synthetic accessions, which this step replaced with 982 real 1001 Genomes
> accessions. The scripts it names (13-16) are not in this repository. For current results, see
> [`manuscript/manuscript.tex`](../manuscript/manuscript.tex).

# Plan: Geographic Mapping, Altitude-Biomarker Analysis, and Manuscript Upgrade

## Summary

Extend the ecotype prediction analysis by mapping the 982 real 1001 Genomes accessions onto actual geographic coordinates with elevation data, testing whether altitude correlates with predicted spaceflight response scores, photorespiration biomarker expression, and plant mineral ion content GWAS gene overlap. Then perform a full manuscript revision incorporating all new results, figures, and discussion.

## Data Sources (all verified accessible)

1. **Coordinates**: `hagax8/arabidopsis_viz` GitHub CSV — 1,135 accessions with real lat/lon, 100% match to our 982 accessions
2. **Elevation**: Open-Meteo Elevation API (`api.open-meteo.com/v1/elevation`) — batch support, free, no key needed
3. **Photorespiration genes**: Curated list of ~15 core genes (GOX1, SHM1/2, PGLP1, GLYK, GGT1/2, CAT2, HPR1, GDC-P, GLU1, GLYDH, GLO1) + all genes annotated to GO:0009853 from our GO data
4. **Ion content GWAS genes**: 10 ion traits in database (Mo98, S34, Co59, As75, Se82, Zn66, At2, M130T666, M172T666, M216T665)
5. **DE results**: 534,547 rows, 53 photorespiration gene entries significant (padj<0.05)
6. **Existing results**: `ecotype_prediction_results.csv` (982 rows), `trait_cluster_analysis.csv`, `trait_category_overlap.csv`

## Implementation Steps

### Step 1: Data Acquisition & Merge (`code/13_geo_altitude_data.py`)

- Download GitHub CSV with coordinates for 1,135 accessions
- Merge with existing `ecotype_prediction_results.csv` on accession ID → 982 rows with real lat/lon
- Fetch elevation for all 982 accessions via Open-Meteo batch API (100 coords per request, ~10 requests, ~15 seconds)
- Save merged dataset as `results/tables/ecotype_prediction_geo.csv` (982 rows with: accession_id, name, country, group, latitude, longitude, elevation, predicted_response_score, genetic_diversity, etc.)
- Compile photorespiration gene set: curated list + GO:0009853 annotations from `data/go/` directory
- Query GWAS database for ion content trait genes (10 traits)
- Save photorespiration gene list and ion content gene list

### Step 2: Geographic Visualizations (`code/14_geo_visualizations.py`)

**Static (SVG/PNG for manuscript):**
- Fig 14: World map (cartopy) with accessions as dots colored by predicted score (colorbar: Very Low → Very High). Point size by genetic diversity. Title: "Predicted Spaceflight Response Scores Across 982 Arabidopsis Accessions"
- Fig 15: Europe-focused map (cartopy, zoomed to 35-70°N, -10-40°E) showing the dense European accessions with score coloring
- Fig 16: Altitude vs predicted response score scatter plot with regression line, Spearman rho annotation, colored by population group

**Interactive (Plotly HTML/JSON for dashboard):**
- `docs/assets/data/fig_geo_world.json` — Plotly scatter_mapbox with hover: name, country, elevation, score, group
- `docs/assets/data/fig_geo_europe.json` — same, Europe zoomed

### Step 3: Statistical Analysis (`code/15_altitude_biomarker_stats.py`)

**Correlations (Spearman, BH-FDR corrected):**
1. Altitude vs predicted response score
2. Altitude vs genetic diversity (n_het + n_alt_hom)
3. Latitude vs predicted response score
4. Altitude vs mean log2FC of DE photorespiration genes (per ecotype if data available, or aggregate)
5. Altitude vs ion content GWAS gene count (per accession from genotype data)

**Multiple Regression (OLS):**
1. `predicted_response_score ~ elevation + latitude + C(population_group)`
2. `genetic_diversity ~ elevation + latitude + C(population_group)`
3. If per-ecotype photorespiration expression available: `mean_photorespiration_log2FC ~ elevation + latitude`

**Ion content analysis:**
- For each of the 10 ion content traits, get GWAS genes and test overlap with spaceflight DEGs (Fisher exact)
- Test whether accessions from high-altitude locations carry more ion content GWAS alleles at the 31 SNP positions
- Correlation: altitude vs count of ion-content-associated alleles per accession

**Output:** `results/tables/altitude_biomarker_correlations.csv`, `results/tables/ion_content_overlap.csv`

### Step 4: Additional Visualizations (`code/16_biomarker_visualizations.py`)

- Fig 17: Multi-panel figure: (A) Altitude vs score with regression, (B) Altitude vs photorespiration gene expression, (C) Altitude vs ion content allele count, (D) Correlation heatmap of all variables
- Fig 18: Ion content trait enrichment bar chart (10 traits, fold enrichment with spaceflight DEGs)
- Fig 19: Photorespiration gene expression heatmap (curated genes × studies/ecotypes, log2FC colored)

### Step 5: Manuscript Full Revision (`manuscript/latex/manuscript.tex`)

**Abstract:** Add mention of geographic analysis, altitude associations, photorespiration and ion content findings, updated accession count (982 real vs 2029 synthetic).

**Results — new subsections:**
1. "Geographic patterns reveal altitude-associated variation in predicted spaceflight response" — world/Europe maps, altitude-score correlation, regression results
2. "Altitude links to photorespiration biomarker expression and ion homeostasis genes" — photorespiration DE genes, ion content GWAS overlap, altitude-biomarker correlations

**Results — update existing:**
- "Ecotype prediction" section: replace synthetic accession language with real 1001 Genomes data (982 accessions, 45 countries, 10 population groups, real genotypes at 31 SNPs)

**Discussion — new paragraphs:**
- Altitude as a proxy for multiple stressors (UV, temperature, CO2 partial pressure) and its connection to spaceflight stress (radiation, altered gas exchange)
- Photorespiration gene dysregulation in spaceflight: mechanistic link between altered atmospheric composition in space and the photorespiratory pathway
- Ion homeostasis as the highest-enrichment category: connection between microgravity-induced fluid dynamics changes and terrestrial ion adaptation alleles
- Geographic structuring of predicted responses: implications for selecting accessions for future spaceflight experiments

**Methods — new subsections:**
- "Geographic coordinate and elevation data"
- "Altitude-biomarker association analysis"
- Update "Ecotype prediction" methods with real VCF data source

**Figures:** Add Figs 14-19 to manuscript, update figure references throughout.

### Step 6: Copy All Outputs to `/mnt/results/`

- Tables: `ecotype_prediction_geo.csv`, `altitude_biomarker_correlations.csv`, `ion_content_overlap.csv`
- Figures: `fig14` through `fig19` (SVG + PNG)
- Updated manuscript: `manuscript.tex`, compiled PDF if possible
- Updated dashboard JSON files

## Compute Estimate

- Elevation API: ~10 batch requests, <30 seconds
- Cartopy maps: CPU-only, <2 minutes per figure
- Statistical analysis: trivial (982 rows), <1 minute
- Manuscript LaTeX compilation: <1 minute if pdflatex available
- **Total estimated runtime: ~15-20 minutes**
- **Machine: worker-0 (default, autoscale) — no additional resources needed**

## Assumptions

1. Cartopy is installed or installable via `uv pip install cartopy`
2. Plotly is installed or installable via `uv pip install plotly`
3. GO:0009853 annotation file exists in `data/go/` directory (will check and download from GO website if missing)
4. Per-ecotype photorespiration expression: if DE results don't have enough ecotype-level resolution for regression, we'll use aggregate scores and correlation only
5. LaTeX compiler (`pdflatex`) may not be available in sandbox — if not, we'll update the .tex file and note that compilation should be done locally
6. All 982 accessions will get elevation data (Open-Meteo has global coverage)
