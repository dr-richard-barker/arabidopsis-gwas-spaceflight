# Provenance of results files

Which code produced each committed table, figure and data file. Written 14 Oct 2026 after an
`/abai check`. Scripts are in `code/` and run from the repository root.

Files marked **Biomni session; code not retained** were produced by scripts 01-19, which ran in a
Biomni workspace and are not in this repository. Their values are reported as they are, but they
cannot be regenerated or checked against code until those scripts are recovered.

## Tables

| File | Produced by |
|---|---|
| `tables/altitude_biomarker_correlations.csv` | **Biomni session; code not retained** (scripts 01-19) |
| `tables/altitude_regression_results.csv` | **Biomni session; code not retained** (scripts 01-19) |
| `tables/de_summary_by_study.csv` | **Biomni session; code not retained** (scripts 01-19) |
| `tables/ecotype_prediction_geo.csv` | **Biomni session; code not retained** (scripts 01-19) |
| `tables/ecotype_prediction_metrics.json` | **Biomni session; code not retained** (scripts 01-19). Its `n_countries: 45` disagrees with the 44 countries in `ecotype_prediction_geo.csv`. |
| `tables/enrichment_tests.csv` | **Biomni session; code not retained** (scripts 01-19) |
| `tables/go_enrichment_ora.csv` | **Biomni session; code not retained** (scripts 01-19) |
| `tables/gwas_de_overlap_genes.csv` | **Biomni session; code not retained** (scripts 01-19) |
| `tables/gwas_spaceflight_integration.csv` | **Biomni session; code not retained** (scripts 01-19) |
| `tables/ion_content_locus_sharing.csv` | **Biomni session; code not retained** (scripts 01-19) |
| `tables/ion_content_overlap.csv` | **Biomni session; code not retained** (scripts 01-19) |
| `tables/ionome_enrichment.csv` | `30_ionome_enrichment.py` |
| `tables/kegg_enrichment.csv` | **Biomni session; code not retained** (scripts 01-19) |
| `tables/meta_analysis_consensus_genes.csv` | **Biomni session; code not retained** (scripts 01-19) |
| `tables/meta_analysis_results.csv` | **Biomni session; code not retained** (scripts 01-19) |
| `tables/ml_feature_importance.csv` | **Biomni session; code not retained** (scripts 01-19) |
| `tables/ml_gene_classifier_results.csv` | **Biomni session; code not retained** (scripts 01-19) |
| `tables/ml_top_genes_for_ecotype.csv` | **Biomni session; code not retained** (scripts 01-19) |
| `tables/photorespiration_de_genes.csv` | **Biomni session; code not retained** (scripts 01-19) |
| `tables/photorespiration_ecotype_expression.csv` | **Biomni session; code not retained** (scripts 01-19) |
| `tables/power_analysis_results.csv` | `20_extended_analyses.py` |
| `tables/top10_flowering_times.csv` | `32_arapheno_flowering_times.py` |
| `tables/top10_phenotypic_profiles.csv` | `21_plot_top10_map.py` |
| `tables/top10_predicted_ecotypes.csv` | **Biomni session; code not retained** (scripts 01-19). Its unsourced `supporting_evidence` column was removed by hand on 14 Oct 2026; all other fields are unchanged. |
| `tables/trait_broad_category_summary.csv` | execution_trace notebook (cells 11-14, Biomni session 28 Aug 2026) |
| `tables/trait_category_overlap.csv` | execution_trace notebook (cells 11-14, Biomni session 28 Aug 2026) |
| `tables/trait_category_summary.csv` | execution_trace notebook (cells 11-14, Biomni session 28 Aug 2026) |
| `tables/trait_cluster_analysis.csv` | **Biomni session; code not retained** (scripts 01-19) |
| `tables/trait_enrichment.csv` | **Biomni session; code not retained** (scripts 01-19) |

## Figures

| Figure | Produced by |
|---|---|
| `figures/fig1_gwas_manhattan` | **Biomni session; code not retained** (scripts 01-19) |
| `figures/fig2_meta_analysis` | **Biomni session; code not retained** (scripts 01-19) |
| `figures/fig3_integration` | **Biomni session; code not retained** (scripts 01-19) |
| `figures/fig4_ml_results` | `29_fig4_ml_results.py` |
| `figures/fig5_knowledge_graph` | **Biomni session; code not retained** (scripts 01-19) |
| `figures/fig07_trait_network` | **Biomni session; code not retained** (scripts 01-19) |
| `figures/fig08_category_spaceflight_overlap` | execution_trace notebook (cells 11-14, Biomni session 28 Aug 2026) |
| `figures/fig09_category_gene_sharing` | **Biomni session; code not retained** (scripts 01-19) |
| `figures/fig10_trait_category_counts` | **Biomni session; code not retained** (scripts 01-19) |
| `figures/fig11_trait_dendrogram` | **Biomni session; code not retained** (scripts 01-19) |
| `figures/fig12_broad_category_overview` | execution_trace notebook (cells 11-14, Biomni session 28 Aug 2026) |
| `figures/fig13_ecotype_prediction_real` | `28_fig13_ecotype_prediction.py` |
| `figures/fig14_world_map` | **Biomni session; code not retained** (scripts 01-19) |
| `figures/fig15_europe_map` | `33_fig15_europe_map.py` |
| `figures/fig16_altitude_score` | **Biomni session; code not retained** (scripts 01-19) |
| `figures/fig17_biomarker_multipanel` | **Biomni session; code not retained** (scripts 01-19) |
| `figures/fig18_ionome_enrichment` | `30_ionome_enrichment.py` |
| `figures/fig19_photorespiration_heatmap` | **Biomni session; code not retained** (scripts 01-19) |
| `figures/fig22_clinal_latitude_gradient` | `20_extended_analyses.py` |
| `figures/fig23_chr4_aragwas_window` | `31_chr4_aragwas_window.py` |
| `figures/fig24_europe_score_heatmap` | `20_extended_analyses.py` |
| `figures/fig26_score_geographic_surface` | `20_extended_analyses.py` |
| `figures/fig27_population_violin` | `20_extended_analyses.py` |
| `figures/fig28_top10_ecotypes_map` | `26_atlas_political_fig28.py` (panel B, cropped 14 Oct 2026; see its header). `21_plot_top10_map.py` also draws a draft. |
| `figures/fig29_score_bar_plot` | `27_generate_dashboard_static_charts.py` |
| `figures/fig30_latitude_elevation_scatter` | `27_generate_dashboard_static_charts.py` |

The PNG copies in `docs/assets/` come from the same source as the figure in `figures/`: scripts write both, and figures without a script were copied.

## Other data

| File | Produced by |
|---|---|
| `data/aragwas/chr4_window_associations.csv` | `31_chr4_aragwas_window.py` |
| `data/aragwas/chr4_window_genes.csv` | `31_chr4_aragwas_window.py` |
| `data/knowledge_graph/kg_edges.csv` | **Biomni session; code not retained** (scripts 01-19) |
| `data/knowledge_graph/kg_nodes.csv` | **Biomni session; code not retained** (scripts 01-19) |
| `data/knowledge_graph/knowledge_graph_dashboard.json` | **Biomni session; code not retained** (scripts 01-19) |
| `docs/assets/fig_geo_europe.json` | **Biomni session; code not retained** (scripts 01-19) |
| `docs/assets/fig_geo_world.json` | **Biomni session; code not retained** (scripts 01-19) |
