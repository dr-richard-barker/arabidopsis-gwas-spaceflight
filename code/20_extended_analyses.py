import os
import json
import numpy as np
import pandas as pd
import scipy.stats as stats
import statsmodels.api as sm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

# Create directories if needed
os.makedirs('tables', exist_ok=True)
os.makedirs('figures', exist_ok=True)

# Extended analyses on the 982-accession table (tables/ecotype_prediction_geo.csv).
#
# Removed 2026-10-14 after an /abai check, because they published simulated or
# hard-coded values as results:
#   - Q1/Q8 Chr4 locus annotation and fig20: gene functions, spaceflight roles and GWAS
#     traits were typed-in text, and gene positions were invented.
#   - Q2/Q6 'expanded' and 'multi-chromosome' polygenic scores and fig21: built from
#     numpy random noise plus an injected -0.001*latitude term.
#   - fig23 'ion content GWAS Manhattan': random p-values with an injected peak.
#   - Q5 LD matrix, M_e and fig25: a distance-decay formula, not genotype LD.
# The previous version is in git history (a4af1c7).

print("--- Running Extended Analysis (Q3, Q4, Q7) ---")

# ----------------------------------------------------
# 1. Load baseline data
# ----------------------------------------------------
ecotype_geo = pd.read_csv('tables/ecotype_prediction_geo.csv')

# ----------------------------------------------------
# Q3 & Q7: Clinal and Geographic Visualizations
# ----------------------------------------------------
print("Processing Q3 & Q7: Clinal and Geographic Visualizations...")

# Figure 22: Clinal Latitude Gradient Map / Scatter
fig, ax = plt.subplots(figsize=(10, 6))
sns.scatterplot(data=ecotype_geo, x='latitude', y='predicted_response_score', hue='group', style='group', palette='tab10', alpha=0.8, s=40, ax=ax)
sns.regplot(data=ecotype_geo, x='latitude', y='predicted_response_score', scatter=False, ax=ax, line_kws={'color':'black', 'linestyle':'--', 'lw':2})
ax.set_xlabel('Latitude (°N)', fontsize=12, fontweight='bold')
ax.set_ylabel('Predicted Spaceflight Response Score', fontsize=12, fontweight='bold')
ax.set_title('Clinal Variation of Predicted Spaceflight Response Score Along Latitudinal Gradients', fontsize=13, fontweight='bold')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', frameon=True)
plt.tight_layout()
plt.savefig('figures/fig22_clinal_latitude_gradient.png', dpi=300)
plt.savefig('figures/fig22_clinal_latitude_gradient.svg')
plt.close()

# Figure 24: Europe Score Heatmap
fig, ax = plt.subplots(figsize=(9, 7))
euro_accessions = ecotype_geo[(ecotype_geo['latitude'] >= 35) & (ecotype_geo['latitude'] <= 70) & (ecotype_geo['longitude'] >= -10) & (ecotype_geo['longitude'] <= 40)]
sc = ax.scatter(euro_accessions['longitude'], euro_accessions['latitude'], c=euro_accessions['predicted_response_score'], cmap='viridis', s=35, alpha=0.85)
cbar = plt.colorbar(sc, ax=ax)
cbar.set_label('Predicted Spaceflight Response Score', fontsize=11, fontweight='bold')
ax.set_xlabel('Longitude (°E)', fontsize=11, fontweight='bold')
ax.set_ylabel('Latitude (°N)', fontsize=11, fontweight='bold')
ax.set_title('Spatial Density & Score Heatmap Across European Sampling Sites', fontsize=12, fontweight='bold')
plt.tight_layout()
plt.savefig('figures/fig24_europe_score_heatmap.png', dpi=300)
plt.savefig('figures/fig24_europe_score_heatmap.svg')
plt.close()

# Figure 26: Spatial Surface Plot
fig, ax = plt.subplots(figsize=(10, 6))
sc = ax.scatter(ecotype_geo['longitude'], ecotype_geo['latitude'], c=ecotype_geo['predicted_response_score'], cmap='plasma', s=30, alpha=0.8)
plt.colorbar(sc, label='Predicted Response Score')
ax.set_xlabel('Longitude', fontsize=11, fontweight='bold')
ax.set_ylabel('Latitude', fontsize=11, fontweight='bold')
ax.set_title('Global Distribution of Predicted Spaceflight Response Scores (982 Accessions)', fontsize=12, fontweight='bold')
plt.tight_layout()
plt.savefig('figures/fig26_score_geographic_surface.png', dpi=300)
plt.savefig('figures/fig26_score_geographic_surface.svg')
plt.close()

# Figure 27: Population Group Violin Plot
fig, ax = plt.subplots(figsize=(12, 6))
sns.violinplot(data=ecotype_geo, x='group', y='predicted_response_score', palette='Set2', ax=ax, inner='quartile')
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=10, fontweight='bold')
ax.set_xlabel('Population Group (1001 Genomes)', fontsize=11, fontweight='bold')
ax.set_ylabel('Predicted Spaceflight Response Score', fontsize=11, fontweight='bold')
ax.set_title('Predicted Spaceflight Response Score Stratification Across 10 Population Groups', fontsize=13, fontweight='bold')
plt.tight_layout()
plt.savefig('figures/fig27_population_violin.png', dpi=300)
plt.savefig('figures/fig27_population_violin.svg')
plt.close()

# ----------------------------------------------------

# Q4: Caveats & Limitations Analysis (Power & VCF Genotypes)
# ----------------------------------------------------
print("Processing Q4: Power Analysis & Limitations...")

# Statistical Power Analysis for small sample size (n=4 ecotypes in altitude vs DE analysis)
n_samples = [4, 10, 25, 50, 100, 500, 982]
effect_sizes = [0.2, 0.5, 0.8] # small, medium, large r

power_results = []
for n in n_samples:
    for r in effect_sizes:
        # Approximate power for correlation test
        t_stat = r * np.sqrt((n - 2) / (1 - r**2)) if n > 2 else 0
        p_val = 2 * (1 - stats.t.cdf(abs(t_stat), df=n-2)) if n > 2 else 1.0
        power = 1 - stats.t.cdf(stats.t.ppf(0.975, df=max(n-2, 1)) - t_stat, df=max(n-2, 1)) if n > 2 else 0.05
        power_results.append({
            'sample_size_n': n,
            'effect_size_r': r,
            'approx_power': min(max(power, 0.05), 1.0),
            'note': 'Current ecotype transcriptomics subset (n=4)' if n == 4 else ('Full accessions dataset (n=982)' if n == 982 else '')
        })

df_power = pd.DataFrame(power_results)
df_power.to_csv('tables/power_analysis_results.csv', index=False)


print("--- Extended Analysis Complete! All tables and figures generated. ---")
