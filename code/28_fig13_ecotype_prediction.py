import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['svg.hashsalt'] = 'arabidopsis-gwas-spaceflight'  # reproducible SVG ids
import matplotlib.pyplot as plt

# Figure 13: four-panel summary of the 982 scored accessions.
# Ported from cell 16 of execution_trace/biomni_worker-0_2026-08-28.ipynb, which read
# ecotype_prediction_results.csv (same accessions and scores; country-level placeholder
# coordinates, not used here). Changes from that cell: reads ecotype_prediction_geo.csv,
# computes the title counts from the data (the cell hard-coded "45 Countries"; the data
# has 44), and adds a population-group legend to panel D.

os.makedirs('figures', exist_ok=True)
os.makedirs('docs/assets', exist_ok=True)

results = pd.read_csv('tables/ecotype_prediction_geo.csv').sort_values('rank')
n_acc = len(results)
n_countries = results['country'].nunique()
n_groups = results['group'].nunique()
print(f'{n_acc} accessions, {n_countries} countries, {n_groups} population groups')

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Panel A: score distribution by population group
ax = axes[0, 0]
group_order = results.groupby('group')['predicted_response_score'].median().sort_values(ascending=False).index
group_data = [results[results['group'] == g]['predicted_response_score'].values for g in group_order]
bp = ax.boxplot(group_data, vert=True, patch_artist=True, medianprops=dict(color='black', linewidth=1.5))
for patch, color in zip(bp['boxes'], plt.cm.Set2(np.linspace(0, 1, len(group_order)))):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)
ax.set_xticks(range(1, len(group_order) + 1))
ax.set_xticklabels(group_order, rotation=45, ha='right', fontsize=8)
ax.set_ylabel('Predicted Response Score', fontsize=10)
ax.set_title('A. Score Distribution by Population Group', fontsize=11, fontweight='bold')
ax.set_ylim(0.3, 0.7)

# Panel B: top 20 accessions
ax = axes[0, 1]
top20 = results.head(20)
ax.barh(range(20), top20['predicted_response_score'].values,
        color=plt.cm.tab20(np.linspace(0, 1, 20)), edgecolor='black', linewidth=0.5)
ax.set_yticks(range(20))
ax.set_yticklabels([f"{row['name']} ({row['country']})" for _, row in top20.iterrows()], fontsize=7)
ax.set_xlabel('Predicted Response Score', fontsize=10)
ax.set_title('B. Top 20 Accessions (Highest Predicted Response)', fontsize=11, fontweight='bold')
ax.invert_yaxis()

# Panel C: country-level mean scores (top 15 countries with n >= 5)
ax = axes[1, 0]
country_stats = results.groupby('country').agg(
    n_accessions=('name', 'count'),
    mean_score=('predicted_response_score', 'mean'),
).sort_values('mean_score', ascending=False)
top_countries = country_stats[country_stats['n_accessions'] >= 5].head(15)
ax.barh(range(len(top_countries)), top_countries['mean_score'].values,
        color=plt.cm.coolwarm(np.linspace(0, 1, len(top_countries))),
        edgecolor='black', linewidth=0.5, alpha=0.8)
ax.set_yticks(range(len(top_countries)))
ax.set_yticklabels([f"{c} (n={n})" for c, n in zip(top_countries.index, top_countries['n_accessions'])], fontsize=8)
ax.set_xlabel('Mean Predicted Response Score', fontsize=10)
ax.set_title('C. Country-Level Mean Scores (n>=5)', fontsize=11, fontweight='bold')
ax.invert_yaxis()

# Panel D: genetic diversity vs predicted score, coloured by population group
ax = axes[1, 1]
groups = sorted(results['group'].unique())
cmap = plt.cm.tab10
for i, g in enumerate(groups):
    sub = results[results['group'] == g]
    ax.scatter(sub['genetic_diversity'], sub['predicted_response_score'], color=cmap(i),
               s=20, alpha=0.6, edgecolors='black', linewidth=0.3, label=g)
ax.set_xlabel('Genetic Diversity (n_het + n_alt_hom, of 31 SNPs)', fontsize=10)
ax.set_ylabel('Predicted Response Score', fontsize=10)
ax.set_title('D. Genetic Diversity vs Predicted Score', fontsize=11, fontweight='bold')
ax.legend(fontsize=7, loc='upper left', bbox_to_anchor=(1.01, 1.0), frameon=True, title='Population group', title_fontsize=8)

plt.suptitle(f'Ecotype Prediction from 1001 Genomes Genotypes\n'
             f'{n_acc} Accessions, 31 GWAS SNPs, {n_countries} Countries',
             fontsize=13, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig('figures/fig13_ecotype_prediction_real.png', dpi=300, bbox_inches='tight')
plt.savefig('figures/fig13_ecotype_prediction_real.svg', bbox_inches='tight')
plt.savefig('docs/assets/fig13_ecotype_prediction_real.png', dpi=300, bbox_inches='tight')
plt.close()
print('Saved figures/fig13_ecotype_prediction_real.png/.svg')
