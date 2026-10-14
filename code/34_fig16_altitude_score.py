import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['svg.hashsalt'] = 'arabidopsis-gwas-spaceflight'  # reproducible SVG ids
import matplotlib.pyplot as plt
from scipy import stats

# Figure 16: elevation and predicted spaceflight response score.
# The original figure's script (14/16, not in this repository) labelled the "Very High" category
# as "High" in panel B. This redraws both panels from tables/ecotype_prediction_geo.csv.

os.makedirs('figures', exist_ok=True)
os.makedirs('docs/assets', exist_ok=True)

eco = pd.read_csv('tables/ecotype_prediction_geo.csv').dropna(subset=['elevation_m'])
x, y = eco['elevation_m'], eco['predicted_response_score']
rho, p_rho = stats.spearmanr(x, y)
r, p_r = stats.pearsonr(x, y)
slope, intercept = np.polyfit(x, y, 1)
print(f'n={len(eco)}  Spearman rho={rho:.3f} (p={p_rho:.2e})  Pearson r={r:.3f} (p={p_r:.2e})')

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6), gridspec_kw={'width_ratios': [1.4, 1]})
for i, g in enumerate(sorted(eco['group'].unique())):
    sub = eco[eco['group'] == g]
    ax1.scatter(sub['elevation_m'], sub['predicted_response_score'], s=14, alpha=0.6, color=plt.cm.tab10(i), label=g)
xs = np.linspace(x.min(), x.max(), 100)
ax1.plot(xs, slope * xs + intercept, 'k--', linewidth=1.5, label='Linear fit')
ax1.set_xlabel('Elevation (m)', fontsize=11)
ax1.set_ylabel('Predicted Response Score', fontsize=11)
ax1.set_title(f'A. Elevation vs Predicted Score (n = {len(eco)})\n'
              f'Spearman rho = {rho:.3f} (p = {p_rho:.2e}); Pearson r = {r:.3f} (p = {p_r:.2e})', fontsize=11)
ax1.legend(fontsize=7, loc='upper left', bbox_to_anchor=(1.01, 1.0), title='Population group', title_fontsize=8)

order = [c for c in ['Very Low', 'Low', 'Moderate', 'High', 'Very High'] if (eco['predicted_category'] == c).any()]
data = [eco.loc[eco['predicted_category'] == c, 'elevation_m'] for c in order]
bp = ax2.boxplot(data, patch_artist=True, medianprops=dict(color='black'))
for patch, color in zip(bp['boxes'], plt.cm.YlOrRd(np.linspace(0.2, 0.8, len(order)))):
    patch.set_facecolor(color)
ax2.set_xticks(range(1, len(order) + 1))
ax2.set_xticklabels([f'{c}\n(n={len(d)})' for c, d in zip(order, data)], fontsize=9)
ax2.set_xlabel('Predicted Response Category', fontsize=11)
ax2.set_ylabel('Elevation (m)', fontsize=11)
ax2.set_title('B. Elevation by Predicted Category', fontsize=11)

plt.tight_layout()
plt.savefig('figures/fig16_altitude_score.png', dpi=300, bbox_inches='tight')
plt.savefig('figures/fig16_altitude_score.svg', bbox_inches='tight')
plt.savefig('docs/assets/fig16_altitude_score.png', dpi=300, bbox_inches='tight')
plt.close()
print('Saved figures/fig16_altitude_score.png/.svg')
