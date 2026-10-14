import os
import pandas as pd
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['svg.hashsalt'] = 'arabidopsis-gwas-spaceflight'  # reproducible SVG ids
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Figure 15: European accessions by predicted category, coloured by predicted score.
# The original figure's script (14_geo_visualizations.py) is not in this repository, and its
# colour bar (0-1) did not match its markers. This redraws it from tables/ecotype_prediction_geo.csv
# using the manuscript's European box (35-72 N, 15 W-45 E) and Natural Earth 50 m borders.

os.makedirs('figures', exist_ok=True)
os.makedirs('docs/assets', exist_ok=True)

LAT, LON = (35, 72), (-15, 45)
eco = pd.read_csv('tables/ecotype_prediction_geo.csv').dropna(subset=['latitude', 'longitude'])
eur = eco[eco['latitude'].between(*LAT) & eco['longitude'].between(*LON)]
print(f'{len(eur)} European accessions')

proj = ccrs.PlateCarree()
fig = plt.figure(figsize=(11, 9))
ax = fig.add_subplot(1, 1, 1, projection=proj)
ax.set_extent([LON[0], LON[1], LAT[0], LAT[1]], crs=proj)
ax.add_feature(cfeature.OCEAN.with_scale('50m'), facecolor='#e0f2fe', zorder=0)
ax.add_feature(cfeature.LAND.with_scale('50m'), facecolor='#f8fafc', zorder=0)
ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.6, edgecolor='#475569')
ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.4, edgecolor='#94a3b8')

markers = {'Very Low': 'o', 'Low': 's', 'Moderate': '^', 'High': 'D', 'Very High': '*'}
vmin, vmax = eco['predicted_response_score'].min(), eco['predicted_response_score'].max()
sc = None
for cat, mk in markers.items():
    sub = eur[eur['predicted_category'] == cat]
    if sub.empty:
        continue
    sc = ax.scatter(sub['longitude'], sub['latitude'], c=sub['predicted_response_score'], cmap='viridis',
                    vmin=vmin, vmax=vmax, marker=mk, s=70 if cat == 'Very High' else 32,
                    edgecolors='black', linewidths=0.4, alpha=0.85, transform=proj, zorder=3,
                    label=f'{cat} (n={len(sub)})')
cbar = fig.colorbar(sc, ax=ax, orientation='horizontal', shrink=0.6, pad=0.06)
cbar.set_label('Predicted spaceflight response score (colour)', fontsize=11)
handles = [Line2D([], [], marker=markers[c], linestyle='', markersize=11 if c == 'Very High' else 7,
                  markerfacecolor='#94a3b8', markeredgecolor='black', label=f'{c} (n={(eur["predicted_category"] == c).sum()})')
           for c in markers if (eur['predicted_category'] == c).any()]
ax.legend(handles=handles, title='Predicted category (marker)', loc='upper right', fontsize=9, title_fontsize=9, framealpha=0.95)
gl = ax.gridlines(draw_labels=True, linestyle=':', alpha=0.5)
gl.top_labels = False
gl.right_labels = False
ax.set_title(f'European Accessions: Predicted Spaceflight Response\n({len(eur)} accessions, 35-72 N, 15 W-45 E)',
             fontsize=12, fontweight='bold')

plt.savefig('figures/fig15_europe_map.png', dpi=300, bbox_inches='tight')
plt.savefig('figures/fig15_europe_map.svg', bbox_inches='tight')
plt.savefig('docs/assets/fig15_europe_map.png', dpi=300, bbox_inches='tight')
plt.close()
print('Saved figures/fig15_europe_map.png/.svg')
