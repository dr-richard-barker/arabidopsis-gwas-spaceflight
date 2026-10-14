import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

os.makedirs('figures', exist_ok=True)
os.makedirs('tables', exist_ok=True)

df_top10 = pd.read_csv('tables/top10_predicted_ecotypes.csv')

# Flowering times are the 1001 Genomes FT10/FT16 values in AraPheno (phenotypes 261 and 262),
# read 2026-10-13. Only accessions whose earlier profile quoted a value are filled; the rest are blank.
# Habitat, ionome and morphology descriptions were removed on 2026-10-14: they had no source, and
# two of the six quoted flowering times were wrong (Istisu-9 FT10, Xan-3 FT16).
pheno_data = [
    {'rank': 1, 'name': 'Istisu-9', 'id': 9099, 'stock_id': 'CS76953', 'flowering_time': 'FT10 100.25 d'},
    {'rank': 2, 'name': 'Xan-3', 'id': 9067, 'stock_id': 'CS78860', 'flowering_time': 'FT16 130 d'},
    {'rank': 3, 'name': 'Xan-5', 'id': 9069, 'stock_id': 'CS78861', 'flowering_time': ''},
    {'rank': 4, 'name': 'Lerik2-3', 'id': 9081, 'stock_id': 'CS77026', 'flowering_time': 'FT10 91.25 d'},
    {'rank': 5, 'name': 'Anz-0', 'id': 9759, 'stock_id': 'CS76439', 'flowering_time': 'FT16 58.5 d'},
    {'rank': 6, 'name': 'Oy-0', 'id': 7288, 'stock_id': 'CS77156', 'flowering_time': 'FT16 70.25 d'},
    {'rank': 7, 'name': 'Ty-1', 'id': 5784, 'stock_id': 'CS78790', 'flowering_time': ''},
    {'rank': 8, 'name': 'Mc-1', 'id': 5757, 'stock_id': 'CS78785', 'flowering_time': 'FT10 87.75 d'},
    {'rank': 9, 'name': 'Sr:3', 'id': 6086, 'stock_id': 'CS77267', 'flowering_time': ''},
    {'rank': 10, 'name': 'Monte-1', 'id': 9966, 'stock_id': 'CS76361', 'flowering_time': ''}
]

df_pheno = pd.DataFrame(pheno_data)
df_merged = pd.merge(df_top10, df_pheno[['rank', 'flowering_time']], on='rank')
df_merged.to_csv('tables/top10_phenotypic_profiles.csv', index=False)
print("Saved tables/top10_phenotypic_profiles.csv")

# ----------------------------------------------------
# Plot Figure 28: Top 10 Ecotypes Geographic Map
# ----------------------------------------------------
fig, ax = plt.subplots(figsize=(12, 8))

# Define map extent over Western Eurasia / Europe / Caucasus
lats = df_top10['latitude'].values
lons = df_top10['longitude'].values
names = df_top10['name'].values
ranks = df_top10['rank'].values
countries = df_top10['country'].values
elevs = df_top10['elevation_m'].values
groups = df_top10['group'].values

# Custom palette for population groups
group_colors = {
    'italy_balkan_caucasus': '#d62728',
    'admixed': '#1f77b4',
    'south_sweden': '#2ca02c'
}

# Plot points
for i, row in df_top10.iterrows():
    c = group_colors.get(row['group'], '#9467bd')
    ax.scatter(row['longitude'], row['latitude'], s=220, color=c, zorder=5, edgecolors='black', linewidth=1.5)
    
    # Label with rank and name
    label_text = f"#{row['rank']} {row['name']}\n({row['country']}, {int(row['elevation_m'])}m)"
    
    # Adjust annotation offsets to avoid collision
    x_offset = 1.2
    y_offset = 0.5
    if row['name'] in ['Xan-3', 'Xan-5']:
        y_offset = -1.2
    elif row['name'] in ['Istisu-9', 'Lerik2-3']:
        y_offset = 1.2
    elif row['name'] == 'Sr:3':
        x_offset = -5.5
        
    ax.annotate(label_text, (row['longitude'], row['latitude']),
                xytext=(row['longitude'] + x_offset, row['latitude'] + y_offset),
                fontsize=9.5, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.85, edgecolor=c),
                arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.1', color=c, lw=1.2))

# Background grid and borders setup
ax.set_xlim(-15, 60)
ax.set_ylim(32, 68)
ax.set_xlabel('Longitude (°E)', fontsize=12, fontweight='bold')
ax.set_ylabel('Latitude (°N)', fontsize=12, fontweight='bold')
ax.set_title('Geographic Origin of the Top 10 Prioritized Spaceflight-Responsive Accessions\n(Annotated by Rank, Country, Elevation, and Population Group)', fontsize=13, fontweight='bold')

# Add legend for population groups
for grp, color in group_colors.items():
    ax.scatter([], [], color=color, s=150, edgecolors='black', label=grp.replace('_', ' ').title())

ax.legend(title='Population Group', loc='upper left', frameon=True, fontsize=10, title_fontsize=11)
ax.grid(True, linestyle=':', alpha=0.6)

plt.tight_layout()
plt.savefig('figures/fig28_top10_ecotypes_map.png', dpi=300)
plt.savefig('figures/fig28_top10_ecotypes_map.svg')
plt.close()

print("Saved figures/fig28_top10_ecotypes_map.png")
