import os
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

os.makedirs('figures', exist_ok=True)
os.makedirs('docs/assets', exist_ok=True)

df = pd.read_csv('tables/top10_predicted_ecotypes.csv')

# 1. Bar Plot: Score & SNP Burden
fig, ax = plt.subplots(figsize=(10, 5))
colors = ['#2b6cb0' if s > 0.65 else '#319795' for s in df['predicted_response_score']]
bars = ax.bar([f"#{r} {n}" for r, n in zip(df['rank'], df['name'])], df['predicted_response_score'], color=colors, edgecolor='black', linewidth=1)

ax.set_ylim(0.635, 0.658)
ax.set_ylabel('Predicted Response Score', fontsize=11, fontweight='bold')
ax.set_xlabel('Accession Rank & Name', fontsize=11, fontweight='bold')
ax.set_title('Predicted Spaceflight Response Score across Top 10 Accessions', fontsize=12, fontweight='bold')
plt.xticks(rotation=25, ha='right', fontweight='bold')

for bar, pct in zip(bars, df['response_percentile']):
    yval = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2.0, yval + 0.0005, f"{pct:.1f}%", ha='center', va='bottom', fontsize=9, fontweight='bold')

ax.grid(True, linestyle=':', alpha=0.5, axis='y')
plt.tight_layout()
plt.savefig('figures/fig29_score_bar_plot.png', dpi=300)
plt.savefig('docs/assets/fig29_score_bar_plot.png', dpi=300)
plt.close()

# 2. Scatter Plot: Latitude vs Elevation
fig, ax = plt.subplots(figsize=(10, 5))
scatter = ax.scatter(df['latitude'], df['elevation_m'], c=df['predicted_response_score'], cmap='viridis', s=200, edgecolors='black', linewidth=1.5, zorder=5)
cbar = plt.colorbar(scatter, ax=ax)
cbar.set_label('Predicted Response Score', fontweight='bold')

label_offsets = {
    'Istisu-9': (-1.5, 25), 'Xan-3': (-1.5, -35), 'Xan-5': (0.8, -35),
    'Lerik2-3': (0.8, 25), 'Anz-0': (0.8, -35), 'Oy-0': (-1.5, 25),
    'Ty-1': (-1.5, 25), 'Mc-1': (-1.5, -35), 'Sr:3': (0.8, 25), 'Monte-1': (-1.5, -35)
}

for i, row in df.iterrows():
    xo, yo = label_offsets.get(row['name'], (0.5, 10))
    ax.annotate(f"#{row['rank']} {row['name']}", (row['latitude'], row['elevation_m']),
                xytext=(row['latitude'] + xo, row['elevation_m'] + yo),
                fontsize=9, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.85, edgecolor='gray'),
                arrowprops=dict(arrowstyle='->', lw=1, color='gray'), zorder=10)

ax.set_xlabel('Latitude (°N)', fontsize=11, fontweight='bold')
ax.set_ylabel('Elevation (m)', fontsize=11, fontweight='bold')
ax.set_title('Latitudinal (°N) & Altitudinal (m) Distribution of Top 10 Accessions', fontsize=12, fontweight='bold')
ax.grid(True, linestyle=':', alpha=0.5)

plt.tight_layout()
plt.savefig('figures/fig30_latitude_elevation_scatter.png', dpi=300)
plt.savefig('docs/assets/fig30_latitude_elevation_scatter.png', dpi=300)
plt.close()

print("Successfully generated static dashboard chart fallbacks!")
