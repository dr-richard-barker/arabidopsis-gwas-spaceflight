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
fig, ax = plt.subplots(figsize=(11, 6))
scatter = ax.scatter(df['latitude'], df['elevation_m'], c=df['predicted_response_score'], cmap='viridis', s=200, edgecolors='black', linewidth=1.5, zorder=5)
cbar = plt.colorbar(scatter, ax=ax)
cbar.set_label('Predicted Response Score', fontweight='bold')

# Accessions at identical coordinates (Xan-3 and Xan-5) share one label.
# Label positions are fixed in data coordinates so that no two labels overlap.
label_pos = {
    'Monte-1': (43.0, 570), 'Lerik2-3': (34.8, 540), 'Istisu-9': (43.0, 200),
    'Xan-3': (43.0, -70), 'Anz-0': (34.8, -150),
    'Mc-1': (51.0, 690), 'Ty-1': (52.5, 300), 'Sr:3': (56.0, -150), 'Oy-0': (61.0, 230),
}
for (lat, elev), grp in df.groupby(['latitude', 'elevation_m'], sort=False):
    label = ' / '.join(f"#{r} {n}" for r, n in zip(grp['rank'], grp['name']))
    tx, ty = label_pos[grp['name'].iloc[0]]
    ax.annotate(label, (lat, elev), xytext=(tx, ty), ha='center', va='center',
                fontsize=9, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.9, edgecolor='gray'),
                arrowprops=dict(arrowstyle='->', lw=1, color='gray'), zorder=10)

ax.set_xlim(32, 64)
ax.set_ylim(-230, 780)
ax.set_xlabel('Latitude (°N)', fontsize=11, fontweight='bold')
ax.set_ylabel('Elevation (m)', fontsize=11, fontweight='bold')
ax.set_title('Latitude and Elevation of the Top 10 Accessions', fontsize=12, fontweight='bold')
ax.grid(True, linestyle=':', alpha=0.5)

plt.tight_layout()
plt.savefig('figures/fig30_latitude_elevation_scatter.png', dpi=300)
plt.savefig('docs/assets/fig30_latitude_elevation_scatter.png', dpi=300)
plt.close()

print("Successfully generated static dashboard chart fallbacks!")
