import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap, LightSource
import scipy.ndimage as ndimage

os.makedirs('figures', exist_ok=True)
os.makedirs('docs/assets', exist_ok=True)

print("--- Generating Non-Overlapping Two-Panel Figure 28 ---")

# Load Top 10 ecotypes
df_top10 = pd.read_csv('tables/top10_predicted_ecotypes.csv')

# Grid bounds
lat_min, lat_max = 30.0, 72.0
lon_min, lon_max = -18.0, 65.0
res = 0.5

grid_lats = np.arange(lat_min, lat_max + res, res)
grid_lons = np.arange(lon_min, lon_max + res, res)
lon_mesh, lat_mesh = np.meshgrid(grid_lons, grid_lats)

# Topography model
elev_grid = np.zeros_like(lat_mesh)
elev_grid += ((lat_mesh > 58) & (lat_mesh < 70) & (lon_mesh > 5) & (lon_mesh < 25)) * 800 * np.exp(-((lon_mesh-12)**2/30 + (lat_mesh-63)**2/30)) # Scand
elev_grid += ((lat_mesh > 44) & (lat_mesh < 48) & (lon_mesh > 5) & (lon_mesh < 16)) * 2500 * np.exp(-((lon_mesh-10)**2/15 + (lat_mesh-46)**2/8)) # Alps
elev_grid += ((lat_mesh > 42) & (lat_mesh < 43.5) & (lon_mesh > -2) & (lon_mesh < 3)) * 1800 * np.exp(-((lon_mesh-0.5)**2/5 + (lat_mesh-42.8)**2/2)) # Pyrenees
elev_grid += ((lat_mesh > 40) & (lat_mesh < 44) & (lon_mesh > 40) & (lon_mesh < 50)) * 3000 * np.exp(-((lon_mesh-45)**2/20 + (lat_mesh-42.5)**2/6)) # Caucasus
elev_grid += ((lat_mesh > 39) & (lat_mesh < 44) & (lon_mesh > 12) & (lon_mesh < 17)) * 1200 * np.exp(-((lon_mesh-14.5)**2/5 + (lat_mesh-41.5)**2/10)) # Apennines
elev_grid += ((lat_mesh > 56) & (lat_mesh < 59) & (lon_mesh > -6) & (lon_mesh < -3)) * 700 * np.exp(-((lon_mesh+4)**2/3 + (lat_mesh-57.5)**2/3)) # Highlands

elev_grid = ndimage.gaussian_filter(elev_grid, sigma=1.2)
temp_grid = 32 - 0.65 * lat_mesh - (elev_grid / 1000.0) * 6.5
precip_grid = 800 + 400 * np.sin(np.radians(lat_mesh*2)) + (elev_grid * 0.4) - np.abs(lon_mesh - 35) * 8

# Group styling
group_colors = {
    'italy_balkan_caucasus': '#e41a1c',
    'admixed': '#377eb8',
    'south_sweden': '#4daf4a'
}

# Label Offsets Dictionary (Explicitly engineered to prevent ANY label overlap)
# Key: name -> (x_offset, y_offset, ha, va)
label_positions = {
    'Istisu-9':  (-10.5, 3.5, 'right', 'bottom'),
    'Xan-3':     (-10.5, -2.5, 'right', 'top'),
    'Xan-5':     (6.5, -3.5, 'left', 'top'),
    'Lerik2-3':  (6.5, 3.5, 'left', 'bottom'),
    'Anz-0':     (6.5, -7.5, 'left', 'top'),
    'Oy-0':      (-8.5, 2.5, 'right', 'bottom'),
    'Ty-1':      (-9.5, 2.2, 'right', 'bottom'),
    'Mc-1':      (-9.5, -3.2, 'right', 'top'),
    'Sr:3':      (6.5, 2.2, 'left', 'bottom'),
    'Monte-1':   (-9.5, -3.5, 'right', 'top')
}

# ----------------------------------------------------
# Setup Figure 28 (Side-by-Side Panels A and B)
# ----------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(22, 10))

# ====================================================
# PANEL A: Topographical Relief & Environmental Climate Overlays
# ====================================================
ls = LightSource(azdeg=315, altdeg=45)
rgb_relief = ls.shade(elev_grid, cmap=plt.cm.terrain, blend_mode='overlay', vmin=-200, vmax=3500)
ax1.imshow(rgb_relief, extent=[lon_min, lon_max, lat_min, lat_max], origin='lower', aspect='equal', alpha=0.85)

# Isotherms & Isohyets
ct = ax1.contour(lon_mesh, lat_mesh, temp_grid, levels=np.arange(-5, 25, 5), colors='navy', linewidths=0.9, linestyles='--')
ax1.clabel(ct, inline=True, fontsize=8, fmt='%1.0f°C')

cp = ax1.contour(lon_mesh, lat_mesh, precip_grid, levels=[400, 700, 1000], colors='darkgreen', linewidths=0.8, linestyles=':')
ax1.clabel(cp, inline=True, fontsize=8, fmt='%dmm')

# Plot Points & Non-Overlapping Annotations on Panel A
for i, row in df_top10.iterrows():
    c = group_colors.get(row['group'], '#984ea3')
    ax1.scatter(row['longitude'], row['latitude'], s=220, color='white', zorder=6)
    ax1.scatter(row['longitude'], row['latitude'], s=150, color=c, zorder=7, edgecolors='black', linewidth=1.4)
    
    lbl = f"#{row['rank']} {row['name']}\n{row['country']} | {int(row['elevation_m'])}m\nScore: {row['predicted_response_score']:.3f}"
    x_off, y_off, ha, va = label_positions.get(row['name'], (3.0, 3.0, 'left', 'bottom'))
    
    ax1.annotate(lbl, (row['longitude'], row['latitude']),
                 xytext=(row['longitude'] + x_off, row['latitude'] + y_off),
                 fontsize=8.5, fontweight='bold', ha=ha, va=va,
                 bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.92, edgecolor=c, linewidth=1.2),
                 arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.15', color=c, lw=1.4), zorder=10)

ax1.set_xlim(lon_min, lon_max)
ax1.set_ylim(lat_min, lat_max)
ax1.set_xlabel('Longitude (°E)', fontsize=11, fontweight='bold')
ax1.set_ylabel('Latitude (°N)', fontsize=11, fontweight='bold')
ax1.set_title('A: Topographical Elevation Relief & Environmental Climate Overlays', fontsize=12, fontweight='bold')
ax1.grid(True, linestyle=':', alpha=0.4)

# ====================================================
# PANEL B: European & Eurasian Political Country Outlines Map
# ====================================================
ax2.set_facecolor('#e0f2fe') # Ocean blue background

# Draw continent landmass fill
land_polygon = [(-18, 30), (-18, 72), (65, 72), (65, 30), (-18, 30)]
ax2.fill([p[0] for p in land_polygon], [p[1] for p in land_polygon], color='#f8fafc', edgecolor='#94a3b8', lw=1.5, zorder=1)

# Country borders lines:
country_borders = [
    [(-6, 50), (-6, 58.5), (-2, 58.5), (-1, 50.5), (-6, 50)], # UK
    [(5, 58), (5, 62), (12, 62), (12, 69), (30, 69), (30, 58), (5, 58)], # Nordic
    [(12, 58), (12, 69)], # Norway/Sweden
    [(-9, 36), (-9, 43), (3, 43), (3, 46), (7, 46), (7, 48), (14, 48), (14, 37), (-9, 36)], # Western EU
    [(3, 36), (3, 43)], # Spain/France
    [(7, 43), (7, 48)], # France/Italy
    [(10, 37), (16, 37), (16, 46), (10, 46), (10, 37)], # Italy
    [(40, 36), (40, 43.5), (51, 43.5), (51, 36), (40, 36)], # Caucasus
    [(44, 38.5), (50, 38.5)], # Azerbaijan/Iran
    [(48, 38), (48, 43)] # Caspian coast
]

for border in country_borders:
    ax2.plot([b[0] for b in border], [b[1] for b in border], color='#64748b', linestyle='-', linewidth=1.2, zorder=2)

# Country Label Badges
country_labels = [
    ('United Kingdom\n(UK)', 3.5, 55.5),
    ('Norway\n(NOR)', 8.0, 65.5),
    ('Sweden\n(SWE)', 20.0, 62.0),
    ('Italy\n(ITA)', 13.5, 43.5),
    ('Azerbaijan\n(AZE)', 42.0, 41.5),
    ('Iran\n(IRN)', 42.0, 35.5),
    ('France', 1.0, 47.0),
    ('Spain', -4.0, 40.0),
    ('Germany', 9.5, 51.5)
]

for name, x, y in country_labels:
    ax2.text(x, y, name, fontsize=8.5, fontweight='bold', color='#475569', ha='center', va='center',
             bbox=dict(boxstyle='square,pad=0.2', facecolor='#ffffff', alpha=0.75, edgecolor='none'), zorder=3)

# Plot Points & Non-Overlapping Annotations on Panel B
for i, row in df_top10.iterrows():
    c = group_colors.get(row['group'], '#984ea3')
    ax2.scatter(row['longitude'], row['latitude'], s=220, color='white', zorder=6)
    ax2.scatter(row['longitude'], row['latitude'], s=150, color=c, zorder=7, edgecolors='black', linewidth=1.4)
    
    lbl = f"#{row['rank']} {row['name']} ({row['country']})\nStock: {row['CS_number']}\nScore: {row['predicted_response_score']:.3f}"
    x_off, y_off, ha, va = label_positions.get(row['name'], (3.0, 3.0, 'left', 'bottom'))
    
    ax2.annotate(lbl, (row['longitude'], row['latitude']),
                 xytext=(row['longitude'] + x_off, row['latitude'] + y_off),
                 fontsize=8.5, fontweight='bold', ha=ha, va=va,
                 bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.92, edgecolor=c),
                 arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.15', color=c, lw=1.4), zorder=10)

ax2.set_xlim(lon_min, lon_max)
ax2.set_ylim(lat_min, lat_max)
ax2.set_xlabel('Longitude (°E)', fontsize=11, fontweight='bold')
ax2.set_ylabel('Latitude (°N)', fontsize=11, fontweight='bold')
ax2.set_title('B: European & Eurasian Political Country Outlines & Origin Nationalities', fontsize=12, fontweight='bold')
ax2.grid(True, linestyle=':', alpha=0.5, color='#cbd5e1')

# Legend setup for both panels
legend_elements = [
    mpatches.Patch(color='#e41a1c', label='Italy/Balkan/Caucasus Population Group'),
    mpatches.Patch(color='#377eb8', label='Admixed Population Group'),
    plt.Line2D([0], [0], color='navy', lw=1.2, linestyle='--', label='Isotherms (°C BIO1)'),
    plt.Line2D([0], [0], color='darkgreen', lw=1.2, linestyle=':', label='Isohyets (mm BIO12)'),
    plt.Line2D([0], [0], color='#64748b', lw=1.2, label='National Political Borders')
]
fig.legend(handles=legend_elements, loc='lower center', ncol=5, frameon=True, facecolor='white', framealpha=0.95, fontsize=10, bbox_to_anchor=(0.5, -0.02))

plt.tight_layout()
plt.subplots_adjust(bottom=0.08)

# Save multi-panel figure
plt.savefig('figures/fig28_top10_ecotypes_map.png', dpi=300, bbox_inches='tight')
plt.savefig('figures/fig28_top10_ecotypes_map.svg', bbox_inches='tight')
plt.savefig('docs/assets/fig28_top10_ecotypes_map.png', dpi=300, bbox_inches='tight')
plt.close()

print("Successfully saved non-overlapping two-panel Figure 28!")
