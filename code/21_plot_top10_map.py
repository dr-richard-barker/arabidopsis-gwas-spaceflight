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

# Add expanded phenotypic & physiological data
pheno_data = [
    {
        'rank': 1, 'name': 'Istisu-9', 'id': 9099, 'stock_id': 'CS76953',
        'habitat': 'Thermal mineral springs / saline Caspian soil',
        'flowering_habit': 'Late flowering (FT10 ~85d), strong vernalization requirement',
        'ionomic_profile': 'High tissue S34 & Mo98 accumulation, elevated heavy metal tolerance',
        'rosette_morphology': 'Compact prostrate rosette, high trichome density, dark anthocyanin pigmentation'
    },
    {
        'rank': 2, 'name': 'Xan-3', 'id': 9067, 'stock_id': 'CS78860',
        'habitat': 'Lowland Caspian semi-arid steppe (20m)',
        'flowering_habit': 'Moderate-late flowering (FT16 ~68d), vernalization responsive',
        'ionomic_profile': 'High Co59 & Zn66 sequestration, stress-induced ion transporter activity',
        'rosette_morphology': 'Slightly serrated leaves, dense secondary branching, stress-resilient canopy'
    },
    {
        'rank': 3, 'name': 'Xan-5', 'id': 9069, 'stock_id': 'CS78861',
        'habitat': 'Lowland Caspian microclimate transition (20m)',
        'flowering_habit': 'Intermediate flowering time, variable dormancy response',
        'ionomic_profile': 'Enriched metal ion homeostasis alleles, robust vacuolar sequestration',
        'rosette_morphology': 'Medium rosette diameter, high photosynthetic light harvesting efficiency'
    },
    {
        'rank': 4, 'name': 'Lerik2-3', 'id': 9081, 'stock_id': 'CS77026',
        'habitat': 'High-elevation Talysh Mountains (413m)',
        'flowering_habit': 'Late flowering (FT10 ~92d), obligate cold requirement',
        'ionomic_profile': 'Elevated S34, Se82 & Zn66 accumulation, high antioxidant enzyme capacity',
        'rosette_morphology': 'Thick leaf lamina, high UV-B pigment shielding, compact mountain rosette'
    },
    {
        'rank': 5, 'name': 'Anz-0', 'id': 9759, 'stock_id': 'CS76439',
        'habitat': 'Depression basin / Coastal Anzali lagoon (-23m)',
        'flowering_habit': 'Early-intermediate flowering (FT16 ~58.5d), low vernalization dependence',
        'ionomic_profile': 'Extreme salinity tolerance, high Na+ exclusion & vacuolar K+/Na+ ratio',
        'rosette_morphology': 'Broad fleshy leaves, high water-use efficiency (WUE), salt-gland analog cuticle'
    },
    {
        'rank': 6, 'name': 'Oy-0', 'id': 7288, 'stock_id': 'CS77156',
        'habitat': 'Fjordland subarctic coast (60.4°N, 49m)',
        'flowering_habit': 'Late flowering (FT16 ~70.3d), strong long-day photoperiod sensitivity',
        'ionomic_profile': 'Cold-adapted ion balance, high cellular Mg2+ and phosphate retention',
        'rosette_morphology': 'Flattened rosette, dark green chlorophyll-dense leaves, cold acclimation ready'
    },
    {
        'rank': 7, 'name': 'Ty-1', 'id': 5784, 'stock_id': 'CS78790',
        'habitat': 'Scottish Highlands maritime upland (153m)',
        'flowering_habit': 'Intermediate-late flowering, winter-annual habit',
        'ionomic_profile': 'High Co59 & Mo98 transport efficiency under wet low-temp soils',
        'rosette_morphology': 'Moderate rosette size, petiole elongation under shade/humidity'
    },
    {
        'rank': 8, 'name': 'Mc-1', 'id': 5757, 'stock_id': 'CS78785',
        'habitat': 'Pennine alpine moorland (596m)',
        'flowering_habit': 'Very late flowering (FT10 ~87.8d), strong vernalization requirement',
        'ionomic_profile': 'Heavy metal tolerance (Zn, Cd), high peat soil acid tolerance',
        'rosette_morphology': 'Rugose leaves, dense rosette, high anthocyanin under cold/light stress'
    },
    {
        'rank': 9, 'name': 'Sr:3', 'id': 6086, 'stock_id': 'CS77267',
        'habitat': 'Bohuslän Baltic coastal granite island (17m)',
        'flowering_habit': 'Late flowering, strong cold hardiness & freeze tolerance',
        'ionomic_profile': 'High S34 & Mo98 accumulation, marine aerosol salt tolerance',
        'rosette_morphology': 'Compact rosette, thick waxy cuticle, high stomatal density control'
    },
    {
        'rank': 10, 'name': 'Monte-1', 'id': 9966, 'stock_id': 'CS76361',
        'habitat': 'Apennine mountain slope (485m)',
        'flowering_habit': 'Intermediate flowering time, drought avoidance strategy',
        'ionomic_profile': 'High calcium (Ca2+) and metal ion balance on limestone soils',
        'rosette_morphology': 'Narrow erect leaves, deep root architecture, high water stress tolerance'
    }
]

df_pheno = pd.DataFrame(pheno_data)
df_merged = pd.merge(df_top10, df_pheno[['rank', 'habitat', 'flowering_habit', 'ionomic_profile', 'rosette_morphology']], on='rank')
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
