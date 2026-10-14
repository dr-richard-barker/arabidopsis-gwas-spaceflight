import os
import sys
import json
import urllib.request
import pandas as pd

# Flowering times of the top 10 accessions from AraPheno (Seren et al. 2017):
# the 1001 Genomes FT10 and FT16 phenotypes (days to flowering at 10 and 16 degC),
# AraPheno phenotypes 261 and 262. Writes tables/top10_flowering_times.csv, which
# 21_plot_top10_map.py reads. A failed fetch stops the script: there is no fallback.

PHENOTYPES = {'FT10': 261, 'FT16': 262}
OUT = 'tables/top10_flowering_times.csv'

os.makedirs('tables', exist_ok=True)
top10 = pd.read_csv('tables/top10_predicted_ecotypes.csv')

values = {}
for name, pk in PHENOTYPES.items():
    url = f'https://arapheno.1001genomes.org/rest/phenotype/{pk}/values.json'
    req = urllib.request.Request(url, headers={'User-Agent': 'arabidopsis-gwas-spaceflight/1.0'})
    with urllib.request.urlopen(req, timeout=90) as r:
        rows = json.load(r)
    if not rows or rows[0].get('phenotype_name') != name:
        sys.exit(f'AraPheno phenotype {pk}: unexpected response for {name}')
    values[name] = {int(x['accession_id']): float(x['phenotype_value']) for x in rows}
    print(f'{name} (AraPheno {pk}): {len(values[name])} accessions')

out = top10[['rank', 'name', 'accession_id', 'CS_number']].copy()
for name in PHENOTYPES:
    out[f'{name}_days'] = out['accession_id'].map(values[name])
out['source'] = 'AraPheno phenotypes 261 (FT10) and 262 (FT16), 1001 Genomes; blank = not phenotyped'
out.to_csv(OUT, index=False)
print(out.to_string(index=False))
