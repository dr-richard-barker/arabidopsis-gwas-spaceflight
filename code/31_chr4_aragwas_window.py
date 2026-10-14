import os
import sys
import json
import time
import urllib.request
import pandas as pd
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['svg.hashsalt'] = 'arabidopsis-gwas-spaceflight'  # reproducible SVG ids
import matplotlib.pyplot as plt

# Figure 23: the chromosome 4 window around the 29 scoring SNPs, from the AraGWAS catalog.
#
# Every association AraGWAS stores on chr4 is fetched for the three LC-MS metabolite-feature
# studies whose hits define the 29 SNPs, and for the cobalt, molybdenum and arsenic ionome
# studies (both the Atwell 2010 and the larger ionome panels), then restricted to
# 1.20-1.32 Mb. Genes in the window come from the AraGWAS genes endpoint.
# This replaces an earlier fig23 whose p-values were drawn with np.random.
#
# The fetched tables are cached in data/aragwas/; run with --refresh to fetch again.
# A failed or incomplete fetch stops the script: there is no fallback.

API = 'https://aragwas.1001genomes.org/api'
WINDOW = (1_200_000, 1_320_000)
STUDIES = {
    144: 'M216T665', 146: 'M172T666', 145: 'M130T666',   # LC-MS metabolite features
    22: 'Co59', 279: 'Co59', 73: 'Mo98', 286: 'Mo98', 67: 'As75', 283: 'As75',  # ICP-MS ionome
}
METABOLITE = {'M216T665', 'M172T666', 'M130T666'}
ASSOC_CSV = 'data/aragwas/chr4_window_associations.csv'
GENES_CSV = 'data/aragwas/chr4_window_genes.csv'

os.makedirs('data/aragwas', exist_ok=True)
os.makedirs('figures', exist_ok=True)
os.makedirs('docs/assets', exist_ok=True)


def get(url):
    req = urllib.request.Request(url, headers={'User-Agent': 'arabidopsis-gwas-spaceflight/1.0'})
    with urllib.request.urlopen(req, timeout=90) as r:
        return json.load(r)


def fetch():
    rows = []
    for pk, trait in STUDIES.items():
        # This endpoint ignores `page` (every page repeats page 1) but honours `offset`.
        offset, got, count = 0, [], None
        while True:
            d = get(f'{API}/studies/{pk}/associations/?chr=chr4&limit=500&offset={offset}')
            count = d['count']
            if not d['results']:
                break
            got += d['results']
            offset += 500
            time.sleep(0.3)
        keys = {(r['snp']['position'], r['score']) for r in got}
        if len(got) != count or len(keys) != count:
            sys.exit(f'study {pk}: fetched {len(got)} ({len(keys)} unique) of {count} associations')
        for r in got:
            assert r['study']['id'] == pk and r['snp']['chr'] == 'chr4'
            rows.append(dict(study_pk=pk, trait=trait, phenotype=r['study']['phenotype']['name'],
                             position=r['snp']['position'], score=r['score'], maf=r['maf'],
                             over_bonferroni=r['overBonferroni'],
                             bonferroni_threshold=next(t['value'] for t in r['study']['thresholds'] if t['name'] == 'bonferroni_threshold05')))
        print(f'study {pk} ({trait}): {count} chr4 associations')
    a = pd.DataFrame(rows)
    assert (a['phenotype'] == a['trait']).all()
    a = a[(a['position'] >= WINDOW[0]) & (a['position'] <= WINDOW[1])].sort_values(['trait', 'position'])
    a.to_csv(ASSOC_CSV, index=False)
    genes = get(f'{API}/genes/?chr=chr4&start={WINDOW[0]}&end={WINDOW[1]}')
    g = pd.DataFrame([dict(gene=x['name'], start=x['positions']['gte'], end=x['positions']['lte'],
                           description=(x.get('isoforms') or [{}])[0].get('shortDescription', ''))
                      for x in genes])
    g.sort_values('start').to_csv(GENES_CSV, index=False)


if '--refresh' in sys.argv or not (os.path.exists(ASSOC_CSV) and os.path.exists(GENES_CSV)):
    fetch()
assoc = pd.read_csv(ASSOC_CSV)
genes = pd.read_csv(GENES_CSV)
meta = json.load(open('tables/ecotype_prediction_metrics.json'))
scoring = sorted(int(p.split(':')[1]) for p in meta['snp_positions'] if p.startswith('4:'))
assert len(scoring) == 29, len(scoring)
print(f'{len(assoc)} associations in the window; {len(genes)} genes; {len(scoring)} scoring SNPs on chr4')
print(assoc.groupby('trait')['score'].agg(['count', 'max']).round(1).to_string())

fig, (ax, axg) = plt.subplots(2, 1, figsize=(12, 7), sharex=True, gridspec_kw={'height_ratios': [4, 1]})
colors = {'M216T665': '#E69F00', 'M172T666': '#D55E00', 'M130T666': '#CC79A7',
          'Co59': '#0072B2', 'Mo98': '#56B4E9', 'As75': '#009E73'}
# Metabolite features first; the much sparser ionome associations are drawn on top.
for trait in ['M130T666', 'M172T666', 'M216T665', 'Co59', 'Mo98', 'As75']:
    sub = assoc[assoc['trait'] == trait]
    ionome = trait not in METABOLITE
    kind = 'ionome' if ionome else 'LC-MS metabolite feature'
    ax.scatter(sub['position'] / 1e6, sub['score'], s=40 if ionome else 12, alpha=0.95 if ionome else 0.5,
               color=colors[trait], edgecolors='black' if ionome else 'none', linewidths=0.6,
               marker='D' if ionome else 'o', zorder=4 if ionome else 2, label=f'{trait} ({kind}, n={len(sub)})')
for p in scoring:
    ax.axvline(p / 1e6, color='grey', alpha=0.15, linewidth=0.8, zorder=0)
lo, hi = assoc['bonferroni_threshold'].min(), assoc['bonferroni_threshold'].max()
ax.axhline(lo, color='black', linestyle=':', linewidth=1,
           label=f'Bonferroni threshold (-log10 p {lo:.1f}-{hi:.1f}, varies by study)')
ax.set_ylabel('-log10(p) in AraGWAS', fontsize=11)
ax.set_title('Chromosome 4, 1.20-1.32 Mb: AraGWAS associations by trait\n'
             '(grey lines: the 29 chromosome-4 SNPs used for ecotype scoring)', fontsize=12)
ax.legend(fontsize=8, loc='upper right')

for i, r in enumerate(genes.itertuples()):
    y = i % 3
    axg.plot([r.start / 1e6, r.end / 1e6], [y, y], color='#334155', linewidth=4, solid_capstyle='butt')
    if r.end >= scoring[0] - 1000 and r.start <= scoring[-1] + 1000:
        axg.text((r.start + r.end) / 2e6, y + 0.35, r.gene, fontsize=7, ha='center', va='bottom')
axg.set_ylim(-0.6, 3.2)
axg.set_yticks([])
axg.set_ylabel('Genes', fontsize=10)
axg.set_xlabel('Chromosome 4 position (Mb)', fontsize=11)
axg.set_xlim(WINDOW[0] / 1e6, WINDOW[1] / 1e6)

plt.tight_layout()
plt.savefig('figures/fig23_chr4_aragwas_window.png', dpi=300, bbox_inches='tight')
plt.savefig('figures/fig23_chr4_aragwas_window.svg', bbox_inches='tight')
plt.savefig('docs/assets/fig23_chr4_aragwas_window.png', dpi=300, bbox_inches='tight')
plt.close()
print('Saved figures/fig23_chr4_aragwas_window.png/.svg')
