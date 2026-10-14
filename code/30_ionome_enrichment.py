import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['svg.hashsalt'] = 'arabidopsis-gwas-spaceflight'  # reproducible SVG ids
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact

# Enrichment of leaf-ionome GWAS genes among spaceflight DEGs, recomputed from committed tables.
#
# The original analysis (script 15, not in this repository) pooled ten AraGWAS traits as
# "ion content". Four of them are not ion traits: M216T665, M172T666 and M130T666 are LC-MS
# leaf metabolite features, and At2 is a bacterial-resistance trait (Atwell et al. 2010).
# This script reproduces its per-trait tests, checks them against tables/ion_content_overlap.csv,
# and then pools only the six ICP-MS ionome traits.
#
# Method, reproduced exactly: a gene is a spaceflight DEG if padj < 0.05 in at least one
# comparison (n_sig > 0); the universe is the 52,920-gene DE background; one-sided Fisher exact.

os.makedirs('tables', exist_ok=True)
os.makedirs('figures', exist_ok=True)
os.makedirs('docs/assets', exist_ok=True)

N_UNIVERSE = 52920
IONOME = {'Co59': 'Cobalt', 'As75': 'Arsenic', 'Mo98': 'Molybdenum', 'S34': 'Sulfur', 'Zn66': 'Zinc', 'Se82': 'Selenium'}
METABOLITE = {'M216T665': 'm/z 216, RT 665 s', 'M172T666': 'm/z 172, RT 666 s', 'M130T666': 'm/z 130, RT 666 s'}
OTHER = {'At2': 'bacterial resistance'}

g = pd.read_csv('tables/gwas_spaceflight_integration.csv', low_memory=False)
ref = pd.read_csv('tables/ion_content_overlap.csv').set_index('ion_trait')
traits_per_gene = g['gwas_traits'].fillna('').str.split(';').apply(lambda L: {t.strip() for t in L if t.strip()})
is_de = (g['n_sig'] > 0).to_numpy()
N_DE = int(round(ref['expected_overlap'].iloc[0] / ref['n_gwas_genes'].iloc[0] * N_UNIVERSE))
assert N_DE == 23013, N_DE


def test(k, n):
    table = [[k, n - k], [N_DE - k, N_UNIVERSE - N_DE - (n - k)]]
    odds, p = fisher_exact(table, alternative='greater')
    return odds, p, k / (n * N_DE / N_UNIVERSE)


rows = []
missing_non_de = {}
for trait, label in {**IONOME, **METABOLITE, **OTHER}.items():
    mask = traits_per_gene.apply(lambda s: trait in s).to_numpy()
    n, k = int(mask.sum()), int((mask & is_de).sum())
    n_ref, k_ref = int(ref.loc[trait, 'n_gwas_genes']), int(ref.loc[trait, 'n_de_overlap'])
    # The integration table lacks a gene that the original count included (Se82: 12 vs 13).
    # DE counts agree, so the missing gene was not DE; carry it as a non-DE gene.
    assert k == k_ref, (trait, k, k_ref)
    assert n <= n_ref, (trait, n, n_ref)
    missing_non_de[trait] = n_ref - n
    odds, p, fe = test(k_ref, n_ref)
    assert np.isclose(p, ref.loc[trait, 'pvalue_fisher'], rtol=1e-6), (trait, p)
    kind = 'ionome (ICP-MS)' if trait in IONOME else ('LC-MS metabolite feature' if trait in METABOLITE else 'bacterial resistance')
    rows.append(dict(trait=trait, description=label, trait_type=kind, n_gwas_genes=n_ref, n_de_overlap=k_ref,
                     fold_enrichment=fe, odds_ratio=odds, pvalue_fisher_one_sided=p))

def bh(p):
    p = np.asarray(p, dtype=float)
    order = np.argsort(p)
    scaled = p[order] * len(p) / np.arange(1, len(p) + 1)
    q = np.minimum.accumulate(scaled[::-1])[::-1]
    out = np.empty(len(p))
    out[order] = np.minimum(q, 1.0)
    return out


df = pd.DataFrame(rows)
# The committed table applied BH across eleven tests (the ten traits plus its pooled row);
# reproducing it checks bh().
ref_p = np.r_[df['pvalue_fisher_one_sided'].to_numpy(), ref.loc['ALL_ION_CONTENT', 'pvalue_fisher']]
ref_q = np.r_[ref.loc[df['trait'], 'pvalue_fdr'].to_numpy(), ref.loc['ALL_ION_CONTENT', 'pvalue_fdr']]
assert np.allclose(bh(ref_p), ref_q, rtol=1e-6)
ion = df['trait_type'] == 'ionome (ICP-MS)'
df.loc[ion, 'pvalue_fdr_ionome'] = bh(df.loc[ion, 'pvalue_fisher_one_sided'])

pool_mask = traits_per_gene.apply(lambda s: bool(s & set(IONOME))).to_numpy()
n_pool = int(pool_mask.sum()) + sum(missing_non_de[t] for t in IONOME)
k_pool = int((pool_mask & is_de).sum())
odds, p, fe = test(k_pool, n_pool)
df = pd.concat([df, pd.DataFrame([dict(trait='POOLED_6_IONOME', description='union of the six ionome traits',
                                       trait_type='ionome (ICP-MS)', n_gwas_genes=n_pool, n_de_overlap=k_pool,
                                       fold_enrichment=fe, odds_ratio=odds, pvalue_fisher_one_sided=p)])], ignore_index=True)
df.to_csv('tables/ionome_enrichment.csv', index=False)
print(df.to_string(index=False))
print(f'\nPooled six ionome traits: {k_pool}/{n_pool} DE, FE {fe:.2f}, OR {odds:.2f}, one-sided p {p:.2e}')

# Figure 18
share = pd.read_csv('tables/ion_content_locus_sharing.csv')
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6), gridspec_kw={'width_ratios': [1.25, 1]})
type_color = {'ionome (ICP-MS)': '#0072B2', 'LC-MS metabolite feature': '#E69F00', 'bacterial resistance': '#999999'}
per = df[df['trait'] != 'POOLED_6_IONOME'].iloc[::-1]
labels = [f"{r.trait} ({r.description})" for r in per.itertuples()]
ax1.barh(range(len(per)), per['fold_enrichment'], color=[type_color[t] for t in per['trait_type']], edgecolor='black', linewidth=0.5)
for i, r in enumerate(per.itertuples()):
    star = ''
    if r.trait_type == 'ionome (ICP-MS)' and r.pvalue_fdr_ionome < 0.05:
        star = ' *'
    ax1.text(r.fold_enrichment + 0.03, i, f"{r.n_de_overlap}/{r.n_gwas_genes}{star}", va='center', fontsize=9)
ax1.axvline(1.0, color='grey', linestyle='--', linewidth=1)
ax1.set_yticks(range(len(per)))
ax1.set_yticklabels(labels, fontsize=9)
ax1.set_xlim(0, 2.75)
ax1.set_xlabel('Fold enrichment among spaceflight DEGs (DE genes / expected)', fontsize=10)
ax1.set_title(f'A. Trait GWAS genes among spaceflight DEGs\nPooled six ionome traits: {k_pool}/{n_pool} DE, '
              f'FE {fe:.2f}, one-sided p = {p:.1e}', fontsize=11)
handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in type_color.values()]
ax1.legend(handles, list(type_color.keys()), loc='upper center', bbox_to_anchor=(0.5, -0.12), ncol=3, fontsize=8,
           title='Trait type  (* FDR < 0.05 among the six ionome traits)', title_fontsize=8)

share = share.iloc[::-1]
share_labels = [f"{t} ({IONOME.get(t) or METABOLITE.get(t)})" for t in share['ion_trait']]
share_colors = ['#0072B2' if t in IONOME else '#E69F00' for t in share['ion_trait']]
ax2.barh(range(len(share)), share['n_shared_with_spaceflight'], color=share_colors, edgecolor='black', linewidth=0.5)
for i, v in enumerate(share['n_shared_with_spaceflight']):
    ax2.text(v + 0.4, i, f'{v}/31', va='center', fontsize=9)
ax2.set_yticks(range(len(share)))
ax2.set_yticklabels(share_labels, fontsize=9)
ax2.set_xlim(0, 35)
ax2.set_xlabel('Scoring SNP positions shared with the trait (of 31)', fontsize=10)
ax2.set_title('B. Which traits share the 31 scoring SNPs\n(ionome traits share none)', fontsize=11)

plt.tight_layout()
plt.savefig('figures/fig18_ionome_enrichment.png', dpi=300, bbox_inches='tight')
plt.savefig('figures/fig18_ionome_enrichment.svg', bbox_inches='tight')
plt.savefig('docs/assets/fig18_ionome_enrichment.png', dpi=300, bbox_inches='tight')
plt.close()
print('Saved figures/fig18_ionome_enrichment.png/.svg and tables/ionome_enrichment.csv')
