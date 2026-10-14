import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['svg.hashsalt'] = 'arabidopsis-gwas-spaceflight'  # reproducible SVG ids
import matplotlib.pyplot as plt
from scipy.stats import rankdata

# Figure 4, rebuilt from committed tables (the original plotting code is not in this repository).
#   A. Random forest feature importances: tables/ml_feature_importance.csv.
#   B. ROC curve of rf_probability against the consensus label: tables/ml_gene_classifier_results.csv.
#      The features include the Fisher and meta-analysis statistics that define the label, so a
#      near-perfect AUC is expected and is not evidence of independent predictive performance.
#      The table repeats a gene once per gene-symbol alias (64,212 rows, 52,901 genes); the curve
#      uses all rows, as the original figure did, and the AUC on one row per gene is also printed.
#   C. Distribution of the 982 predicted ecotype scores: tables/ecotype_prediction_geo.csv.
#      This replaces the earlier panel C, which plotted about 2,000 synthetic accessions.

os.makedirs('figures', exist_ok=True)
os.makedirs('docs/assets', exist_ok=True)


def roc_points(y, s):
    order = np.argsort(-s, kind='mergesort')
    y, s = y[order], s[order]
    distinct = np.r_[np.where(np.diff(s))[0], len(s) - 1]
    tps = np.cumsum(y)[distinct]
    fps = (distinct + 1) - tps
    return np.r_[0, fps / fps[-1]], np.r_[0, tps / tps[-1]]


def auc_rank(y, s):
    r = rankdata(s)
    n1 = y.sum()
    n0 = len(y) - n1
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


def average_precision(y, s):
    y = y[np.argsort(-s, kind='mergesort')]
    tp = np.cumsum(y)
    return float((tp / np.arange(1, len(y) + 1) * y).sum() / y.sum())


fi = pd.read_csv('tables/ml_feature_importance.csv').head(15)
ml = pd.read_csv('tables/ml_gene_classifier_results.csv', low_memory=False)
y = (ml['consensus_class'] == 'consensus_significant').to_numpy().astype(int)
s = ml['rf_probability'].to_numpy()
auc = auc_rank(y, s)
ap = average_precision(y, s)
ml1 = ml.drop_duplicates('gene_id')
auc_dedup = auc_rank((ml1['consensus_class'] == 'consensus_significant').to_numpy().astype(int), ml1['rf_probability'].to_numpy())
print(f'AUC {auc:.4f}, AP {ap:.4f} on {len(ml)} rows; AUC {auc_dedup:.4f} on {len(ml1)} unique genes')

eco = pd.read_csv('tables/ecotype_prediction_geo.csv')

fig, axes = plt.subplots(1, 3, figsize=(18, 5))

ax = axes[0]
ax.barh(range(len(fi)), fi['importance'].values, color='#0072B2', alpha=0.85)
ax.set_yticks(range(len(fi)))
ax.set_yticklabels(fi['feature'], fontsize=9)
ax.invert_yaxis()
ax.set_xlabel('Importance', fontsize=11)
ax.set_title('A. RF Feature Importance', fontsize=12)

ax = axes[1]
fpr, tpr = roc_points(y, s)
ax.plot(fpr, tpr, color='#D55E00', lw=2, label=f'AUC = {auc:.3f}')
ax.plot([0, 1], [0, 1], '--', color='grey')
ax.set_xlabel('False Positive Rate', fontsize=11)
ax.set_ylabel('True Positive Rate', fontsize=11)
ax.set_title('B. RF ROC Curve (label-defining features included)', fontsize=12)
ax.legend(loc='lower right')

ax = axes[2]
ax.hist(eco['predicted_response_score'], bins=40, color='#009E73', alpha=0.8, edgecolor='white')
ax.set_xlabel('Predicted Response Score', fontsize=11)
ax.set_ylabel('N accessions', fontsize=11)
ax.set_title(f'C. Ecotype Score Distribution ({len(eco)} accessions)', fontsize=12)

plt.tight_layout()
plt.savefig('figures/fig4_ml_results.png', dpi=300, bbox_inches='tight')
plt.savefig('figures/fig4_ml_results.svg', bbox_inches='tight')
plt.savefig('docs/assets/fig4_ml_results.png', dpi=300, bbox_inches='tight')
plt.close()
print('Saved figures/fig4_ml_results.png/.svg')
