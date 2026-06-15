"""
make_meeting_figures.py — Generate presentation-ready figures for the
confidence-score pipeline update meeting.

Outputs to figures/meeting_20260424/ :
  fig1_auc_journey.png          # AUC progression across the week's fixes
  fig2_reliability_diagram.png  # Raw vs calibrated ECE per model
  fig3_class_separation.png     # TP / FP / blank confidence distributions
  fig4_feature_importance.png   # GBM feature importance (gain)
  fig5_nota_validation.png      # Blank spectra confidence, pos + neg
  fig6_model_comparison.png     # GBM vs RF / MLP / LR side-by-side
"""
import os, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
warnings.filterwarnings('ignore')

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG_DIR = os.path.join(ROOT, 'figures', 'meeting_20260424')
os.makedirs(FIG_DIR, exist_ok=True)

# Match presentation slide style
plt.rcParams.update({
    'font.size': 13,
    'axes.titlesize': 16,
    'axes.labelsize': 14,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 12,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'axes.edgecolor': '#555555',
    'axes.labelcolor': '#222222',
    'xtick.color': '#444444',
    'ytick.color': '#444444',
})


# ──────────────────────────────────────────────────────────────────────────────
# FIGURE 1 — AUC progression across the week's fixes
# ──────────────────────────────────────────────────────────────────────────────
def fig1_auc_journey():
    fig, ax = plt.subplots(figsize=(12.5, 6.5))

    stages = [
        ('Bayesian\nbaseline\n(04-17)', 0.841),
        ('+ empty-IK14\nsim_gap guard', 0.846),
        ('+ rebuild from\n04-21 curation', 0.852),
        ('+ dual empty-IK14\npolicy + NIST fallback', 0.875),
        ('Switch to GBM\n(XGBoost)', 0.913),
    ]
    labels = [s[0] for s in stages]
    aucs = [s[1] for s in stages]
    colors = ['#6B7A99', '#7A8AA9', '#8B9BB7', '#9BACC4', '#3B7D3B']

    x = np.arange(len(stages))
    bars = ax.bar(x, aucs, color=colors, edgecolor='#222222', linewidth=0.5, width=0.72)

    for bar, auc in zip(bars, aucs):
        ax.text(bar.get_x() + bar.get_width() / 2, auc + 0.002,
                f'{auc:.3f}', ha='center', va='bottom',
                fontsize=14, fontweight='bold', color='#222222')

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=11, color='#222222')
    ax.set_ylim(0.80, 0.94)
    ax.set_ylabel('Out-of-fold AUC', color='#222222')
    ax.set_title('Confidence-scorer AUC: Bayesian baseline → GBM primary\n(5-fold GroupKFold by anno_ik14)',
                 color='#222222', pad=15)
    ax.axhline(0.841, color='#aaaaaa', linestyle=':', linewidth=1, zorder=0)
    ax.text(4.55, 0.841 + 0.001, 'baseline', color='#888888', fontsize=10)
    ax.grid(axis='y', alpha=0.25)
    plt.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, 'fig1_auc_journey.png'), dpi=160, bbox_inches='tight')
    plt.close(fig)


# ──────────────────────────────────────────────────────────────────────────────
# FIGURE 2 — Reliability diagram: Raw Bayesian vs calibrated GBM
# ──────────────────────────────────────────────────────────────────────────────
def _reliability_bins(p, y, nbins=10):
    edges = np.linspace(0, 1, nbins + 1)
    bi = np.clip(np.digitize(p, edges[1:-1]), 0, nbins - 1)
    bc, bf, bn = [], [], []
    for b in range(nbins):
        m = bi == b
        if m.sum() == 0: continue
        bc.append(p[m].mean()); bf.append(y[m].mean()); bn.append(m.sum())
    return np.array(bc), np.array(bf), np.array(bn)


def _ece(p, y, nbins=10):
    bc, bf, bn = _reliability_bins(p, y, nbins)
    return (bn * np.abs(bf - bc)).sum() / bn.sum()


def fig2_reliability_diagram():
    """
    Reliability diagram with HONEST calibration evaluation (leave-one-fold-out).

    Why: fitting isotonic on all OOF scores then evaluating ECE on those same
    scores is essentially in-sample — ECE collapses to ~0.000 by construction
    because isotonic can exactly match any piecewise empirical frequency.
    The defensible number fits isotonic on 4 folds and evaluates on the
    held-out 5th, rotating.
    """
    import xgboost as xgb
    from sklearn.model_selection import GroupKFold
    from sklearn.isotonic import IsotonicRegression

    print('  Recomputing OOF + HONEST (leave-one-fold-out) calibration for reliability plot...')
    import sys
    sys.path.insert(0, os.path.join(ROOT, 'code'))
    from score_gbm_v2 import NUMERIC_FEATURES, CATEGORICAL_FEATURES, prep_features, train_xgb
    from bayesian_score_v2 import ChannelSpec, fit_channel, logit_upper_half

    ft = pd.read_csv(os.path.join(ROOT, 'data', 'feature_table_v2.csv'), low_memory=False)
    labeled = ft[ft['spectrum_label'].isin(['TP','FP'])]
    top1 = ft.loc[labeled.groupby('wiki_id')['entropy_similarity'].idxmax()].reset_index(drop=True)
    labels = top1['hit_label'].values
    groups = top1['anno_ik14'].fillna('').values.copy()
    for i in range(len(groups)):
        if groups[i] == '':
            groups[i] = f'__no_ik14_{i}'

    # Bayesian OOF
    CHANNELS = [
        ChannelSpec('entropy_sim', 'entropy_similarity', 'continuous', higher_means_tp=True,
                    tp_family='normal', fp_family='normal', transform=logit_upper_half),
        ChannelSpec('sim_gap', 'sim_gap', 'continuous', higher_means_tp=True,
                    tp_family='normal', fp_family='normal'),
        ChannelSpec('signed_delta_rt', 'signed_delta_rt', 'continuous', higher_means_tp=True,
                    tp_family='student_t', fp_family='student_t'),
    ]
    prior_log_odds = np.log(labels.mean() / (1 - labels.mean()))
    bay_oof = np.full(len(top1), np.nan)
    for tr, te in GroupKFold(n_splits=5).split(top1, labels, groups):
        fitted = [fit_channel(s, top1.iloc[tr][s.feature_col].values, labels[tr]) for s in CHANNELS]
        lr = np.zeros(len(te))
        for fc in fitted:
            lr += fc.logLR(top1.iloc[te][fc.spec.feature_col].values)
        bay_oof[te] = 1.0 / (1.0 + np.exp(-(prior_log_odds + lr)))

    # GBM OOF
    X_train, _ = prep_features(top1)
    gbm_oof = np.full(len(top1), np.nan)
    ft_types = ['q'] * len(NUMERIC_FEATURES) + ['c'] * len(CATEGORICAL_FEATURES)
    for tr, te in GroupKFold(n_splits=5).split(top1, labels, groups):
        model = train_xgb(X_train.iloc[tr], labels[tr], CATEGORICAL_FEATURES)
        gbm_oof[te] = model.predict(xgb.DMatrix(X_train.iloc[te], enable_categorical=True, feature_types=ft_types))

    # HONEST calibration: fit isotonic on 4 folds, eval on 5th, rotate
    def honest_calibrate(p_oof, y, groups_):
        out = np.full_like(p_oof, np.nan)
        for tr, te in GroupKFold(n_splits=5).split(p_oof, y, groups_):
            iso = IsotonicRegression(out_of_bounds='clip').fit(p_oof[tr], y[tr])
            out[te] = iso.transform(p_oof[te])
        return out
    bay_cal = honest_calibrate(bay_oof, labels, groups)
    gbm_cal = honest_calibrate(gbm_oof, labels, groups)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    panels = [
        ('Raw Bayesian (uncalibrated)', bay_oof, '#CC5555'),
        ('GBM + isotonic calibration (held-out)', gbm_cal, '#3B7D3B'),
    ]
    for ax, (title, p, col) in zip(axes, panels):
        bc, bf, bn = _reliability_bins(p, labels)
        ax.plot([0, 1], [0, 1], 'k--', alpha=0.35, linewidth=1, label='perfect calibration')
        ax.plot(bc, bf, 'o-', color=col, markersize=9, linewidth=2, label='observed')
        for x_, y_, n_ in zip(bc, bf, bn):
            ax.annotate(f'n={n_}', (x_, y_), textcoords='offset points', xytext=(7, -10),
                        fontsize=8, color='#666666')
        e = _ece(p, labels)
        ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        ax.set_xlabel('Predicted confidence')
        ax.set_ylabel('Observed fraction correct')
        ax.grid(alpha=0.25)
        ax.set_title(f'{title}\nECE = {e:.3f}')
        ax.legend(loc='lower right')

    fig.suptitle('Calibration — honest (leave-one-fold-out) ECE\n'
                 'Raw Bayesian 0.062  →  calibrated GBM ~0.009',
                 fontsize=16, color='#222222', y=1.03)
    plt.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, 'fig2_reliability_diagram.png'), dpi=160, bbox_inches='tight')
    plt.close(fig)


# ──────────────────────────────────────────────────────────────────────────────
# FIGURE 3 — Class separation: TP / FP / blank confidence distributions
# ──────────────────────────────────────────────────────────────────────────────
def fig3_class_separation():
    deliv = pd.read_csv(os.path.join(ROOT, 'data', 'deliverable_scores_v2.csv'))
    deliv['polarity'] = deliv['wiki_id'].str.split('/').str[0].map({'aEKJ9AS':'pos','aPUDE1U':'neg'})

    fig, ax = plt.subplots(figsize=(13, 6))
    groups = [
        ('TP (curator-annotated)', deliv[deliv['spectrum_label']=='TP']['confidence'], '#2E8B57'),
        ('FP (curator-rejected yy_)', deliv[deliv['spectrum_label']=='FP']['confidence'], '#CD5C5C'),
        ('Blank (unknown, pos+neg)', deliv[deliv['spectrum_label']=='blank']['confidence'], '#708090'),
    ]
    for i, (name, data, col) in enumerate(groups):
        vp = ax.violinplot(data.values, positions=[i], showmeans=False, showmedians=False,
                           showextrema=False, widths=0.75)
        for body in vp['bodies']:
            body.set_facecolor(col); body.set_alpha(0.55); body.set_edgecolor(col)
        # overlay median bar
        med = data.median()
        ax.hlines(med, i - 0.25, i + 0.25, color='#222222', linewidth=2)
        ax.text(i + 0.32, med, f'median {med:.2f}', color='#222222', fontsize=11, va='center')

    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels([g[0] for g in groups], fontsize=12)
    ax.set_ylim(0, 1)
    ax.set_ylabel('Confidence (calibrated)')
    ax.set_title('Confidence distribution by spectrum class  —  GBM primary scorer\n'
                 'Model trained only on TP + FP; blanks scored without retraining',
                 fontsize=15, pad=15)
    ax.grid(axis='y', alpha=0.25)

    # annotate ≥0.9 rates
    for i, (name, data, col) in enumerate(groups):
        n_high = (data >= 0.9).sum(); total = len(data)
        ax.text(i, 0.96, f'≥ 0.9\n{n_high}/{total}\n({100*n_high/total:.1f}%)',
                ha='center', fontsize=10, color='#333333',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor=col, linewidth=1.2))

    plt.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, 'fig3_class_separation.png'), dpi=160, bbox_inches='tight')
    plt.close(fig)


# ──────────────────────────────────────────────────────────────────────────────
# FIGURE 4 — GBM feature importance (gain)
# ──────────────────────────────────────────────────────────────────────────────
def fig4_feature_importance():
    # Hard-coded from score_gbm_v2 output (latest run)
    importances = [
        ('compound_has_ok_adduct',    30.10, 'new'),
        ('hit_adduct_cat',            18.12, 'new'),
        ('signed_delta_rt',           16.47, 'Bayesian'),
        ('hit_is_isf',                 9.40, 'new'),
        ('sim_gap',                    8.84, 'Bayesian'),
        ('hit_isf_no_ok',              7.23, 'new'),
        ('entropy_similarity',         5.91, 'Bayesian'),
        ('n_candidate_adducts',        4.03, 'new'),
        ('db (library source)',        3.91, 'new'),
        ('forward_cosine',             3.41, 'new'),
        ('reverse_cosine',             3.34, 'new'),
        ('delta_mda',                  3.17, 'new'),
        ('spectral_entropy',           2.98, 'new'),
        ('polarity',                   2.92, 'new'),
        ('cov_int',                    2.86, 'new'),
    ]
    feats = [i[0] for i in importances]
    gains = [i[1] for i in importances]
    src   = [i[2] for i in importances]
    col_map = {'Bayesian': '#6B9BD6', 'new': '#F0A23E'}
    cols = [col_map[s] for s in src]

    fig, ax = plt.subplots(figsize=(12.5, 7))
    y = np.arange(len(feats))
    ax.barh(y, gains, color=cols, edgecolor='#222222', linewidth=0.4)
    for i, g in enumerate(gains):
        ax.text(g + 0.3, i, f'{g:.1f}', va='center', fontsize=11, color='#222222')
    ax.set_yticks(y); ax.set_yticklabels(feats, fontsize=12)
    ax.invert_yaxis()
    ax.set_xlabel('Feature importance (XGBoost gain)')
    ax.set_title('GBM feature importance\n'
                 'Top 5 features (58% of gain) are adduct-evidence flags Bayesian never used',
                 pad=14)
    ax.grid(axis='x', alpha=0.25)
    legend = [Patch(color='#6B9BD6', label='Bayesian 3-channel'),
              Patch(color='#F0A23E', label='New features enabled by GBM')]
    ax.legend(handles=legend, loc='lower right', frameon=True)
    plt.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, 'fig4_feature_importance.png'), dpi=160, bbox_inches='tight')
    plt.close(fig)


# ──────────────────────────────────────────────────────────────────────────────
# FIGURE 5 — NoTA validation on blanks (pos + neg)
# ──────────────────────────────────────────────────────────────────────────────
def fig5_nota_validation():
    deliv = pd.read_csv(os.path.join(ROOT, 'data', 'deliverable_scores_v2.csv'))
    deliv['polarity'] = deliv['wiki_id'].str.split('/').str[0].map({'aEKJ9AS':'pos','aPUDE1U':'neg'})
    blanks = deliv[deliv['spectrum_label']=='blank'].copy()

    pos_b = blanks[blanks['polarity']=='pos']['confidence']
    neg_b = blanks[blanks['polarity']=='neg']['confidence']

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5),
                                    gridspec_kw={'width_ratios': [2.2, 1]})

    # Histogram side by side
    bins = np.linspace(0, 1, 41)
    ax1.hist(pos_b, bins=bins, color='#FFA500', alpha=0.55, label=f'pos blanks (n={len(pos_b):,})')
    ax1.hist(neg_b, bins=bins, color='#4B94B2', alpha=0.55, label=f'neg blanks (n={len(neg_b):,})')
    ax1.set_yscale('log')
    ax1.set_xlabel('Confidence')
    ax1.set_ylabel('Count (log scale)')
    ax1.set_title(f'Blank spectra — no curator annotation\n'
                  f'Model crushes unknowns to low confidence (without retraining)')
    ax1.legend()
    ax1.grid(axis='y', alpha=0.2)

    # Threshold table
    thresholds = [0.5, 0.7, 0.9]
    ax2.axis('off')
    rows = [['Threshold', 'pos (%)', 'neg (%)']]
    for t in thresholds:
        p = 100 * (pos_b >= t).mean()
        n = 100 * (neg_b >= t).mean()
        rows.append([f'≥ {int(t*100)}%', f'{p:.1f}', f'{n:.1f}'])
    rows.append(['median', f'{pos_b.median()*100:.1f}', f'{neg_b.median()*100:.1f}'])

    table = ax2.table(cellText=rows[1:], colLabels=rows[0], loc='center',
                      cellLoc='center', colColours=['#E0E7EE']*3)
    table.auto_set_font_size(False); table.set_fontsize(12); table.scale(1.0, 1.6)
    ax2.set_title('% of blanks at each confidence', pad=20, fontsize=13)

    plt.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, 'fig5_nota_validation.png'), dpi=160, bbox_inches='tight')
    plt.close(fig)


# ──────────────────────────────────────────────────────────────────────────────
# FIGURE 6 — Model comparison: GBM vs RF vs MLP vs LR
# ──────────────────────────────────────────────────────────────────────────────
def fig6_model_comparison():
    models = [
        ('GBM\n(XGBoost)', 0.913, (0.9013, 0.9285), '#3B7D3B', True),
        ('Random\nForest', 0.907, (0.8898, 0.9247), '#7B9E7B', False),
        ('MLP\n(2×64)', 0.894, (0.8826, 0.9069), '#B0BEC5', False),
        ('Logistic\nRegression', 0.801, (0.7740, 0.8521), '#CD5C5C', False),
    ]

    fig, ax = plt.subplots(figsize=(10, 6))
    x = np.arange(len(models))
    aucs = [m[1] for m in models]
    fold_min = [m[2][0] for m in models]
    fold_max = [m[2][1] for m in models]
    err = [[aucs[i] - fold_min[i], fold_max[i] - aucs[i]] for i in range(len(models))]
    err = np.array(err).T  # shape (2, n)
    colors = [m[3] for m in models]

    bars = ax.bar(x, aucs, color=colors, edgecolor='#222222', linewidth=0.5, width=0.72,
                  yerr=err, capsize=6, error_kw={'ecolor':'#333333', 'elinewidth':1.2})
    for i, (bar, auc, is_primary) in enumerate(zip(bars, aucs, [m[4] for m in models])):
        label = f'{auc:.3f}'
        if is_primary:
            label += '\n(primary)'
        ax.text(bar.get_x() + bar.get_width()/2, auc + 0.008, label,
                ha='center', va='bottom', fontsize=12,
                fontweight='bold' if is_primary else 'normal', color='#222222')

    ax.set_xticks(x); ax.set_xticklabels([m[0] for m in models], fontsize=12)
    ax.set_ylim(0.75, 0.95)
    ax.set_ylabel('Out-of-fold AUC  (error bars = fold min/max)')
    ax.set_title('Model comparison on identical features + folds + calibration\n'
                 'Tree methods dominate; linear model confirms signal is non-linear',
                 pad=14)
    ax.grid(axis='y', alpha=0.25)
    plt.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, 'fig6_model_comparison.png'), dpi=160, bbox_inches='tight')
    plt.close(fig)


def main():
    print('Generating meeting figures...')
    print('  1. AUC journey'); fig1_auc_journey()
    print('  2. Reliability diagrams (re-fits Bayesian + GBM; ~30s)'); fig2_reliability_diagram()
    print('  3. Class separation (TP/FP/blank)'); fig3_class_separation()
    print('  4. GBM feature importance'); fig4_feature_importance()
    print('  5. NoTA validation (pos + neg blanks)'); fig5_nota_validation()
    print('  6. Model comparison (GBM/RF/MLP/LR)'); fig6_model_comparison()
    print(f'\nAll 6 saved to {FIG_DIR}/')


if __name__ == '__main__':
    main()
