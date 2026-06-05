"""bench_relational_frozen.py — HONEST per-fold-graph test of the relational features.

The shipped bench (bench_relational.py) merged a relational_features.csv precomputed over
ALL labeled bins, then GroupKFold'd. That lets a test bin's different-compound neighbours
come from its own held-out fold (their labels aren't available at production scoring time).
This script rebuilds the graph from the TRAIN fold only, inside each fold — the same thing
score_gbm_v2.py now does. It answers: does the +0.0070 survive production-faithful scoring?

Reports OOF AUC base vs base+rel, overall and on the 16% confusable band, no bootstrap.
"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd, xgboost as xgb
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold
from bench_harness import BASE_NUMERIC, CATEGORICAL, XGB_PARAMS, NUM_BOOST_ROUND
from relational_graph import Graph, REL_COLS

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FT = os.path.join(ROOT, 'data', 'feature_table_v2.csv')
QP = os.path.join(ROOT, 'data', 'query_peaks_cache_v2.json')
BIN_COLS = ['wiki_id', 'precursor_mz', 'polarity', 'anno_ik14', 'anno_name_lower', 'hit_label']


def prep(df, num_feats, cat_maps=None):
    X = pd.DataFrame(index=df.index)
    for c in num_feats:
        X[c] = pd.to_numeric(df[c], errors='coerce') if c in df else np.nan
    fit = cat_maps is None
    if fit: cat_maps = {}
    for c in CATEGORICAL:
        s = df[c].astype(str).fillna('missing') if c in df else pd.Series('missing', index=df.index)
        if fit:
            cat_maps[c] = {v: i for i, v in enumerate(sorted(s.unique()))}
        mx = max(cat_maps[c].values()) + 1 if cat_maps[c] else 0
        X[c] = s.map(lambda v: cat_maps[c].get(v, mx)).astype('int32')
    return X[num_feats + CATEGORICAL], cat_maps


def train_pred(Xtr, ytr, Xte, n_num):
    ft_types = ['q'] * n_num + ['c'] * len(CATEGORICAL)
    dtr = xgb.DMatrix(Xtr, label=ytr, enable_categorical=True, feature_types=ft_types)
    dte = xgb.DMatrix(Xte, enable_categorical=True, feature_types=ft_types)
    return xgb.train(XGB_PARAMS, dtr, NUM_BOOST_ROUND).predict(dte)


def main():
    ft = pd.read_csv(FT, low_memory=False)
    ft = ft[ft['hit_ik14'].fillna('').ne('')].reset_index(drop=True)
    lab = ft[ft['spectrum_label'].isin(['TP', 'FP']) & ft['anno_ik14'].fillna('').ne('')]
    t = ft.loc[lab.groupby('wiki_id')['entropy_similarity'].idxmax()].reset_index(drop=True)
    y = t['hit_label'].values
    print(f'labeled top-1 bins: {len(t):,}  prior={y.mean():.3f}')

    peaks = {k: np.asarray(v, np.float64) for k, v in json.load(open(QP)).items() if v}
    groups = t['anno_ik14'].fillna('').values.copy()
    for i in range(len(groups)):
        if groups[i] == '': groups[i] = f'__no_{i}'

    Xbase, _ = prep(t, BASE_NUMERIC)
    oof_b = np.full(len(t), np.nan)
    oof_r = np.full(len(t), np.nan)
    rel_nn = np.zeros(len(t))
    for tr, te in GroupKFold(5).split(t, y, groups):
        # base arm
        oof_b[te] = train_pred(Xbase.iloc[tr], y[tr], Xbase.iloc[te], len(BASE_NUMERIC))
        # +rel arm with TRAIN-ONLY graph (honest)
        g = Graph(t.iloc[tr][BIN_COLS], peaks)
        ftr = g.features(t.iloc[tr][BIN_COLS]).set_index('wiki_id')
        fte = g.features(t.iloc[te][BIN_COLS]).set_index('wiki_id')
        Xtr = Xbase.iloc[tr].copy(); Xte = Xbase.iloc[te].copy()
        for c in REL_COLS:
            Xtr.insert(len(BASE_NUMERIC) + REL_COLS.index(c), c,
                       pd.Series(t.iloc[tr]['wiki_id'].values).map(ftr[c]).fillna(0.0).values)
            Xte.insert(len(BASE_NUMERIC) + REL_COLS.index(c), c,
                       pd.Series(t.iloc[te]['wiki_id'].values).map(fte[c]).fillna(0.0).values)
        rel_nn[te] = Xte['rel_nn_sim'].values
        oof_r[te] = train_pred(Xtr, y[tr], Xte, len(BASE_NUMERIC) + len(REL_COLS))

    auc_b = roc_auc_score(y, oof_b)
    auc_r = roc_auc_score(y, oof_r)
    print(f'\nOOF AUC  base            : {auc_b:.4f}')
    print(f'OOF AUC  base + rel (honest): {auc_r:.4f}   Δ = {auc_r-auc_b:+.4f}')

    band = rel_nn >= 0.70
    print(f'\nconfusable band (rel_nn_sim>=0.70): {band.mean()*100:.1f}% of bins ({band.sum()})')
    if band.sum() > 30:
        ab_b = roc_auc_score(y[band], oof_b[band]); ab_r = roc_auc_score(y[band], oof_r[band])
        print(f'  band AUC base {ab_b:.4f}  +rel {ab_r:.4f}   Δ = {ab_r-ab_b:+.4f}')
        nb = ~band
        nb_b = roc_auc_score(y[nb], oof_b[nb]); nb_r = roc_auc_score(y[nb], oof_r[nb])
        print(f'  rest AUC base {nb_b:.4f}  +rel {nb_r:.4f}   Δ = {nb_r-nb_b:+.4f}')


if __name__ == '__main__':
    main()
