"""bench_rt_confirmation.py — P1: same-compound, same-RT, clean-adduct confirmation feature.

Oliver 2026-06-09: an ISF hit (e.g. 2-aminoadipic acid [M+H-H2O]+) should NOT be penalized
when the same compound is confirmed at the same RT with a clean adduct (its [M+H] parent
co-eluting). Reproduction confirmed the effect (50% vs 98% same compound) AND that the naive
"ok adduct anywhere" fix is a trap (35,735 rows, 4% TP). The fix needs the RT gate.

FEATURE (LABEL-FREE — uses only observed adduct/RT/esim, never hit_label, so no label leak):
  rt_confirm      1 if the same hit_ik14 appears in ANOTHER bin within ±RT_WIN s (observed RT)
                  with a clean 'ok' adduct AND entropy_similarity >= ESIM_MIN (label-free quality)
  rt_confirm_esim max entropy_similarity of such a confirming sibling (0 if none) — graded version

CIRCULARITY GUARD (pre-registered, decide BEFORE looking at AUC):
  Oliver labels by RT-disagreement-within-compound (yy_ FP = RT-disagreement). A same-RT
  confirmation feature is the COMPLEMENT of that rule, so an overall AUC lift would be largely
  circular (re-reading the label). Therefore we DO NOT justify this by overall ΔAUC. We justify it as:
    C1 RESCUE: the genuine confirmed-ISF cases (2-aminoadipic type) get rt_confirm=1.
    C2 SPECIFICITY: the same-RT gate filters the junk — the 4%-TP "ok-adduct-anywhere" ISF
       population is mostly NOT confirmed at the same RT (confirmed rate << anywhere rate).
    C3 LEAVE-OUT-THE-RULE: ΔAUC on the NON-yy_ slice (labels not set by the RT rule) is >= 0,
       i.e. the feature does not HURT non-RT-determined cases. Overall ΔAUC reported but
       labelled "rule-encoding" not "generalization".
  SHIP if C1 and C2 hold and C3 >= 0. Frame to Oliver as "implements your confirmed-ISF rule".
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd, xgboost as xgb
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold
from bench_harness import BASE_NUMERIC, CATEGORICAL, XGB_PARAMS, NUM_BOOST_ROUND

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FT = os.path.join(ROOT, 'data', 'feature_table_v2.csv')
RT_WIN = 10.0      # seconds — same-compound co-elution window
ESIM_MIN = 0.50    # label-free quality gate on the confirming sibling
NEW = ['rt_confirm', 'rt_confirm_esim']


def compute_confirmation(ft):
    """Label-free: for each candidate row, is its compound confirmed by a co-eluting clean-adduct sibling bin?"""
    ft = ft.copy()
    ft['rt_obs'] = pd.to_numeric(ft['rt_obs'], errors='coerce')
    ft['entropy_similarity'] = pd.to_numeric(ft['entropy_similarity'], errors='coerce')
    # "clean confirmers": rows that are an ok adduct with a good spectral match
    conf = ft[(ft['hit_adduct_cat'] == 'ok') & (ft['entropy_similarity'] >= ESIM_MIN) & ft['rt_obs'].notna()]
    # per-compound bin-level confirmers: (wiki_id, rt_obs, best esim) — one per bin
    cb = (conf.groupby(['hit_ik14', 'wiki_id'])
              .agg(rt=('rt_obs', 'first'), esim=('entropy_similarity', 'max')).reset_index())
    by_ik = {ik: g[['wiki_id', 'rt', 'esim']].values for ik, g in cb.groupby('hit_ik14')}

    rc = np.zeros(len(ft)); rce = np.zeros(len(ft))
    ik_arr = ft['hit_ik14'].values; wid_arr = ft['wiki_id'].values; rt_arr = ft['rt_obs'].values
    for i in range(len(ft)):
        ik = ik_arr[i]; rt = rt_arr[i]
        if ik not in by_ik or not np.isfinite(rt):
            continue
        best = 0.0
        for w2, rt2, es2 in by_ik[ik]:
            if w2 == wid_arr[i]:           # exclude the bin's own confirmers (self)
                continue
            if abs(rt - rt2) <= RT_WIN:    # co-elution at the SAME observed RT
                if es2 > best:
                    best = es2
        if best > 0:
            rc[i] = 1; rce[i] = best
    ft['rt_confirm'] = rc; ft['rt_confirm_esim'] = rce
    return ft


def main():
    ft = pd.read_csv(FT, low_memory=False)
    ft = ft[ft['hit_ik14'].fillna('').ne('')].reset_index(drop=True)
    ft = compute_confirmation(ft)

    # ---- the canonical case ----
    r = ft[(ft.wiki_id == 'aEKJ9AS/3034') & (ft.hit_ik14 == 'OYIFNHCXNCRBQI')]
    if len(r):
        print(f'C1 RESCUE — 2-aminoadipic [M+H-H2O]+ (bin /3034): rt_confirm={int(r.iloc[0].rt_confirm)} '
              f'(esim of co-eluting clean parent {r.iloc[0].rt_confirm_esim:.3f})')

    # ---- C2 specificity: ISF-no-ok population, confirmed at SAME RT vs "ok-adduct anywhere" ----
    isf = ft[ft.hit_isf_no_ok == 1]
    xbin_ok = ft.assign(_ok=(ft.hit_adduct_cat == 'ok').astype(int)).groupby('hit_ik14')['_ok'].max()
    anywhere = isf['hit_ik14'].map(xbin_ok).fillna(0).astype(int)
    print(f'\nC2 SPECIFICITY — ISF-no-ok candidate rows: {len(isf):,}')
    print(f'   ok-adduct ANYWHERE (the naive trap): {int((anywhere==1).sum()):,}')
    print(f'   confirmed at SAME RT (rt_confirm=1) : {int(isf.rt_confirm.sum()):,}  <- the surgical set')
    lab_isf = isf[isf.spectrum_label.isin(['TP', 'FP'])]
    for nm, msk in [('confirmed@RT', lab_isf.rt_confirm == 1), ('NOT confirmed', lab_isf.rt_confirm == 0)]:
        s = lab_isf[msk]
        if len(s):
            print(f'   labeled ISF-no-ok {nm:13s}: n={len(s):5d}  TP rate {s.hit_label.mean()*100:4.0f}%')

    # ---- OOF AUC (GroupKFold) base vs +confirmation, overall + leave-out-the-rule (non-yy_) ----
    lab = ft[ft.spectrum_label.isin(['TP', 'FP']) & ft.anno_ik14.fillna('').ne('')]
    t = ft.loc[lab.groupby('wiki_id').entropy_similarity.idxmax()].reset_index(drop=True)
    y = t.hit_label.values
    groups = t.anno_ik14.fillna('').values.copy()
    for i in range(len(groups)):
        if groups[i] == '': groups[i] = f'__no_{i}'
    is_yy = t['name'].fillna('').str.lower().str.startswith('yy_').values   # RT-rule-determined FP

    def prep(df, num, cm=None):
        X = pd.DataFrame(index=df.index); fit = cm is None; cm = cm or {}
        for c in num: X[c] = pd.to_numeric(df[c], errors='coerce') if c in df else np.nan
        for c in CATEGORICAL:
            s = df[c].astype(str).fillna('missing') if c in df else pd.Series('missing', index=df.index)
            if fit: cm[c] = {v: i for i, v in enumerate(sorted(s.unique()))}
            mx = max(cm[c].values())+1 if cm[c] else 0
            X[c] = s.map(lambda v: cm[c].get(v, mx)).astype('int32')
        return X[num+CATEGORICAL], cm

    def oof(num):
        o = np.full(len(t), np.nan)
        for tr, te in GroupKFold(5).split(t, y, groups):
            Xtr, cm = prep(t.iloc[tr], num); Xte, _ = prep(t.iloc[te], num, cm)
            ft_types = ['q']*len(num)+['c']*len(CATEGORICAL)
            m = xgb.train(XGB_PARAMS, xgb.DMatrix(Xtr, label=y[tr], enable_categorical=True, feature_types=ft_types), NUM_BOOST_ROUND)
            o[te] = m.predict(xgb.DMatrix(Xte, enable_categorical=True, feature_types=ft_types))
        return o
    ob = oof(BASE_NUMERIC); oc = oof(BASE_NUMERIC + NEW)
    nonrt = ~is_yy   # TP + non-yy_ FP = labels NOT set by the RT-disagreement rule
    print(f'\nOOF AUC  base {roc_auc_score(y,ob):.4f}  +confirm {roc_auc_score(y,oc):.4f}  '
          f'Δ {roc_auc_score(y,oc)-roc_auc_score(y,ob):+.4f}  (overall — rule-encoding, not a generalization claim)')
    print(f'C3 LEAVE-OUT-THE-RULE — non-yy_ slice (n={int(nonrt.sum())}, FP={int((nonrt&(y==0)).sum())}): '
          f'base {roc_auc_score(y[nonrt],ob[nonrt]):.4f}  +confirm {roc_auc_score(y[nonrt],oc[nonrt]):.4f}  '
          f'Δ {roc_auc_score(y[nonrt],oc[nonrt])-roc_auc_score(y[nonrt],ob[nonrt]):+.4f}  (must be >= 0)')


if __name__ == '__main__':
    main()
