"""bench_name_rt_propagate.py — RT-anchor label propagation (KG prototype layer 2).

Paradigm 1 of the knowledge-graph pivot (project_knowledge_graph_pivot): propagate
curated evidence along same-identity edges to the UNCURATED bins that layer 1
(bench_name_rt_graph.py) structurally could not reach.

Idea
----
Curated TRUE-POSITIVE bins establish, per compound, a reference elution time
(anchor RT). An uncurated bin that *claims* that compound is then scored by whether
it elutes where the compound is KNOWN to elute — independent of the GBM's own
confidence, which is selection-biased low on unannotated bins (only 245/9,175
uncurated bins clear conf 0.5; project_selection_bias_blanks).

  anchor_rt(X)  = median observed RT of curated TP bins assigned to compound X
  rt_agreement  = exp(-(rt_bin - anchor_rt)^2 / (2 sigma^2)),  sigma = RT_SIGMA

Two propagated outputs (both go to Oliver, neither is a calibrated score):
  MISSED ANNOTATION  — uncurated bin, anchor exists, rt_agreement high AND decent MS2
                       => elutes at the known RT with matching spectrum but never
                          annotated. Candidate missed annotation (NOT confirmed — could
                          be isomer/ISF; for curator review). Cf. feedback_annotated_only
                          (unannotated assumed FP) — this challenges that for the subset
                          that looks real.
  CONTRADICTED       — uncurated bin, anchor exists, rt_agreement low => claims a known
                       compound at the wrong RT. Likely FP / ISF / misassignment.

Validation (leave-one-out on curated, NOT circular)
---------------------------------------------------
For curated bins whose compound has >=1 OTHER curated TP bin, recompute the anchor
from the peers (excluding self) and ask whether rt_agreement-to-peers separates this
bin's own hit_label (TP vs FP). This tests the propagation signal directly. It shares
the within-name RT mechanism's known ~1/5 yy_ coverage ceiling (project_knowledge_graph_pivot,
project_yy_formation) — report AUC descriptively, do not over-claim.

Outputs
  data/name_rt_propagated.csv          — per-bin propagated score (all bins w/ an anchor)
  data/missed_annotation_candidates.csv — uncurated, high agreement + MS2, ranked
  data/contradicted_uncurated.csv       — uncurated, anchor exists, wrong RT
"""
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parent.parent
GRAPH = ROOT / 'data' / 'bench_name_rt_graph.csv'  # layer-1 artifact (composable)
OUT_PROP = ROOT / 'data' / 'name_rt_propagated.csv'
OUT_MISSED = ROOT / 'data' / 'missed_annotation_candidates.csv'
OUT_CONTRA = ROOT / 'data' / 'contradicted_uncurated.csv'

RT_SIGMA = 20.0          # s — agreement kernel width (within ~1 lab peak tolerance)
AGREE_HI = 0.5           # rt_agreement >= 0.5  => within ~24 s of anchor
AGREE_LO = 0.1           # rt_agreement <  0.1  => >~43 s off anchor (contradiction)
MS2_OK = 0.5             # entropy_similarity floor for a "missed annotation" candidate


def rt_agreement(rt, anchor):
    d = rt - anchor
    return np.exp(-(d * d) / (2.0 * RT_SIGMA * RT_SIGMA))


def main():
    print(f'Loading layer-1 graph {GRAPH.name}')
    df = pd.read_csv(GRAPH)
    df['is_curated'] = df['is_curated'].astype(bool)

    # ---- Anchors: median RT of curated TRUE-POSITIVE bins per compound ----------
    tp = df[df['is_curated'] & (df['hit_label'] == 1)].copy()
    anchor = (tp.groupby('norm_name')
                .agg(anchor_rt=('rt', 'median'),
                     anchor_n=('rt', 'size'),
                     anchor_rt_std=('rt', 'std'))
                .reset_index())
    anchor['anchor_rt_std'] = anchor['anchor_rt_std'].fillna(0.0)
    print(f'  compounds with a curated-TP anchor: {len(anchor):,}')

    # Multi-peak compounds: anchor (single median) is unreliable when curated TP bins
    # themselves span a wide RT — flag so we don't over-trust agreement there.
    anchor['anchor_multimodal'] = anchor['anchor_rt_std'] > RT_SIGMA

    df = df.merge(anchor, on='norm_name', how='left')
    df['anchor_multimodal'] = df['anchor_multimodal'].fillna(False).astype(bool)
    has_anchor = df['anchor_rt'].notna()
    df['rt_agreement'] = np.where(has_anchor, rt_agreement(df['rt'], df['anchor_rt']), np.nan)
    df['delta_to_anchor'] = np.where(has_anchor, df['rt'] - df['anchor_rt'], np.nan)

    # ---- Coverage: how far does propagation reach into the uncurated set? -------
    unc = df[~df['is_curated']]
    print('\n=== Propagation reach ===')
    print(f'  uncurated bins:                 {len(unc):,}')
    print(f'  uncurated WITH an anchor:       {int(unc["anchor_rt"].notna().sum()):,}  '
          f'(claim a compound the curator established)')

    # ---- MISSED-ANNOTATION candidates -------------------------------------------
    missed = df[
        (~df['is_curated'])
        & df['anchor_rt'].notna()
        & (~df['anchor_multimodal'])
        & (df['rt_agreement'] >= AGREE_HI)
        & (df['entropy_similarity'] >= MS2_OK)
    ].copy()
    missed['propagated_score'] = missed['rt_agreement'] * missed['entropy_similarity']
    missed = missed.sort_values('propagated_score', ascending=False)

    # ---- CONTRADICTED uncurated -------------------------------------------------
    contra = df[
        (~df['is_curated'])
        & df['anchor_rt'].notna()
        & (~df['anchor_multimodal'])
        & (df['rt_agreement'] < AGREE_LO)
        & (df['entropy_similarity'] >= MS2_OK)
    ].copy().sort_values('entropy_similarity', ascending=False)

    print('\n=== Propagated findings (uncurated) ===')
    print(f'  MISSED-ANNOTATION candidates (right RT + MS2>={MS2_OK}): {len(missed):,}')
    print(f'  CONTRADICTED (claims known compound, wrong RT):         {len(contra):,}')

    # ---- Leave-one-out validation on curated bins -------------------------------
    # Recompute each curated bin's anchor from its TP PEERS (exclude self), then test
    # whether rt_agreement-to-peers separates its own TP/FP label.
    print('\n=== Leave-one-out validation (curated bins w/ >=1 TP peer) ===')
    cur = df[df['is_curated'] & df['hit_label'].isin([0, 1])].copy()
    # peer anchor = median RT of *other* curated TP bins of same compound
    tp_sum = tp.groupby('norm_name')['rt'].agg(['sum', 'size']).rename(columns={'sum': 'tp_sum', 'size': 'tp_n'})
    cur = cur.merge(tp_sum, on='norm_name', how='left')
    is_self_tp = (cur['hit_label'] == 1).astype(float)
    peer_n = cur['tp_n'].fillna(0) - is_self_tp
    peer_sum = cur['tp_sum'].fillna(0.0) - is_self_tp * cur['rt']
    cur['peer_anchor'] = np.where(peer_n > 0, peer_sum / peer_n.replace(0, np.nan), np.nan)
    loo = cur[cur['peer_anchor'].notna()].copy()
    loo['rt_agreement_peer'] = rt_agreement(loo['rt'], loo['peer_anchor'])
    n_tp = int((loo['hit_label'] == 1).sum())
    n_fp = int((loo['hit_label'] == 0).sum())
    if n_tp >= 5 and n_fp >= 5:
        auc = roc_auc_score(loo['hit_label'].values, loo['rt_agreement_peer'].values)
        print(f'  n={len(loo):,} (TP={n_tp:,}, FP={n_fp:,}) | '
              f'rt_agreement-to-peers AUC for TP/FP = {auc:.3f}')
        print('  (descriptive: same within-name RT mechanism, ~1/5 yy_ coverage ceiling)')
    else:
        print(f'  too few peer-anchored curated bins to validate (TP={n_tp}, FP={n_fp})')

    # ---- Write artifacts --------------------------------------------------------
    keep = ['wiki_id', 'splash', 'polarity', 'rt', 'is_curated', 'name', 'norm_name',
            'curator_name', 'hit_ik14', 'adduct', 'hit_label', 'gbm_cal',
            'entropy_similarity', 'anchor_rt', 'anchor_n', 'anchor_multimodal',
            'delta_to_anchor', 'rt_agreement']
    keep = [c for c in keep if c in df.columns]
    df[df['anchor_rt'].notna()][keep].to_csv(OUT_PROP, index=False)
    print(f'\nWrote {OUT_PROP} ({int(df["anchor_rt"].notna().sum()):,} rows)')

    miss_cols = keep + ['propagated_score']
    missed[[c for c in miss_cols if c in missed.columns]].to_csv(OUT_MISSED, index=False)
    print(f'Wrote {OUT_MISSED} ({len(missed):,} rows)')
    contra[keep].to_csv(OUT_CONTRA, index=False)
    print(f'Wrote {OUT_CONTRA} ({len(contra):,} rows)')

    # ---- Preview ----------------------------------------------------------------
    print('\n=== Top missed-annotation candidates ===')
    for _, b in missed.head(8).iterrows():
        print(f'   wiki_id={str(b["wiki_id"]):<16} rt={b["rt"]:7.1f}s anchor={b["anchor_rt"]:7.1f}s '
              f'(n={int(b["anchor_n"])}) agree={b["rt_agreement"]:.2f} esim={b["entropy_similarity"]:.2f} '
              f'gbm={b["gbm_cal"]:.2f} pol={b["polarity"]} "{str(b["name"])[:28]}"')

    print('\n=== Top contradicted uncurated (claims known compound, wrong RT) ===')
    for _, b in contra.head(6).iterrows():
        print(f'   wiki_id={str(b["wiki_id"]):<16} rt={b["rt"]:7.1f}s anchor={b["anchor_rt"]:7.1f}s '
              f'Δ={b["delta_to_anchor"]:+7.1f}s esim={b["entropy_similarity"]:.2f} '
              f'pol={b["polarity"]} "{str(b["name"])[:28]}"')


if __name__ == '__main__':
    main()
