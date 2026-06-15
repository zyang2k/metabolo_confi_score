"""bench_joint_inference_v1.py — Phase 1.0 of cross-bin joint inference.

Sibling repulsion only, implemented as a greedy global assignment with hard
uniqueness constraint on (hit_ik14, adduct).

Each (hit_ik14, adduct) pair can be the top-1 of at most one bin. Iterate
(bin, candidate) rows by gbm_cal descending; first row claims the bin AND the
pair. Subsequent rows for the same pair must fall back to whichever next-best
candidate the bin still has available.

This is the λ_sib = ∞ limit of the soft-constraint ILP. If it shows signal
(curator-visible flips, Brier improvement on labeled set), Phase 1.1 moves to
soft penalty + PuLP. If it doesn't, the entire cross-bin direction is closed
cheaply.

Outputs:
    data/bench_joint_inference_v1.csv   — joint top-1 per bin
    data/bench_joint_inference_v1_flips.csv — bins where joint differs from GBM top-1
    Console: flip count, fallback-label retention, Min-108 robustness, AUC delta
"""
import os
from pathlib import Path
import numpy as np
import pandas as pd

from sklearn.metrics import roc_auc_score, brier_score_loss

ROOT = Path(__file__).resolve().parent.parent
CAND = ROOT / 'data' / 'candidate_scores_v2.csv'
DEL  = ROOT / 'data' / 'deliverable_scores_v2.csv'
MIN  = ROOT / 'benchmark' / '_internal' / 'min_verified_entries.csv'
OUT_JOINT = ROOT / 'data' / 'bench_joint_inference_v1.csv'
OUT_FLIPS = ROOT / 'data' / 'bench_joint_inference_v1_flips.csv'


def greedy_assign(cand: pd.DataFrame) -> pd.DataFrame:
    """Greedy global assignment with sibling repulsion as a hard constraint.

    For each (bin, candidate) row in descending order of gbm_cal:
      - claim the bin if it doesn't already have an assignment, AND
      - the (hit_ik14, adduct) pair isn't already claimed by another bin.

    Bins where every candidate is blocked by a higher-scoring sibling end up
    unassigned (NaN). These are rare in practice but worth surfacing.
    """
    df = cand.sort_values('gbm_cal', ascending=False, kind='mergesort').reset_index(drop=True)
    bin_taken: set = set()
    pair_taken: set = set()
    chosen_idx: list = []

    for i, row in df.iterrows():
        b = row['wiki_id']
        if b in bin_taken:
            continue
        ik = row['hit_ik14']
        ad = row['adduct']
        ik_valid = pd.notna(ik) and ik != ''
        if ik_valid and (ik, ad) in pair_taken:
            continue
        bin_taken.add(b)
        if ik_valid:
            pair_taken.add((ik, ad))
        chosen_idx.append(i)

    return df.loc[chosen_idx].copy()


def report_metrics(tag: str, df: pd.DataFrame) -> dict:
    labelled = df[df['hit_label'].notna()].copy()
    y = labelled['hit_label'].astype(int).values
    s = labelled['gbm_cal'].astype(float).values
    if len(set(y)) > 1:
        auc = roc_auc_score(y, s)
        brier = brier_score_loss(y, s)
    else:
        auc = float('nan'); brier = float('nan')
    tp = int((y == 1).sum()); fp = int((y == 0).sum())
    fdr_at_09 = float((s[y == 0] >= 0.9).sum()) / max(1, int((s >= 0.9).sum()))
    print(f'  [{tag}] n={len(df):,}  labelled={len(labelled):,}  TP={tp:,}  FP={fp:,}  '
          f'AUC={auc:.4f}  Brier={brier:.4f}  FDR@0.9={fdr_at_09:.4f}')
    return {'tag': tag, 'n': len(df), 'labelled': len(labelled), 'auc': auc,
            'brier': brier, 'fdr_at_09': fdr_at_09}


def main():
    print(f'Loading {CAND}')
    cand = pd.read_csv(CAND)
    cand = cand[cand['gbm_cal'].notna()].copy()
    print(f'  {len(cand):,} rows, {cand["wiki_id"].nunique():,} bins, '
          f'mean {len(cand)/cand["wiki_id"].nunique():.1f} cand/bin')

    print('\n-- baseline: GBM top-1 (no joint constraint) --')
    gbm_top1 = cand.loc[cand.groupby('wiki_id')['gbm_cal'].idxmax()].copy()
    base_metrics = report_metrics('gbm_top1', gbm_top1)

    # How many sibling collisions does the unconstrained top-1 have?
    sib = gbm_top1[gbm_top1['hit_ik14'].notna() & (gbm_top1['hit_ik14'] != '')]
    pair_counts = sib.groupby(['hit_ik14', 'adduct']).size()
    n_collisions = int((pair_counts >= 2).sum())
    n_rows_in_collision = int(pair_counts[pair_counts >= 2].sum())
    print(f'  sibling collisions (same hit_ik14+adduct in top-1 of >=2 bins): '
          f'{n_collisions:,} pairs, {n_rows_in_collision:,} rows '
          f'({100*n_rows_in_collision/len(gbm_top1):.1f}% of bins)')

    print('\n-- joint: greedy sibling repulsion (lambda_sib = inf) --')
    joint_top1 = greedy_assign(cand)
    joint_metrics = report_metrics('joint_top1', joint_top1)

    # Verify zero collisions remain
    j_sib = joint_top1[joint_top1['hit_ik14'].notna() & (joint_top1['hit_ik14'] != '')]
    j_pair_counts = j_sib.groupby(['hit_ik14', 'adduct']).size()
    assert (j_pair_counts >= 2).sum() == 0, 'sibling constraint violated'
    print(f'  sibling collisions after joint: 0 (by construction)')

    # Compare GBM top-1 vs joint top-1 per bin
    g = gbm_top1.set_index('wiki_id')
    j = joint_top1.set_index('wiki_id')

    cmp = g[['hit_ik14', 'adduct', 'name', 'gbm_cal', 'hit_label']].rename(
        columns={'hit_ik14': 'gbm_ik14', 'adduct': 'gbm_adduct',
                 'name': 'gbm_name', 'gbm_cal': 'gbm_score',
                 'hit_label': 'gbm_label'}
    ).join(
        j[['hit_ik14', 'adduct', 'name', 'gbm_cal', 'hit_label']].rename(
            columns={'hit_ik14': 'joint_ik14', 'adduct': 'joint_adduct',
                     'name': 'joint_name', 'gbm_cal': 'joint_score',
                     'hit_label': 'joint_label'}
        )
    )

    flipped = cmp[
        (cmp['gbm_ik14'].fillna('') != cmp['joint_ik14'].fillna('')) |
        (cmp['gbm_adduct'].fillna('') != cmp['joint_adduct'].fillna(''))
    ].copy()

    print(f'\n-- flips: bins where joint changes top-1 --')
    print(f'  total flips: {len(flipped):,}  ({100*len(flipped)/len(cmp):.1f}% of bins)')

    # Of the flips, what does the label transition look like?
    both_lab = flipped[flipped['gbm_label'].notna() & flipped['joint_label'].notna()].copy()
    if len(both_lab):
        tp_to_tp = ((both_lab['gbm_label'] == 1) & (both_lab['joint_label'] == 1)).sum()
        tp_to_fp = ((both_lab['gbm_label'] == 1) & (both_lab['joint_label'] == 0)).sum()
        fp_to_tp = ((both_lab['gbm_label'] == 0) & (both_lab['joint_label'] == 1)).sum()
        fp_to_fp = ((both_lab['gbm_label'] == 0) & (both_lab['joint_label'] == 0)).sum()
        print(f'  flips where both old & new have labels (n={len(both_lab):,}):')
        print(f'    TP -> TP: {tp_to_tp:,}    (lateral move, both correct)')
        print(f'    TP -> FP: {tp_to_fp:,}    (constraint HURT — replaced TP with FP)')
        print(f'    FP -> TP: {fp_to_tp:,}    (constraint HELPED — replaced FP with TP)')
        print(f'    FP -> FP: {fp_to_fp:,}    (lateral move, both wrong)')

    score_drop = (flipped['gbm_score'] - flipped['joint_score']).describe()
    print(f'  score drop on flips (gbm_score - joint_score): '
          f'median={score_drop["50%"]:.3f}  mean={score_drop["mean"]:.3f}  max={score_drop["max"]:.3f}')

    # Min-108 verification
    if MIN.exists():
        mn = pd.read_csv(MIN)
        print(f'\n-- Min-{len(mn)} verified set --')
        if 'wiki_id' in mn.columns:
            mn_joint = j.loc[j.index.intersection(mn['wiki_id'])].copy()
            mn_gbm   = g.loc[g.index.intersection(mn['wiki_id'])].copy()
            # Match by IK14 if available in the verified file, else by name
            if 'verified_ik14' in mn.columns:
                truth = mn.set_index('wiki_id')['verified_ik14']
                gbm_correct = (mn_gbm['hit_ik14'] == truth.loc[mn_gbm.index]).sum()
                joint_correct = (mn_joint['hit_ik14'] == truth.loc[mn_joint.index]).sum()
                print(f'  GBM correct:   {gbm_correct}/{len(mn_gbm)}')
                print(f'  joint correct: {joint_correct}/{len(mn_joint)}')
            else:
                print(f'  (no verified_ik14 column in {MIN.name}; manual inspection needed)')
                print(f'  columns: {list(mn.columns)}')

    # Save outputs
    joint_top1.reset_index(drop=True).to_csv(OUT_JOINT, index=False)
    flipped.reset_index().to_csv(OUT_FLIPS, index=False)
    print(f'\nWrote:\n  {OUT_JOINT}\n  {OUT_FLIPS}')


if __name__ == '__main__':
    main()
