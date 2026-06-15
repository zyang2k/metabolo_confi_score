"""bench_joint_inference_flags.py — Sibling-conflict curator triage flags.

Recasts joint inference v1 from "re-scorer" to "label-quality detector / triage."

For each (hit_ik14, adduct) pair appearing as GBM top-1 in ≥2 bins, computes:
  - the bins in the collision group
  - the observed-RT range across those bins
  - severity tier (RT spread + group size + dual-annotation name pattern)

Outputs:
  data/bench_joint_inference_flags.csv          — per-bin flag columns on all bins
  data/oliver_review_sibling_conflicts.csv      — sorted shortlist for curator review

Severity tiers (highest takes precedence):
  3 — RT range > 60s   (very strong incompatibility under fixed column/method)
  2 — RT range > 30s   (genuine incompatibility)
  1 — RT range ≤ 30s   (probably same chromatographic peak — peak-pick redundancy)
  0 — not in a collision group

Severity boosts (additive):
  +1 if compound name has underscore pattern (dual annotation, e.g. "Methylhistamine_Hexanal")
  +1 if collision group has ≥5 bins (broader curator-attention payoff)
"""
import os
import re
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
CAND = ROOT / 'data' / 'candidate_scores_v2.csv'
SPEC_NEG = ROOT / 'data' / 'Orbitrap_HILIC_negESI_curated_042126.csv'
SPEC_POS = ROOT / 'data' / 'Orbitrap_HILIC_posESI_curated_042126.csv'
OUT_FLAGS = ROOT / 'data' / 'bench_joint_inference_flags.csv'
OUT_REVIEW = ROOT / 'data' / 'oliver_review_sibling_conflicts.csv'
OUT_GROUPS = ROOT / 'data' / 'oliver_review_groups_summary.csv'

RT_TIGHT = 30.0  # seconds
RT_STRONG = 60.0  # seconds

# Dual-annotation pattern: two capitalized words joined by underscore. Tight
# heuristic — matches "Methylhistamine_Hexanal", not "Benzoic_acid" or "2_5H_furanone".
DUAL_NAME_RE = re.compile(r'^[A-Z][a-z]+[A-Za-z0-9-]*_[A-Z][a-z]+[A-Za-z0-9-]*')

# Confidence threshold for "actually worth curator attention" — a collision group
# is only interesting if the model is confident in at least one of its bins.
MIN_CONF_FOR_REVIEW = 0.5


def has_dual_name(name) -> bool:
    if not isinstance(name, str):
        return False
    head = name.split(' ')[0]
    if not head or head.startswith('yy_'):
        return False
    return bool(DUAL_NAME_RE.match(head))


def severity(rt_range: float, group_size: int, max_conf: float, dual_name: bool) -> int:
    """Severity = chemistry-incompatibility tier × confidence × group attention.

    Returns a sortable int where higher = more curator-worthy.
    """
    if rt_range > RT_STRONG:
        rt_tier = 3
    elif rt_range > RT_TIGHT:
        rt_tier = 2
    elif group_size >= 2:
        rt_tier = 1
    else:
        rt_tier = 0
    # confidence tier — only severity if at least one bin in group is confident
    if max_conf >= 0.8:
        conf_tier = 3
    elif max_conf >= 0.5:
        conf_tier = 2
    elif max_conf >= 0.3:
        conf_tier = 1
    else:
        conf_tier = 0
    boost = (1 if dual_name else 0) + (1 if group_size >= 5 else 0)
    return rt_tier * 100 + conf_tier * 10 + boost


def main():
    print(f'Loading {CAND}')
    cand = pd.read_csv(CAND)
    cand = cand[cand['gbm_cal'].notna()].copy()

    # GBM top-1 + GBM 2nd-best per bin
    cand_sorted = cand.sort_values(['wiki_id', 'gbm_cal'], ascending=[True, False])
    cand_sorted['rank_in_bin'] = cand_sorted.groupby('wiki_id').cumcount()
    top1 = cand_sorted[cand_sorted['rank_in_bin'] == 0].drop(columns=['rank_in_bin']).copy()
    top2 = cand_sorted[cand_sorted['rank_in_bin'] == 1][
        ['wiki_id', 'hit_ik14', 'name', 'adduct', 'gbm_cal']
    ].rename(columns={
        'hit_ik14': 'alt2_ik14', 'name': 'alt2_name',
        'adduct': 'alt2_adduct', 'gbm_cal': 'alt2_gbm_cal',
    }).copy()
    top1 = top1.merge(top2, on='wiki_id', how='left')
    print(f'  top-1 rows: {len(top1):,}  (with alt2 attached)')

    # Observed RT + SPLASH + curator-annotation status per bin from the curated tables.
    # SPLASH is the stable identifier (wiki_id drifts on BinBase scan-order shifts —
    # see reference_splash_vs_wiki_id memory).
    print(f'Loading spectra metadata from {SPEC_NEG.name} + {SPEC_POS.name}')
    spec_cols = ['wiki_id', 'rt', 'raw_splash', 'is_manual_annotated']
    neg = pd.read_csv(SPEC_NEG, usecols=spec_cols, low_memory=False)
    neg['polarity'] = 'neg'
    pos = pd.read_csv(SPEC_POS, usecols=spec_cols, low_memory=False)
    pos['polarity'] = 'pos'
    spec = pd.concat([neg, pos], ignore_index=True).drop_duplicates(subset='wiki_id', keep='first')
    spec = spec.rename(columns={'raw_splash': 'splash'})
    top1 = top1.merge(
        spec[['wiki_id', 'rt', 'polarity', 'splash', 'is_manual_annotated']],
        on='wiki_id', how='left',
    )

    # Collision groups: (hit_ik14, adduct) appearing as top-1 in ≥2 bins.
    # Restrict to CURATOR-ANNOTATED bins (is_manual_annotated=True). Unannotated
    # bins have no SMILES/InChI from the curator — they don't belong in any
    # group_size or rt_range that ends up in a curator-review artifact. The
    # pipeline still scores them for production deliverable purposes, but they
    # shouldn't bias the triage stats here.
    valid_ik = top1[
        top1['hit_ik14'].notna()
        & (top1['hit_ik14'] != '')
        & top1['is_manual_annotated'].fillna(False)
    ].copy()
    print(f'  curator-annotated top-1 rows (basis for collision groups): {len(valid_ik):,}')
    grp = valid_ik.groupby(['hit_ik14', 'adduct'])
    group_stats = grp.agg(
        group_size=('wiki_id', 'size'),
        rt_min=('rt', 'min'),
        rt_max=('rt', 'max'),
        rt_median=('rt', 'median'),
        rt_std=('rt', 'std'),
        max_group_conf=('gbm_cal', 'max'),
    ).reset_index()
    group_stats['rt_range'] = group_stats['rt_max'] - group_stats['rt_min']
    collision = group_stats[group_stats['group_size'] >= 2].copy()

    print(f'\nCollision groups (same hit_ik14+adduct as top-1 in ≥2 bins): {len(collision):,}')
    print(f'  RT range > 30s: {(collision["rt_range"] > RT_TIGHT).sum():,}')
    print(f'  RT range > 60s: {(collision["rt_range"] > RT_STRONG).sum():,}')

    # Attach group stats back to per-bin rows
    flags = top1.merge(
        collision[['hit_ik14', 'adduct', 'group_size', 'rt_range', 'rt_median', 'max_group_conf']],
        on=['hit_ik14', 'adduct'], how='left',
    )
    flags['in_collision'] = flags['group_size'].notna()
    flags['group_size'] = flags['group_size'].fillna(1).astype(int)
    flags['rt_range'] = flags['rt_range'].fillna(0.0)
    flags['max_group_conf'] = flags['max_group_conf'].fillna(0.0)

    flags['flag_sibling_conflict'] = flags['in_collision']
    flags['flag_rt_incompatible'] = flags['rt_range'] > RT_TIGHT
    flags['flag_rt_strongly_incompatible'] = flags['rt_range'] > RT_STRONG
    flags['flag_confident_in_group'] = flags['max_group_conf'] >= MIN_CONF_FOR_REVIEW
    flags['flag_dual_name'] = flags['name'].apply(has_dual_name)
    flags['severity'] = flags.apply(
        lambda r: severity(
            r['rt_range'], int(r['group_size']),
            float(r['max_group_conf']), bool(r['flag_dual_name']),
        ),
        axis=1,
    )

    n_total = len(flags)
    print(f'\nFlag summary across {n_total:,} bins:')
    print(f'  in any collision:          {flags["flag_sibling_conflict"].sum():,}  '
          f'({100*flags["flag_sibling_conflict"].mean():.1f}%)')
    print(f'  RT-incompatible (>30s):    {flags["flag_rt_incompatible"].sum():,}  '
          f'({100*flags["flag_rt_incompatible"].mean():.1f}%)')
    print(f'  Strongly incompat (>60s):  {flags["flag_rt_strongly_incompatible"].sum():,}  '
          f'({100*flags["flag_rt_strongly_incompatible"].mean():.1f}%)')
    print(f'  Dual-name annotation:      {flags["flag_dual_name"].sum():,}  '
          f'({100*flags["flag_dual_name"].mean():.1f}%)')

    # Save full per-bin flag artifact (includes SPLASH for stable cross-reference)
    flag_cols = [
        'wiki_id', 'splash', 'polarity', 'rt', 'is_manual_annotated',
        'name', 'adduct', 'hit_ik14', 'hit_label',
        'gbm_cal', 'gbm_raw',
        'alt2_name', 'alt2_ik14', 'alt2_adduct', 'alt2_gbm_cal',
        'group_size', 'rt_range', 'max_group_conf',
        'flag_sibling_conflict', 'flag_rt_incompatible',
        'flag_rt_strongly_incompatible', 'flag_confident_in_group',
        'flag_dual_name',
        'severity',
    ]
    flags[flag_cols].to_csv(OUT_FLAGS, index=False)
    print(f'\nWrote {OUT_FLAGS} ({len(flags):,} rows)')

    # Oliver review shortlist — only flagged rows worth curator attention.
    # Three filters now:
    #   1. RT-incompatible group (>30s spread under fixed column/method)
    #   2. THIS BIN itself has gbm_cal >= MIN_CONF_FOR_REVIEW (not just group max —
    #      otherwise unannotated bins with conf=0 sneak in via collision-group lookup)
    #   3. THIS BIN was curator-annotated (is_manual_annotated=True). Unannotated
    #      bins (raw_name "unknown_*", is_manual_annotated=False) get top-1 from
    #      candidate generator but Oliver never said anything about them — they
    #      don't belong in a curator-review list at all.
    review = flags[
        flags['flag_rt_incompatible']
        & (flags['gbm_cal'] >= MIN_CONF_FOR_REVIEW)
        & flags['is_manual_annotated'].fillna(False)
    ].copy()
    review['collision_group_id'] = (
        review['hit_ik14'].fillna('') + '|' + review['adduct'].fillna('')
    )
    review = review.sort_values(
        ['severity', 'collision_group_id', 'gbm_cal'],
        ascending=[False, True, False],
    )
    review_cols = [
        'severity', 'collision_group_id', 'group_size', 'rt_range', 'max_group_conf',
        'splash', 'wiki_id', 'polarity', 'rt', 'name', 'adduct', 'hit_ik14',
        'gbm_cal', 'hit_label',
        'alt2_name', 'alt2_ik14', 'alt2_adduct', 'alt2_gbm_cal',
        'flag_rt_strongly_incompatible', 'flag_dual_name',
    ]
    review[review_cols].to_csv(OUT_REVIEW, index=False)
    print(f'Wrote {OUT_REVIEW} ({len(review):,} rows)')

    # One-row-per-collision-group summary — skim-friendly for curator.
    # Includes top-bin confidence, label profile, and a hint about what joint
    # inference would do.
    groups_view = review.copy()
    groups_view['n_tp'] = groups_view['hit_label'].fillna(-1).eq(1).astype(int)
    groups_view['n_fp'] = groups_view['hit_label'].fillna(-1).eq(0).astype(int)
    grp_summary = groups_view.groupby('collision_group_id').agg(
        severity=('severity', 'max'),
        name=('name', 'first'),
        hit_ik14=('hit_ik14', 'first'),
        adduct=('adduct', 'first'),
        n_bins=('wiki_id', 'size'),
        rt_range=('rt_range', 'first'),
        max_conf=('gbm_cal', 'max'),
        min_conf=('gbm_cal', 'min'),
        n_TP=('n_tp', 'sum'),
        n_FP=('n_fp', 'sum'),
        polarities=('polarity', lambda s: ','.join(sorted(set(s.dropna())))),
        flag_dual_name=('flag_dual_name', 'first'),
    ).reset_index().sort_values('severity', ascending=False)
    grp_summary.to_csv(OUT_GROUPS, index=False)
    print(f'Wrote {OUT_GROUPS} ({len(grp_summary):,} groups)')

    # Console: top 5 collision groups by severity, fully shown
    print('\n=== Top severity collision groups (preview) ===')
    top_groups = (
        review.drop_duplicates('collision_group_id')
        .sort_values(['severity', 'group_size'], ascending=[False, False])
        .head(5)
    )
    for _, g in top_groups.iterrows():
        members = review[review['collision_group_id'] == g['collision_group_id']]
        print(f'\nseverity={g["severity"]}  group={g["collision_group_id"][:30]}  '
              f'n={int(g["group_size"])}  rt_range={g["rt_range"]:.1f}s  '
              f'name="{g["name"]}"  pol={g["polarity"]}')
        for _, m in members.head(6).iterrows():
            print(f'   wiki_id={m["wiki_id"]:<18}  rt={m["rt"]:6.1f}s  conf={m["gbm_cal"]:.3f}  '
                  f'label={m["hit_label"]}  '
                  f'alt2={str(m["alt2_name"])[:30]} (conf={m["alt2_gbm_cal"]:.3f})')
        if len(members) > 6:
            print(f'   ... +{len(members)-6} more in group')


if __name__ == '__main__':
    main()
