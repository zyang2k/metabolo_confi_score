"""Bench MolRex-derived features against baseline GBM.

Three arms (uses code/bench_harness.py primitives):
  baseline                                         current production
  +struct_sim_gap                                  cleanest single feature
  +struct_sim_gap +struct_logit +struct_centroid_dist   kitchen sink

Plus a stratified diagnostic: report each arm's AUC on rows where
signed_delta_rt is available vs. missing — tells us whether MolRex
features pay rent specifically on the null-RT slice.

Side table:  data/molrex_features.csv  (built by code/build_molrex_features.py)
Join keys:   (wiki_id, hit_ik14)
"""
from __future__ import annotations

import os

import pandas as pd
from sklearn.metrics import roc_auc_score

# --- sibling-import path shim (code/ root) ---
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))

from bench_harness import (
    FEATURE_TABLE,
    build_top1,
    make_groups,
    print_importance,
    run_arm,
    summarize,
)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXTRA_TABLE = os.path.join(ROOT, "data", "molrex_features.csv")
OUT_PATH = os.path.join(ROOT, "data", "bench_molrex_summary.csv")

EXTRA_COLS = ["struct_sim_gap", "struct_logit", "struct_centroid_dist"]


def main() -> None:
    print(f"Loading {FEATURE_TABLE}")
    ft = pd.read_csv(FEATURE_TABLE, low_memory=False)
    print(f"  {len(ft):,} rows, {ft.wiki_id.nunique():,} unique spectra")

    print(f"Loading side table {EXTRA_TABLE}")
    side = pd.read_csv(EXTRA_TABLE)
    keep = ["wiki_id", "hit_ik14"] + EXTRA_COLS
    side = side[[c for c in keep if c in side.columns]]
    print(f"  side rows: {len(side):,}")

    ft = ft.merge(side, on=["wiki_id", "hit_ik14"], how="left")
    for c in EXTRA_COLS:
        cov = ft[c].notna().mean() * 100
        print(f"  merged {c:24s}: {cov:5.2f}% overall coverage")

    top1 = build_top1(ft)
    print(f"  {len(top1):,} top-1 rows; prior_TP={top1.hit_label.mean():.3f}")
    for c in EXTRA_COLS:
        n_nan = top1[c].isna().sum()
        print(f"  in top-1: {c:24s} NaN: {n_nan:>4d} "
              f"({n_nan / len(top1) * 100:5.1f}%)")

    groups = make_groups(top1)

    print("\n=== ARM A: baseline (current production features) ===")
    arm_a = run_arm(top1, [], groups, "baseline")

    print("\n=== ARM B: +struct_sim_gap (RT-independent, intra-bin) ===")
    arm_b = run_arm(top1, ["struct_sim_gap"], groups, "+sim_gap")

    print("\n=== ARM C: +sim_gap +logit +centroid (kitchen sink) ===")
    arm_c = run_arm(top1, EXTRA_COLS, groups, "+all_three")

    y = top1["hit_label"].values
    summarize([arm_a, arm_b, arm_c], y, OUT_PATH)
    print_importance([arm_a, arm_b, arm_c], EXTRA_COLS, topk=18)

    # --- Stratified-by-RT-availability diagnostic ---
    has_rt = top1["signed_delta_rt"].notna().values
    n_with, n_without = int(has_rt.sum()), int((~has_rt).sum())
    print("\n=== Stratified by signed_delta_rt availability ===")
    print(f"  with RT   : {n_with:>5,d} ({n_with / len(top1):.1%})   "
          f"TP prior={y[has_rt].mean():.3f}")
    print(f"  without RT: {n_without:>5,d} ({n_without / len(top1):.1%})   "
          f"TP prior={y[~has_rt].mean():.3f}")
    print(f"  {'arm':<24s}  {'AUC|RT':>8s}  {'AUC|noRT':>9s}  "
          f"{'Δ_RT':>7s}  {'Δ_noRT':>8s}")
    base_auc_with = roc_auc_score(y[has_rt], arm_a.oof[has_rt])
    base_auc_without = (roc_auc_score(y[~has_rt], arm_a.oof[~has_rt])
                        if n_without >= 30 else float("nan"))
    for arm in (arm_a, arm_b, arm_c):
        auc_with = roc_auc_score(y[has_rt], arm.oof[has_rt])
        auc_without = (roc_auc_score(y[~has_rt], arm.oof[~has_rt])
                       if n_without >= 30 else float("nan"))
        d_with = auc_with - base_auc_with
        d_without = (auc_without - base_auc_without
                     if n_without >= 30 else float("nan"))
        print(f"  {arm.name:<24s}  {auc_with:>8.4f}  {auc_without:>9.4f}  "
              f"{d_with:+8.4f}  {d_without:+8.4f}")


if __name__ == "__main__":
    main()
