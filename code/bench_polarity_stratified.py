"""Polarity-stratified OOF AUC for the production GBM (with struct_logit).

Confirms the single-model design (polarity as feature) serves both
polarities without an asymmetric weakness.
"""
from __future__ import annotations

import os

import numpy as np
import pandas as pd
from sklearn.isotonic import IsotonicRegression
from sklearn.metrics import brier_score_loss, roc_auc_score

from bench_harness import (
    FEATURE_TABLE,
    build_top1,
    make_groups,
    run_arm,
)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _ece(p, y, nbins=10) -> float:
    edges = np.linspace(0, 1, nbins + 1)
    bi = np.clip(np.digitize(p, edges[1:-1]), 0, nbins - 1)
    total, n = 0.0, 0
    for b in range(nbins):
        m = bi == b
        if m.sum() == 0:
            continue
        total += m.sum() * abs(p[m].mean() - y[m].mean())
        n += m.sum()
    return total / max(n, 1)


def main() -> None:
    print(f"Loading {FEATURE_TABLE}")
    ft = pd.read_csv(FEATURE_TABLE, low_memory=False)
    side = pd.read_csv(os.path.join(ROOT, "data", "molrex_features.csv"))
    side["wiki_id"] = side.wiki_id.astype(str)
    ft["wiki_id"] = ft.wiki_id.astype(str)
    ft = ft.merge(
        side[["wiki_id", "hit_ik14", "struct_logit"]],
        on=["wiki_id", "hit_ik14"], how="left",
    )

    top1 = build_top1(ft).reset_index(drop=True)
    groups = make_groups(top1)
    print(f"  top-1 rows: {len(top1):,}   "
          f"polarity counts: {top1.polarity.value_counts().to_dict()}")

    print("\nRunning OOF with current production features + struct_logit ...")
    arm = run_arm(top1, ["struct_logit"], groups, "production")
    y = top1["hit_label"].values
    oof = arm.oof

    # Global isotonic on full OOF (matches bench_harness convention)
    iso = IsotonicRegression(out_of_bounds="clip")
    iso.fit(oof, y)
    p_cal = iso.transform(oof)

    print("\n" + "=" * 70)
    print("POLARITY-STRATIFIED METRICS")
    print("=" * 70)
    rows = []
    for pol_val, pol_name in [(0, "neg"), (1, "pos")]:
        mask = (top1.polarity == pol_val).values
        n = int(mask.sum())
        n_tp = int((mask & (y == 1)).sum())
        n_fp = int((mask & (y == 0)).sum())
        if n < 30:
            print(f"  polarity={pol_name}: skipped (n={n} too small)")
            continue
        auc = roc_auc_score(y[mask], oof[mask])
        brier_raw = brier_score_loss(y[mask], oof[mask])
        brier_cal = brier_score_loss(y[mask], p_cal[mask])
        ece_raw = _ece(oof[mask], y[mask])
        ece_cal = _ece(p_cal[mask], y[mask])
        rows.append({
            "polarity": pol_name, "n": n, "n_TP": n_tp, "n_FP": n_fp,
            "prior_TP": y[mask].mean(),
            "AUC": auc, "Brier_raw": brier_raw, "Brier_cal": brier_cal,
            "ECE_raw": ece_raw, "ECE_cal": ece_cal,
        })

    # Overall row
    auc_all = roc_auc_score(y, oof)
    brier_raw_all = brier_score_loss(y, oof)
    brier_cal_all = brier_score_loss(y, p_cal)
    rows.append({
        "polarity": "all", "n": len(y),
        "n_TP": int((y == 1).sum()), "n_FP": int((y == 0).sum()),
        "prior_TP": y.mean(),
        "AUC": auc_all, "Brier_raw": brier_raw_all, "Brier_cal": brier_cal_all,
        "ECE_raw": _ece(oof, y), "ECE_cal": _ece(p_cal, y),
    })

    df = pd.DataFrame(rows)
    cols = ["polarity", "n", "n_TP", "n_FP", "prior_TP", "AUC",
            "Brier_raw", "Brier_cal", "ECE_raw", "ECE_cal"]
    print(df[cols].to_string(index=False,
                             float_format=lambda x: f"{x:.4f}"))

    out = os.path.join(ROOT, "data", "bench_polarity_summary.csv")
    df.to_csv(out, index=False)
    print(f"\nWrote {out}")

    print("\n--- Asymmetry check ---")
    pos = df[df.polarity == "pos"].iloc[0]
    neg = df[df.polarity == "neg"].iloc[0]
    print(f"  pos − neg:  ΔAUC = {pos.AUC - neg.AUC:+.4f}   "
          f"ΔBrier_cal = {pos.Brier_cal - neg.Brier_cal:+.4f}")


if __name__ == "__main__":
    main()
