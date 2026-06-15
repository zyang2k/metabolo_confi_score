"""4-arm isolation bench on the FULL (production-training) slice.

Same arms as bench_molrex_isolate.py but no golden filter:
  baseline
  + struct_logit only
  + struct_centroid_dist only
  + both

Used to check whether struct_logit's golden-slice lift (+0.0022 AUC,
-0.0016 Brier_cal) survives on the noisier full set the production GBM
actually trains on. If yes, struct_logit is ship-ready on the production
training distribution. If no, the golden lift was a slice artifact.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from bench_harness import (
    FEATURE_TABLE,
    build_top1,
    make_groups,
    print_importance,
    run_arm,
    summarize,
)

ROOT = Path(__file__).resolve().parent.parent
SIDE = ROOT / "data" / "molrex_features.csv"
OUT = ROOT / "data" / "bench_molrex_isolate_full_summary.csv"

LOGIT = "struct_logit"
CENT = "struct_centroid_dist"


def main() -> None:
    print(f"Loading {FEATURE_TABLE}")
    ft = pd.read_csv(FEATURE_TABLE, low_memory=False)
    print(f"  {len(ft):,} rows, {ft.wiki_id.nunique():,} unique spectra")

    top1 = build_top1(ft).reset_index(drop=True)
    n_tp = (top1.spectrum_label == "TP").sum()
    n_fp = (top1.spectrum_label == "FP").sum()
    print(f"\nFull slice: {len(top1):,} bins "
          f"({n_tp:,} TP / {n_fp:,} FP)   "
          f"prior_TP_row={top1.hit_label.mean():.3f}")

    side = pd.read_csv(SIDE)
    side["wiki_id"] = side.wiki_id.astype(str)
    top1["wiki_id"] = top1.wiki_id.astype(str)
    top1 = top1.merge(
        side[["wiki_id", "hit_ik14", LOGIT, CENT]],
        on=["wiki_id", "hit_ik14"], how="left",
    )
    print(f"  {LOGIT} NaN: {top1[LOGIT].isna().sum()}   "
          f"{CENT} NaN: {top1[CENT].isna().sum()}")

    groups = make_groups(top1)
    y = top1.hit_label.values

    arms = []
    print("\n=== ARM A: baseline ===")
    arms.append(run_arm(top1, [], groups, "baseline"))

    print(f"\n=== ARM B: +{LOGIT} only ===")
    arms.append(run_arm(top1, [LOGIT], groups, f"+{LOGIT}"))

    print(f"\n=== ARM C: +{CENT} only ===")
    arms.append(run_arm(top1, [CENT], groups, f"+{CENT}"))

    print(f"\n=== ARM D: +{LOGIT} + {CENT} ===")
    arms.append(run_arm(top1, [LOGIT, CENT], groups, "+both"))

    summarize(arms, y, str(OUT))
    print_importance(arms, [LOGIT, CENT], topk=18)

    base, b, c, d = (a.auc for a in arms)
    print("\n" + "=" * 60)
    print(f"ISOLATION READOUT (FULL slice, n={len(top1):,})")
    print("=" * 60)
    print(f"  baseline                       : AUC={base:.4f}")
    print(f"  +{LOGIT:24s}      : AUC={b:.4f}   Δ={b-base:+.4f}")
    print(f"  +{CENT:24s}      : AUC={c:.4f}   Δ={c-base:+.4f}")
    print(f"  +both                          : AUC={d:.4f}   Δ={d-base:+.4f}")
    print(f"  redundancy: (B+C)-base - (D-base) = "
          f"{(b - base) + (c - base) - (d - base):+.4f}")


if __name__ == "__main__":
    main()
