"""Isolate which across-set MolRex feature carries the +0.002 AUC lift.

Four arms on the golden slice (same filter as bench_molrex_golden.py):
  baseline
  + struct_logit         (OOF logistic on raw emb -> label predictor)
  + struct_centroid_dist (distance to training-fold TP centroid)
  + both

struct_sim_gap is intentionally excluded — three independent confirmations
that within-spectrum candidate similarity is dead (Tanimoto-BP 2026-04-28
ΔAUC −0.0035; struct_sim_gap full-slice +0.0006; struct_sim_gap
golden-slice −0.0017).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
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

ROOT = Path(__file__).resolve().parent.parent
OH_PATH = ROOT / "data" / "orbitrap_hits_v2.csv"
SIDE = ROOT / "data" / "molrex_features.csv"
OUT = ROOT / "data" / "bench_molrex_isolate_summary.csv"

RT_TOLERANCE_S = 10.0
LOGIT = "struct_logit"
CENT = "struct_centroid_dist"


def pull_measured_rt(oh_path: Path) -> pd.DataFrame:
    chunks = []
    for chunk in pd.read_csv(
        oh_path,
        usecols=["wiki_id", "predicted_rt_hilic", "delta_predicted_rt"],
        dtype={"wiki_id": "string"},
        chunksize=200_000,
    ):
        chunk = chunk.dropna(subset=["predicted_rt_hilic",
                                     "delta_predicted_rt"])
        chunk["measured_rt"] = (chunk.predicted_rt_hilic
                                + chunk.delta_predicted_rt)
        chunks.append(chunk[["wiki_id", "measured_rt"]]
                      .groupby("wiki_id", as_index=False).first())
    return (pd.concat(chunks, ignore_index=True)
              .groupby("wiki_id", as_index=False).first())


def build_golden_slice() -> pd.DataFrame:
    print(f"Loading {FEATURE_TABLE}")
    ft = pd.read_csv(FEATURE_TABLE, low_memory=False)
    print(f"  {len(ft):,} rows, {ft.wiki_id.nunique():,} unique spectra")

    print("\nPulling measured RT ...")
    rt_df = pull_measured_rt(OH_PATH)
    ft["wiki_id"] = ft.wiki_id.astype(str)
    rt_df["wiki_id"] = rt_df.wiki_id.astype(str)
    ft = ft.merge(rt_df, on="wiki_id", how="left")

    top1 = build_top1(ft).reset_index(drop=True)
    tp = top1[top1.spectrum_label == "TP"].copy()

    agg = (tp.dropna(subset=["measured_rt", "anno_ik14"])
             .groupby("anno_ik14")
             .agg(n_bins=("wiki_id", "nunique"),
                  n_adducts=("adduct", "nunique"),
                  rt_min=("measured_rt", "min"),
                  rt_max=("measured_rt", "max"))
             .reset_index())
    agg["rt_spread"] = agg.rt_max - agg.rt_min
    golden_iks = set(
        agg.query("n_bins>=2 and n_adducts>=2 and rt_spread<=@RT_TOLERANCE_S")
           .anno_ik14
    )
    golden_tp_bins = set(tp[tp.anno_ik14.isin(golden_iks)].wiki_id)
    fp_bins = set(top1[top1.spectrum_label == "FP"].wiki_id)
    keep = golden_tp_bins | fp_bins
    golden = top1[top1.wiki_id.isin(keep)].reset_index(drop=True)

    side = pd.read_csv(SIDE)
    side["wiki_id"] = side.wiki_id.astype(str)
    golden = golden.merge(
        side[["wiki_id", "hit_ik14", LOGIT, CENT]],
        on=["wiki_id", "hit_ik14"], how="left",
    )

    n_tp = (golden.spectrum_label == "TP").sum()
    n_fp = (golden.spectrum_label == "FP").sum()
    print(f"\nGolden slice: {len(golden):,} bins "
          f"({n_tp:,} TP / {n_fp:,} FP)   "
          f"prior_TP_row={golden.hit_label.mean():.3f}")
    print(f"  {LOGIT} NaN: {golden[LOGIT].isna().sum()}   "
          f"{CENT} NaN: {golden[CENT].isna().sum()}")
    return golden


def main() -> None:
    golden = build_golden_slice()
    groups = make_groups(golden)
    y = golden.hit_label.values

    arms = []
    print("\n=== ARM A: baseline ===")
    arms.append(run_arm(golden, [], groups, "baseline"))

    print(f"\n=== ARM B: +{LOGIT} only ===")
    arms.append(run_arm(golden, [LOGIT], groups, f"+{LOGIT}"))

    print(f"\n=== ARM C: +{CENT} only ===")
    arms.append(run_arm(golden, [CENT], groups, f"+{CENT}"))

    print(f"\n=== ARM D: +{LOGIT} + {CENT} ===")
    arms.append(run_arm(golden, [LOGIT, CENT], groups, "+both"))

    summarize(arms, y, str(OUT))
    print_importance(arms, [LOGIT, CENT], topk=18)

    # ---- isolation summary ----
    base, b, c, d = (a.auc for a in arms)
    print("\n" + "=" * 60)
    print("ISOLATION READOUT")
    print("=" * 60)
    print(f"  baseline                : AUC={base:.4f}")
    print(f"  +{LOGIT:24s}: AUC={b:.4f}   Δ={b-base:+.4f}")
    print(f"  +{CENT:24s}: AUC={c:.4f}   Δ={c-base:+.4f}")
    print(f"  +both                          : AUC={d:.4f}   Δ={d-base:+.4f}")
    print(f"  redundancy check: (B+C)-base - (D-base) = "
          f"{(b - base) + (c - base) - (d - base):+.4f}")
    print("    positive => features carry overlapping signal (redundant)")
    print("    near 0   => features are additive")
    print("    negative => features are synergistic (rare)")


if __name__ == "__main__":
    main()
