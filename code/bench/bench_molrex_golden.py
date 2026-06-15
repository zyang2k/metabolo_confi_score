"""Re-run MolRex 3-arm bench on the strict multi-adduct-same-RT golden TP slice.

Golden TP criterion:
  spectrum_label == 'TP'  AND
  the bin's anno_ik14 appears in ≥2 distinct bins  AND
  those bins span ≥2 distinct adducts  AND
  the measured-RT spread across those bins ≤ RT_TOLERANCE_S seconds

FPs are kept as-is (the negatives are already curator-rejected).

Reports golden-slice 3-arm AUC + Brier + ECE, plus pre-filter inventory.
"""
from __future__ import annotations

import os
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
OUT = ROOT / "data" / "bench_molrex_golden_summary.csv"

RT_TOLERANCE_S = 10.0
EXTRA_COLS = ["struct_sim_gap", "struct_logit", "struct_centroid_dist"]


def pull_measured_rt(oh_path: Path) -> pd.DataFrame:
    """One measured_rt per wiki_id via predicted_rt_hilic + delta_predicted_rt."""
    chunks = []
    for chunk in pd.read_csv(
        oh_path,
        usecols=["wiki_id", "predicted_rt_hilic", "delta_predicted_rt"],
        dtype={"wiki_id": "string"},
        chunksize=200_000,
    ):
        chunk = chunk.dropna(subset=["predicted_rt_hilic", "delta_predicted_rt"])
        chunk["measured_rt"] = (chunk.predicted_rt_hilic
                                + chunk.delta_predicted_rt)
        chunks.append(chunk[["wiki_id", "measured_rt"]]
                      .groupby("wiki_id", as_index=False).first())
    df = pd.concat(chunks, ignore_index=True)
    return df.groupby("wiki_id", as_index=False).first()


def main() -> None:
    print(f"Loading {FEATURE_TABLE}")
    ft = pd.read_csv(FEATURE_TABLE, low_memory=False)
    print(f"  {len(ft):,} rows, {ft.wiki_id.nunique():,} unique spectra")

    print(f"\nPulling measured RT from {OH_PATH} ...")
    rt_df = pull_measured_rt(OH_PATH)
    print(f"  measured_rt for {len(rt_df):,} bins")

    ft["wiki_id"] = ft.wiki_id.astype(str)
    rt_df["wiki_id"] = rt_df.wiki_id.astype(str)
    ft = ft.merge(rt_df, on="wiki_id", how="left")
    print(f"  measured_rt coverage in ft: "
          f"{ft.measured_rt.notna().mean():.1%}")

    # top-1 across full ft, same as production bench
    top1 = build_top1(ft).reset_index(drop=True)
    print(f"\nTop-1 bins (TP+FP): {len(top1):,}   "
          f"prior_TP={top1.hit_label.mean():.3f}")

    # ---- Golden criterion ----
    tp = top1[top1.spectrum_label == "TP"].copy()
    print(f"\nTP top-1 candidates: {len(tp):,}")
    print(f"  TPs with measured_rt: {tp.measured_rt.notna().sum():,}")

    agg = (tp.dropna(subset=["measured_rt", "anno_ik14"])
             .groupby("anno_ik14")
             .agg(n_bins=("wiki_id", "nunique"),
                  n_adducts=("adduct", "nunique"),
                  rt_min=("measured_rt", "min"),
                  rt_max=("measured_rt", "max"))
             .reset_index())
    agg["rt_spread"] = agg.rt_max - agg.rt_min
    print(f"\nUnique anno_ik14 in TPs with RT: {len(agg):,}")
    print(f"  with n_bins>=2:                "
          f"{(agg.n_bins>=2).sum():,}")
    print(f"  with n_bins>=2 & n_adducts>=2: "
          f"{((agg.n_bins>=2) & (agg.n_adducts>=2)).sum():,}")
    print(f"  + rt_spread<={RT_TOLERANCE_S}s:       "
          f"{((agg.n_bins>=2) & (agg.n_adducts>=2) & (agg.rt_spread<=RT_TOLERANCE_S)).sum():,}")

    golden_iks = set(
        agg.query("n_bins>=2 and n_adducts>=2 and rt_spread<=@RT_TOLERANCE_S")
           .anno_ik14
    )
    print(f"\nGolden IK14s: {len(golden_iks):,}")

    golden_tp_bins = set(tp[tp.anno_ik14.isin(golden_iks)].wiki_id)
    fp_bins = set(top1[top1.spectrum_label == "FP"].wiki_id)
    keep_bins = golden_tp_bins | fp_bins
    golden = top1[top1.wiki_id.isin(keep_bins)].reset_index(drop=True)
    n_tp = (golden.spectrum_label == "TP").sum()
    n_fp = (golden.spectrum_label == "FP").sum()
    print(f"\nGolden bench scope: {len(golden):,} bins  "
          f"({n_tp:,} TP / {n_fp:,} FP)   "
          f"prior_TP_row={golden.hit_label.mean():.3f}")
    print(f"  golden TP rows hit_label=1: "
          f"{(golden[golden.spectrum_label=='TP'].hit_label==1).sum():,}")
    print(f"  golden TP rows hit_label=0: "
          f"{(golden[golden.spectrum_label=='TP'].hit_label==0).sum():,}")

    # ---- Merge MolRex features ----
    side = pd.read_csv(SIDE)
    side = side[["wiki_id", "hit_ik14"] + EXTRA_COLS]
    side["wiki_id"] = side.wiki_id.astype(str)
    golden = golden.merge(side, on=["wiki_id", "hit_ik14"], how="left")
    for c in EXTRA_COLS:
        n_nan = golden[c].isna().sum()
        print(f"  {c:24s} NaN in golden: {n_nan} "
              f"({n_nan / len(golden) * 100:.1f}%)")

    # ---- 3-arm bench ----
    groups = make_groups(golden)
    print("\n=== ARM A: baseline (golden slice) ===")
    arm_a = run_arm(golden, [], groups, "baseline_golden")
    print("\n=== ARM B: +struct_sim_gap (golden slice) ===")
    arm_b = run_arm(golden, ["struct_sim_gap"], groups, "+sim_gap_golden")
    print("\n=== ARM C: +all three (golden slice) ===")
    arm_c = run_arm(golden, EXTRA_COLS, groups, "+all_three_golden")

    y = golden.hit_label.values
    summarize([arm_a, arm_b, arm_c], y, str(OUT))
    print_importance([arm_a, arm_b, arm_c], EXTRA_COLS, topk=18)


if __name__ == "__main__":
    main()
