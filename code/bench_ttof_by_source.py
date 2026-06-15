"""Stratify TTOF bench AUC by library source (db column).

If TTOF AUC is dramatically different across library sources, that's
evidence the cross-platform transfer issue is dominated by
query-vs-library platform mismatch (not features alone).

Re-uses build_ttof_table from bench_ttof_golden.py to produce predictions,
then groups by `db` on the top-1 labeled rows.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xgboost as xgb
from sklearn.metrics import roc_auc_score

sys.path.insert(0, str(Path(__file__).resolve().parent))
from bench_harness import (
    BASE_NUMERIC,
    CATEGORICAL,
    FEATURE_TABLE,
    build_top1,
    prep,
    train,
)
from bench_ttof_golden import (
    LIB_PEAKS,
    TTOF_NEG_HITS, TTOF_NEG_PEAKS, TTOF_NEG_SPEC,
    TTOF_POS_HITS, TTOF_POS_PEAKS, TTOF_POS_SPEC,
    build_ttof_table,
    load_adduct_taxonomy,
)

ROOT = Path(__file__).resolve().parent.parent
MOLREX_FEATURES = ROOT / "data" / "molrex_features.csv"
OUT = ROOT / "data" / "bench_ttof_by_source_summary.csv"


def main() -> None:
    adduct_tax = load_adduct_taxonomy()
    print("Loading library peaks ...")
    with open(LIB_PEAKS) as f:
        lib_peaks = json.load(f)
    print(f"  library peaks: {len(lib_peaks):,} entries")

    ttof_neg = build_ttof_table(TTOF_NEG_HITS, TTOF_NEG_SPEC, TTOF_NEG_PEAKS,
                                lib_peaks, 0, adduct_tax)
    ttof_pos = build_ttof_table(TTOF_POS_HITS, TTOF_POS_SPEC, TTOF_POS_PEAKS,
                                lib_peaks, 1, adduct_tax)
    ttof = pd.concat([ttof_neg, ttof_pos], ignore_index=True)
    ttof_top1 = build_top1(ttof).reset_index(drop=True)

    # Train Orbitrap model
    orb = pd.read_csv(FEATURE_TABLE, low_memory=False)
    side = pd.read_csv(MOLREX_FEATURES)
    side["wiki_id"] = side.wiki_id.astype(str)
    orb["wiki_id"] = orb.wiki_id.astype(str)
    orb = orb.merge(
        side[["wiki_id", "hit_ik14", "struct_logit"]],
        on=["wiki_id", "hit_ik14"], how="left",
    )
    orb_top1 = build_top1(orb).reset_index(drop=True)

    print("\nTraining Orbitrap model ...")
    X_tr = prep(orb_top1, ["struct_logit"])
    y_tr = orb_top1["hit_label"].values
    n_num = len(BASE_NUMERIC) + 1
    n_cat = len(CATEGORICAL)
    model = train(X_tr, y_tr, n_num, n_cat)

    X_ttof = prep(ttof_top1, ["struct_logit"])
    dttof = xgb.DMatrix(X_ttof, enable_categorical=True,
                        feature_types=["q"] * n_num + ["c"] * n_cat)
    ttof_top1["pred"] = model.predict(dttof)

    # Stratify by db on labeled subset only
    print("\n" + "=" * 70)
    print("AUC stratified by library db (TTOF labeled bins)")
    print("=" * 70)
    rows = []
    for pol_val, pol_name in [(0, "neg"), (1, "pos")]:
        sub = ttof_top1[(ttof_top1.polarity == pol_val)
                        & ttof_top1.spectrum_label.isin(["TP", "FP"])]
        print(f"\n--- TTOF {pol_name} (n={len(sub):,}) ---")
        for db_val, grp in sub.groupby("db"):
            y = grp.hit_label.values
            p = grp.pred.values
            n = len(grp)
            n_tp = int((y == 1).sum())
            n_fp = int((y == 0).sum())
            if n_tp < 5 or n_fp < 5:
                print(f"  {db_val:18s}  n={n:>4d} TP={n_tp:>3d} FP={n_fp:>3d}  "
                      f"AUC=  (too few)")
                continue
            auc = roc_auc_score(y, p)
            print(f"  {db_val:18s}  n={n:>4d} TP={n_tp:>3d} FP={n_fp:>3d}  "
                  f"AUC={auc:.4f}  TP_med={np.median(p[y==1]):.3f}  "
                  f"FP_med={np.median(p[y==0]):.3f}")
            rows.append(dict(polarity=pol_name, db=db_val, n=n,
                             n_TP=n_tp, n_FP=n_fp, AUC=auc,
                             TP_median=float(np.median(p[y == 1])),
                             FP_median=float(np.median(p[y == 0]))))

    df = pd.DataFrame(rows)
    df.to_csv(OUT, index=False)
    print(f"\nWrote {OUT}")

    # Orbitrap comparison
    print("\n--- Orbitrap stratified (in-domain) for reference ---")
    # Use simple model prediction (not OOF) for comparison
    dorb = xgb.DMatrix(X_tr, enable_categorical=True,
                       feature_types=["q"] * n_num + ["c"] * n_cat)
    orb_top1["pred_in_sample"] = model.predict(dorb)
    sub = orb_top1[orb_top1.spectrum_label.isin(["TP", "FP"])]
    for db_val, grp in sub.groupby("db"):
        y = grp.hit_label.values
        p = grp["pred_in_sample"].values
        n = len(grp)
        n_tp = int((y == 1).sum())
        n_fp = int((y == 0).sum())
        if n_tp < 5 or n_fp < 5:
            continue
        auc = roc_auc_score(y, p)
        print(f"  {db_val:18s}  n={n:>5d} TP={n_tp:>4d} FP={n_fp:>4d}  "
              f"AUC={auc:.4f}  (in-sample, optimistic)")


if __name__ == "__main__":
    main()
