"""Train-scope architectural experiment: does training on noisy labels hurt?

2x2 cells, all evaluated on identical golden test folds (so AUC numbers
are apples-to-apples). For each test fold the test compounds (anno_ik14)
are excluded from BOTH possible training sets.

  cell                       train scope                test scope
  TF_B  (train-full, baseline)        full set     golden held-out fold
  TG_B  (train-golden, baseline)      golden       golden held-out fold
  TF_L  (train-full, +struct_logit)   full set     golden held-out fold
  TG_L  (train-golden, +struct_logit) golden       golden held-out fold

Comparisons of interest:
  TG_B vs TF_B: does training on noisy labels hurt the BASELINE model?
  TG_L vs TF_L: does training on noisy labels hurt struct_logit's lift?
  TF_L vs TF_B: struct_logit's lift under current production design
  TG_L vs TG_B: struct_logit's lift under golden-only training
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.isotonic import IsotonicRegression
from sklearn.metrics import brier_score_loss, roc_auc_score
from sklearn.model_selection import GroupKFold

from bench_harness import (
    FEATURE_TABLE,
    build_top1,
    make_groups,
    prep,
    train,
)
import xgboost as xgb

ROOT = Path(__file__).resolve().parent.parent
OH_PATH = ROOT / "data" / "orbitrap_hits_v2.csv"
SIDE = ROOT / "data" / "molrex_features.csv"
OUT = ROOT / "data" / "bench_train_scope_summary.csv"

RT_TOLERANCE_S = 10.0
LOGIT = "struct_logit"


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


def identify_golden_iks(top1: pd.DataFrame) -> set:
    tp = top1[top1.spectrum_label == "TP"]
    agg = (tp.dropna(subset=["measured_rt", "anno_ik14"])
             .groupby("anno_ik14")
             .agg(n_bins=("wiki_id", "nunique"),
                  n_adducts=("adduct", "nunique"),
                  rt_min=("measured_rt", "min"),
                  rt_max=("measured_rt", "max"))
             .reset_index())
    agg["rt_spread"] = agg.rt_max - agg.rt_min
    return set(
        agg.query("n_bins>=2 and n_adducts>=2 and rt_spread<=@RT_TOLERANCE_S")
           .anno_ik14
    )


def ece(p, y, nbins=10) -> float:
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


def run_cell(top1: pd.DataFrame, golden_mask: np.ndarray,
             groups: np.ndarray, extra: list, train_scope: str,
             cell_name: str) -> dict:
    """One cell: train on `train_scope` slice, evaluate on golden test folds.

    For each fold, test = test_idx ∩ golden_mask. Train rows are restricted
    by train_scope ('full' = all train_idx rows, 'golden' = train_idx rows
    whose anno_ik14 ∈ golden_iks).
    """
    from bench_harness import BASE_NUMERIC, CATEGORICAL
    X = prep(top1, extra)
    y = top1.hit_label.values
    n_num = len(BASE_NUMERIC) + len(extra)
    n_cat = len(CATEGORICAL)

    gkf = GroupKFold(n_splits=5)
    oof = np.full(len(top1), np.nan)
    fold_aucs = []
    for fold, (tr, te) in enumerate(gkf.split(top1, y, groups)):
        # restrict test to golden compounds only
        te_g = te[golden_mask[te]]
        # restrict train as per scope
        if train_scope == "full":
            tr_use = tr
        elif train_scope == "golden":
            tr_use = tr[golden_mask[tr]]
        else:
            raise ValueError(train_scope)
        if len(tr_use) < 50 or len(te_g) < 20:
            print(f"  [{cell_name}] fold {fold}: skipped "
                  f"(n_tr={len(tr_use)}, n_te_g={len(te_g)})")
            continue
        m = train(X.iloc[tr_use], y[tr_use], n_num, n_cat)
        dte = xgb.DMatrix(X.iloc[te_g], enable_categorical=True,
                          feature_types=["q"] * n_num + ["c"] * n_cat)
        oof[te_g] = m.predict(dte)
        a = roc_auc_score(y[te_g], oof[te_g])
        fold_aucs.append(a)
        print(f"  [{cell_name}] fold {fold}: n_tr={len(tr_use):,}  "
              f"n_te_g={len(te_g):,}  AUC={a:.4f}")

    valid = ~np.isnan(oof) & golden_mask
    p = oof[valid]
    y_g = y[valid]
    auc = roc_auc_score(y_g, p)
    brier_raw = brier_score_loss(y_g, p)
    iso = IsotonicRegression(out_of_bounds="clip")
    iso.fit(p, y_g)
    pc = iso.transform(p)
    brier_cal = brier_score_loss(y_g, pc)
    ece_cal = ece(pc, y_g)
    print(f"  [{cell_name}] golden OOF AUC={auc:.4f}  "
          f"Brier_cal={brier_cal:.4f}  ECE_cal={ece_cal:.4f}  "
          f"per-fold=[{', '.join(f'{a:.4f}' for a in fold_aucs)}]")
    return dict(cell=cell_name, train_scope=train_scope,
                features=",".join(extra) or "baseline",
                auc=auc, brier_cal=brier_cal, ece_cal=ece_cal,
                n_test=int(valid.sum()),
                fold_aucs=fold_aucs)


def main() -> None:
    print(f"Loading {FEATURE_TABLE}")
    ft = pd.read_csv(FEATURE_TABLE, low_memory=False)
    rt_df = pull_measured_rt(OH_PATH)
    ft["wiki_id"] = ft.wiki_id.astype(str)
    rt_df["wiki_id"] = rt_df.wiki_id.astype(str)
    ft = ft.merge(rt_df, on="wiki_id", how="left")

    top1 = build_top1(ft).reset_index(drop=True)
    golden_iks = identify_golden_iks(top1)
    print(f"\nGolden IK14s: {len(golden_iks):,}")
    golden_mask = top1.anno_ik14.isin(golden_iks).values
    print(f"Golden rows in top-1: {golden_mask.sum():,}  "
          f"(TP only)   FP rows: {(top1.spectrum_label=='FP').sum():,}")
    # FPs are kept in TEST regardless (they're the negatives we want to discriminate)
    golden_test_mask = golden_mask | (top1.spectrum_label == "FP").values
    n_tp_test = (golden_mask & (top1.hit_label == 1).values).sum()
    n_fp_test = (top1.spectrum_label == "FP").sum()
    print(f"Golden test scope: TP={n_tp_test:,} + FP={n_fp_test:,} = "
          f"{golden_test_mask.sum():,} rows")

    side = pd.read_csv(SIDE)
    side["wiki_id"] = side.wiki_id.astype(str)
    top1 = top1.merge(
        side[["wiki_id", "hit_ik14", LOGIT]],
        on=["wiki_id", "hit_ik14"], how="left",
    )

    groups = make_groups(top1)

    print("\n" + "=" * 60)
    print("Cell TF_B: train=FULL, features=baseline")
    print("=" * 60)
    tf_b = run_cell(top1, golden_test_mask, groups, [],
                    "full", "TF_B")

    print("\n" + "=" * 60)
    print("Cell TG_B: train=GOLDEN, features=baseline")
    print("=" * 60)
    tg_b = run_cell(top1, golden_test_mask, groups, [],
                    "golden", "TG_B")

    print("\n" + "=" * 60)
    print(f"Cell TF_L: train=FULL, features=baseline+{LOGIT}")
    print("=" * 60)
    tf_l = run_cell(top1, golden_test_mask, groups, [LOGIT],
                    "full", "TF_L")

    print("\n" + "=" * 60)
    print(f"Cell TG_L: train=GOLDEN, features=baseline+{LOGIT}")
    print("=" * 60)
    tg_l = run_cell(top1, golden_test_mask, groups, [LOGIT],
                    "golden", "TG_L")

    rows = [tf_b, tg_b, tf_l, tg_l]
    summary = pd.DataFrame([
        {k: v for k, v in r.items() if k != "fold_aucs"} for r in rows
    ])
    summary.to_csv(OUT, index=False)
    print("\n" + "=" * 70)
    print("FINAL SUMMARY")
    print("=" * 70)
    cols = ["cell", "train_scope", "features", "auc", "brier_cal", "ece_cal", "n_test"]
    print(summary[cols].to_string(index=False,
                                   float_format=lambda x: f"{x:.4f}"))

    print("\nKey comparisons:")
    print(f"  TG_B vs TF_B  (golden vs full train, baseline)        : "
          f"ΔAUC = {tg_b['auc']-tf_b['auc']:+.4f}   "
          f"ΔBrier = {tg_b['brier_cal']-tf_b['brier_cal']:+.4f}")
    print(f"  TG_L vs TF_L  (golden vs full train, +{LOGIT})      : "
          f"ΔAUC = {tg_l['auc']-tf_l['auc']:+.4f}   "
          f"ΔBrier = {tg_l['brier_cal']-tf_l['brier_cal']:+.4f}")
    print(f"  TF_L vs TF_B  (struct_logit lift under FULL train)   : "
          f"ΔAUC = {tf_l['auc']-tf_b['auc']:+.4f}   "
          f"ΔBrier = {tf_l['brier_cal']-tf_b['brier_cal']:+.4f}")
    print(f"  TG_L vs TG_B  (struct_logit lift under GOLDEN train) : "
          f"ΔAUC = {tg_l['auc']-tg_b['auc']:+.4f}   "
          f"ΔBrier = {tg_l['brier_cal']-tg_b['brier_cal']:+.4f}")


if __name__ == "__main__":
    main()
