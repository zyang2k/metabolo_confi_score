"""Coverage sanity check: does TG_L (golden-train + struct_logit) collapse
on non-golden compounds, or does it stay behaviorally sane?

Generates OOF predictions for both TF_B (current production design: full
train, baseline features) and TG_L (proposed: golden train, +struct_logit),
across the SAME 5-fold GroupKFold splits, predicting on the FULL test
fold (not restricted to golden). Then splits predictions by golden vs
non-golden bins and reports:

  1. Score-distribution shape (histograms by spectrum_label x slice)
  2. Agreement rate (fraction TG_L>0.5 on non-golden TPs vs TF_B)
  3. Worst-disagreement cases (|p_TG_L - p_TF_B| > 0.3, top 15)
  4. Calibration on non-golden subset (Brier_cal, ECE_cal)

Outputs CSV: data/pred_coverage_check.csv  [wiki_id, hit_ik14, name,
adduct, anno_name_lower, spectrum_label, hit_label, is_golden, p_TF_B,
p_TG_L, diff]
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import xgboost as xgb
from sklearn.isotonic import IsotonicRegression
from sklearn.metrics import brier_score_loss, roc_auc_score
from sklearn.model_selection import GroupKFold

# --- sibling-import path shim (code/ root) ---
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))

from bench_harness import (
    BASE_NUMERIC,
    CATEGORICAL,
    FEATURE_TABLE,
    build_top1,
    make_groups,
    prep,
    train,
)

ROOT = Path(__file__).resolve().parent.parent
OH_PATH = ROOT / "data" / "orbitrap_hits_v2.csv"
SIDE = ROOT / "data" / "molrex_features.csv"
OUT_PATH = ROOT / "data" / "pred_coverage_check.csv"

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


def oof_predict(top1: pd.DataFrame, golden_mask: np.ndarray,
                groups: np.ndarray, extra: list, train_scope: str) -> np.ndarray:
    """Compute OOF predictions on the FULL test fold each iteration."""
    X = prep(top1, extra)
    y = top1.hit_label.values
    n_num = len(BASE_NUMERIC) + len(extra)
    n_cat = len(CATEGORICAL)
    oof = np.full(len(top1), np.nan)
    gkf = GroupKFold(n_splits=5)
    for fold, (tr, te) in enumerate(gkf.split(top1, y, groups)):
        tr_use = tr if train_scope == "full" else tr[golden_mask[tr]]
        m = train(X.iloc[tr_use], y[tr_use], n_num, n_cat)
        dte = xgb.DMatrix(X.iloc[te], enable_categorical=True,
                          feature_types=["q"] * n_num + ["c"] * n_cat)
        oof[te] = m.predict(dte)
    return oof


def isotonize(oof: np.ndarray, y: np.ndarray, fit_mask: np.ndarray) -> np.ndarray:
    """Fit isotonic on fit_mask rows, apply everywhere."""
    iso = IsotonicRegression(out_of_bounds="clip")
    iso.fit(oof[fit_mask], y[fit_mask])
    return iso.transform(oof)


def main() -> None:
    print(f"Loading {FEATURE_TABLE}")
    ft = pd.read_csv(FEATURE_TABLE, low_memory=False)
    rt_df = pull_measured_rt(OH_PATH)
    ft["wiki_id"] = ft.wiki_id.astype(str)
    rt_df["wiki_id"] = rt_df.wiki_id.astype(str)
    ft = ft.merge(rt_df, on="wiki_id", how="left")

    top1 = build_top1(ft).reset_index(drop=True)
    golden_iks = identify_golden_iks(top1)
    print(f"Golden IK14s: {len(golden_iks):,}")
    golden_mask = top1.anno_ik14.isin(golden_iks).values
    print(f"Golden TP rows: {(golden_mask & (top1.hit_label==1).values).sum():,}")
    print(f"Non-golden TP rows: "
          f"{(~golden_mask & (top1.spectrum_label=='TP').values).sum():,}")
    print(f"FP rows (kept in both train scopes): "
          f"{(top1.spectrum_label=='FP').sum():,}")

    side = pd.read_csv(SIDE)
    side["wiki_id"] = side.wiki_id.astype(str)
    top1 = top1.merge(
        side[["wiki_id", "hit_ik14", LOGIT]],
        on=["wiki_id", "hit_ik14"], how="left",
    )

    # FPs participate in both train scopes; golden_mask is TP-side only.
    # Define `train_eligible_mask` for the GOLDEN train scope: golden-TP or any FP.
    train_eligible = golden_mask | (top1.spectrum_label == "FP").values

    groups = make_groups(top1)
    y = top1.hit_label.values

    print("\n--- Generating OOF predictions ---")
    print("TF_B (train=full, features=baseline) ...")
    oof_tfb = oof_predict(top1, np.ones(len(top1), dtype=bool), groups, [], "full")
    print("TG_L (train=golden+all FPs, features=baseline+struct_logit) ...")
    oof_tgl = oof_predict(top1, train_eligible, groups, [LOGIT], "golden")

    # Calibrate using golden test slice (where labels are clean)
    eval_golden = golden_mask | (top1.spectrum_label == "FP").values
    p_tfb = isotonize(oof_tfb, y, eval_golden)
    p_tgl = isotonize(oof_tgl, y, eval_golden)

    # ---- Diagnostics ----
    is_golden = eval_golden
    is_nongolden_tp = (~golden_mask) & (top1.spectrum_label == "TP").values
    is_fp = (top1.spectrum_label == "FP").values

    print("\n" + "=" * 70)
    print("DIAGNOSTIC 1 — Score-distribution shape")
    print("=" * 70)
    for label, mask in [
        ("golden TPs        ", is_golden & (y == 1)),
        ("golden FPs        ", is_golden & (y == 0)),
        ("non-golden TPs    ", is_nongolden_tp & (y == 1)),
        ("non-golden hit=0  ", is_nongolden_tp & (y == 0)),
    ]:
        n = int(mask.sum())
        if n == 0:
            continue
        tfb_mean, tfb_med = p_tfb[mask].mean(), np.median(p_tfb[mask])
        tgl_mean, tgl_med = p_tgl[mask].mean(), np.median(p_tgl[mask])
        tfb_std, tgl_std = p_tfb[mask].std(), p_tgl[mask].std()
        print(f"\n{label}  (n={n:,})")
        print(f"  TF_B: mean={tfb_mean:.3f}  median={tfb_med:.3f}  std={tfb_std:.3f}")
        print(f"  TG_L: mean={tgl_mean:.3f}  median={tgl_med:.3f}  std={tgl_std:.3f}")
        # coarse histogram
        hist_tfb = np.histogram(p_tfb[mask], bins=[0, 0.2, 0.4, 0.6, 0.8, 1.0])[0]
        hist_tgl = np.histogram(p_tgl[mask], bins=[0, 0.2, 0.4, 0.6, 0.8, 1.0])[0]
        print(f"  bins  [.0–.2 .2–.4 .4–.6 .6–.8 .8–1.0]")
        print(f"  TF_B   {hist_tfb}")
        print(f"  TG_L   {hist_tgl}")

    print("\n" + "=" * 70)
    print("DIAGNOSTIC 2 — Agreement rate")
    print("=" * 70)
    for label, mask in [
        ("non-golden TPs", is_nongolden_tp & (y == 1)),
        ("FPs           ", is_fp),
    ]:
        n = mask.sum()
        if n == 0:
            continue
        tfb_pos = (p_tfb[mask] > 0.5).sum()
        tgl_pos = (p_tgl[mask] > 0.5).sum()
        agree = ((p_tfb[mask] > 0.5) == (p_tgl[mask] > 0.5)).sum()
        print(f"  {label}  (n={n:,})  "
              f"TF_B p>0.5: {tfb_pos:,} ({tfb_pos/n:.1%})   "
              f"TG_L p>0.5: {tgl_pos:,} ({tgl_pos/n:.1%})   "
              f"agree: {agree:,} ({agree/n:.1%})")

    print("\n" + "=" * 70)
    print("DIAGNOSTIC 3 — Calibration on non-golden")
    print("=" * 70)
    nongolden_mask = is_nongolden_tp | is_fp
    yng = y[nongolden_mask]
    print(f"  scope: non-golden TPs ({is_nongolden_tp.sum():,}) "
          f"+ FPs ({is_fp.sum():,}) = {nongolden_mask.sum():,}")
    auc_tfb_ng = roc_auc_score(yng, p_tfb[nongolden_mask])
    auc_tgl_ng = roc_auc_score(yng, p_tgl[nongolden_mask])
    bri_tfb_ng = brier_score_loss(yng, p_tfb[nongolden_mask])
    bri_tgl_ng = brier_score_loss(yng, p_tgl[nongolden_mask])
    ece_tfb_ng = ece(p_tfb[nongolden_mask], yng)
    ece_tgl_ng = ece(p_tgl[nongolden_mask], yng)
    print(f"  TF_B: AUC={auc_tfb_ng:.4f}  Brier={bri_tfb_ng:.4f}  ECE={ece_tfb_ng:.4f}")
    print(f"  TG_L: AUC={auc_tgl_ng:.4f}  Brier={bri_tgl_ng:.4f}  ECE={ece_tgl_ng:.4f}")
    print(f"  Δ:    AUC={auc_tgl_ng-auc_tfb_ng:+.4f}  "
          f"Brier={bri_tgl_ng-bri_tfb_ng:+.4f}  "
          f"ECE={ece_tgl_ng-ece_tfb_ng:+.4f}")
    print("  (caveat: non-golden labels are noisy, so these numbers reflect "
          "label noise too)")

    print("\n" + "=" * 70)
    print("DIAGNOSTIC 4 — Worst disagreements (top 15)")
    print("=" * 70)
    disagree = pd.DataFrame({
        "wiki_id": top1.wiki_id,
        "name": top1["name"],
        "anno_name_lower": top1.anno_name_lower,
        "adduct": top1.adduct,
        "spectrum_label": top1.spectrum_label,
        "hit_label": top1.hit_label,
        "is_golden": is_golden,
        "p_TF_B": p_tfb,
        "p_TG_L": p_tgl,
        "diff": p_tgl - p_tfb,
    })
    disagree["abs_diff"] = disagree["diff"].abs()
    big = disagree[(~disagree.is_golden) & (disagree.abs_diff > 0.3)] \
        .sort_values("abs_diff", ascending=False).head(15)
    print(big[["spectrum_label", "hit_label", "name", "anno_name_lower",
              "adduct", "p_TF_B", "p_TG_L", "diff"]]
          .to_string(index=False, float_format=lambda x: f"{x:.3f}"))

    disagree.to_csv(OUT_PATH, index=False)
    print(f"\nWrote {OUT_PATH}  ({len(disagree):,} rows)")


if __name__ == "__main__":
    main()
