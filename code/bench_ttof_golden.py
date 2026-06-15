"""Cross-platform TTOF golden bench.

Builds a minimal TTOF feature table from library_hits + spectrum metadata,
identifies golden TPs (compound has ≥2 adducts at same RT) and yy_ FPs,
scores with the current Orbitrap-trained GBM (incl. struct_logit), reports
per-polarity AUC.

Features computed (most of production NUMERIC_FEATURES):
  entropy_similarity, sim_gap, signed_delta_rt, delta_mda,
  forward_cosine, reverse_cosine, cov_count, cov_int, spectral_entropy,
  n_candidates, n_candidate_adducts, compound_has_ok_adduct,
  hit_is_isf, hit_is_dubious, hit_isf_no_ok, struct_logit (from cache)

MS2 cosines computed via build_features_v2.compute_ms2_scores on labeled-bin
candidates only (not all 737K rows) for speed.

Categoricals: hit_adduct_cat, db, polarity (all computable).
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem.inchi import InchiToInchiKey, MolToInchi
from sklearn.metrics import brier_score_loss, roc_auc_score

RDLogger.DisableLog("rdApp.*")

from bench_harness import (
    BASE_NUMERIC,
    CATEGORICAL,
    FEATURE_TABLE,
    XGB_PARAMS,
    NUM_BOOST_ROUND,
    build_top1,
    prep,
    train,
)
from build_features_v2 import compute_ms2_scores
import json
import xgboost as xgb

ROOT = Path(__file__).resolve().parent.parent
TTOF_NEG_HITS = ROOT / "data" / "library_hits" / "ttof_hilic_neg_masswiki_hits.csv"
TTOF_POS_HITS = ROOT / "data" / "library_hits" / "ttof_hilic_pos_masswiki_hits.csv"
TTOF_NEG_SPEC = ROOT / "data" / "TTOF_HILIC_negESI_uncurated_041326.csv"
TTOF_POS_SPEC = ROOT / "data" / "TTOF_HILIC_posESI_uncurated_041326.csv"
TTOF_NEG_PEAKS = ROOT / "data" / "ttof_neg_query_peaks_cache.json"
TTOF_POS_PEAKS = ROOT / "data" / "ttof_pos_query_peaks_cache.json"
LIB_PEAKS = ROOT / "data" / "library_peaks_cache.json"
ADDUCT_TAX = ROOT / "data" / "adduct_taxonomy_oliver.csv"
EMB_NPZ = ROOT / "data" / "molrex_embeddings.npz"
MOLREX_FEATURES = ROOT / "data" / "molrex_features.csv"
OUT_PATH = ROOT / "data" / "bench_ttof_golden_summary.csv"

RT_TOL_S = 10.0


def get_ik14(smi: str) -> str:
    if not isinstance(smi, str) or not smi.strip():
        return ""
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return ""
    ik = InchiToInchiKey(MolToInchi(mol))
    return ik[:14] if ik else ""


def load_adduct_taxonomy() -> dict:
    tax = pd.read_csv(ADDUCT_TAX)
    return dict(zip(tax["adduct"].str.strip(), tax["category"]))


def build_ttof_table(hits_path: Path, spec_path: Path, peaks_path: Path,
                     lib_peaks: dict, polarity: int,
                     adduct_tax: dict) -> pd.DataFrame:
    print(f"\n--- Building TTOF polarity={polarity} ---")
    hits = pd.read_csv(hits_path, low_memory=False)
    spec = pd.read_csv(spec_path, low_memory=False)
    print(f"  hits: {len(hits):,}   spec rows: {len(spec):,}")

    # Add precursor_mz, annotation_ik14, annotation status, spectral_entropy from spec
    spec["annotation_ik14"] = spec["annotation-smiles"].apply(get_ik14)
    spec["is_yy"] = spec["name"].astype(str).str.startswith("yy_")
    spec_min = spec[["wiki_id", "precursor_mz", "rt", "entropy",
                     "annotation_ik14", "annotation-adduct",
                     "name", "is_yy"]].copy()
    spec_min = spec_min.rename(columns={
        "rt": "measured_rt",
        "annotation-adduct": "annotation_adduct",
        "entropy": "spectral_entropy",
    })

    df = hits.merge(spec_min, on="wiki_id", how="left")
    df["polarity"] = polarity

    # Compute hit_ik14 from candidate's smiles
    print(f"  computing hit_ik14 ...")
    df["hit_ik14"] = df["smiles"].apply(get_ik14)
    df = df[df.hit_ik14.ne("")].reset_index(drop=True)

    # Drop candidates lacking RT predictions so top-1 selection picks among
    # rows the production model can actually use. Without this, TTOF top-1
    # by entropy_similarity systematically lands on candidates with no
    # MolRex RT coverage (production model gain rank #4 feature missing).
    before = len(df)
    df = df[df.delta_predicted_rt.notna()].reset_index(drop=True)
    print(f"  dropped {before - len(df):,} candidates without delta_predicted_rt; "
          f"kept {len(df):,}")

    # delta_mda = |precursor_mz - lib_precursor_mz| * 1000
    df["delta_mda"] = (df.precursor_mz - df.lib_precursor_mz).abs() * 1000.0
    # signed_delta_rt is precomputed in library_hits as delta_predicted_rt.
    # Matches the build_features_v2.py convention for Orbitrap.
    df["signed_delta_rt"] = df["delta_predicted_rt"]

    # Per-bin features
    per_bin = df.groupby("wiki_id").agg(
        n_candidates=("hit_ik14", "count"),
        n_candidate_adducts=("adduct", "nunique"),
    ).reset_index()
    df = df.merge(per_bin, on="wiki_id", how="left")

    # sim_gap: top - 2nd top entropy per bin
    df = df.sort_values(["wiki_id", "entropy_similarity"], ascending=[True, False])
    df["entropy_rank"] = df.groupby("wiki_id").cumcount()
    top2 = df[df.entropy_rank < 2][["wiki_id", "entropy_rank", "entropy_similarity"]]
    pv = top2.pivot(index="wiki_id", columns="entropy_rank",
                    values="entropy_similarity").reset_index()
    pv.columns = ["wiki_id", "_top1_esim", "_top2_esim"]
    pv["sim_gap"] = pv._top1_esim - pv._top2_esim.fillna(0)
    df = df.merge(pv[["wiki_id", "sim_gap"]], on="wiki_id", how="left")

    # Adduct categorization
    df["hit_adduct_cat"] = df.adduct.astype(str).str.strip().map(adduct_tax).fillna("unknown")
    df["compound_has_ok_adduct"] = (df.hit_adduct_cat == "ok").astype(int)
    df["hit_is_isf"] = (df.hit_adduct_cat == "isf").astype(int)
    df["hit_is_dubious"] = (df.hit_adduct_cat == "dubious").astype(int)
    df["hit_isf_no_ok"] = ((df.hit_is_isf == 1)
                          & (df.compound_has_ok_adduct == 0)).astype(int)

    # spectral_entropy was carried from spec (per-spectrum scalar).
    df["spectral_entropy"] = pd.to_numeric(df["spectral_entropy"], errors="coerce")

    # MS2 cosines: only compute on labeled bins (TP/FP) to save time.
    # See build_features_v2.compute_ms2_scores for definitions.
    print(f"  loading TTOF query peaks: {peaks_path}")
    with open(peaks_path) as f:
        q_peaks = json.load(f)
    print(f"  query peaks: {len(q_peaks):,}")

    # Initialize MS2 cosine columns; will fill labeled rows below.
    df["forward_cosine"] = np.nan
    df["reverse_cosine"] = np.nan
    df["cov_count"] = np.nan
    df["cov_int"] = np.nan

    # struct_logit from cache (OOF labeled set covers ~30K IK14s)
    # For TTOF candidates whose IK14 is NOT in cache, we'd need to compute
    # the inference logistic. For this bench, take the cached values where
    # available and NaN elsewhere.
    molrex = pd.read_csv(MOLREX_FEATURES)
    ik14_to_logit = dict(zip(molrex["hit_ik14"], molrex["struct_logit"]))
    df["struct_logit"] = df["hit_ik14"].map(ik14_to_logit)

    # Label assembly
    # golden TP IK14 set
    annot_subset = spec_min.dropna(subset=["annotation_ik14",
                                            "annotation_adduct",
                                            "measured_rt"])
    annot_subset = annot_subset[annot_subset.annotation_ik14.ne("")]
    agg = annot_subset.groupby("annotation_ik14").agg(
        n_bins=("wiki_id", "nunique"),
        n_adducts=("annotation_adduct", "nunique"),
        rt_min=("measured_rt", "min"),
        rt_max=("measured_rt", "max"),
    ).reset_index()
    agg["rt_spread"] = agg.rt_max - agg.rt_min
    golden_iks = set(
        agg.query("n_bins>=2 and n_adducts>=2 and rt_spread<=@RT_TOL_S")
           .annotation_ik14
    )
    print(f"  golden IK14s (compound-level): {len(golden_iks):,}")

    # Per bin: spectrum_label + hit_label
    df["spectrum_label"] = "unlabeled"
    df.loc[df.annotation_ik14.isin(golden_iks)
           & df.annotation_ik14.ne(""), "spectrum_label"] = "TP"
    df.loc[df.is_yy, "spectrum_label"] = "FP"
    df["hit_label"] = (df.hit_ik14 == df.annotation_ik14).astype(int)
    df.loc[df.spectrum_label == "FP", "hit_label"] = 0  # FP bins all hit_label=0

    # anno_ik14 for harness compat
    df["anno_ik14"] = df.annotation_ik14
    df["db"] = df["db"].astype(str)

    n_tp_bins = df[df.spectrum_label == "TP"].wiki_id.nunique()
    n_fp_bins = df[df.spectrum_label == "FP"].wiki_id.nunique()
    print(f"  TP bins: {n_tp_bins:,}   FP bins: {n_fp_bins:,}")

    # Compute MS2 cosines for labeled-bin candidates only
    labeled_mask = df.spectrum_label.isin(["TP", "FP"])
    n_labeled = int(labeled_mask.sum())
    print(f"  computing MS2 cosines for {n_labeled:,} labeled candidates ...")
    fc = df["forward_cosine"].values.copy()
    rc = df["reverse_cosine"].values.copy()
    cc = df["cov_count"].values.copy()
    ci = df["cov_int"].values.copy()
    wids = df["wiki_id"].values
    lids = df["library_wiki_id"].values
    n_hit = 0
    for idx in np.where(labeled_mask.values)[0]:
        q = q_peaks.get(str(wids[idx]))
        l = lib_peaks.get(str(lids[idx])) if isinstance(lids[idx], str) else None
        if q is None or l is None:
            continue
        f, r, c, i = compute_ms2_scores(q, l)
        fc[idx], rc[idx], cc[idx], ci[idx] = f, r, c, i
        n_hit += 1
    df["forward_cosine"] = fc
    df["reverse_cosine"] = rc
    df["cov_count"] = cc
    df["cov_int"] = ci
    print(f"    populated MS2 features for {n_hit:,}/{n_labeled:,} "
          f"({n_hit/max(n_labeled,1):.1%})")
    return df


def main() -> None:
    adduct_tax = load_adduct_taxonomy()
    print("Loading library peaks cache (335 MB) ...")
    with open(LIB_PEAKS) as f:
        lib_peaks = json.load(f)
    print(f"  library peaks: {len(lib_peaks):,} entries")
    print("\nLoading TTOF data ...")
    ttof_neg = build_ttof_table(TTOF_NEG_HITS, TTOF_NEG_SPEC, TTOF_NEG_PEAKS,
                                lib_peaks, 0, adduct_tax)
    ttof_pos = build_ttof_table(TTOF_POS_HITS, TTOF_POS_SPEC, TTOF_POS_PEAKS,
                                lib_peaks, 1, adduct_tax)
    ttof = pd.concat([ttof_neg, ttof_pos], ignore_index=True)
    print(f"\nTTOF combined: {len(ttof):,} candidate rows, "
          f"{ttof.wiki_id.nunique():,} bins")

    ttof_top1 = build_top1(ttof).reset_index(drop=True)
    n_tp = (ttof_top1.spectrum_label == "TP").sum()
    n_fp = (ttof_top1.spectrum_label == "FP").sum()
    print(f"TTOF top-1 labeled: TP={n_tp:,}  FP={n_fp:,}  "
          f"prior_TP_row={ttof_top1.hit_label.mean():.3f}")

    print("\nCoverage of feature columns on TTOF top-1:")
    for c in BASE_NUMERIC + ["struct_logit"]:
        if c in ttof_top1.columns:
            cov = ttof_top1[c].notna().mean() * 100
            print(f"  {c:24s}: {cov:5.1f}%")
        else:
            print(f"  {c:24s}: MISSING (all NaN)")
            ttof_top1[c] = np.nan

    # Train production model on Orbitrap labeled data, predict on TTOF
    print("\nLoading Orbitrap feature table for training ...")
    orb = pd.read_csv(FEATURE_TABLE, low_memory=False)
    side = pd.read_csv(MOLREX_FEATURES)
    side["wiki_id"] = side.wiki_id.astype(str)
    orb["wiki_id"] = orb.wiki_id.astype(str)
    orb = orb.merge(
        side[["wiki_id", "hit_ik14", "struct_logit"]],
        on=["wiki_id", "hit_ik14"], how="left",
    )
    orb_top1 = build_top1(orb).reset_index(drop=True)
    print(f"  Orbitrap top-1 train: {len(orb_top1):,}")

    X_train = prep(orb_top1, ["struct_logit"])
    y_train = orb_top1["hit_label"].values
    n_num = len(BASE_NUMERIC) + 1  # +struct_logit
    n_cat = len(CATEGORICAL)
    model = train(X_train, y_train, n_num, n_cat)

    X_ttof = prep(ttof_top1, ["struct_logit"])
    dttof = xgb.DMatrix(X_ttof, enable_categorical=True,
                        feature_types=["q"] * n_num + ["c"] * n_cat)
    ttof_top1["pred"] = model.predict(dttof)

    # Per-polarity AUC
    print("\n" + "=" * 60)
    print("TTOF cross-platform AUC (Orbitrap-trained model)")
    print("=" * 60)
    rows = []
    for pol_val, pol_name in [(0, "neg"), (1, "pos")]:
        m = ttof_top1.polarity == pol_val
        m_lab = m & ttof_top1.spectrum_label.isin(["TP", "FP"])
        y = ttof_top1.loc[m_lab, "hit_label"].values
        p = ttof_top1.loc[m_lab, "pred"].values
        if len(y) < 30 or len(np.unique(y)) < 2:
            print(f"  {pol_name}: skipped (n={len(y)}, labels={np.unique(y)})")
            continue
        auc = roc_auc_score(y, p)
        brier_raw = brier_score_loss(y, p)
        n_tp = int((y == 1).sum())
        n_fp = int((y == 0).sum())
        # TP score distribution
        tp_mask = y == 1
        fp_mask = y == 0
        print(f"\n  TTOF {pol_name}:")
        print(f"    labeled rows: {len(y):,}  (TP={n_tp:,}  FP={n_fp:,})")
        print(f"    AUC = {auc:.4f}   Brier_raw = {brier_raw:.4f}")
        print(f"    TP scores: median={np.median(p[tp_mask]):.3f}  "
              f"mean={p[tp_mask].mean():.3f}   "
              f"≥0.5: {(p[tp_mask]>=0.5).sum()/n_tp:.1%}   "
              f"≥0.7: {(p[tp_mask]>=0.7).sum()/n_tp:.1%}")
        print(f"    FP scores: median={np.median(p[fp_mask]):.3f}  "
              f"mean={p[fp_mask].mean():.3f}   "
              f"≥0.5: {(p[fp_mask]>=0.5).sum()/n_fp:.1%}   "
              f"≥0.7: {(p[fp_mask]>=0.7).sum()/n_fp:.1%}")
        rows.append(dict(polarity=pol_name, n=len(y), n_TP=n_tp, n_FP=n_fp,
                         AUC=auc, Brier_raw=brier_raw))

    # Combined
    m_lab = ttof_top1.spectrum_label.isin(["TP", "FP"])
    y_all = ttof_top1.loc[m_lab, "hit_label"].values
    p_all = ttof_top1.loc[m_lab, "pred"].values
    auc_all = roc_auc_score(y_all, p_all)
    rows.append(dict(polarity="all", n=len(y_all),
                     n_TP=int((y_all == 1).sum()),
                     n_FP=int((y_all == 0).sum()),
                     AUC=auc_all,
                     Brier_raw=brier_score_loss(y_all, p_all)))

    df = pd.DataFrame(rows)
    df.to_csv(OUT_PATH, index=False)
    print(f"\nWrote {OUT_PATH}")
    print("\n" + df.to_string(index=False, float_format=lambda x: f"{x:.4f}"))

    # Comparison to Orbitrap (today, with struct_logit)
    print("\n--- Comparison vs Orbitrap (today's production with struct_logit) ---")
    print("  Orbitrap neg: AUC 0.9083  (1,931 bins)")
    print("  Orbitrap pos: AUC 0.9176  (3,707 bins)")
    print("  Orbitrap all: AUC 0.9145  (5,638 bins)")
    print("\n  April-era TTOF (auto-labels, pre-fix features):")
    print("    TTOF neg: 0.838    TTOF pos: 0.850")


if __name__ == "__main__":
    main()
