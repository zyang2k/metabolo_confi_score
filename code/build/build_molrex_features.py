"""Stage 3: derive scalar MolRex features for the GBM bench.

Three features per top-1 row:
  struct_sim_gap       cos(emb[top1], emb[rank2])  within bin  (NO OOF needed)
  struct_logit         5-fold OOF logistic on raw 512-d -> logit P(label=1)
  struct_centroid_dist 1 - cos(emb, TP-centroid_train_fold)   per-fold

Output:
  data/molrex_features.csv   [wiki_id, hit_ik14, struct_sim_gap,
                              struct_logit, struct_centroid_dist]

Diagnostics printed at end:
  - per-feature TP vs FP mean/std (does it discriminate?)
  - univariate ROC AUC per feature (single-feature label predictor)
  - Pearson correlation with existing top GBM features
    (entropy_similarity, sim_gap, signed_delta_rt)
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import logit
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold

EMB_NPZ = Path("data/molrex_embeddings.npz")
FT_PATH = Path("data/feature_table_v2.csv")
OUT_CSV = Path("data/molrex_features.csv")

EPS = 1e-12


def _cos(a: np.ndarray, b: np.ndarray) -> float:
    return float(a @ b / (np.linalg.norm(a) * np.linalg.norm(b) + EPS))


def main() -> None:
    print(f"loading {EMB_NPZ} ...")
    d = np.load(EMB_NPZ, allow_pickle=True)
    ik14_arr = d["ik14"]
    emb_mat = d["embedding"]  # (N, 512) float32
    ik14_to_idx: dict[str, int] = {ik: i for i, ik in enumerate(ik14_arr)}
    print(f"  {len(ik14_arr):,} IK14s, emb {emb_mat.shape}")

    print(f"\nloading {FT_PATH} ...")
    ft = pd.read_csv(
        FT_PATH,
        usecols=["wiki_id", "hit_ik14", "anno_ik14", "entropy_similarity",
                 "spectrum_label", "hit_label", "rank",
                 "sim_gap", "signed_delta_rt"],
    )
    labeled_mask = ft.spectrum_label.isin(["TP", "FP"]) & \
                   ft.anno_ik14.fillna("").ne("")
    ft_labeled = ft[labeled_mask].copy()
    print(f"  labeled rows: {len(ft_labeled):,}  "
          f"({ft_labeled.wiki_id.nunique():,} bins)")

    # ---- top-1 per bin (matches score_gbm_v2.py) ----
    top1_idx = ft_labeled.groupby("wiki_id")["entropy_similarity"].idxmax()
    top1 = ft.loc[top1_idx].reset_index(drop=True)
    labels = top1["hit_label"].values.astype(int)
    print(f"  top1_train rows: {len(top1):,}   prior={labels.mean():.3f}")

    # ---- struct_sim_gap: per-bin cos(emb_top1, emb_rank2) ----
    print("\ncomputing struct_sim_gap ...")
    sim_gap_struct = np.full(len(top1), np.nan, dtype=np.float32)
    n_single = 0
    n_missing = 0
    # Pre-sort all labeled candidates by entropy_similarity descending within bin
    ft_sorted = ft_labeled.sort_values(
        ["wiki_id", "entropy_similarity"], ascending=[True, False]
    )
    grouped = ft_sorted.groupby("wiki_id", sort=False)
    rank2_lookup: dict[str, str | None] = {}
    for wiki_id, grp in grouped:
        if len(grp) < 2:
            rank2_lookup[wiki_id] = None
            continue
        rank2_lookup[wiki_id] = grp.iloc[1]["hit_ik14"]

    for i, row in top1.iterrows():
        wiki_id = row["wiki_id"]
        top1_ik = row["hit_ik14"]
        rank2_ik = rank2_lookup.get(wiki_id)
        if rank2_ik is None:
            n_single += 1
            continue
        t_idx = ik14_to_idx.get(top1_ik)
        r_idx = ik14_to_idx.get(rank2_ik)
        if t_idx is None or r_idx is None:
            n_missing += 1
            continue
        sim_gap_struct[i] = _cos(emb_mat[t_idx], emb_mat[r_idx])
    n_valid = int(np.isfinite(sim_gap_struct).sum())
    print(f"  valid: {n_valid:,}/{len(top1):,}  "
          f"(single-candidate bins: {n_single:,}, "
          f"missing IK14: {n_missing:,})")

    # ---- assemble emb_top1 matrix ----
    top1_emb = np.zeros((len(top1), emb_mat.shape[1]), dtype=np.float32)
    has_emb = np.zeros(len(top1), dtype=bool)
    for i, ik in enumerate(top1["hit_ik14"]):
        idx = ik14_to_idx.get(ik)
        if idx is not None:
            top1_emb[i] = emb_mat[idx]
            has_emb[i] = True
    print(f"\ntop1 rows with cached embedding: "
          f"{int(has_emb.sum()):,}/{len(top1):,}")

    # ---- struct_logit + struct_centroid_dist: 5-fold OOF ----
    print("\nfitting 5-fold OOF (same GroupKFold(anno_ik14) as GBM) ...")
    groups = top1["anno_ik14"].fillna("").values.astype(object).copy()
    for i in range(len(groups)):
        if groups[i] == "":
            groups[i] = f"__no_ik14_{i}"
    gkf = GroupKFold(n_splits=5)

    struct_logit = np.full(len(top1), np.nan, dtype=np.float32)
    struct_cent = np.full(len(top1), np.nan, dtype=np.float32)
    fold_aucs_logit = []
    for fold, (tr, te) in enumerate(gkf.split(top1, labels, groups)):
        tr_use = tr[has_emb[tr]]
        te_use = te[has_emb[te]]
        if len(tr_use) == 0 or len(te_use) == 0:
            print(f"  fold {fold}: skipped (no embedded rows)")
            continue
        X_tr, y_tr = top1_emb[tr_use], labels[tr_use]
        X_te = top1_emb[te_use]

        # struct_logit: L2-regularized logistic on raw 512-d
        clf = LogisticRegression(C=1.0, max_iter=2000, solver="lbfgs")
        clf.fit(X_tr, y_tr)
        p_te = clf.predict_proba(X_te)[:, 1]
        p_te = np.clip(p_te, 1e-6, 1 - 1e-6)
        struct_logit[te_use] = logit(p_te).astype(np.float32)
        fold_auc = roc_auc_score(labels[te_use], p_te)
        fold_aucs_logit.append(fold_auc)

        # struct_centroid_dist: 1 - cos(test_emb, TP_train_centroid)
        tp_mask = y_tr == 1
        if tp_mask.sum() < 2:
            print(f"  fold {fold}: <2 TP in train, skipping centroid")
        else:
            centroid = X_tr[tp_mask].mean(axis=0)
            cent_norm = np.linalg.norm(centroid) + EPS
            te_norms = np.linalg.norm(X_te, axis=1) + EPS
            cos_te = (X_te @ centroid) / (te_norms * cent_norm)
            struct_cent[te_use] = (1.0 - cos_te).astype(np.float32)

        print(f"  fold {fold}: n_tr={len(tr_use):,} n_te={len(te_use):,}  "
              f"logit_AUC={fold_auc:.4f}")

    n_valid_logit = int(np.isfinite(struct_logit).sum())
    n_valid_cent = int(np.isfinite(struct_cent).sum())
    print(f"\nstruct_logit valid: {n_valid_logit:,}   "
          f"struct_centroid_dist valid: {n_valid_cent:,}")
    print(f"per-fold logit AUC: "
          f"mean={np.mean(fold_aucs_logit):.4f}  "
          f"[{min(fold_aucs_logit):.4f}, {max(fold_aucs_logit):.4f}]")

    # ---- assemble + save ----
    out = pd.DataFrame({
        "wiki_id": top1["wiki_id"].values,
        "hit_ik14": top1["hit_ik14"].values,
        "struct_sim_gap": sim_gap_struct,
        "struct_logit": struct_logit,
        "struct_centroid_dist": struct_cent,
    })
    out.to_csv(OUT_CSV, index=False)
    print(f"\nwrote {OUT_CSV}  ({len(out):,} rows)")

    # =========================== DIAGNOSTICS ===========================
    print("\n" + "=" * 60)
    print("DIAGNOSTICS")
    print("=" * 60)

    for col in ["struct_sim_gap", "struct_logit", "struct_centroid_dist"]:
        v = out[col].values
        ok = np.isfinite(v)
        v_tp = v[ok & (labels == 1)]
        v_fp = v[ok & (labels == 0)]
        if len(v_tp) and len(v_fp):
            uni_auc = roc_auc_score(labels[ok], v[ok])
        else:
            uni_auc = np.nan
        print(f"\n{col}:")
        print(f"  n_valid: {ok.sum():,}/{len(v):,}")
        print(f"  TP (n={len(v_tp):,}):  "
              f"mean={np.nanmean(v_tp):+.4f}  std={np.nanstd(v_tp):.4f}")
        print(f"  FP (n={len(v_fp):,}):  "
              f"mean={np.nanmean(v_fp):+.4f}  std={np.nanstd(v_fp):.4f}")
        print(f"  univariate AUC (single-feature label predictor): "
              f"{uni_auc:.4f}")

    print("\nPearson corr with existing top GBM features:")
    top1_corr = top1[["entropy_similarity", "sim_gap", "signed_delta_rt"]].copy()
    for col in ["struct_sim_gap", "struct_logit", "struct_centroid_dist"]:
        top1_corr[col] = out[col].values
    corr_block = top1_corr.corr().loc[
        ["struct_sim_gap", "struct_logit", "struct_centroid_dist"],
        ["entropy_similarity", "sim_gap", "signed_delta_rt"],
    ]
    print(corr_block.round(3).to_string())


if __name__ == "__main__":
    main()
