"""Train and freeze the inference logistic for struct_logit deployment.

This is the artifact you'd hand Fanzhou's team for MassWiki integration.
A single L2-logistic regression on all labeled rows (no folds — folds were
only needed for OOF training-time features) mapping the 512-d MolRex
embedding to P(hit_label=1). The output logit is the same scalar
`struct_logit` the production GBM consumes.

Outputs:
  data/molrex_logistic_inference.npz    weights + bias + metadata
  data/molrex_logistic_inference.json   human-readable metadata + recipe

Validation: correlation of inference predictions vs OOF predictions on
labeled rows. Should be high (~0.95+); they're fit on the same target,
just final model uses all data instead of held-out folds.
"""
from __future__ import annotations

import json
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import logit
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parent.parent
EMB_NPZ = ROOT / "data" / "molrex_embeddings.npz"
FT_PATH = ROOT / "data" / "feature_table_v2.csv"
OOF_PATH = ROOT / "data" / "molrex_features.csv"
OUT_NPZ = ROOT / "data" / "molrex_logistic_inference.npz"
OUT_JSON = ROOT / "data" / "molrex_logistic_inference.json"

C_VALUE = 1.0  # L2 reg strength; matches OOF training in build_molrex_features.py


def main() -> None:
    print(f"Loading {EMB_NPZ}")
    d = np.load(EMB_NPZ, allow_pickle=True)
    ik14 = d["ik14"]
    emb = d["embedding"].astype(np.float32)
    ik14_to_idx = {k: i for i, k in enumerate(ik14)}
    print(f"  embeddings: {emb.shape}")

    print(f"\nLoading labeled set from {FT_PATH}")
    ft = pd.read_csv(
        FT_PATH,
        usecols=["wiki_id", "hit_ik14", "spectrum_label", "hit_label",
                 "anno_ik14", "entropy_similarity"],
    )
    labeled = ft[ft.spectrum_label.isin(["TP", "FP"])
                 & ft.anno_ik14.fillna("").ne("")]
    top1_idx = labeled.groupby("wiki_id")["entropy_similarity"].idxmax()
    top1 = ft.loc[top1_idx].reset_index(drop=True)
    print(f"  top-1 labeled bins: {len(top1):,}   "
          f"prior TP: {top1.hit_label.mean():.3f}")

    # Build X, y for labeled rows that have embeddings
    has_emb = top1.hit_ik14.map(lambda k: k in ik14_to_idx).values
    print(f"  rows with cached embedding: {int(has_emb.sum()):,}/{len(top1):,}")
    use = top1[has_emb].reset_index(drop=True)
    X = np.stack([emb[ik14_to_idx[k]] for k in use.hit_ik14]).astype(np.float32)
    y = use.hit_label.values.astype(int)
    print(f"  X={X.shape}  prior_TP={y.mean():.3f}")

    # Fit single inference model on ALL labeled rows
    print(f"\nFitting L2 logistic (C={C_VALUE}) on all labeled rows ...")
    clf = LogisticRegression(C=C_VALUE, max_iter=2000, solver="lbfgs")
    clf.fit(X, y)
    w = clf.coef_[0].astype(np.float32)
    b = float(clf.intercept_[0])
    print(f"  trained: w shape={w.shape}  b={b:+.4f}  "
          f"||w||={np.linalg.norm(w):.4f}")

    # In-sample sanity
    p_in = clf.predict_proba(X)[:, 1]
    auc_in = roc_auc_score(y, p_in)
    print(f"  in-sample AUC (optimistic): {auc_in:.4f}")

    # Validate against OOF predictions (the per-fold model from training time)
    oof_df = pd.read_csv(OOF_PATH, usecols=["wiki_id", "hit_ik14",
                                             "struct_logit"])
    use_w_oof = use.merge(oof_df, on=["wiki_id", "hit_ik14"], how="left")
    valid = use_w_oof.struct_logit.notna().values
    if valid.sum() > 0:
        p_oof = use_w_oof.struct_logit.values[valid]
        p_inf_logit = logit(np.clip(p_in[valid], 1e-6, 1 - 1e-6))
        corr = float(np.corrcoef(p_oof, p_inf_logit)[0, 1])
        print(f"\n  Correlation(OOF struct_logit vs inference logit): "
              f"r={corr:.4f}  n={int(valid.sum()):,}")
        print(f"  Mean diff (inference - OOF): "
              f"{(p_inf_logit - p_oof).mean():+.4f}  "
              f"std={(p_inf_logit - p_oof).std():.4f}")

    # Save artifact
    np.savez_compressed(
        OUT_NPZ,
        weights=w, bias=np.float32(b),
        feature_dim=np.int32(w.shape[0]),
        training_corpus=np.array(
            f"Orbitrap HILIC labeled top-1 n={len(use):,} as of "
            f"{date.today().isoformat()}", dtype=object),
        prior_TP=np.float32(y.mean()),
        in_sample_auc=np.float32(auc_in),
        regularization_C=np.float32(C_VALUE),
    )
    print(f"\nWrote {OUT_NPZ}  "
          f"({OUT_NPZ.stat().st_size / 1024:.1f} KB)")

    meta = {
        "name": "molrex_struct_logit_inference",
        "version": date.today().isoformat(),
        "encoder": "Graphormer3D (Fanzhou Kong, MolRex), d_model=512",
        "encoder_checkpoint": "data/encoder_best.pt",
        "head_type": "L2-logistic regression on raw 512-d emb",
        "head_recipe": (
            "logit_struct = w @ emb + b   (then can be used directly as a "
            "feature, or sigmoid'd to get P)"
        ),
        "feature_dim": int(w.shape[0]),
        "bias": float(b),
        "weights_norm": float(np.linalg.norm(w)),
        "regularization": f"L2, C={C_VALUE} (sklearn default)",
        "training_corpus": (
            f"Orbitrap HILIC labeled top-1 (n={len(use)}, prior_TP={y.mean():.3f}), "
            f"trained {date.today().isoformat()}"
        ),
        "in_sample_auc": float(auc_in),
        "oof_correlation": float(corr) if valid.sum() > 0 else None,
        "consumer": "production GBM feature `struct_logit` in code/score_gbm_v2.py",
        "to_apply": (
            "1) compute emb = embed_smiles(encoder_checkpoint, smiles)  "
            "2) struct_logit = float(weights @ emb + bias)  "
            "3) emit struct_logit alongside predicted_rt_hilic in MassWiki "
            "candidate response"
        ),
    }
    with open(OUT_JSON, "w") as f:
        json.dump(meta, f, indent=2)
    print(f"Wrote {OUT_JSON}")


if __name__ == "__main__":
    main()
