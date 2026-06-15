"""Probe: how much of HILIC retention time does the MolRex emb encode?

Two ridge regressions on the TP slice, 5-fold CV:

  Probe A: emb (512) -> measured_rt                    (Fanzhou-style)
  Probe B: emb (512) -> signed_delta_rt (per spectrum) (direct circularity test)

A high R² for A means the encoder encodes HILIC RT (Fanzhou seminar
reported 0.874 on MassWiki HILIC). A high R² for B means structure alone
predicts whether the RT predictor will be wrong on a given candidate —
which is the direct mechanism by which struct_logit could recapitulate
signed_delta_rt as a label proxy.

measured_rt is reconstructed as predicted_rt_hilic + delta_predicted_rt.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error, r2_score
from sklearn.model_selection import KFold

EMB_PATH = "data/molrex_embeddings.npz"
FT_PATH = "data/feature_table_v2.csv"
OH_PATH = "data/orbitrap_hits_v2.csv"


def cv_ridge(X: np.ndarray, y: np.ndarray, name: str, alpha: float = 1.0,
             n_splits: int = 5, seed: int = 42) -> float:
    print(f"\n--- Probe: {name} ---")
    print(f"  X={X.shape}  y range=[{y.min():.3f}, {y.max():.3f}]   "
          f"std={y.std():.3f}")
    kf = KFold(n_splits=n_splits, shuffle=True, random_state=seed)
    oof = np.zeros_like(y, dtype=np.float64)
    fold_r2 = []
    for fold, (tr, te) in enumerate(kf.split(X)):
        m = Ridge(alpha=alpha)
        m.fit(X[tr], y[tr])
        oof[te] = m.predict(X[te])
        r2 = r2_score(y[te], oof[te])
        fold_r2.append(r2)
        print(f"  fold {fold}: n_tr={len(tr):,} n_te={len(te):,}   "
              f"R²={r2:+.4f}")
    r2_oof = r2_score(y, oof)
    rmse = float(np.sqrt(((y - oof) ** 2).mean()))
    mae = mean_absolute_error(y, oof)
    print(f"  OOF R²={r2_oof:+.4f}   RMSE={rmse:.3f}   MAE={mae:.3f}   "
          f"folds=[{min(fold_r2):+.4f}, {max(fold_r2):+.4f}]")
    return r2_oof


def main() -> None:
    print(f"loading {EMB_PATH} ...")
    d = np.load(EMB_PATH, allow_pickle=True)
    ik14_arr = d["ik14"]
    emb = d["embedding"]
    ik14_to_idx = {k: i for i, k in enumerate(ik14_arr)}

    print(f"\nloading TP join keys from {FT_PATH} ...")
    ft = pd.read_csv(
        FT_PATH,
        usecols=["wiki_id", "library_wiki_id", "hit_ik14",
                 "hit_label", "spectrum_label", "signed_delta_rt"],
        dtype={"library_wiki_id": "string", "wiki_id": "string"},
    )
    tp = ft[(ft.spectrum_label == "TP") & (ft.hit_label == 1)][
        ["wiki_id", "library_wiki_id", "hit_ik14", "signed_delta_rt"]
    ].copy()
    print(f"  TP rows: {len(tp):,}   "
          f"signed_delta_rt non-null: {tp.signed_delta_rt.notna().sum():,}")

    print(f"\nstreaming RT cols from {OH_PATH} ...")
    chunks = []
    for chunk in pd.read_csv(
        OH_PATH,
        usecols=["wiki_id", "library_wiki_id",
                 "predicted_rt_hilic", "delta_predicted_rt"],
        dtype={"library_wiki_id": "string", "wiki_id": "string"},
        chunksize=200_000,
    ):
        m = chunk.merge(tp, on=["wiki_id", "library_wiki_id"], how="inner")
        chunks.append(m)
    rt_df = pd.concat(chunks, ignore_index=True)
    rt_df["measured_rt"] = (rt_df.predicted_rt_hilic
                            + rt_df.delta_predicted_rt)
    print(f"  joined rows: {len(rt_df):,}")
    print(f"  predicted_rt_hilic non-null: {rt_df.predicted_rt_hilic.notna().sum():,}")
    print(f"  delta_predicted_rt non-null: {rt_df.delta_predicted_rt.notna().sum():,}")
    print(f"  measured_rt non-null      : {rt_df.measured_rt.notna().sum():,}")

    # ============ Probe A: emb -> measured_rt (Fanzhou-style) ============
    df_a = (rt_df.dropna(subset=["measured_rt"])
                .groupby("hit_ik14")["measured_rt"]
                .median()
                .reset_index())
    df_a = df_a[df_a.hit_ik14.isin(ik14_to_idx)].reset_index(drop=True)
    print(f"\nProbe A scope: {len(df_a):,} IK14s with both emb + measured_rt")
    X_a = np.stack([emb[ik14_to_idx[k]] for k in df_a.hit_ik14]).astype(np.float32)
    y_a = df_a.measured_rt.values.astype(np.float64)
    r2_a = cv_ridge(X_a, y_a, "emb -> measured_rt")

    # =============== Probe B: emb -> signed_delta_rt ===============
    df_b = (tp.dropna(subset=["signed_delta_rt"])
              .groupby("hit_ik14")["signed_delta_rt"]
              .median()
              .reset_index())
    df_b = df_b[df_b.hit_ik14.isin(ik14_to_idx)].reset_index(drop=True)
    print(f"\nProbe B scope: {len(df_b):,} IK14s with both emb + signed_delta_rt")
    X_b = np.stack([emb[ik14_to_idx[k]] for k in df_b.hit_ik14]).astype(np.float32)
    y_b = df_b.signed_delta_rt.values.astype(np.float64)
    r2_b = cv_ridge(X_b, y_b, "emb -> signed_delta_rt")

    # =============== Interpretation ===============
    print("\n" + "=" * 70)
    print("INTERPRETATION")
    print("=" * 70)
    print(f"Fanzhou seminar baseline (MassWiki HILIC): R²=0.874, RMSE=16.91 s")
    print(f"Our Probe A (Orbitrap HILIC):              R²={r2_a:+.4f}")
    print(f"Our Probe B (emb -> signed_delta_rt):      R²={r2_b:+.4f}")
    print()
    if r2_a > 0.5:
        print(f"[A] R²={r2_a:.3f} -> encoder encodes HILIC RT. struct_logit's")
        print("    predictive signal is *potentially* circular with signed_delta_rt.")
    else:
        print(f"[A] R²={r2_a:.3f} -> encoder does NOT strongly encode HILIC RT.")
        print("    struct_logit signal is more credible as structure-only.")
    print()
    if r2_b > 0.2:
        print(f"[B] R²={r2_b:.3f} -> structure alone meaningfully predicts")
        print("    signed_delta_rt. Direct evidence of circularity for struct_logit.")
        print("    Stage 4: prioritize struct_sim_gap (RT-free) over struct_logit.")
    else:
        print(f"[B] R²={r2_b:.3f} -> structure alone does NOT strongly predict")
        print("    signed_delta_rt. struct_logit's signal is mostly NOT a")
        print("    structure-to-RT-mismatch shortcut — credible to bench.")


if __name__ == "__main__":
    main()
