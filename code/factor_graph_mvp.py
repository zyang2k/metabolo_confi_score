"""
factor_graph_mvp.py — Factor graph confidence scoring with structural competition.

Builds on top of existing GBM scores (candidate_scores_gbm.csv).
Adds one new factor: structural competition between candidates in the same
spectrum. Candidates that are structurally similar suppress each other's
beliefs via loopy belief propagation.

Three factors per spectrum:
  1. Evidence factor  : GBM calibrated score (existing, unchanged)
  2. Competition factor: Tanimoto-weighted suppression from structural neighbors
  3. NoTA factor      : fixed prior for the no-target hypothesis

Output: data/candidate_scores_fg.csv — same schema as candidate_scores_gbm.csv
        plus fg_cal (factor graph posterior) column.

Usage:
    python code/factor_graph_mvp.py
"""

import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from sklearn.metrics import roc_auc_score
from sklearn.isotonic import IsotonicRegression


# ── Parameters ────────────────────────────────────────────────────────────────

ALPHA      = 1.0   # competition suppression strength — tune this
N_ITER     = 5     # BP iterations; converges fast for K=10-20
NOTA_PRIOR = 0.10  # prior probability mass held for NoTA
TOP_K      = 20    # cap candidates per spectrum (for speed)


# ── Fingerprint helpers ────────────────────────────────────────────────────────

def get_fp(smiles: str):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
    except Exception:
        return None


def tanimoto_matrix(smiles_list: list) -> np.ndarray:
    """Return K×K pairwise Tanimoto matrix. Diagonal = 0 (no self-competition)."""
    fps = [get_fp(s) for s in smiles_list]
    K = len(fps)
    T = np.zeros((K, K))
    for i in range(K):
        for j in range(i + 1, K):
            if fps[i] is not None and fps[j] is not None:
                sim = DataStructs.TanimotoSimilarity(fps[i], fps[j])
                T[i, j] = sim
                T[j, i] = sim
    return T


# ── Belief propagation ─────────────────────────────────────────────────────────

def run_bp(evidence: np.ndarray,
           tanimoto: np.ndarray,
           alpha: float = ALPHA,
           n_iter: int = N_ITER,
           nota_prior: float = NOTA_PRIOR) -> np.ndarray:
    """
    Loopy belief propagation over a single spectrum's candidate graph.

    evidence  : (K,) array of GBM calibrated scores in [0,1]
    tanimoto  : (K,K) pairwise structural similarity matrix
    returns   : (K,) posterior beliefs, normalized to sum ≤ (1 - nota_prior)

    Competition message from j to k:
        m_{j→k} = exp(-alpha * tanimoto[k,j] * belief[j])

    Each candidate's new belief:
        belief[k] = evidence[k] * prod_j m_{j→k}
                  = evidence[k] * exp(-alpha * sum_j tanimoto[k,j] * belief[j])

    This suppresses k proportionally to how similar and confident its
    competitors are. Structurally isolated candidates are unaffected.
    """
    K = len(evidence)
    beliefs = evidence.copy()

    for _ in range(n_iter):
        # Competition pressure on each candidate from all others
        pressure = tanimoto @ beliefs          # (K,) weighted sum of competitor beliefs
        new_beliefs = evidence * np.exp(-alpha * pressure)
        new_beliefs = np.clip(new_beliefs, 0.0, 1.0)

        # Damp updates for stability
        beliefs = 0.7 * new_beliefs + 0.3 * beliefs

    # Normalize: reserve nota_prior mass for NoTA
    total = beliefs.sum()
    if total > (1.0 - nota_prior):
        beliefs = beliefs * (1.0 - nota_prior) / total

    return beliefs


# ── Per-spectrum scoring ───────────────────────────────────────────────────────

def score_spectrum(group: pd.DataFrame) -> pd.DataFrame:
    """Apply factor graph BP to one spectrum's candidate block."""
    group = group.copy()

    # Cap at TOP_K by GBM score for speed
    if len(group) > TOP_K:
        group = group.nlargest(TOP_K, "gbm_cal")

    K = len(group)

    if K == 1:
        # No competition possible — pass through GBM score
        group["fg_cal"] = group["gbm_cal"].values
        return group

    smiles_list = group["smiles"].fillna("").tolist()
    evidence    = group["gbm_cal"].values.copy()

    T = tanimoto_matrix(smiles_list)
    posteriors = run_bp(evidence, T)

    group["fg_cal"] = posteriors
    return group


# ── Calibration ───────────────────────────────────────────────────────────────

def calibrate(scores: np.ndarray, labels: np.ndarray) -> np.ndarray:
    """Isotonic regression calibration on labeled rows."""
    mask = labels >= 0
    ir = IsotonicRegression(out_of_bounds="clip")
    ir.fit(scores[mask], labels[mask])
    return ir.transform(scores)


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    print("Loading data...")
    scores = pd.read_csv("data/candidate_scores_gbm.csv")
    hits   = pd.read_csv("data/orbitrap_hits_v2.csv", low_memory=False)

    smiles_map = (hits[["library_wiki_id", "smiles"]]
                  .dropna(subset=["smiles"])
                  .drop_duplicates("library_wiki_id")
                  .set_index("library_wiki_id")["smiles"])

    scores["smiles"] = scores["library_wiki_id"].map(smiles_map)

    print(f"Spectra: {scores['wiki_id'].nunique():,}  |  Candidates: {len(scores):,}")
    print(f"SMILES coverage: {scores['smiles'].notna().mean():.1%}")
    print()

    print("Running belief propagation...")
    results = (scores
               .groupby("wiki_id", group_keys=False)
               .apply(score_spectrum))

    # Calibrate fg_cal using labeled rows (same isotonic approach as GBM)
    labeled = results["hit_label"].isin([0, 1])
    results.loc[labeled, "fg_cal"] = calibrate(
        results.loc[labeled, "fg_cal"].values,
        results.loc[labeled, "hit_label"].values
    )

    # ── Evaluation ────────────────────────────────────────────────────────────
    eval_df = results[results["hit_label"].isin([0, 1])].copy()

    auc_gbm = roc_auc_score(eval_df["hit_label"], eval_df["gbm_cal"])
    auc_fg  = roc_auc_score(eval_df["hit_label"], eval_df["fg_cal"])

    print(f"GBM baseline AUC : {auc_gbm:.4f}")
    print(f"Factor graph AUC : {auc_fg:.4f}")
    print(f"Delta            : {auc_fg - auc_gbm:+.4f}")
    print()

    # Hard FP analysis: does competition suppress them?
    hard_fp = eval_df[(eval_df["hit_label"] == 0) & (eval_df["gbm_cal"] >= 0.7)]
    print(f"Hard FPs (gbm_cal ≥ 0.7): {len(hard_fp)}")
    print(f"  GBM mean score : {hard_fp['gbm_cal'].mean():.3f}")
    print(f"  FG  mean score : {hard_fp['fg_cal'].mean():.3f}")
    print(f"  Suppressed below 0.7: {(hard_fp['fg_cal'] < 0.7).sum()} "
          f"({100*(hard_fp['fg_cal'] < 0.7).mean():.1f}%)")
    print()

    # Near-identical pairs: does competition specifically suppress them?
    print("Competition effect by Tanimoto tier (hard FPs only):")
    print("  (requires joining back to Tanimoto — see figures below)")
    print()

    out_path = "data/candidate_scores_fg.csv"
    results.drop(columns=["smiles"]).to_csv(out_path, index=False)
    print(f"Saved → {out_path}")


if __name__ == "__main__":
    import os
    os.chdir("/Users/ellayoung/Desktop/metabolo_confi_score")
    main()
