"""
proposal_mvp.py
===============
Implements the Bayesian compound ID confidence scoring model from proposal_mvp.

Model (per proposal):
    P(compound_i | D) ∝ P(compound_i | context)
                       · P(MS2, H | compound_i)
                       · P(MS1    | compound_i)
                       · P(RT     | compound_i)   [drops out — no lib RT available]

Terms
-----
  Prior:     P(compound_i | context)  = LC-BinBase confirmed occurrence counts
             Normalized to sum to (1 - P_novel) within each candidate set.
  MS2:       λ(H)·sᵢ/Σsⱼ + (1-λ(H))·1/N   where λ(H) = 1 - exp(-α·H)
             LR_MS2 = p_ms2_i * N   (vs. uniform null 1/N)
  MS1:       N(Δppm; 0, σ_M)
             LR_MS1 = N(Δppm; 0, σ_M) / N(0; 0, σ_broad)
             where σ_broad = 50 ppm models the "novel" null mass distribution
  RT:        N(ΔRT; 0, σ_RT)  — not used here (lib_rt unavailable in dataset)
  Novel:     unnorm_novel = P_novel  (reference = 1, other terms divide by LR_novel)

Posterior:
  unnorm_i   = prior_i * LR_MS2_i * LR_MS1_i
  unnorm_nov = P_novel
  post_i     = unnorm_i / (Σ unnorm_j + unnorm_nov)
  P_novel    = unnorm_nov / (Σ unnorm_j + unnorm_nov)

Parameters fitted from confirmed (is_manual_annotated=True) LC-BinBase records:
  σ_M   — mass accuracy spread (trimmed std of Δppm from name-matched confirmed hits)
  α     — spectral complexity weight (maximize log-posterior on confirmed set)
  σ_broad — broad null, fixed at 50 ppm (instrument mass range)

Output (CSV):
  out_proposal/assertions.csv   — per-(spectrum, candidate) posteriors
  out_proposal/top_calls.csv    — best call per spectrum, or abstain if P_novel > threshold
  out_proposal/reliability.png  — calibration diagram vs. baselines
"""

import os, sys
import numpy as np
import pandas as pd
from scipy.stats import norm
from scipy.optimize import minimize_scalar
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ─── Paths ───────────────────────────────────────────────────────────────────
ROOT      = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR  = os.path.join(ROOT, "data")
OUT_DIR   = os.path.join(ROOT, "out_proposal")
os.makedirs(OUT_DIR, exist_ok=True)

SPECTRA_PATH  = os.path.join(DATA_DIR, "ttof+neg+hilic.csv")
HITS_PATH     = os.path.join(DATA_DIR, "hilic_masswiki_reference_hits.csv")
BINBASE_PATH  = os.path.join(DATA_DIR, "hilic_all.csv")

# ─── 1. Load & join ──────────────────────────────────────────────────────────
print("Loading data…")
spectra = pd.read_csv(SPECTRA_PATH, low_memory=False)
hits    = pd.read_csv(HITS_PATH,    low_memory=False)

# Rename hits columns to avoid collision after merge
hits = hits.rename(columns={
    "precursor_mz": "lib_precursor_mz",
    "rt":           "lib_rt",
    "name":         "lib_name",
})

# Keep relevant spectrum columns
spec_keep = ["wiki_id", "rt", "precursor_mz", "entropy",
             "is_manual_annotated", "name", "annotation-name"]
spec_sub  = spectra[spec_keep].copy()
spec_sub  = spec_sub.rename(columns={"name": "obs_name",
                                      "annotation-name": "annot_name"})

# Merge: every library hit gets the observed spectrum row
joint = hits.merge(spec_sub, on="wiki_id", how="inner")

# Computed channels
joint["delta_ppm"] = (
    (joint["precursor_mz"] - joint["lib_precursor_mz"])
    / joint["lib_precursor_mz"] * 1e6
)
joint["delta_rt"] = np.where(
    joint["lib_rt"].notna(),
    joint["rt"] - joint["lib_rt"],
    np.nan,
)

# De-duplicate by (wiki_id, lib_name): keep max entropy_similarity per compound per spectrum
joint = (joint
         .sort_values("entropy_similarity", ascending=False)
         .drop_duplicates(subset=["wiki_id", "lib_name"])
         .reset_index(drop=True))

print(f"  Spectra: {spec_sub['wiki_id'].nunique():,}  |  "
      f"Spectra with hits: {joint['wiki_id'].nunique():,}  |  "
      f"Total candidate rows: {len(joint):,}")


# ─── 2. Ground-truth labels (name match for confirmed spectra) ───────────────
joint["correct"] = (
    joint["lib_name"].str.strip().str.lower()
    == joint["obs_name"].str.strip().str.lower()
)
joint["confirmed"] = joint["is_manual_annotated"].fillna(False).astype(bool)


# ─── 3. Fit σ_M ──────────────────────────────────────────────────────────────
print("Fitting σ_M…")
# Use confirmed + name-matched correct hits
correct_hits = joint[joint["confirmed"] & joint["correct"]]
ppm          = correct_hits["delta_ppm"].dropna()
# Robust: MAD-based sigma, then fallback to trimmed std
ppm_trimmed  = ppm[ppm.abs() <= 15]
sigma_M      = float(ppm_trimmed.std(ddof=1)) if len(ppm_trimmed) > 10 else 5.0
print(f"  σ_M = {sigma_M:.3f} ppm  (n={len(ppm_trimmed)} trimmed correct hits)")


# ─── 4. Biological prior from LC-BinBase ─────────────────────────────────────
print("Building biological prior from LC-BinBase confirmed records…")
ha = pd.read_csv(BINBASE_PATH, usecols=["name", "target_type"], low_memory=False)
confirmed_ha = ha[
    (ha["target_type"] == "CONFIRMED")
    & (~ha["name"].str.startswith("zz ", na=False))
    & (~ha["name"].str.startswith("yy ", na=False))
].copy()
prior_counts = confirmed_ha["name"].value_counts()
prior_total  = prior_counts.sum()
EPSILON      = 1e-4   # floor for unseen compounds (relative, before per-spectrum renorm)
SIGMA_BROAD  = 50.0   # ppm: null mass distribution for novel hypothesis
P_MS1_NULL   = norm.pdf(0.0, 0.0, SIGMA_BROAD)  # = 0.00798

def get_prior(name: str) -> float:
    c = prior_counts.get(name, 0)
    return max(c / prior_total, EPSILON) if prior_total > 0 else EPSILON

print(f"  {len(prior_counts):,} unique confirmed compounds; "
      f"total observations: {prior_total:,}")


def log_likelihood_alpha(alpha: float, df_conf: pd.DataFrame,
                          sigma_M_: float, prior_fn,
                          p_novel_: float) -> float:
    """
    For each confirmed spectrum i, compute log P(correct_i | D_i; α).
    Uses same LR formulation as score_spectrum_group.
    Minimize the negative.
    """
    total_ll = 0.0
    for wid, g in df_conf.groupby("wiki_id"):
        if not g["correct"].any():
            continue

        H   = float(g["entropy"].iloc[0])
        lam = 1.0 - np.exp(-alpha * H)
        N   = len(g)
        sum_s = g["entropy_similarity"].sum()

        if sum_s > 1e-12:
            p_ms2 = lam * g["entropy_similarity"].values / sum_s + (1.0 - lam) / N
        else:
            p_ms2 = np.full(N, 1.0 / N)
        LR_ms2 = p_ms2 * N

        p_ms1  = norm.pdf(g["delta_ppm"].values, loc=0.0, scale=sigma_M_)
        LR_ms1 = p_ms1 / P_MS1_NULL

        raw_priors = np.array([prior_fn(n) for n in g["lib_name"].values])
        sum_raw    = raw_priors.sum()
        priors = (raw_priors / sum_raw * (1.0 - p_novel_)
                  if sum_raw > 1e-300
                  else np.full(N, (1.0 - p_novel_) / N))

        unnorm        = priors * LR_ms2 * LR_ms1
        unnorm_novel  = p_novel_
        total         = unnorm.sum() + unnorm_novel
        if total < 1e-300:
            continue

        unnorm_correct = unnorm[g["correct"].values]
        ll = np.log(unnorm_correct.sum() / total + 1e-300)
        total_ll += ll

    return -total_ll


# ─── 6. Estimate P(novel) ────────────────────────────────────────────────────
# Fraction of all spectra with no library hits
all_wids     = set(spectra["wiki_id"].astype(str))
hit_wids     = set(joint["wiki_id"].astype(str))
p_novel_est  = len(all_wids - hit_wids) / len(all_wids) if all_wids else 0.1
print(f"P(novel) estimate: {p_novel_est:.3f}  "
      f"({len(all_wids - hit_wids):,} spectra with no hits out of {len(all_wids):,})")


# ─── 5. Fit α (spectral complexity weight) ──────────────────────────────────
print("\nFitting α…")

# Only use confirmed spectra that have at least one correct hit
df_conf_fit = joint[joint["confirmed"] & joint["wiki_id"].isin(
    joint[joint["confirmed"] & joint["correct"]]["wiki_id"]
)].copy()

print(f"  Fitting on {df_conf_fit['wiki_id'].nunique():,} confirmed spectra "
      f"with correct hits…")

result = minimize_scalar(
    lambda a: log_likelihood_alpha(a, df_conf_fit, sigma_M, get_prior, p_novel_est),
    bounds=(0.01, 10.0), method="bounded",
    options={"xatol": 1e-4, "maxiter": 300}
)
alpha = float(result.x)
print(f"  α = {alpha:.4f}  (neg-LL = {result.fun:.2f})")


# ─── 7. Score all spectra ────────────────────────────────────────────────────
print("\nScoring all spectra…")

def score_spectrum_group(g: pd.DataFrame,
                          sigma_M_: float,
                          alpha_: float,
                          prior_fn,
                          p_novel_: float,
                          sigma_RT_: float = None) -> pd.DataFrame:
    """
    Compute full posterior for all candidates of a single spectrum.

    Uses likelihood ratios (LR) vs. the novel null hypothesis so that
    unnorm_i and unnorm_novel are on a comparable probability scale:

        unnorm_i   = prior_i * LR_MS2_i * LR_MS1_i
        unnorm_nov = p_novel_  (reference term)

    prior_i is normalized to sum to (1 - p_novel_) within the candidate set
    so that the prior alone recovers the specified P_novel in absence of evidence.
    """
    H   = float(g["entropy"].iloc[0])
    lam = 1.0 - np.exp(-alpha_ * H)
    N   = len(g)
    sum_s = g["entropy_similarity"].sum()

    # ── MS2: complexity-weighted probability, then LR vs uniform null ──────
    if sum_s > 1e-12:
        p_ms2 = lam * g["entropy_similarity"].values / sum_s + (1.0 - lam) / N
    else:
        p_ms2 = np.full(N, 1.0 / N)
    LR_ms2 = p_ms2 * N          # ratio to uniform null (1/N)

    # ── MS1: Gaussian on Δppm, LR vs broad null ────────────────────────────
    p_ms1  = norm.pdf(g["delta_ppm"].values, loc=0.0, scale=sigma_M_)
    LR_ms1 = p_ms1 / P_MS1_NULL  # ratio to σ_broad=50 ppm null

    # ── RT: drops out (no lib_rt in dataset); placeholder for future use ───
    LR_rt = np.ones(N)
    if sigma_RT_ is not None and g["delta_rt"].notna().any():
        p_rt_null = norm.pdf(0.0, 0.0, 30.0)  # broad RT null (~30 s)
        p_rt  = np.where(
            g["delta_rt"].notna(),
            norm.pdf(g["delta_rt"].values, loc=0.0, scale=sigma_RT_),
            norm.pdf(0.0, loc=0.0, scale=sigma_RT_),   # use peak when missing
        )
        LR_rt = p_rt / p_rt_null

    # ── Biological prior, normalized within candidate set ─────────────────
    raw_priors = np.array([prior_fn(n) for n in g["lib_name"].values])
    sum_raw    = raw_priors.sum()
    if sum_raw > 1e-300:
        priors = raw_priors / sum_raw * (1.0 - p_novel_)
    else:
        priors = np.full(N, (1.0 - p_novel_) / N)

    # ── Unnormalized posterior ─────────────────────────────────────────────
    unnorm       = priors * LR_ms2 * LR_ms1 * LR_rt
    unnorm_novel = p_novel_   # reference; LR_novel ≡ 1

    total = unnorm.sum() + unnorm_novel
    if total < 1e-300:
        post    = np.full(N, 1.0 / N)
        p_nova_ = p_novel_
    else:
        post    = unnorm / total
        p_nova_ = unnorm_novel / total

    out = g.copy()
    out["LR_ms2"]   = LR_ms2
    out["LR_ms1"]   = LR_ms1
    out["p_ms2"]    = p_ms2
    out["p_ms1"]    = p_ms1
    out["prior"]    = priors
    out["post"]     = post
    out["P_novel"]  = p_nova_
    out["PEP"]      = 1.0 - out["post"]
    return out


scored_parts = []
for wid, g in joint.groupby("wiki_id"):
    scored_parts.append(
        score_spectrum_group(g, sigma_M, alpha, get_prior, p_novel_est)
    )
scored = pd.concat(scored_parts, ignore_index=True)

# Mark top candidate per spectrum
scored["is_top"] = False
idx_top = scored.groupby("wiki_id")["post"].idxmax()
scored.loc[idx_top, "is_top"] = True

print(f"  Scored {scored['wiki_id'].nunique():,} spectra, "
      f"{len(scored):,} candidate rows")


# ─── 8. Top calls ────────────────────────────────────────────────────────────
NOTA_THRESHOLD = 0.6

def make_top_calls(df: pd.DataFrame, nota_thresh: float = NOTA_THRESHOLD) -> pd.DataFrame:
    rows = []
    for wid, g in df.groupby("wiki_id"):
        p_nova_ = float(g["P_novel"].iloc[0])
        abstain  = p_nova_ >= nota_thresh
        if abstain:
            rows.append({"wiki_id": wid, "lib_name_top": None, "post_top": np.nan,
                          "PEP_top": np.nan, "P_novel": p_nova_,
                          "abstain": True, "reason": "P_novel_high"})
        else:
            row = g.loc[g["post"].idxmax()]
            rows.append({
                "wiki_id":      wid,
                "lib_name_top": row["lib_name"],
                "post_top":     float(row["post"]),
                "PEP_top":      float(row["PEP"]),
                "P_novel":      p_nova_,
                "abstain":      False,
                "reason":       "",
            })
    return pd.DataFrame(rows)

top_calls = make_top_calls(scored)
print(f"  Top calls: {(~top_calls['abstain']).sum():,} called, "
      f"{top_calls['abstain'].sum():,} abstained")


# ─── 9. Assertions table ─────────────────────────────────────────────────────
assertions_cols = [
    "wiki_id", "lib_name", "entropy_similarity", "delta_ppm",
    "LR_ms2", "LR_ms1", "p_ms2", "p_ms1", "prior",
    "post", "P_novel", "PEP",
    "is_top", "correct", "confirmed", "entropy",
]
assertions = scored[[c for c in assertions_cols if c in scored.columns]].copy()


# ─── 10. Calibration / reliability diagram ──────────────────────────────────
print("\nBuilding calibration plot…")

# Use confirmed spectra top calls.
# Note: correct_top only meaningful for spectra where the correct compound is
# in the library (358 of 585 confirmed spectra). Evaluate on that subset.
conf_top = scored[scored["confirmed"] & scored["is_top"]].copy()
conf_top["correct_top"] = (
    conf_top["lib_name"].str.strip().str.lower()
    == conf_top["obs_name"].str.strip().str.lower()
)

# Restrict calibration to spectra that HAVE a correct library hit
# (avoids contaminating the diagram with impossible-to-be-correct cases)
has_correct_hit = set(joint[joint["confirmed"] & joint["correct"]]["wiki_id"])
conf_top_eval   = conf_top[conf_top["wiki_id"].isin(has_correct_hit)].copy()

N_BINS = 8
bins   = np.linspace(0, 1, N_BINS + 1)

def calibration_curve(posts, corrects, bins_):
    frac_pos, mean_pred, n_in_bin = [], [], []
    for lo, hi in zip(bins_[:-1], bins_[1:]):
        mask = (posts >= lo) & (posts < hi)
        n    = mask.sum()
        n_in_bin.append(n)
        if n < 3:
            frac_pos.append(np.nan)
            mean_pred.append(0.5 * (lo + hi))
        else:
            frac_pos.append(float(corrects[mask].mean()))
            mean_pred.append(float(posts[mask].mean()))
    return np.array(frac_pos), np.array(mean_pred), np.array(n_in_bin)

# Proposal posteriors
frac_prop, pred_prop, n_prop = calibration_curve(
    conf_top_eval["post"].values,
    conf_top_eval["correct_top"].values,
    bins,
)

# Baseline 1: raw entropy_similarity as pseudo-probability
sim_vals = conf_top_eval["entropy_similarity"].values
frac_sim, pred_sim, _ = calibration_curve(
    sim_vals,
    conf_top_eval["correct_top"].values,
    bins,
)

# Baseline 2: 1/N uniform
set_sizes = scored[scored["confirmed"]].groupby("wiki_id").size()
conf_top_eval["uniform_post"] = conf_top_eval["wiki_id"].map(
    lambda w: 1.0 / set_sizes.get(w, 1)
)
frac_unif, pred_unif, _ = calibration_curve(
    conf_top_eval["uniform_post"].values,
    conf_top_eval["correct_top"].values,
    bins,
)

fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# — Reliability diagram —
ax = axes[0]
ax.plot([0, 1], [0, 1], "k--", lw=1, label="Perfect calibration")
valid = ~np.isnan(frac_prop)
ax.plot(pred_prop[valid], frac_prop[valid], "o-", color="steelblue",
        lw=2, ms=7, label="Proposal model")
valid_s = ~np.isnan(frac_sim)
ax.plot(pred_sim[valid_s], frac_sim[valid_s], "s--", color="tomato",
        lw=1.5, ms=5, label="Entropy sim (raw)")
valid_u = ~np.isnan(frac_unif)
ax.plot(pred_unif[valid_u], frac_unif[valid_u], "^:", color="gray",
        lw=1.5, ms=5, label="Uniform 1/N")
ax.set_xlabel("Mean predicted probability")
ax.set_ylabel("Fraction correct")
ax.set_title(
    f"Reliability diagram\n(confirmed spectra with correct library hit, n={len(conf_top_eval)})"
)
ax.legend(fontsize=8)
ax.set_xlim(0, 1); ax.set_ylim(0, 1)

# — Posterior histogram: correct vs incorrect —
ax2 = axes[1]
ax2.hist(conf_top_eval.loc[ conf_top_eval["correct_top"], "post"],
         bins=20, alpha=0.65, color="steelblue", label="Correct")
ax2.hist(conf_top_eval.loc[~conf_top_eval["correct_top"], "post"],
         bins=20, alpha=0.65, color="tomato", label="Incorrect")
ax2.set_xlabel("Posterior probability")
ax2.set_ylabel("Count")
ax2.set_title("Posterior: correct vs incorrect top calls")
ax2.legend(fontsize=8)

# — λ(H) complexity weighting illustration —
ax3 = axes[2]
H_range = np.linspace(0, 4.5, 200)
for a_val in [0.5, 1.0, 2.0, alpha]:
    lam = 1 - np.exp(-a_val * H_range)
    label = f"α={a_val:.1f}" if a_val != alpha else f"α={alpha:.2f} (fitted)"
    ax3.plot(H_range, lam, label=label)
ax3.set_xlabel("Spectral entropy H")
ax3.set_ylabel("λ(H) = 1 − exp(−αH)")
ax3.set_title("MS2 complexity weighting λ(H)")
ax3.legend(fontsize=8)
ax3.axvline(x=1.0, color="gray", linestyle=":", lw=1)
ax3.text(1.05, 0.05, "H=1\n(few peaks)", fontsize=7, color="gray")

fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "reliability.png"), dpi=150)
print(f"  Saved reliability.png")

# ─── 11. Summary stats ───────────────────────────────────────────────────────
n_correct_top = conf_top["correct_top"].sum()
n_conf_total  = len(conf_top)
n_conf_eval   = len(conf_top_eval)
n_correct_eval = conf_top_eval["correct_top"].sum()
acc_proposal_all  = n_correct_top / n_conf_total if n_conf_total else 0
acc_proposal_eval = n_correct_eval / n_conf_eval  if n_conf_eval  else 0

# Baselines on evaluation set (confirmed with correct library hit)
top_sim_eval = scored[scored["confirmed"] & scored["wiki_id"].isin(has_correct_hit)].loc[
    scored[scored["confirmed"] & scored["wiki_id"].isin(has_correct_hit)]
    .groupby("wiki_id")["entropy_similarity"].idxmax()
].copy()
top_sim_eval["correct_top"] = (
    top_sim_eval["lib_name"].str.strip().str.lower()
    == top_sim_eval["obs_name"].str.strip().str.lower()
)
acc_baseline = top_sim_eval["correct_top"].mean()

# Library coverage
n_conf_total_all  = spectra["is_manual_annotated"].sum()
lib_coverage = len(has_correct_hit) / int(n_conf_total_all) if n_conf_total_all else 0

print(f"\n── Results ──────────────────────────────────────────────────────")
print(f"  Fitted params:   σ_M = {sigma_M:.3f} ppm   α = {alpha:.2f}   σ_broad = {SIGMA_BROAD:.0f} ppm")
print(f"  P(novel) = {p_novel_est:.3f}  (spectra with no library hits)")
print(f"")
print(f"  Library coverage:  {len(has_correct_hit)}/{int(n_conf_total_all)} "
      f"({lib_coverage:.1%}) confirmed spectra have correct compound in library")
print(f"")
print(f"  Top-call accuracy (all confirmed, n={n_conf_total}):")
print(f"    Proposal model:      {acc_proposal_all:.3f}  ({n_correct_top}/{n_conf_total})")
print(f"")
print(f"  Top-call accuracy (confirmed with library hit, n={n_conf_eval}):")
print(f"    Proposal model:      {acc_proposal_eval:.3f}  ({n_correct_eval}/{n_conf_eval})")
print(f"    Baseline (top sim):  {acc_baseline:.3f}")
print(f"")

# ─── 12. Save outputs ────────────────────────────────────────────────────────
assertions.to_csv(os.path.join(OUT_DIR, "assertions.csv"), index=False)
top_calls.to_csv(os.path.join(OUT_DIR, "top_calls.csv"), index=False)
print(f"  Saved assertions.csv and top_calls.csv  →  {OUT_DIR}")
print("Done.")
