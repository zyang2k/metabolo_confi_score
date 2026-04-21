# gof_channels.py
# Goodness-of-fit diagnostics for Δppm and entropy_similarity channels
import sys, os
sys.path.append("/Users/ellayoung/Desktop/metabolo_confi_score/code/0918mvp")  # add your code folder to Python path

import os, numpy as np, pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import kstest, t, laplace, norm, probplot

# --- Your modules (same as in the MVP) ---
# Uses the same APIs you already use to normalize and score
from mvp_ingest import load_spectra, load_hits, build_confusable_sets  # :contentReference[oaicite:6]{index=6}
from mvp_score  import fit_channel_params, transform_entropy_similarity   # :contentReference[oaicite:7]{index=7}
from mvp_ingest import load_spectra, load_hits, join_hits_to_spectra, build_confusable_sets

# ---------------------------
# 0) CONFIG: set your paths
# ---------------------------
SPECTRA_CSV = "/Users/ellayoung/Desktop/metabolo_confi_score/data/ttof+neg+hilic.csv"
HITS_CSV    = "/Users/ellayoung/Desktop/metabolo_confi_score/data/hilic_masswiki_reference_hits.csv"
from pathlib import Path

OUTDIR = Path("/Users/ellayoung/Desktop/metabolo_confi_score/out_gof")
OUTDIR.mkdir(parents=True, exist_ok=True)

MS2_MIN     = 0.60     # same gate you used when forming confusable sets  :contentReference[oaicite:8]{index=8}

# optional gates for sanity / stability
PPM_CLIP_ABS = 20.0      # ignore extreme ppm outliers (abs <= 20)

os.makedirs(OUTDIR, exist_ok=True)
# gof_channels.py
# Goodness-of-fit diagnostics for delta_ppm and entropy_similarity channels
# using manual ground-truth labels from spec (annotation-name/adduct + is_manual_annotated)

import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt

from scipy.stats import kstest, t, laplace, norm, probplot
from scipy.optimize import minimize

# ---------------------------
# 0) CONFIG: set your paths
# ---------------------------
SPECTRA_CSV = "/Users/ellayoung/Desktop/metabolo_confi_score/data/ttof+neg+hilic.csv"
HITS_CSV    = "/Users/ellayoung/Desktop/metabolo_confi_score/data/hilic_masswiki_reference_hits.csv"
OUTDIR      = "/Users/ellayoung/Desktop/metabolo_confi_score/out_gof"
MS2_MIN     = 0.60     # same gate you used when forming confusable sets  :contentReference[oaicite:8]{index=8}
PPM_MAX     = 20.0     # optional, clip extreme ppm for stability (diagnostics still OK)

os.makedirs(OUTDIR, exist_ok=True)
# ---------------------------
# 1) Load data
# ---------------------------
spec = pd.read_csv(SPECTRA_CSV)
hits = pd.read_csv(HITS_CSV)

required_spec_cols = {
    "wiki_id", "name", "adduct", "precursor_mz",
    "annotation-name", "annotation-adduct", "is_manual_annotated"
}
missing_spec = required_spec_cols - set(spec.columns)
if missing_spec:
    raise RuntimeError(f"Missing columns in spec: {missing_spec}")

required_hit_cols = {
    "wiki_id", "db", "id", "name", "adduct", "precursor_mz",
    "entropy_similarity", "library_type", "rt", "ri", "rank"
}
missing_hits = required_hit_cols - set(hits.columns)
if missing_hits:
    # we allow a subset: at minimum need wiki_id, name, adduct, precursor_mz, entropy_similarity
    minimal = {"wiki_id","name","adduct","precursor_mz","entropy_similarity"}
    if minimal - set(hits.columns):
        raise RuntimeError(f"Missing required columns in hits: {minimal - set(hits.columns)}")

# ---------------------------
# 2) Join and compute delta_ppm
# ---------------------------
# Keep only spectra with manual annotations (ground truth)
spec_annot = spec.loc[spec["is_manual_annotated"].astype(bool)].copy()

def _norm_name(s: pd.Series) -> pd.Series:
    return s.fillna("").str.strip().str.lower()

def _norm_adduct(s: pd.Series) -> pd.Series:
    return s.fillna("").str.replace(r"\s+","", regex=True).str.upper()

# Normalize manual labels
spec_annot["manual_name_norm"]   = _norm_name(spec_annot["annotation-name"])
spec_annot["manual_adduct_norm"] = _norm_adduct(spec_annot["annotation-adduct"])

# Join: each hit for an annotated spectrum
lab = hits.merge(
    spec_annot[["wiki_id","precursor_mz","manual_name_norm","manual_adduct_norm"]],
    on="wiki_id", how="inner", validate="m:1"
).rename(columns={"precursor_mz_x":"hit_precursor_mz",
                  "precursor_mz_y":"spec_precursor_mz"})

# Normalize hit-side fields
lab["hit_name_norm"]   = _norm_name(lab["name"])
lab["hit_adduct_norm"] = _norm_adduct(lab["adduct"])

# Compute delta_ppm relative to the spectrum m/z
lab["delta_ppm"] = 1e6 * (lab["spec_precursor_mz"] - lab["hit_precursor_mz"]) / lab["spec_precursor_mz"]

# Label each hit: TP if hit name+adduct matches the manual name+adduct
lab["label"] = ((lab["hit_name_norm"] == lab["manual_name_norm"]) &
                (lab["hit_adduct_norm"] == lab["manual_adduct_norm"])).astype(int)

# Optional: clip extreme ppm outliers for stability (diagnostics still valid)
if PPM_CLIP_ABS is not None:
    lab = lab.loc[lab["delta_ppm"].abs() <= PPM_CLIP_ABS].copy()

print(f"[JOIN] Annotated spectra: {spec_annot.shape[0]}, labeled hits: {lab.shape[0]}, "
      f"TP={int((lab['label']==1).sum())}, FP={int((lab['label']==0).sum())}")

# ---------------------------
# 3) Transform for MS2 channel (upper-half logit on entropy_similarity)
# ---------------------------
def logit_upper_half(x: pd.Series, eps=1e-6) -> pd.Series:
    """
    Map x in [0.5, 1] to (0,1) via u=(x-0.5)/0.5, clip to (eps,1-eps), then z=log(u/(1-u)).
    """
    x = x.astype(float)
    mask = x >= 0.5
    u = (x[mask] - 0.5) / 0.5
    u = u.clip(eps, 1-eps)
    z = np.log(u / (1.0 - u))
    out = pd.Series(np.nan, index=x.index, dtype=float)
    out.loc[mask] = z
    return out

lab["ms2_logit_upper"] = logit_upper_half(lab["entropy_similarity"])

# ---------------------------
# 4) Helpers: QQ, posterior-predictive overlays, AIC/BIC
# ---------------------------
from scipy.stats import probplot, t as dist_t, laplace as dist_laplace, norm as dist_norm

def qq_plot(ax, data, dist, sparams=(), loc=0.0, scale=1.0, title=""):
    """
    Make a QQ-plot by standardizing data with (loc, scale) and passing only shape params to probplot.
    - dist: a scipy.stats distribution object (e.g., dist_t, dist_laplace, dist_norm)
    - sparams: tuple of shape params for that dist (e.g., (df_hat,) for t)
    """
    x = np.asarray(data, dtype=float)
    x = (x - loc) / (scale if scale > 0 else 1.0)  # standardize to the fitted location/scale

    (theo, samp), (slope, intercept, r) = probplot(x, sparams=sparams, dist=dist)

    ax.scatter(theo, samp, s=10, alpha=0.6)
    lo, hi = np.percentile(np.r_[theo, samp], [1, 99])
    ax.plot([lo, hi], [lo, hi], lw=1)
    ax.set_title(title)
    ax.set_xlabel("Theoretical quantiles")
    ax.set_ylabel("Sample quantiles")


def overlay(ax, data, sample_fn, label, bins=60):
    ax.hist(data, bins=bins, density=True, alpha=0.45, label=f"Empirical {label}")
    sim = sample_fn(len(data))
    ax.hist(sim, bins=bins, density=True, histtype="step", linewidth=1.5, label=f"Posterior-pred {label}")
    ax.legend()

def aic_bic_from_logpdf(logpdf_fn, x, k_params):
    if len(x) == 0: return (np.nan, np.nan, np.nan)
    ll = float(np.sum(logpdf_fn(x)))
    k  = int(k_params)
    n  = len(x)
    aic = 2*k - 2*ll
    bic = k*np.log(n) - 2*ll
    return ll, aic, bic

# ---------------------------
# 5) Δppm GOF: TP ~ Student-t, FP ~ Laplace
# ---------------------------
tp_ppm = lab.loc[lab["label"]==1, "delta_ppm"].dropna().to_numpy(float)
fp_ppm = lab.loc[lab["label"]==0, "delta_ppm"].dropna().to_numpy(float)

results = []

# Fit parameters by MLE (let SciPy estimate)
ppm_fig, axes = plt.subplots(2, 2, figsize=(10, 8))

if len(tp_ppm) >= 10:
    # Fit t-distribution (df, loc, scale)
    df_hat, loc_hat, scale_hat = t.fit(tp_ppm, floc=0)  # keep centered at 0
    ks_tp = kstest(tp_ppm, lambda v: t.cdf(v, df=df_hat, loc=0.0, scale=scale_hat))
    results.append(("ppm_tp_ks", ks_tp.statistic, ks_tp.pvalue, len(tp_ppm)))

    qq_plot(
    axes[0,0],
    tp_ppm,
    dist=dist_t,
    sparams=(df_hat,),
    loc=0.0,
    scale=scale_hat,
    title=f"Δppm TP QQ ~ t(df={df_hat:.1f}, scale={scale_hat:.2g})"
)

    overlay(axes[0,1], tp_ppm,
            sample_fn=lambda n: t.rvs(df=df_hat, loc=0.0, scale=scale_hat, size=n, random_state=0),
            label="TP Δppm")
    # AIC/BIC (params: df + scale => k=2; loc fixed)
    logpdf_tp = lambda x: t.logpdf(x, df=df_hat, loc=0.0, scale=scale_hat)
    ll, aic, bic = aic_bic_from_logpdf(logpdf_tp, tp_ppm, k_params=2)
    results.append(("ppm_tp_ll_aic_bic", ll, aic, bic))
else:
    axes[0,0].axis("off"); axes[0,1].axis("off")

if len(fp_ppm) >= 10:
    loc_fp, b_fp = laplace.fit(fp_ppm, floc=0.0)  # center at 0
    ks_fp = kstest(fp_ppm, lambda v: laplace.cdf(v, loc=0.0, scale=b_fp))
    results.append(("ppm_fp_ks", ks_fp.statistic, ks_fp.pvalue, len(fp_ppm)))

    qq_plot(
    axes[1,0],
    fp_ppm,
    dist=dist_laplace,
    sparams=(),
    loc=0.0,
    scale=b_fp,
    title=f"Δppm FP QQ ~ Laplace(b={b_fp:.2g})"
)

    overlay(axes[1,1], fp_ppm,
            sample_fn=lambda n: laplace.rvs(loc=0.0, scale=b_fp, size=n, random_state=0),
            label="FP Δppm")
    logpdf_fp = lambda x: laplace.logpdf(x, loc=0.0, scale=b_fp)
    ll, aic, bic = aic_bic_from_logpdf(logpdf_fp, fp_ppm, k_params=1)  # scale only; loc fixed
    results.append(("ppm_fp_ll_aic_bic", ll, aic, bic))
else:
    axes[1,0].axis("off"); axes[1,1].axis("off")


ppm_fig.tight_layout()
plt.show()
# ppm_fig.savefig(OUTDIR / "gof_delta_ppm.png", dpi=180)

# ---------------------------
# 6) MS2 GOF: use upper-half logit transform; TP/FP ~ Normal
# ---------------------------
tp_ms2_raw = lab.loc[lab["label"]==1, "entropy_similarity"].astype(float)
fp_ms2_raw = lab.loc[lab["label"]==0, "entropy_similarity"].astype(float)

tp_ms2 = logit_upper_half(tp_ms2_raw).dropna().to_numpy()
fp_ms2 = logit_upper_half(fp_ms2_raw).dropna().to_numpy()

ms2_fig, axes = plt.subplots(2, 2, figsize=(10, 8))

if len(tp_ms2) >= 20:
    mu_tp, sd_tp = float(np.mean(tp_ms2)), float(np.std(tp_ms2, ddof=1)) or 1e-6
    ks_tp = kstest(tp_ms2, lambda z: norm.cdf(z, loc=mu_tp, scale=sd_tp))
    results.append(("ms2_tp_ks", ks_tp.statistic, ks_tp.pvalue, len(tp_ms2)))

    qq_plot(
    axes[0,0],  # TP panel in your MS2 figure
    tp_ms2,
    dist=dist_norm,
    sparams=(),
    loc=mu_tp,
    scale=sd_tp,
    title=f"MS2 (logit≥0.5) TP QQ ~ N({mu_tp:.2f},{sd_tp:.2f})"
)

    overlay(axes[0,1], tp_ms2,
            sample_fn=lambda n: norm.rvs(loc=mu_tp, scale=sd_tp, size=n, random_state=0),
            label="TP logit(MS2)")
    logpdf_tp_ms2 = lambda z: norm.logpdf(z, loc=mu_tp, scale=sd_tp)
    ll, aic, bic = aic_bic_from_logpdf(logpdf_tp_ms2, tp_ms2, k_params=2)  # mu, sd
    results.append(("ms2_tp_ll_aic_bic", ll, aic, bic))
else:
    axes[0,0].axis("off"); axes[0,1].axis("off")

if len(fp_ms2) >= 20:
    mu_fp, sd_fp = float(np.mean(fp_ms2)), float(np.std(fp_ms2, ddof=1)) or 1e-6
    ks_fp = kstest(fp_ms2, lambda z: norm.cdf(z, loc=mu_fp, scale=sd_fp))
    results.append(("ms2_fp_ks", ks_fp.statistic, ks_fp.pvalue, len(fp_ms2)))

    qq_plot(
    axes[1,0],  # FP panel in your MS2 figure
    fp_ms2,
    dist=dist_norm,
    sparams=(),
    loc=mu_fp,
    scale=sd_fp,
    title=f"MS2 (logit≥0.5) FP QQ ~ N({mu_fp:.2f},{sd_fp:.2f})"
)
    overlay(axes[1,1], fp_ms2,
            sample_fn=lambda n: norm.rvs(loc=mu_fp, scale=sd_fp, size=n, random_state=0),
            label="FP logit(MS2)")
    logpdf_fp_ms2 = lambda z: norm.logpdf(z, loc=mu_fp, scale=sd_fp)
    ll, aic, bic = aic_bic_from_logpdf(logpdf_fp_ms2, fp_ms2, k_params=2)
    results.append(("ms2_fp_ll_aic_bic", ll, aic, bic))
else:
    axes[1,0].axis("off"); axes[1,1].axis("off")

ms2_fig.tight_layout()
plt.show()
# ms2_fig.savefig(OUTDIR / "gof_ms2_logit.png", dpi=180)

# ---------------------------
# 7) Save summary
# ---------------------------
# df_res = pd.DataFrame(results, columns=["metric","stat_or_ll","pvalue_or_aic","n_or_bic"])
# df_res.to_csv(OUTDIR / "gof_summary.csv", index=False)

# print("Saved outputs:")
# print(" -", OUTDIR / "gof_delta_ppm.png")
# print(" -", OUTDIR / "gof_ms2_logit.png")
# print(" -", OUTDIR / "gof_summary.csv")
