
# mvp_score.py
import pandas as pd
import numpy as np
from math import log, exp, isfinite
from typing import Dict, Any
from scipy.stats import norm, laplace, t

# -------- Module 03 · Feature transforms --------

def _safe_logit_ms2(s: float, eps: float = 1e-6) -> float:
    # entropy_similarity assumed in [0,1]; focus on [0.5,1]
    s = 0.5 + max(0.0, s - 0.5)  # clamp lower half to 0.5 for stability
    x = (s - 0.5) / 0.5
    x = min(max(x, eps), 1.0 - eps)
    return np.log(x / (1 - x))

def transform_entropy_similarity(series: pd.Series) -> pd.Series:
    return series.astype(float).map(_safe_logit_ms2)

# -------- Module 04 · Channel Models (fit or defaults) --------

def fit_channel_params(label_df: pd.DataFrame) -> Dict[str, Dict[str, Any]]:
    params = {}
    # estimate t-scale to match tp variance (df fixed=3 works well)
    df_t = 3.0
    tp_ppm = label_df.loc[label_df["label"]==1, "delta_ppm"].dropna().astype(float)
    var_tp = np.var(tp_ppm, ddof=1) if len(tp_ppm)>5 else 4.6795171539**2  # from your diag
    scale_t = np.sqrt(max(var_tp, 1e-6) * (df_t-2)/df_t)
    b_ppm = float(max(np.mean(np.abs(label_df.loc[label_df["label"]==0, "delta_ppm"])), 1.0)) if (label_df["label"]==0).any() else 3.3925560859
    params["ppm"] = {"model":"t", "df":df_t, "scale":scale_t, "b":b_ppm}


    # zRT: TP ~ N(0,1); FP ~ N(0, sigma_fp^2)
    tp_zrt = label_df.loc[label_df["label"] == 1, "zrt"].dropna().astype(float)
    fp_zrt = label_df.loc[label_df["label"] == 0, "zrt"].dropna().astype(float)
    sigma_fp = float(max(fp_zrt.std(ddof=1), 1.2)) if len(fp_zrt) > 5 else 2.0
    params["zrt"] = {"sigma_tp": 1.0, "sigma_fp": sigma_fp}

    # MS2 entropy similarity (logit-transformed): TP, FP ~ Normal(mu, sigma^2)
    tp_ms2 = transform_entropy_similarity(label_df.loc[label_df["label"] == 1, "entropy_similarity"].dropna())
    fp_ms2 = transform_entropy_similarity(label_df.loc[label_df["label"] == 0, "entropy_similarity"].dropna())
    mu_tp = float(tp_ms2.mean()) if len(tp_ms2) > 5 else 1.2
    sd_tp = float(max(tp_ms2.std(ddof=1), 0.5)) if len(tp_ms2) > 5 else 0.9
    mu_fp = float(fp_ms2.mean()) if len(fp_ms2) > 5 else -0.2
    sd_fp = float(max(fp_ms2.std(ddof=1), 0.6)) if len(fp_ms2) > 5 else 1.2
    params["ms2"] = {"mu_tp": mu_tp, "sd_tp": sd_tp, "mu_fp": mu_fp, "sd_fp": sd_fp}
    return params

# -------- Module 04b · Per-hit log LR --------

def per_hit_logLR(df: pd.DataFrame, params: Dict[str, Dict[str, Any]]) -> pd.Series:
    out = np.zeros(len(df), dtype=float)

    # ppm channel
    if "delta_ppm" in df.columns:
        x = df["delta_ppm"].astype(float).to_numpy()
        if params["ppm"].get("model") == "t":
            log_num = t.logpdf(x, df=params["ppm"]["df"], loc=0.0, scale=params["ppm"]["scale"])
        else:
            log_num = norm.logpdf(x, loc=0.0, scale=params["ppm"]["sigma"])
        log_den = laplace.logpdf(x, loc=0.0, scale=params["ppm"]["b"])
        out += (log_num - log_den)



    # zRT channel (skip when NaN)
    if "zrt" in df.columns:
        x = df["zrt"].astype(float).to_numpy()
        # TP N(0,1); FP N(0,sigma_fp^2)
        log_num = norm.logpdf(x, loc=0.0, scale=params["zrt"]["sigma_tp"])
        log_den = norm.logpdf(x, loc=0.0, scale=params["zrt"]["sigma_fp"])
        # NaNs -> 0 contribution
        add = (log_num - log_den)
        add[~np.isfinite(x)] = 0.0
        out += add

    # MS2 channel
    if "entropy_similarity" in df.columns:
        x = transform_entropy_similarity(df["entropy_similarity"]).to_numpy()
        log_num = norm.logpdf(x, loc=params["ms2"]["mu_tp"], scale=params["ms2"]["sd_tp"])
        log_den = norm.logpdf(x, loc=params["ms2"]["mu_fp"], scale=params["ms2"]["sd_fp"])
        out += (log_num - log_den)

    return pd.Series(out, index=df.index, name="logLR_hit")

# -------- Module 05 · Aggregate duplicates → structure score --------

def aggregate_to_structure(df: pd.DataFrame, quality_cols=("entropy_similarity",), topK=None, lambda_var: float = 0.0) -> pd.DataFrame:
    # df must contain: wiki_id, library_id, name, logLR_hit
    req = ["wiki_id","library_id","name","logLR_hit"]
    missing = [c for c in req if c not in df.columns]
    if missing:
        raise ValueError(f"aggregate_to_structure: missing {missing}")

    # quality weights α ∝ mean(quality_cols)
    if quality_cols and all(c in df.columns for c in quality_cols):
        q = df[list(quality_cols)].astype(float).mean(axis=1)
        df = df.assign(_q=q)
    else:
        df = df.assign(_q=1.0)

    # Optionally select topK hits by logLR within each (wiki_id, library_id)
    def _agg(g):
        if topK is not None and len(g) > topK:
            g = g.nlargest(topK, "logLR_hit")
        # weights
        w = g["_q"].to_numpy()
        w = np.maximum(w, 1e-6)
        w = w / w.sum()
        # log-sum-exp with weights
        vals = g["logLR_hit"].to_numpy()
        m = np.max(vals)
        lse = m + np.log(np.sum(w * np.exp(vals - m)))
        penalty = 0.0
        if lambda_var > 0.0 and len(vals) > 1:
            penalty = lambda_var * np.var(vals, ddof=0)
        return pd.Series({
            "theta": lse - penalty,
            "member_hits": len(g),
            "quality_mean": float(g["_q"].mean())
        })

    agg = df.groupby(["wiki_id","library_id","name"], as_index=False).apply(_agg).reset_index(drop=True)
    return agg

# -------- Module 05b · Local rank probability --------

def local_rank_probability(struct_scores: pd.DataFrame) -> pd.DataFrame:
    req = ["wiki_id","library_id","theta"]
    if any(c not in struct_scores.columns for c in req):
        raise ValueError("local_rank_probability: missing required columns")
    def _softmax_block(g):
        x = g["theta"].to_numpy()
        m = np.max(x)
        ex = np.exp(x - m)
        s = ex.sum()
        cloc = ex / s
        # delta to runner-up
        order = np.argsort(-x)
        delta = float(x[order[0]] - x[order[1]]) if len(x) > 1 else np.nan
        out = g.copy()
        out["C_loc"] = cloc
        out["delta_theta"] = delta
        out["set_size"] = len(g)
        return out
    return struct_scores.groupby("wiki_id", as_index=False, group_keys=False).apply(_softmax_block)

# -------- Module 06 · Global posterior with NoTA --------

def compute_global_posterior(struct_scores: pd.DataFrame, prior_policy: str = "uniform", lr_nota: float = 0.5, prior_nota: float = 0.5) -> pd.DataFrame:
    if "theta" not in struct_scores.columns:
        raise ValueError("compute_global_posterior: theta missing")
    df = struct_scores.copy()
    df["LR"] = np.exp(df["theta"].astype(float))

    def _block(g):
        n = len(g)
        if prior_policy == "uniform":
            priors = np.full(n, 1.0 / n)
        else:
            priors = np.full(n, 1.0 / n)
        numer = priors * g["LR"].to_numpy()
        p_not_a = prior_nota * lr_nota
        denom = numer.sum() + p_not_a
        post = numer / denom
        p_nota = p_not_a / denom
        idx_top = int(np.argmax(post))
        pep = 1.0 - float(post[idx_top])
        out = g.copy()
        out["post"] = post
        out["P_NoTA"] = p_nota
        out["PEP_local"] = 1.0 - out["post"]  # per-structure
        # annotate top within block
        out["is_top"] = False
        out.iloc[idx_top, out.columns.get_loc("is_top")] = True
        out["abstain"] = p_nota >= 0.6
        out["reason"] = np.where(out["abstain"], "NoTA", "")
        return out

    return df.groupby("wiki_id", as_index=False, group_keys=False).apply(_block)
