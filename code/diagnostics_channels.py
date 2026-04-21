# diagnostics_channels.py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

# ---------- helpers ----------
def _nanstrip(x):
    x = np.asarray(x, float)
    return x[np.isfinite(x)]

def _aic(loglik, k, n):
    return 2*k - 2*loglik

def _loglik_normal(x, loc=0.0, scale=1.0):
    return np.sum(stats.norm.logpdf(x, loc=loc, scale=scale))

def _loglik_laplace(x, loc=0.0, scale=1.0):
    return np.sum(stats.laplace.logpdf(x, loc=loc, scale=scale))

def _loglik_student(x, df=3, loc=0.0, scale=1.0):
    return np.sum(stats.t.logpdf(x, df, loc=loc, scale=scale))

def _qqplot(ax, data, dist='norm', dist_params=()):
    # theoretical quantiles vs sample quantiles
    data = np.sort(_nanstrip(data))
    n = len(data)
    if n < 5:
        ax.text(0.1, 0.5, "Too few points", transform=ax.transAxes)
        return
    probs = (np.arange(1, n+1) - 0.5) / n
    if dist == 'norm':
        theo = stats.norm.ppf(probs, *dist_params)
    elif dist == 'laplace':
        theo = stats.laplace.ppf(probs, *dist_params)
    elif dist == 't':
        df, loc, scale = dist_params
        theo = stats.t.ppf(probs, df, loc, scale)
    else:
        raise ValueError("Unknown dist")
    ax.plot(theo, data, marker='.', linestyle='None')
    # y=x reference
    lo = min(np.min(theo), np.min(data))
    hi = max(np.max(theo), np.max(data))
    ax.plot([lo, hi], [lo, hi], linewidth=1)
    ax.set_xlabel("Theoretical quantiles")
    ax.set_ylabel("Sample quantiles")

def _pit(ax, x, cdf):
    u = cdf(_nanstrip(x))
    ax.hist(u, bins=20, density=True)
    ax.set_xlabel("PIT values (should be uniform)")
    ax.set_ylabel("Density")
    # KS vs Uniform
    stat, p = stats.kstest(u, 'uniform')
    ax.set_title(f"PIT KS p={p:.3g}")

def _overlay_pdf(ax, x, pdfs, labels):
    # data histogram
    x = _nanstrip(x)
    ax.hist(x, bins=50, density=True, alpha=0.5)
    # overlay model pdfs
    xs = np.linspace(np.percentile(x, 1), np.percentile(x, 99), 400)
    for pdf, lab in zip(pdfs, labels):
        ax.plot(xs, pdf(xs), label=lab)
    ax.legend()
    ax.set_xlabel("Value")
    ax.set_ylabel("Density")

# ---------- MS2 transform ----------
def ms2_logit_upper_half(s, eps=1e-6):
    s = np.asarray(s, float)
    s = np.maximum(s, 0.5)
    x = (s - 0.5) / 0.5
    x = np.clip(x, eps, 1 - eps)
    return np.log(x / (1 - x))

# ---------- Channel diagnostics ----------
def diag_delta_ppm(df, label_col='label', value_col='delta_ppm', title="Δppm"):
    """Fit TP: Normal(0,σ), FP: Laplace(0,b); compare with alternatives; make plots."""
    tp = _nanstrip(df.loc[df[label_col]==1, value_col])
    fp = _nanstrip(df.loc[df[label_col]==0, value_col])
    results = {}

    # --- fit params (robust centers fixed at 0 for MVP) ---
    sigma_tp = np.std(tp, ddof=1) if len(tp)>5 else 3.0
    b_fp = np.mean(np.abs(fp)) if len(fp)>5 else 8.0
    # alts:
    # TP Student-t
    df_t = 3
    scale_t = np.sqrt(np.var(tp, ddof=1)* (df_t-2)/df_t) if len(tp)>5 else sigma_tp
    # FP Uniform over observed window
    W = np.percentile(np.abs(fp), 99.5) if len(fp)>20 else max(10.0, np.max(np.abs(fp)) if len(fp)>0 else 10.0)

    # --- log-likelihoods & AIC ---
    ll_tp_norm = _loglik_normal(tp, 0.0, max(sigma_tp, 1e-6))
    ll_tp_t    = _loglik_student(tp, df=df_t, loc=0.0, scale=max(scale_t, 1e-6))
    ll_fp_lap  = _loglik_laplace(fp, 0.0, max(b_fp, 1e-6))
    # uniform loglik
    if len(fp)>0:
        ll_fp_uni = np.sum(np.log(np.where(np.abs(fp)<=W, 1.0/(2*W), 1e-12)))
    else:
        ll_fp_uni = np.nan

    results['TP_Normal_sigma'] = sigma_tp
    results['TP_StudentT_df']  = df_t
    results['TP_StudentT_scale'] = scale_t
    results['FP_Laplace_b'] = b_fp
    results['FP_Uniform_W'] = W

    results['AIC_TP_Normal']  = _aic(ll_tp_norm, k=1, n=len(tp))
    results['AIC_TP_StudentT']= _aic(ll_tp_t,    k=2, n=len(tp))
    results['AIC_FP_Laplace'] = _aic(ll_fp_lap,  k=1, n=len(fp))
    results['AIC_FP_Uniform'] = _aic(ll_fp_uni,  k=1, n=len(fp)) if np.isfinite(ll_fp_uni) else np.nan

    # --- plots ---
    fig, axes = plt.subplots(2, 3, figsize=(12, 7))
    fig.suptitle(f"{title} diagnostics")

    # TP overlays & QQ/PIT
    if len(tp)>5:
        _overlay_pdf(
            axes[0,0], tp,
            pdfs=[
                lambda z: stats.norm.pdf(z, 0.0, max(sigma_tp, 1e-6)),
                lambda z: stats.t.pdf(z, df_t, 0.0, max(scale_t, 1e-6)),
            ],
            labels=[f"TP Normal σ={sigma_tp:.2f}", f"TP Student-t df={df_t}"]
        )
        _qqplot(axes[0,1], tp, dist='norm', dist_params=(0.0, max(sigma_tp, 1e-6)))
        _pit(axes[0,2], tp, cdf=lambda z: stats.norm.cdf(z, 0.0, max(sigma_tp, 1e-6)))
    else:
        axes[0,0].text(0.1,0.5,"Too few TP points"); axes[0,1].axis('off'); axes[0,2].axis('off')

    # FP overlays & QQ/PIT
    if len(fp)>5:
        _overlay_pdf(
            axes[1,0], fp,
            pdfs=[
                lambda z: stats.laplace.pdf(z, 0.0, max(b_fp, 1e-6)),
                lambda z: np.where(np.abs(z)<=W, 1.0/(2*W), 0.0),
            ],
            labels=[f"FP Laplace b={b_fp:.2f}", f"FP Uniform ±{W:.1f}"]
        )
        _qqplot(axes[1,1], fp, dist='laplace', dist_params=(0.0, max(b_fp, 1e-6)))
        _pit(axes[1,2], fp, cdf=lambda z: stats.laplace.cdf(z, 0.0, max(b_fp, 1e-6)))
    else:
        axes[1,0].text(0.1,0.5,"Too few FP points"); axes[1,1].axis('off'); axes[1,2].axis('off')

    plt.tight_layout()
    return results

def diag_ms2(df, label_col='label', value_col='entropy_similarity', title="MS2 (entropy sim)"):
    tp_raw = _nanstrip(df.loc[df[label_col]==1, value_col])
    fp_raw = _nanstrip(df.loc[df[label_col]==0, value_col])
    tp = ms2_logit_upper_half(tp_raw)
    fp = ms2_logit_upper_half(fp_raw)

    # fit normals on transformed scale
    mu_tp, sd_tp = (np.mean(tp), np.std(tp, ddof=1)) if len(tp)>5 else (1.2, 0.9)
    mu_fp, sd_fp = (np.mean(fp), np.std(fp, ddof=1)) if len(fp)>5 else (-0.2, 1.2)

    results = dict(mu_tp=mu_tp, sd_tp=sd_tp, mu_fp=mu_fp, sd_fp=sd_fp)

    # AIC & KS on transformed scale
    ll_tp = _loglik_normal(tp, mu_tp, max(sd_tp, 1e-6))
    ll_fp = _loglik_normal(fp, mu_fp, max(sd_fp, 1e-6))
    results['AIC_TP_Normal'] = _aic(ll_tp, k=2, n=len(tp))
    results['AIC_FP_Normal'] = _aic(ll_fp, k=2, n=len(fp))

    # plots
    fig, axes = plt.subplots(2, 3, figsize=(12, 7))
    fig.suptitle(f"{title} diagnostics (logit upper-half)")

    if len(tp)>5:
        _overlay_pdf(
            axes[0,0], tp,
            pdfs=[lambda z: stats.norm.pdf(z, mu_tp, max(sd_tp, 1e-6))],
            labels=[f"TP Normal μ={mu_tp:.2f}, σ={sd_tp:.2f}"]
        )
        _qqplot(axes[0,1], tp, dist='norm', dist_params=(mu_tp, max(sd_tp, 1e-6)))
        _pit(axes[0,2], tp, cdf=lambda z: stats.norm.cdf(z, mu_tp, max(sd_tp, 1e-6)))
    else:
        axes[0,0].text(0.1,0.5,"Too few TP points"); axes[0,1].axis('off'); axes[0,2].axis('off')

    if len(fp)>5:
        _overlay_pdf(
            axes[1,0], fp,
            pdfs=[lambda z: stats.norm.pdf(z, mu_fp, max(sd_fp, 1e-6))],
            labels=[f"FP Normal μ={mu_fp:.2f}, σ={sd_fp:.2f}"]
        )
        _qqplot(axes[1,1], fp, dist='norm', dist_params=(mu_fp, max(sd_fp, 1e-6)))
        _pit(axes[1,2], fp, cdf=lambda z: stats.norm.cdf(z, mu_fp, max(sd_fp, 1e-6)))
    else:
        axes[1,0].text(0.1,0.5,"Too few FP points"); axes[1,1].axis('off'); axes[1,2].axis('off')

    plt.tight_layout()
    return results

def diag_zrt(df, label_col='label', value_col='zrt', title="zRT"):
    tp = _nanstrip(df.loc[df[label_col]==1, value_col])
    fp = _nanstrip(df.loc[df[label_col]==0, value_col])
    results = {}

    if len(tp)>5:
        sigma_tp = np.std(tp, ddof=1)
        ll_tp = _loglik_normal(tp, 0.0, max(sigma_tp, 1e-6))
        results['TP_sigma'] = sigma_tp
        results['AIC_TP_Normal'] = _aic(ll_tp, k=1, n=len(tp))
    else:
        sigma_tp = 1.0

    if len(fp)>5:
        sigma_fp = np.std(fp, ddof=1)
        ll_fp = _loglik_normal(fp, 0.0, max(sigma_fp, 1e-6))
        results['FP_sigma'] = sigma_fp
        results['AIC_FP_Normal'] = _aic(ll_fp, k=1, n=len(fp))
    else:
        sigma_fp = 2.0

    # plots
    fig, axes = plt.subplots(2, 3, figsize=(12, 7))
    fig.suptitle(f"{title} diagnostics")

    if len(tp)>5:
        _overlay_pdf(
            axes[0,0], tp,
            pdfs=[lambda z: stats.norm.pdf(z, 0.0, max(sigma_tp, 1e-6))],
            labels=[f"TP Normal σ={sigma_tp:.2f}"]
        )
        _qqplot(axes[0,1], tp, dist='norm', dist_params=(0.0, max(sigma_tp, 1e-6)))
        _pit(axes[0,2], tp, cdf=lambda z: stats.norm.cdf(z, 0.0, max(sigma_tp, 1e-6)))
    else:
        axes[0,0].text(0.1,0.5,"Too few TP points"); axes[0,1].axis('off'); axes[0,2].axis('off')

    if len(fp)>5:
        _overlay_pdf(
            axes[1,0], fp,
            pdfs=[lambda z: stats.norm.pdf(z, 0.0, max(sigma_fp, 1e-6))],
            labels=[f"FP Normal σ={sigma_fp:.2f}"]
        )
        _qqplot(axes[1,1], fp, dist='norm', dist_params=(0.0, max(sigma_fp, 1e-6)))
        _pit(axes[1,2], fp, cdf=lambda z: stats.norm.cdf(z, 0.0, max(sigma_fp, 1e-6)))
    else:
        axes[1,0].text(0.1,0.5,"Too few FP points"); axes[1,1].axis('off'); axes[1,2].axis('off')

    plt.tight_layout()
    return results

# ---------- Example usage ----------
if __name__ == "__main__":
    # Expect a DataFrame 'cset' with columns: delta_ppm, entropy_similarity, zrt (optional),
    # and a hit-level label column 'label' (1=TP, 0=FP). Adapt this loader to your environment.
    # For your pipeline you likely do:
    # labeled = cset.assign(label=cset['is_manual_annotated'].astype(int))

    raise SystemExit("Import these functions and call diag_delta_ppm/diag_ms2/diag_zrt on your labeled DataFrame.")
