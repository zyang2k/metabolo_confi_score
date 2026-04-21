"""
bayesian_score_v2.py — Per-channel Bayesian logLR scoring engine.

Scores (spectrum, candidate) pairs. Per-channel TP/FP density estimation →
log likelihood ratio → sum → sigmoid → calibrated posterior.

Usage:
    from bayesian_score_v2 import fit_all_channels, score, calibrate
"""

import numpy as np
import pandas as pd
from scipy import stats
from dataclasses import dataclass, field
from typing import Optional, Callable, List, Dict
from sklearn.isotonic import IsotonicRegression


# ── Channel specification ───────────────────────────────────────────────────

@dataclass
class ChannelSpec:
    """Specification for one scoring channel."""
    name: str
    feature_col: str
    channel_type: str       # 'continuous', 'binary', 'discrete'
    higher_means_tp: bool   # True if higher values → more likely TP
    tp_family: str = ''     # e.g. 'normal', 'lognormal', 'beta' (set after EDA)
    fp_family: str = ''
    transform: Optional[Callable] = None  # applied before fitting
    clamp: tuple = (-5.0, 5.0)


@dataclass
class FittedChannel:
    """Fitted TP and FP distributions for one channel."""
    spec: ChannelSpec
    tp_dist: object = None     # scipy.stats frozen distribution or dict
    fp_dist: object = None
    tp_n: int = 0
    fp_n: int = 0

    def logLR(self, x: np.ndarray) -> np.ndarray:
        """Compute clamped log-likelihood ratio for array of feature values."""
        x = np.asarray(x, dtype=float)
        result = np.zeros_like(x)
        valid = np.isfinite(x)

        if valid.sum() == 0:
            return result

        xv = x[valid]

        if self.spec.transform is not None:
            xv = self.spec.transform(xv)
            finite_after = np.isfinite(xv)
            # If transform produces non-finite, those contribute 0
            tmp = np.zeros_like(xv)
            if finite_after.sum() > 0:
                tmp[finite_after] = self._raw_logLR(xv[finite_after])
            result[valid] = tmp
        else:
            result[valid] = self._raw_logLR(xv)

        # Clamp
        lo, hi = self.spec.clamp
        result = np.clip(result, lo, hi)
        return result

    def _raw_logLR(self, x: np.ndarray) -> np.ndarray:
        """Unclamped logLR dispatch."""
        if self.spec.channel_type == 'continuous':
            return self._continuous_logLR(x)
        elif self.spec.channel_type == 'binary':
            return self._binary_logLR(x)
        elif self.spec.channel_type == 'discrete':
            return self._discrete_logLR(x)
        return np.zeros_like(x)

    def _continuous_logLR(self, x: np.ndarray) -> np.ndarray:
        tp_ll = self.tp_dist.logpdf(x)
        fp_ll = self.fp_dist.logpdf(x)
        return tp_ll - fp_ll

    def _binary_logLR(self, x: np.ndarray) -> np.ndarray:
        # tp_dist and fp_dist are dicts: {0: log_prob, 1: log_prob}
        result = np.zeros_like(x)
        for val in [0, 1]:
            mask = x == val
            if mask.sum() > 0:
                result[mask] = self.tp_dist[val] - self.fp_dist[val]
        return result

    def _discrete_logLR(self, x: np.ndarray) -> np.ndarray:
        # tp_dist and fp_dist are dicts: {value: log_prob}
        result = np.zeros_like(x)
        tp_default = self.tp_dist.get('_default', -10)
        fp_default = self.fp_dist.get('_default', -10)
        for val in set(list(self.tp_dist.keys()) + list(self.fp_dist.keys())):
            if val == '_default':
                continue
            mask = x == val
            if np.any(mask):
                tp_lp = self.tp_dist.get(val, tp_default)
                fp_lp = self.fp_dist.get(val, fp_default)
                result[mask] = tp_lp - fp_lp
        return result


# ── Fitting functions ───────────────────────────────────────────────────────

def _fit_continuous(values: np.ndarray, family: str):
    """Fit a continuous distribution. Returns a scipy frozen distribution."""
    values = values[np.isfinite(values)]
    if len(values) < 5:
        return None

    if family == 'normal':
        mu, sigma = np.mean(values), max(np.std(values, ddof=1), 1e-6)
        return stats.norm(loc=mu, scale=sigma)

    elif family == 'half_normal':
        # Half-normal for non-negative data (like delta_mda)
        sigma = max(np.sqrt(np.mean(values ** 2)), 1e-6)  # MLE for half-normal
        return stats.halfnorm(loc=0, scale=sigma)

    elif family == 'exponential':
        lam = max(np.mean(values), 1e-6)
        return stats.expon(loc=0, scale=lam)

    elif family == 'student_t':
        # Fix df=3 for heavy tails, fit loc and scale
        loc = np.median(values)
        scale = max(np.std(values, ddof=1) * 0.5, 1e-6)
        return stats.t(df=3, loc=loc, scale=scale)

    elif family == 'kde':
        # Gaussian KDE — nonparametric
        try:
            kernel = stats.gaussian_kde(values)
            # Wrap to have logpdf method
            class KDEDist:
                def __init__(self, kde):
                    self._kde = kde
                def logpdf(self, x):
                    pdf = self._kde(x)
                    return np.log(np.maximum(pdf, 1e-300))
            return KDEDist(kernel)
        except Exception:
            return None

    else:
        raise ValueError(f'Unknown family: {family}')


def _fit_binary(values: np.ndarray, alpha: float = 1.0):
    """Fit Bernoulli with Laplace smoothing. Returns {0: log_p, 1: log_p}."""
    values = values[np.isfinite(values)]
    n1 = np.sum(values == 1)
    n0 = np.sum(values == 0)
    n = n0 + n1
    p1 = (n1 + alpha) / (n + 2 * alpha)
    p0 = 1 - p1
    return {0: np.log(max(p0, 1e-10)), 1: np.log(max(p1, 1e-10))}


def _fit_discrete(values: np.ndarray, alpha: float = 1.0, max_val: int = 10):
    """Fit empirical PMF with Laplace smoothing. Returns {val: log_prob}."""
    values = values[np.isfinite(values)].astype(int)
    # Bin values >= max_val into max_val
    values = np.minimum(values, max_val)
    unique_vals = sorted(set(values))
    n = len(values)
    n_categories = len(unique_vals)
    pmf = {}
    for v in unique_vals:
        count = np.sum(values == v)
        pmf[v] = np.log((count + alpha) / (n + n_categories * alpha))
    # Default for unseen values
    pmf['_default'] = np.log(alpha / (n + n_categories * alpha))
    return pmf


def fit_channel(spec: ChannelSpec,
                feature_values: np.ndarray,
                labels: np.ndarray) -> FittedChannel:
    """Fit TP and FP distributions for one channel."""
    tp_mask = labels == 1
    fp_mask = labels == 0
    tp_vals = feature_values[tp_mask]
    fp_vals = feature_values[fp_mask]

    # Apply transform
    if spec.transform is not None:
        tp_vals = spec.transform(tp_vals.copy())
        fp_vals = spec.transform(fp_vals.copy())

    fc = FittedChannel(spec=spec, tp_n=int(tp_mask.sum()), fp_n=int(fp_mask.sum()))

    if spec.channel_type == 'continuous':
        fc.tp_dist = _fit_continuous(tp_vals, spec.tp_family)
        fc.fp_dist = _fit_continuous(fp_vals, spec.fp_family)
        if fc.tp_dist is None or fc.fp_dist is None:
            return None
    elif spec.channel_type == 'binary':
        fc.tp_dist = _fit_binary(tp_vals)
        fc.fp_dist = _fit_binary(fp_vals)
    elif spec.channel_type == 'discrete':
        fc.tp_dist = _fit_discrete(tp_vals)
        fc.fp_dist = _fit_discrete(fp_vals)

    return fc


def fit_all_channels(feature_table: pd.DataFrame,
                     channel_specs: List[ChannelSpec],
                     label_col: str = 'hit_label') -> List[FittedChannel]:
    """Fit all channels. Returns list of FittedChannel objects."""
    labels = feature_table[label_col].values
    fitted = []
    for spec in channel_specs:
        vals = feature_table[spec.feature_col].values
        fc = fit_channel(spec, vals, labels)
        if fc is not None:
            fitted.append(fc)
            print(f'  {spec.name:30s}: fitted (TP n={fc.tp_n}, FP n={fc.fp_n})')
        else:
            print(f'  {spec.name:30s}: SKIPPED (insufficient data)')
    return fitted


# ── Scoring ─────────────────────────────────────────────────────────────────

def compute_logLR_table(feature_table: pd.DataFrame,
                        fitted_channels: List[FittedChannel]) -> pd.DataFrame:
    """Compute per-channel clamped logLR for every row."""
    result = feature_table[['wiki_id', 'hit_ik14', 'hit_label']].copy()

    total = np.zeros(len(feature_table))
    for fc in fitted_channels:
        col = fc.spec.feature_col
        vals = feature_table[col].values
        lr = fc.logLR(vals)
        result[f'logLR_{fc.spec.name}'] = lr
        total += lr

    result['total_logLR'] = total
    return result


def compute_posterior(logLR_table: pd.DataFrame,
                      prior_tp: float) -> pd.DataFrame:
    """Convert total logLR to posterior probability."""
    prior_log_odds = np.log(prior_tp / (1 - prior_tp))
    log_post_odds = prior_log_odds + logLR_table['total_logLR'].values

    posterior = 1.0 / (1.0 + np.exp(-log_post_odds))

    out = logLR_table.copy()
    out['posterior'] = posterior
    out['prior_log_odds'] = prior_log_odds
    return out


# ── Calibration ─────────────────────────────────────────────────────────────

def calibrate_isotonic(posteriors: np.ndarray,
                       labels: np.ndarray) -> IsotonicRegression:
    """Fit isotonic regression on held-out posteriors."""
    iso = IsotonicRegression(y_min=0, y_max=1, out_of_bounds='clip')
    iso.fit(posteriors, labels)
    return iso


# ── Channel definitions ─────────────────────────────────────────────────────

def logit_upper_half(x):
    """Transform [0.5, 1] → bounded logit. Values <0.5 clamped.

    Capped at ±4 to prevent extreme values for entropy_sim near 0.5 or 1.0.
    Without capping, sim=1.0 maps to logit=13 which falls far outside the
    fitted distributions and produces unreliable logLRs.
    """
    x = np.asarray(x, dtype=float)
    x = np.clip(x, 0.5, 1.0 - 1e-6)
    t = (x - 0.5) / 0.5
    t = np.clip(t, 1e-6, 1.0 - 1e-6)
    result = np.log(t / (1 - t))
    return np.clip(result, -4.0, 4.0)


# Default channel specs for hit-level candidate ranking
DEFAULT_CHANNELS = [
    # Primary spectral match
    ChannelSpec('entropy_sim', 'entropy_similarity', 'continuous',
                higher_means_tp=True, tp_family='normal', fp_family='normal',
                transform=logit_upper_half),

    # delta_mda: DROPPED from ranking model. All candidates were found by precursor
    # mass search (~10 ppm window), so delta_mda is small for everyone. It adds noise
    # rather than signal for within-spectrum ranking. Useful for annotation evaluation
    # (Problem B) but not candidate ranking (Problem A).

    # RT prediction error (signed)
    ChannelSpec('signed_delta_rt', 'signed_delta_rt', 'continuous',
                higher_means_tp=True, tp_family='student_t', fp_family='student_t'),

    # Similarity gap (this candidate vs next best)
    ChannelSpec('sim_gap', 'sim_gap', 'continuous',
                higher_means_tp=True, tp_family='normal', fp_family='normal'),

    # Spectral entropy (query spectrum quality — shared across candidates)
    ChannelSpec('spectral_entropy', 'spectral_entropy', 'continuous',
                higher_means_tp=True, tp_family='normal', fp_family='normal'),

    # ISF adduct flag
    ChannelSpec('hit_is_isf', 'hit_is_isf', 'binary',
                higher_means_tp=False),

    # ISF with no ok adduct for this compound in this spectrum
    ChannelSpec('hit_isf_no_ok', 'hit_isf_no_ok', 'binary',
                higher_means_tp=False),

    # Compound has at least one ok adduct in this spectrum
    ChannelSpec('compound_has_ok_adduct', 'compound_has_ok_adduct', 'binary',
                higher_means_tp=True),

    # Number of distinct adducts for this compound in this spectrum
    ChannelSpec('n_candidate_adducts', 'n_candidate_adducts', 'discrete',
                higher_means_tp=True),

    # Polarity
    ChannelSpec('polarity', 'polarity', 'binary',
                higher_means_tp=False),
]


# ── Convenience: full pipeline ──────────────────────────────────────────────

def run_pipeline(feature_table: pd.DataFrame,
                 channel_specs: List[ChannelSpec] = None,
                 prior_tp: float = None,
                 label_col: str = 'hit_label') -> pd.DataFrame:
    """Fit channels, score all rows, compute posteriors."""
    if channel_specs is None:
        channel_specs = DEFAULT_CHANNELS

    if prior_tp is None:
        prior_tp = feature_table[label_col].mean()

    print(f'Fitting channels (prior_tp={prior_tp:.4f})...')
    fitted = fit_all_channels(feature_table, channel_specs, label_col)

    print('Computing logLR...')
    lr_table = compute_logLR_table(feature_table, fitted)

    print('Computing posteriors...')
    result = compute_posterior(lr_table, prior_tp)

    return result, fitted
