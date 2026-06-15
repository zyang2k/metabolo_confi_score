"""simgap_significance.py — ion-counting uncertainty on the candidate similarity gap.

MOTIVATION (Oliver 2026-06-09, methylcytidine bin).
The production `sim_gap` (build_features_v2.py) is a bare point difference in API
entropy_similarity between a spectrum's best and 2nd-best *identifiable* candidate. The
GBM treated a 0.006 gap as decisive and ranked the wrong isomer first, while RT favoured
the other. But a 0.006 gap means something very different on a dense, well-counted MS²
than on a near-single-ion one — a flat floor (e.g. "ignore gaps < 0.05") is too brutal
because it ignores the spectrum.

PRINCIPLE.
Each fragment-ion intensity is proportional to an ion count, so it carries multinomial
counting statistics. Re-draw the query (and library) spectra within that counting noise,
recompute every candidate's entropy similarity, and watch how much the top-vs-runner-up
GAP moves. The gap is "significant" only when it is large relative to its own
counting-induced spread:

        sim_gap_z = gap_mean / sigma_gap          (≈ a per-spectrum z-score)

On the methylcytidine bin this is the whole story: the gap is reproducible in *direction*
but only ~1-2σ of its own uncertainty → MS² cannot separate the isomers here → call them
tied and let RT decide.

LOCAL ENTROPY SIMILARITY.
entropy_similarity in the feature table comes from the MassWiki API as a scalar — it can't
be bootstrapped. We recompute it locally with ms_entropy. At ms2_tolerance=0.02 Da the
local value tracks the API value at corr≈0.97, median |Δ|≈0.012 (validated in
bench_simgap_significance.py). The z-score is scale-free, so the local-vs-API offset does
not matter — what matters is that the *spread* of the gap under resampling is faithful.

ONE FREE PARAMETER: N_eff, the effective total ion count (counting budget). Absolute
calibration of N_eff needs raw ion counts / TIC, which the MassWiki query peaks (relative-
abundance normalised) don't carry — that is the Quentin/Fanzhou ask. Until then we use a
nominal budget and report sensitivity. The PER-SPECTRUM self-adaptation still comes for
free from each spectrum's own peak structure under multinomial resampling: a gap propped
up by one tiny peak is unstable; a gap supported by many concordant peaks is stable. That
is exactly the "dense vs single-ion" distinction.

LABEL-FREE: uses only observed/library peaks, never hit_label / spectrum_label.
"""
from __future__ import annotations

import numpy as np
import ms_entropy as me

# tol that best matches the MassWiki API entropy_similarity (see module docstring).
ENTROPY_TOL_DA = 0.02
# nominal effective ion count (counting budget) — the one free parameter; see docstring.
DEFAULT_N_EFF = 1000
# bootstrap draws per spectrum.
DEFAULT_B = 200
# only the top-K candidates by similarity can plausibly win the top-1-vs-2 gap under
# resampling; deeper candidates never reorder into the top two. Caps compute.
TOPK_COMPETITORS = 8


def local_entropy_sim(q_peaks, l_peaks, tol: float = ENTROPY_TOL_DA) -> float:
    """ms_entropy similarity between two [[mz, intensity], ...] peak lists."""
    if q_peaks is None or l_peaks is None or len(q_peaks) == 0 or len(l_peaks) == 0:
        return np.nan
    q = np.asarray(q_peaks, dtype=np.float64)
    l = np.asarray(l_peaks, dtype=np.float64)
    return float(me.calculate_entropy_similarity(q, l, ms2_tolerance_in_da=tol,
                                                 clean_spectra=True))


def multinomial_resample(peaks, n_eff: int, rng: np.random.Generator):
    """Re-draw a spectrum under ion-counting noise.

    Intensities are treated as relative abundances p_i; we draw n_eff ions ~
    Multinomial(n_eff, p) and return [[mz, count/n_eff], ...] keeping only peaks that drew
    ≥1 ion. m/z values are unchanged (mass axis is not the counting-noisy quantity here).
    """
    arr = np.asarray(peaks, dtype=np.float64)
    mz = arr[:, 0]
    inten = arr[:, 1]
    tot = inten.sum()
    if tot <= 0 or len(mz) == 0:
        return arr
    p = inten / tot
    counts = rng.multinomial(n_eff, p)
    keep = counts > 0
    if not keep.any():
        return arr
    out = np.column_stack([mz[keep], counts[keep].astype(np.float64) / n_eff])
    return out


def gap_significance_for_spectrum(
    query_peaks,
    competitors: list[dict],
    n_eff: int = DEFAULT_N_EFF,
    B: int = DEFAULT_B,
    tol: float = ENTROPY_TOL_DA,
    resample_library: bool = True,
    topk: int = TOPK_COMPETITORS,
    seed: int = 0,
) -> dict:
    """Ion-counting significance of the top-1-vs-2 similarity gap for ONE spectrum.

    `competitors` is a list of dicts, one per *identifiable* candidate (hit_ik14 != ''),
    each with keys: 'library_wiki_id' (str) and 'lib_peaks' ([[mz,int],...] or None).
    Candidates whose library peaks are missing are dropped (cannot be bootstrapped).

    Returns a dict with the point gap, its bootstrap mean/σ, the z-score, the fraction of
    draws in which the point top-1 stays on top (direction stability), and bookkeeping.
    All-NaN if fewer than 2 bootstrappable competitors or no query peaks.
    """
    nan_out = dict(n_competitors=0, gap_point=np.nan, gap_mean=np.nan,
                   gap_sigma=np.nan, sim_gap_z=np.nan, p_top1_stable=np.nan,
                   top1_library_wiki_id=None, sim_top1_local=np.nan, sim_top2_local=np.nan)
    if query_peaks is None or len(query_peaks) == 0:
        return nan_out

    comps = [c for c in competitors if c.get('lib_peaks') is not None
             and len(c['lib_peaks']) > 0]
    if len(comps) < 2:
        return nan_out

    # Point estimate: local entropy sim of each competitor against the *observed* query.
    point = np.array([local_entropy_sim(query_peaks, c['lib_peaks'], tol) for c in comps])
    order = np.argsort(point)[::-1]            # descending
    comps = [comps[i] for i in order]
    point = point[order]
    # Cap to the top-K candidates that could realistically reorder into the top two.
    if len(comps) > topk:
        comps = comps[:topk]
        point = point[:topk]

    gap_point = float(point[0] - point[1])
    top1_lib = comps[0]['library_wiki_id']

    rng = np.random.default_rng(seed)
    lib_arrays = [np.asarray(c['lib_peaks'], dtype=np.float64) for c in comps]
    gaps = np.empty(B, dtype=np.float64)
    top1_stable = 0
    for b in range(B):
        q_b = multinomial_resample(query_peaks, n_eff, rng)
        sims = np.empty(len(comps), dtype=np.float64)
        for j, lib in enumerate(lib_arrays):
            l_b = multinomial_resample(lib, n_eff, rng) if resample_library else lib
            sims[j] = local_entropy_sim(q_b, l_b, tol)
        s_sorted = np.sort(sims)[::-1]
        gaps[b] = s_sorted[0] - s_sorted[1]
        if np.argmax(sims) == 0:               # index 0 == point top-1 (comps reordered)
            top1_stable += 1

    gap_mean = float(np.mean(gaps))
    gap_sigma = float(np.std(gaps, ddof=1))
    # z-score: how many counting-σ the gap sits above zero. σ floored to avoid /0 on
    # degenerate (identical-peak) competitors — those are maximally "tied", z→large only
    # if the gap itself is real, so a tiny floor keeps the ratio honest.
    z = gap_mean / max(gap_sigma, 1e-6)
    return dict(
        n_competitors=len(comps),
        gap_point=gap_point,
        gap_mean=gap_mean,
        gap_sigma=gap_sigma,
        sim_gap_z=float(z),
        p_top1_stable=top1_stable / B,
        top1_library_wiki_id=top1_lib,
        sim_top1_local=float(point[0]),
        sim_top2_local=float(point[1]),
    )


def build_significance_table(
    ft,
    query_cache: dict,
    lib_cache: dict,
    wiki_ids=None,
    n_eff: int = DEFAULT_N_EFF,
    B: int = DEFAULT_B,
    tol: float = ENTROPY_TOL_DA,
    resample_library: bool = True,
    progress_every: int = 500,
):
    """Compute per-spectrum gap significance for every wiki_id (or a given subset).

    Competitors per spectrum = candidate rows with a non-empty hit_ik14 (mirrors the
    sim_gap competitor pool in build_features_v2.py). Returns a DataFrame keyed by
    wiki_id with the significance columns from `gap_significance_for_spectrum`.
    """
    import pandas as pd

    if wiki_ids is None:
        wiki_ids = ft['wiki_id'].unique()
    wid_set = set(wiki_ids)
    sub = ft[ft['wiki_id'].isin(wid_set) & ft['hit_ik14'].fillna('').ne('')]
    # Stable per-wiki_id seed so the bench is reproducible regardless of row order.
    grouped = sub.groupby('wiki_id')
    rows = []
    for i, (wid, g) in enumerate(grouped):
        if progress_every and i % progress_every == 0:
            print(f'  significance: {i}/{len(grouped)} spectra', flush=True)
        q = query_cache.get(wid)
        comps = [
            {'library_wiki_id': lwid, 'lib_peaks': lib_cache.get(str(lwid))}
            for lwid in g['library_wiki_id'].tolist()
        ]
        seed = abs(hash(wid)) % (2 ** 32)
        res = gap_significance_for_spectrum(q, comps, n_eff=n_eff, B=B, tol=tol,
                                            resample_library=resample_library, seed=seed)
        res['wiki_id'] = wid
        rows.append(res)
    return pd.DataFrame(rows)
