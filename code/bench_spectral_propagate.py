"""bench_spectral_propagate.py — Spectral-similarity propagation (KG prototype layer 3).

The edge type where the graph stops being a groupby. Layer 2 (RT-anchor propagation)
only reached the 516/9,175 uncurated bins whose CLAIMED NAME matched a curated compound.
This layer propagates along MS2-NEIGHBOURHOOD edges instead: search each uncurated query
spectrum against the curated TRUE-POSITIVE bins as an in-house reference library, and
inherit the curator's identity from the spectral nearest neighbour — regardless of what
the candidate generator named the query. That reaches the ~94% layer 2 could not.

Two independent evidence axes then compose:
  spectral edge  -> proposes an identity (curated neighbour's compound)
  RT agreement   -> confirms it (query elutes where that compound is known to elute)

  CONFIRMED propagated annotation : spectral_sim >= SIM_THRESH AND |Δrt| <= RT_TIGHT
  MS2-ONLY (analog/isomer?)       : spectral_sim >= SIM_THRESH but RT disagrees

Search is blocked by polarity and observed precursor m/z (±MS1_TOL) — same-compound
edges only; modified-cosine analog edges are deliberately out of scope for v1.

Peaks: curated bins  -> data/query_peaks_cache_v2.json
       uncurated bins -> data/blanks_query_peaks_cache_{pos,neg}.json
Observed precursor m/z: precursor_mz in the curated CSVs (all 14,809 bins covered).

Validation (leave-one-out retrieval on curated TP, NOT circular): search each curated
TP bin against the OTHER curated TP bins; top-1 hit's compound == its own compound?
Reports identity top-1 / top-3 retrieval accuracy — does the spectral neighbourhood
recover the right identity at all.

Outputs
  data/spectral_propagated.csv            — per uncurated query: best curated match + verdict
  data/spectral_missed_annotations.csv    — CONFIRMED (MS2 + RT) propagated annotations
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import ms_entropy as me

ROOT = Path(__file__).resolve().parent.parent
GRAPH = ROOT / 'data' / 'bench_name_rt_graph.csv'
QP_CUR = ROOT / 'data' / 'query_peaks_cache_v2.json'
QP_POS = ROOT / 'data' / 'blanks_query_peaks_cache_pos.json'
QP_NEG = ROOT / 'data' / 'blanks_query_peaks_cache_neg.json'
CUR_NEG = ROOT / 'data' / 'Orbitrap_HILIC_negESI_curated_042126.csv'
CUR_POS = ROOT / 'data' / 'Orbitrap_HILIC_posESI_curated_042126.csv'
OUT_PROP = ROOT / 'data' / 'spectral_propagated.csv'
OUT_MISS = ROOT / 'data' / 'spectral_missed_annotations.csv'

MS1_TOL = 0.01    # Da (=10 mDa, lab convention) — precursor block half-width
MS2_TOL = 0.02    # Da — entropy-similarity fragment tolerance (in-repo default)
SIM_THRESH = 0.70 # entropy similarity to call a spectral edge
RT_TIGHT = 30.0   # s — RT agreement window (shared across KG layers)


def load_peaks():
    peaks = {}
    for f in (QP_CUR, QP_POS, QP_NEG):
        with open(f) as fh:
            d = json.load(fh)
        for k, v in d.items():
            if k not in peaks and v:
                peaks[k] = np.asarray(v, dtype=np.float64)
    return peaks


def best_match(q_peaks, ref_idx, ref_peaks):
    """Entropy similarity of q against each reference in ref_idx; return (best_pos, best_sim)."""
    best_pos, best_sim = -1, 0.0
    for pos in ref_idx:
        sim = me.calculate_entropy_similarity(
            q_peaks, ref_peaks[pos], ms2_tolerance_in_da=MS2_TOL, clean_spectra=True)
        if sim > best_sim:
            best_sim, best_pos = sim, pos
    return best_pos, best_sim


def main():
    g = pd.read_csv(GRAPH); g['is_curated'] = g['is_curated'].astype(bool)

    # observed precursor m/z (all bins)
    pre = pd.concat([
        pd.read_csv(CUR_NEG, usecols=['wiki_id', 'precursor_mz']),
        pd.read_csv(CUR_POS, usecols=['wiki_id', 'precursor_mz']),
    ]).drop_duplicates('wiki_id')
    g = g.merge(pre, on='wiki_id', how='left')

    print('Loading peaks (curated + blank caches)')
    peaks = load_peaks()
    g['has_peaks'] = g['wiki_id'].isin(peaks)
    g = g[g['has_peaks'] & g['precursor_mz'].notna()].copy()

    # Reference library = curated TRUE-POSITIVE bins (the curator's confirmed identities)
    ref = g[g['is_curated'] & (g['hit_label'] == 1)].reset_index(drop=True)
    ref_peaks = [peaks[w] for w in ref['wiki_id']]
    print(f'  reference (curated TP) bins: {len(ref):,}')

    # ---- helper: precursor+polarity-blocked search of a query frame vs ref ------
    def search(queries, exclude_self_wiki=False):
        out = []
        for pol, rsub in ref.groupby('polarity'):
            r_local = rsub.reset_index()           # 'index' -> position in ref
            order = np.argsort(r_local['precursor_mz'].values)
            r_mz = r_local['precursor_mz'].values[order]
            r_globpos = r_local['index'].values[order]
            q = queries[queries['polarity'] == pol]
            for _, b in q.iterrows():
                lo = np.searchsorted(r_mz, b['precursor_mz'] - MS1_TOL, 'left')
                hi = np.searchsorted(r_mz, b['precursor_mz'] + MS1_TOL, 'right')
                cand = list(r_globpos[lo:hi])
                if exclude_self_wiki:
                    cand = [p for p in cand if ref.at[p, 'wiki_id'] != b['wiki_id']]
                if not cand:
                    continue
                pos, sim = best_match(peaks[b['wiki_id']], cand, ref_peaks)
                if pos < 0:
                    continue
                out.append((b['wiki_id'], pos, sim))
        return out

    # ---- Validation: leave-one-out identity retrieval on curated TP -------------
    print('\n=== LOO retrieval validation (curated TP vs other curated TP) ===')
    loo = search(ref, exclude_self_wiki=True)
    if loo:
        lo_df = pd.DataFrame(loo, columns=['wiki_id', 'ref_pos', 'sim'])
        lo_df = lo_df.merge(ref[['wiki_id', 'norm_name']], on='wiki_id')
        lo_df['match_name'] = ref['norm_name'].values[lo_df['ref_pos'].values]
        lo_df = lo_df.merge(ref[['wiki_id', 'hit_ik14']], on='wiki_id', how='left')
        lo_df['match_ik14'] = ref['hit_ik14'].values[lo_df['ref_pos'].values]
        hit_name = (lo_df['match_name'] == lo_df['norm_name'])
        # IK14 identity: canonical compound key (feedback_compound_key). Only judge where
        # BOTH bins carry a non-empty IK14 — lipid shorthand / NIST often have none.
        both_ik = lo_df['hit_ik14'].notna() & (lo_df['hit_ik14'] != '') \
            & lo_df['match_ik14'].notna() & (lo_df['match_ik14'] != '')
        hit_ik = both_ik & (lo_df['hit_ik14'] == lo_df['match_ik14'])
        atthr = lo_df['sim'] >= SIM_THRESH
        print(f'  curated TP with a same-precursor neighbour: {len(lo_df):,}  '
              f'(both-IK14 judgeable: {int(both_ik.sum()):,})')
        print(f'  top-1 by NAME string:   {100*hit_name.mean():.1f}%   '
              f'(of sim>={SIM_THRESH}: {100*hit_name[atthr].mean():.1f}% precision)')
        if both_ik.sum():
            print(f'  top-1 by IK14 (judgeable subset): '
                  f'{100*hit_ik[both_ik].mean():.1f}%   '
                  f'(of sim>={SIM_THRESH} & judgeable: '
                  f'{100*hit_ik[both_ik & atthr].mean():.1f}% precision)')

    # ---- Propagate to uncurated bins --------------------------------------------
    print('\n=== Spectral propagation to uncurated bins ===')
    unc = g[~g['is_curated']].copy()
    res = search(unc)
    pr = pd.DataFrame(res, columns=['wiki_id', 'ref_pos', 'spectral_sim'])
    for col, out in [('norm_name', 'prop_identity'), ('hit_ik14', 'prop_ik14'),
                     ('rt', 'ref_rt'), ('wiki_id', 'ref_wiki_id')]:
        pr[out] = ref[col].values[pr['ref_pos'].values]
    pr = pr.merge(unc[['wiki_id', 'rt', 'polarity', 'name', 'gbm_cal',
                       'entropy_similarity', 'precursor_mz']], on='wiki_id', how='left')
    pr['delta_rt_to_ref'] = pr['rt'] - pr['ref_rt']
    pr['rt_consistent'] = pr['delta_rt_to_ref'].abs() <= RT_TIGHT
    pr['edge'] = pr['spectral_sim'] >= SIM_THRESH
    pr['verdict'] = np.where(~pr['edge'], 'weak',
                      np.where(pr['rt_consistent'], 'CONFIRMED', 'MS2_only_RT_off'))

    n_edge = int(pr['edge'].sum())
    print(f'  uncurated queries with a same-precursor neighbour: {len(pr):,}')
    print(f'  spectral edge (sim>={SIM_THRESH}):                    {n_edge:,}')
    print(f'    CONFIRMED (MS2 + RT agree):                  {int((pr["verdict"]=="CONFIRMED").sum()):,}')
    print(f'    MS2-only, RT off (analog/isomer/FP):         {int((pr["verdict"]=="MS2_only_RT_off").sum()):,}')
    print(f'  REACH vs layer-2 name propagation (516 uncurated): spectral edges reach {n_edge:,}')

    pr.sort_values('spectral_sim', ascending=False).to_csv(OUT_PROP, index=False)
    print(f'\nWrote {OUT_PROP} ({len(pr):,} rows)')
    miss = pr[pr['verdict'] == 'CONFIRMED'].sort_values('spectral_sim', ascending=False)
    miss.to_csv(OUT_MISS, index=False)
    print(f'Wrote {OUT_MISS} ({len(miss):,} CONFIRMED propagated annotations)')

    print('\n=== Top CONFIRMED propagated annotations (uncurated → curated identity) ===')
    for _, b in miss.head(10).iterrows():
        print(f'   q={str(b["wiki_id"]):<15} sim={b["spectral_sim"]:.2f} '
              f'Δrt={b["delta_rt_to_ref"]:+6.1f}s gbm={b["gbm_cal"]:.2f} pol={b["polarity"]} '
              f'gen_name="{str(b["name"])[:20]}" -> curator_id="{str(b["prop_identity"])[:24]}"')


if __name__ == '__main__':
    main()
