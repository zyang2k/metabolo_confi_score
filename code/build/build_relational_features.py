"""build_relational_features.py — relational spectral-network features (GNPS-style).

Turns the molecular-network edge (project_knowledge_graph_pivot, edge 3) into per-bin
*relational* features for the GBM. NOT the bare nearest-neighbour cosine — that is
already captured by entropy_similarity (project_learned_similarity_mvp: Δ +0.0002,
redundant). The new signal is the structure of a bin's spectral neighbourhood.

For each labeled top-1 bin we search its MS2 against the OTHER labeled bins (the curated
reference graph), blocked by polarity + observed precursor m/z (±10 mDa, lab convention).

LEAKAGE GUARD (critical): neighbours sharing this bin's anno_ik14 are EXCLUDED. The GBM
is evaluated with GroupKFold on anno_ik14, so a same-IK14 neighbour sits in the same
held-out fold; using its label would leak the bin's own label via its replicates
(a TP's same-compound siblings are mostly TP). Excluding same anno_ik14 also subsumes
same-batch replicate blocking. Empty-IK14 bins fall back to excluding same anno_name_lower.
Net effect: the features describe the *cross-compound confusability* neighbourhood.

Features (per wiki_id):
  rel_nn_sim         max entropy sim to any different-compound labeled bin (confusability)
  rel_n_confirmed_nbr  # different-compound TP bins with sim >= SIM_THRESH (dense confirmed nbhd)
  rel_nbr_purity     mean(hit_label) over neighbours with sim >= SIM_THRESH (NaN if none)
  rel_claim_corrob   1 if nearest such neighbour's curator identity (anno_ik14) == this
                     bin's CLAIMED hit_ik14 (the spectral nbhd points at the claim)

Output: data/relational_features.csv  (merge on wiki_id)
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import ms_entropy as me

ROOT = Path(__file__).resolve().parent.parent
FT = ROOT / 'data' / 'feature_table_v2.csv'
QP = ROOT / 'data' / 'query_peaks_cache_v2.json'
OUT = ROOT / 'data' / 'relational_features.csv'

MS1_TOL = 0.01      # Da — precursor block half-width (10 mDa, lab convention)
MS2_TOL = 0.02      # Da — entropy fragment tolerance (in-repo default)
SIM_THRESH = 0.70   # entropy similarity to count a spectral edge


def build_top1(ft):
    lab = ft[ft['spectrum_label'].isin(['TP', 'FP'])]
    idx = lab.groupby('wiki_id')['entropy_similarity'].idxmax()
    return ft.loc[idx].reset_index(drop=True)


def main():
    ft = pd.read_csv(FT, low_memory=False)
    t = build_top1(ft)
    t['precursor_mz'] = pd.to_numeric(t['precursor_mz'], errors='coerce')
    t['anno_ik14'] = t['anno_ik14'].fillna('').astype(str)
    t['hit_ik14'] = t['hit_ik14'].fillna('').astype(str)
    t['anno_name_lower'] = t['anno_name_lower'].fillna('').astype(str)
    print(f'labeled top-1 bins: {len(t):,}')

    peaks = {k: np.asarray(v, dtype=np.float64)
             for k, v in json.load(open(QP)).items() if v}
    t = t[t['wiki_id'].isin(peaks) & t['precursor_mz'].notna()].reset_index(drop=True)
    print(f'  with peaks + precursor: {len(t):,}')

    # block key for replicate/leakage exclusion: anno_ik14, else anno_name_lower
    block = np.where(t['anno_ik14'] != '', t['anno_ik14'], 'name::' + t['anno_name_lower'])
    t['block'] = block

    rows = []
    for pol, sub in t.groupby('polarity'):
        sub = sub.reset_index(drop=True)
        mz = sub['precursor_mz'].values
        order = np.argsort(mz)
        mz_sorted = mz[order]
        wid = sub['wiki_id'].values
        blk = sub['block'].values
        lab = sub['hit_label'].values.astype(float)
        nbr_ik = sub['anno_ik14'].values
        for i in range(len(sub)):
            lo = np.searchsorted(mz_sorted, mz[i] - MS1_TOL, 'left')
            hi = np.searchsorted(mz_sorted, mz[i] + MS1_TOL, 'right')
            cand = [order[p] for p in range(lo, hi)
                    if order[p] != i and blk[order[p]] != blk[i]]
            qp = peaks[wid[i]]
            sims = []
            for j in cand:
                s = me.calculate_entropy_similarity(
                    qp, peaks[wid[j]], ms2_tolerance_in_da=MS2_TOL, clean_spectra=True)
                sims.append((s, j))
            claim = sub['hit_ik14'].values[i]
            if not sims:
                rows.append((wid[i], 0.0, 0, np.nan, 0))
                continue
            sims.sort(reverse=True)
            nn_sim, nn_j = sims[0]
            edge = [(s, j) for s, j in sims if s >= SIM_THRESH]
            n_conf = sum(1 for s, j in edge if lab[j] == 1)
            purity = np.mean([lab[j] for s, j in edge]) if edge else np.nan
            corrob = int(claim != '' and nbr_ik[nn_j] == claim and nn_sim >= SIM_THRESH)
            rows.append((wid[i], float(nn_sim), int(n_conf), purity, corrob))

    out = pd.DataFrame(rows, columns=[
        'wiki_id', 'rel_nn_sim', 'rel_n_confirmed_nbr', 'rel_nbr_purity', 'rel_claim_corrob'])
    out.to_csv(OUT, index=False)
    print(f'\nWrote {OUT} ({len(out):,} rows)')
    print(f'  rel_nn_sim       mean {out.rel_nn_sim.mean():.3f}  '
          f'frac>={SIM_THRESH}: {(out.rel_nn_sim>=SIM_THRESH).mean():.3f}')
    print(f'  rel_n_confirmed_nbr  mean {out.rel_n_confirmed_nbr.mean():.2f}  '
          f'max {out.rel_n_confirmed_nbr.max()}')
    print(f'  rel_nbr_purity   non-NaN {out.rel_nbr_purity.notna().mean():.3f}  '
          f'mean {out.rel_nbr_purity.mean():.3f}')
    print(f'  rel_claim_corrob frac=1: {out.rel_claim_corrob.mean():.3f}')


if __name__ == '__main__':
    main()
