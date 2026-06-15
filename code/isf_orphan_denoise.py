"""isf_orphan_denoise.py — reverse/containment filter to denoise bin-unmatched MS/MS.

PURPOSE
  In a regular LCB study, hundreds of MS/MS do not match any Bin. Most are NOT novel
  compounds — they are in-source fragments / adducts / isotopes of compounds that ARE
  already binned, plus some genuine noise. This tool decides, for each such orphan MS/MS,
  whether it is explained as a relational artifact of a CO-ELUTING reference bin
  (-> suppress / link to parent) or has no co-eluting explanation (-> genuine candidate).

WHY REVERSE / CONTAINMENT (not forward cosine)
  An ISF/adduct orphan O is a SUBSET of its richer parent C: C fragments more, so it
  carries O's peaks plus many others. Forward cosine is penalized by C's extra peaks and
  misses the relationship. The one-sided containment(O in C) — "are O's peaks present in C"
  — stays high. This is the NIST/ISFrag reverse-match logic pointed at the orphan->parent pair.

DECISION RULE (per orphan O, against co-eluting reference bins C)
  flag O as relational noise if EXISTS C with:
    (1) co-elution: |rt(O) - rt(C)| <= RT_WIN
    (2) mass relation: prec(C)-prec(O) is a neutral loss (ISF), or |prec(C)-prec(O)| is an
        adduct-pair delta or an isotope spacing (related ion)
    (3) containment(O in C) >= CONT_MIN   [reverse score]
  else O is an unexplained candidate (novel / missed annotation -> worklist).

USAGE
  # Demonstration on curated Orbitrap bins (proxy for the LCB orphan problem):
  python code/isf_orphan_denoise.py --demo
  # Real LCB run (when BinBase/carrot-prod access is available):
  python code/isf_orphan_denoise.py \
      --orphans data/<study>_unmatched_msms.csv \
      --bins    data/<study>_bins.csv \
      --peaks   data/<study>_peaks_cache.json \
      --out     data/<study>_orphan_denoise.csv
  CSV cols required: wiki_id, rt (seconds), precursor_mz. Peaks: JSON keyed by wiki_id ->
  list of [mz, intensity]; or pass --peaks-in-csv with a 'peaks' column of "mz:int;mz:int".

STATUS: kernel validated on binned data (anchor containment 0.976; random-Δm/z control -> 0%).
  Demo discrimination below is the "does it work" evidence on our data. The real LCB orphan
  population is gated on BinBase/carrot-prod access; this script is ready to fire on it.
"""
import os, sys, json, argparse
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TOL = 0.01            # MS2 peak match tolerance (Da)
RT_WIN = 4.0          # co-elution window (seconds)
MTOL = 0.006          # precursor Δm/z tolerance (Da)
CONT_MIN = 0.5        # containment threshold (reverse score)

# Mass relations that mark an orphan as a related ion of a co-eluting bin.
NEUTRAL_LOSSES = {    # parent = orphan + loss   (orphan is the in-source fragment)
    'H2O': 18.0106, '2H2O': 36.0211, 'NH3': 17.0265, 'CO': 27.9949,
    'CO2': 43.9898, 'HCOOH': 46.0055, 'CH2O': 30.0106, 'hexose': 162.0528,
}
ADDUCT_DELTAS = {     # |prec(C) - prec(O)| for same M, different adduct
    'Na-H': 21.9819, 'K-H': 37.9559, 'NH4-H': 17.0265,
}
ISOTOPE = {'13C': 1.00336, '2x13C': 2.00671}


def containment(qF, qP, tol=TOL):
    """Fraction of orphan F's intensity whose peaks are present in reference P (reverse score)."""
    F = np.asarray(qF, float); P = np.asarray(qP, float)
    if F.size == 0 or P.size == 0:
        return np.nan
    pmz = np.sort(P[:, 0]); tot = F[:, 1].sum(); m = 0.0
    for mz, it in F:
        j = np.searchsorted(pmz, mz)
        if (j < len(pmz) and abs(pmz[j] - mz) <= tol) or (j > 0 and abs(pmz[j - 1] - mz) <= tol):
            m += it
    return m / max(1e-12, tot)


CONT_MIN_RELATED = 0.7   # tighter bar for adduct/isotope (no precursor-in-parent signature)


def _relation(d):
    """Given prec(C)-prec(O)=d, return relation_name or None. d>0 means C heavier than O."""
    for nm, L in NEUTRAL_LOSSES.items():
        if abs(d - L) <= MTOL:           # parent heavier by a neutral loss -> O is the in-source fragment
            return 'ISF:' + nm
    for nm, L in ADDUCT_DELTAS.items():
        if abs(abs(d) - L) <= MTOL:
            return 'adduct:' + nm
    for nm, L in ISOTOPE.items():
        if abs(abs(d) - L) <= MTOL:
            return 'isotope:' + nm
    return None


def _peak_present(parent_pk, mz, tol=TOL):
    """Is m/z present as a peak in the parent spectrum? (the strong ISF signature:
    the fragment ion itself appears among the parent's fragments)."""
    P = np.asarray(parent_pk, float)
    if P.size == 0:
        return False
    pmz = np.sort(P[:, 0]); j = np.searchsorted(pmz, mz)
    return (j < len(pmz) and abs(pmz[j] - mz) <= tol) or (j > 0 and abs(pmz[j - 1] - mz) <= tol)


def denoise(orphans, refs, peaks, rt_win=RT_WIN, cont_min=CONT_MIN, require_prec_in_parent=True):
    """orphans, refs: DataFrames with wiki_id, rt, precursor_mz (refs may include orphans).
    peaks: dict wiki_id -> [[mz,int],...]. Returns per-orphan result DataFrame.

    A relation is accepted only if reverse/containment passes AND, for ISF relations, the
    orphan's precursor ion is present as a peak in the parent (the validated strong signature
    that cuts false-suppression of real primaries ~15-22% -> ~6%). Adduct/isotope relations
    have no such signature, so they require a tighter containment (CONT_MIN_RELATED)."""
    refs = refs.dropna(subset=['rt', 'precursor_mz']).copy()
    rt = refs['rt'].values.astype(float); prec = refs['precursor_mz'].values.astype(float)
    wid = refs['wiki_id'].values
    order = np.argsort(rt); rts = rt[order]
    rows = []
    for _, o in orphans.iterrows():
        ow = o['wiki_id']; ort = o['rt']; opr = o['precursor_mz']; opk = peaks.get(ow)
        best = (False, None, None, np.nan)
        if opk is not None and np.isfinite(ort) and np.isfinite(opr):
            lo = np.searchsorted(rts, ort - rt_win); hi = np.searchsorted(rts, ort + rt_win)
            for k in range(lo, hi):
                j = order[k]
                if wid[j] == ow:
                    continue
                rel = _relation(prec[j] - opr)
                if rel is None:
                    continue
                ppk = peaks.get(wid[j])
                cc = containment(opk, ppk)
                if cc is None or np.isnan(cc):
                    continue
                if rel.startswith('ISF'):
                    ok = cc >= cont_min and (not require_prec_in_parent or _peak_present(ppk, opr))
                else:                                   # adduct / isotope: tighter bar
                    ok = cc >= CONT_MIN_RELATED
                if ok and (best[3] != best[3] or cc > best[3]):
                    best = (True, wid[j], rel, cc)
        rows.append((ow, *best))
    return pd.DataFrame(rows, columns=['wiki_id', 'explained', 'parent', 'relation', 'containment'])


# ---------------------------------------------------------------------------
def _adduct_class(a):
    s = '' if not isinstance(a, str) else a
    n = s.strip().strip('[]').rstrip('+-12 ')
    if any(t in n for t in ['-H2O', '+H-', '-H-', '-CH', '-CO', '-NH', '-C']):
        return 'isf'
    if n in ('M+H', 'M-H', 'M+Na', 'M+K', 'M+NH4', 'M+Cl', 'M+FA-H', 'M+HCOO', 'M-H2O') or n.startswith('M+H') or n.startswith('M-H'):
        return 'primary'
    return 'other'


def run_demo():
    """Proxy demo on curated Orbitrap bins. Each bin is treated as an incoming MS/MS and
    tested against the other co-eluting bins. We then check the filter's discrimination:
    it should FLAG in-source-fragment bins (artifacts to suppress) and SPARE clean primary
    compounds (real, must keep). Adduct string = independent eval label (not used by filter)."""
    cols = ['wiki_id', 'name', 'adduct', 'rt', 'precursor_mz']
    pos = pd.read_csv(os.path.join(ROOT, 'data', 'Orbitrap_HILIC_posESI_curated_042126.csv'), usecols=cols)
    neg = pd.read_csv(os.path.join(ROOT, 'data', 'Orbitrap_HILIC_negESI_curated_042126.csv'), usecols=cols)
    cur = pd.concat([pos, neg], ignore_index=True)
    peaks = json.load(open(os.path.join(ROOT, 'data', 'query_peaks_cache_v2.json')))
    B = cur[cur.wiki_id.isin(peaks)].reset_index(drop=True)
    print(f'demo bins with peaks: {len(B)}')

    res = denoise(B, B, peaks)
    B = B.merge(res, on='wiki_id', how='left')
    B['adduct_class'] = B['adduct'].map(_adduct_class)
    B['is_yy'] = B['name'].fillna('').str.lower().str.startswith('yy_')

    print('\n=== DOES IT WORK? flag (=explained-by-co-eluting-parent) rate by bin type ===')
    for lab, mask in [
        ('IN-SOURCE-FRAGMENT adduct  (artifact -> SHOULD flag)', B.adduct_class == 'isf'),
        ('curator-rejected  yy_      (artifact -> SHOULD flag)', B.is_yy),
        ('clean PRIMARY adduct       (real     -> should NOT)', (B.adduct_class == 'primary') & (~B.is_yy)),
    ]:
        s = B[mask]
        if len(s):
            print('  %-52s n=%5d  flagged %5.1f%%' % (lab, len(s), 100 * s.explained.mean()))

    isf = B[(B.adduct_class == 'isf')]
    prim = B[(B.adduct_class == 'primary') & (~B.is_yy)]
    print('\n  -> suppression of ISF artifacts:        %.1f%%' % (100 * isf.explained.mean()))
    print('  -> FALSE-suppression of real primaries: %.1f%%  (want low)' % (100 * prim.explained.mean()))
    print('     discrimination gap: %.1f points' % (100 * (isf.explained.mean() - prim.explained.mean())))
    print('\n  relation types among flagged bins:')
    print('   ' + B[B.explained].relation.value_counts().head(10).to_string().replace('\n', '\n   '))

    tot = len(B); expl = int(B.explained.sum())
    print(f'\n  If this were a study\'s orphan pool: {expl}/{tot} ({100*expl/tot:.0f}%) would be explained away as')
    print('  relational noise (ISF/adduct/isotope of a co-eluting bin); the rest -> novel-candidate worklist.')
    out = os.path.join(ROOT, 'data', 'demo_orphan_denoise.csv')
    B[['wiki_id', 'name', 'adduct', 'rt', 'precursor_mz', 'explained', 'parent', 'relation', 'containment']].to_csv(out, index=False)
    print(f'\n  wrote per-orphan table: {out}')


def run_real(args):
    orphans = pd.read_csv(args.orphans)
    bins = pd.read_csv(args.bins)
    if args.peaks_in_csv:
        peaks = {}
        for df in (orphans, bins):
            for _, r in df.iterrows():
                peaks[r['wiki_id']] = [[float(x.split(':')[0]), float(x.split(':')[1])]
                                       for x in str(r['peaks']).split(';') if ':' in x]
    else:
        peaks = json.load(open(args.peaks))
    # reference set = the study's bins; orphans = the unmatched MS/MS
    res = denoise(orphans, bins, peaks)
    res = res.merge(orphans[['wiki_id', 'rt', 'precursor_mz']], on='wiki_id', how='left')
    n = len(res); e = int(res.explained.sum())
    print(f'orphans: {n} | explained as relational noise: {e} ({100*e/max(1,n):.0f}%) | candidates: {n-e}')
    print(res.relation.value_counts().head(10).to_string())
    res.to_csv(args.out, index=False)
    print('wrote', args.out)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--demo', action='store_true')
    ap.add_argument('--orphans'); ap.add_argument('--bins')
    ap.add_argument('--peaks'); ap.add_argument('--peaks-in-csv', action='store_true')
    ap.add_argument('--out', default=os.path.join(ROOT, 'data', 'orphan_denoise.csv'))
    a = ap.parse_args()
    if a.demo or not a.orphans:
        run_demo()
    else:
        run_real(a)
