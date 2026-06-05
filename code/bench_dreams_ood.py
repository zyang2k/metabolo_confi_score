"""bench_dreams_ood.py — DreaMS embedding as an OOD / abstention GATE (NOT a GBM feature).

THESIS. The GBM ceiling (~0.92) is a labels/Q1↔Q2 problem; a better representation can't
lift ranking. But confidence has a second, label-free half: *knowing when the score itself
is untrustworthy*. entropy_similarity is pairwise (no "where does this spectrum sit"); a
DreaMS embedding is a COORDINATE, so we can measure distance-from-the-trained-manifold and
abstain when a confident GBM call lands on an unfamiliar spectrum.

This does NOT touch ranking and uses NO curator labels to build the manifold — so it dodges
the Q1/Q2 mismatch by construction. It plugs into the existing NoTA/abstention machinery.

OOD distance = mean cosine distance to the k nearest Orbitrap-labeled embeddings, computed
OUT-OF-FOLD (5-fold) so a spectrum is never its own neighbour. The density model is fit
UNSUPERVISED (TP+FP together) — it never sees the TP/FP label, so TEST 1 is not circular.

WHY BLANKS, NOT FP: prep showed the CALIBRATED model makes ~0 confident FP on the labeled
slice (1/1/0 at conf>=.7/.8/.9 of 1,139 FP) — selection bias made numerical, the labeled
slice is exactly where the model is reliable. So a TP-vs-FP gate has nothing to catch here.
The confident mistakes live in the UNLABELED population, where the only proxy-negatives are
the 108 BLANKS the pipeline scored as confident annotations (confidence_pct 74.5-97.7).

PRE-REGISTERED TESTS (decide before looking):
  TEST 1 (primary, reframed) — does OOD distance flag the should-have-abstained high-confidence
    BLANKS as out-of-manifold, while leaving genuine confident TP in-distribution? TP confidence
    is restricted to the blank confidence span so the score can't trivially separate. Report
    AUC(OOD→blank) vs AUC((1-confidence)→blank) baseline, + risk-coverage (retained blank-rate).
  TEST 2 (control, real distribution shift) — TTOF (different platform) must be measurably MORE
    OOD than Orbitrap. If not, the embedding doesn't capture platform shift (contradicts the
    06-05 probe) → debug.

KILL CRITERIA (gate is dead, park DreaMS):
  - TEST 1 AUC(OOD→blank) <= 0.60, OR
  - OOD does not beat the (1-confidence) baseline by > +0.05 (the score already abstains as well).
CAVEAT: GBM trains on hand-features; OOD is in DreaMS space — proxy for "spectrally unfamiliar
→ score less reliable", not the GBM's own input manifold. Blanks are proxy-negatives (curated
as noise), not adjudicated FP; N=108.

STAGES:
  --prep     build labeled Orbitrap (TP+FP) + TTOF eval set → MGF + meta            (3.9 env)
  --embed    dreams_embeddings(MGF) → data/dreams_ood/emb.npy                        (DreaMS env)
  --analyze  OOD distance + TEST 1/2                                                 (3.9 env)
"""
from __future__ import annotations
import json, csv, argparse
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
OUT  = ROOT/'data'/'dreams_ood'
FT   = ROOT/'data'/'feature_table_v2.csv'
DELIV = ROOT/'data'/'deliverable_scores_v2.csv'
QP_ORBI = ROOT/'data'/'query_peaks_cache_v2.json'
BLANK_HI = ROOT/'data'/'blank_high_confidence.csv'        # 108 blanks the pipeline scored as confident annotations
BLANK_PK_POS = ROOT/'data'/'blanks_query_peaks_cache_pos.json'
BLANK_PK_NEG = ROOT/'data'/'blanks_query_peaks_cache_neg.json'
TTOF_NEG_SPEC = ROOT/'data'/'TTOF_HILIC_negESI_uncurated_041326.csv'
TTOF_POS_SPEC = ROOT/'data'/'TTOF_HILIC_posESI_uncurated_041326.csv'
TTOF_NEG_PK = ROOT/'data'/'ttof_neg_query_peaks_cache.json'
TTOF_POS_PK = ROOT/'data'/'ttof_pos_query_peaks_cache.json'

MGF  = OUT/'eval.mgf'
META = OUT/'meta.csv'
EMB  = OUT/'emb.npy'
THRESH = 0.90      # "confident" band for TEST 1
K_NN   = 5         # neighbours for OOD distance
TTOF_CAP = 1500    # cap TTOF rows (control only; keep embed run modest)


def fnum(x, d=None):
    try: return float(x)
    except (TypeError, ValueError): return d


def prep():
    OUT.mkdir(parents=True, exist_ok=True)
    qp = {k: v for k, v in json.load(open(QP_ORBI)).items() if v}

    # wiki_id -> precursor_mz, polarity  (deliverable lacks precursor_mz)
    feat = {}
    for r in csv.DictReader(open(FT)):
        w = r['wiki_id']
        if w in feat: continue
        feat[w] = (fnum(r.get('precursor_mz')), r.get('polarity'))

    rows, peaks_out = [], {}
    # Orbitrap labeled eval set: TP+FP that have a confidence AND peaks AND precursor
    n_lab = {'TP': 0, 'FP': 0}
    for r in csv.DictReader(open(DELIV)):
        lab = r['spectrum_label']
        if lab not in ('TP', 'FP'): continue
        w = r['wiki_id']
        if w not in qp or w not in feat: continue
        mz, pol = feat[w]
        if mz is None: continue
        conf = fnum(r.get('confidence')); craw = fnum(r.get('confidence_raw'))
        esd = fnum(r.get('ensemble_sd'))
        sid = f'orbi|{lab}|{w}'
        rows.append({'id': sid, 'platform': 'orbi', 'label': lab, 'wiki_id': w,
                     'precursor_mz': mz, 'polarity': pol,
                     'confidence': conf, 'confidence_raw': craw, 'ensemble_sd': esd})
        peaks_out[sid] = qp[w]
        n_lab[lab] += 1
    print(f'Orbitrap eval: TP {n_lab["TP"]}  FP {n_lab["FP"]}')

    # SHOULD-HAVE-ABSTAINED set: blanks the pipeline scored as confident annotations.
    # These are the proxy-negatives the OOD gate must flag (no curator FP labels exist here).
    bpk = {**{k: v for k, v in json.load(open(BLANK_PK_POS)).items() if v},
           **{k: v for k, v in json.load(open(BLANK_PK_NEG)).items() if v}}
    n_blank = 0
    for r in csv.DictReader(open(BLANK_HI)):
        w = r['wiki_id']
        if w not in bpk: continue
        mz = fnum(r.get('obs_precursor_mz'))
        cp = fnum(r.get('confidence_pct'))
        if mz is None or cp is None: continue
        sid = f'blank|blank|{w}'
        rows.append({'id': sid, 'platform': 'blank', 'label': 'blank', 'wiki_id': w,
                     'precursor_mz': mz, 'polarity': r.get('polarity'),
                     'confidence': cp/100.0, 'confidence_raw': '', 'ensemble_sd': ''})
        peaks_out[sid] = bpk[w]
        n_blank += 1
    print(f'High-confidence blanks (should-abstain): {n_blank}')

    # TTOF control (annotated, capped) — used only as the OOD-positive control in TEST 2
    n_ttof = 0
    for spec_f, pk_f, pol in [(TTOF_NEG_SPEC, TTOF_NEG_PK, 'neg'), (TTOF_POS_SPEC, TTOF_POS_PK, 'pos')]:
        pk = {k: v for k, v in json.load(open(pk_f)).items() if v}
        for r in csv.DictReader(open(spec_f)):
            if n_ttof >= TTOF_CAP: break
            w = r.get('wiki_id')
            if not w or w not in pk: continue
            mz = fnum(r.get('precursor_mz'))
            if mz is None: continue
            sid = f'ttof|na|{w}'
            if sid in peaks_out: continue
            rows.append({'id': sid, 'platform': 'ttof', 'label': 'na', 'wiki_id': w,
                         'precursor_mz': mz, 'polarity': pol,
                         'confidence': '', 'confidence_raw': '', 'ensemble_sd': ''})
            peaks_out[sid] = pk[w]
            n_ttof += 1
    print(f'TTOF control: {n_ttof}')

    with open(MGF, 'w') as f:
        for r in rows:
            f.write('BEGIN IONS\n')
            f.write(f'TITLE={r["id"]}\n')
            f.write(f'PEPMASS={r["precursor_mz"]}\n')
            f.write('CHARGE=1-\n' if r['polarity'] in ('neg', '0') else 'CHARGE=1+\n')
            for mz, it in peaks_out[r['id']]:
                f.write(f'{mz} {it}\n')
            f.write('END IONS\n')
    cols = ['id', 'platform', 'label', 'wiki_id', 'precursor_mz', 'polarity',
            'confidence', 'confidence_raw', 'ensemble_sd']
    with open(META, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(rows)
    print(f'\nWrote {MGF}\n      {META}  ({len(rows)} spectra)')
    print(f'\nNext: in a DreaMS env →  python code/bench_dreams_ood.py --embed')


def embed():
    # uv editable-install .pth hook doesn't always register; dreams is a flat pkg at DreaMS/.
    import sys
    repo = ROOT/'DreaMS'
    if repo.is_dir() and str(repo) not in sys.path:
        sys.path.insert(0, str(repo))
    from dreams.api import dreams_embeddings
    embs = np.asarray(dreams_embeddings(str(MGF)))
    np.save(EMB, embs)
    print(f'Wrote {EMB}  shape={embs.shape}')


def _auc(scores, pos):
    """AUC of `scores` predicting boolean `pos` (Mann-Whitney U / rank)."""
    s = np.asarray(scores, float); y = np.asarray(pos, bool)
    if y.sum() == 0 or (~y).sum() == 0: return float('nan')
    order = np.argsort(s); ranks = np.empty(len(s)); ranks[order] = np.arange(1, len(s)+1)
    return (ranks[y].sum() - y.sum()*(y.sum()+1)/2) / (y.sum()*(~y).sum())


def _ood_distance(query, ref, k=K_NN, exclude_self=False):
    """Mean cosine distance from each query row to its k nearest ref rows."""
    qn = query / (np.linalg.norm(query, axis=1, keepdims=True) + 1e-9)
    rn = ref / (np.linalg.norm(ref, axis=1, keepdims=True) + 1e-9)
    sim = qn @ rn.T                      # cosine similarity
    dist = 1.0 - sim
    if exclude_self:                     # query is subset of ref aligned by row → mask diagonal
        np.fill_diagonal(dist, np.inf)
    part = np.partition(dist, k, axis=1)[:, :k]
    return part.mean(axis=1)


def analyze():
    meta = list(csv.DictReader(open(META)))
    emb = np.load(EMB)
    assert len(meta) == len(emb), f'meta {len(meta)} vs emb {len(emb)}'
    plat = np.array([m['platform'] for m in meta])
    lab  = np.array([m['label'] for m in meta])
    oi = np.where(plat == 'orbi')[0]
    bi = np.where(plat == 'blank')[0]
    ti = np.where(plat == 'ttof')[0]
    O = emb[oi]
    conf_o = np.array([fnum(meta[i]['confidence']) for i in oi], float)
    is_tp = (lab[oi] == 'TP')
    conf_b = np.array([fnum(meta[i]['confidence']) for i in bi], float)
    print(f'orbi {len(oi)} (TP {int(is_tp.sum())}, FP {int((~is_tp).sum())})  '
          f'blank(should-abstain) {len(bi)}  ttof {len(ti)}')

    # Manifold = Orbitrap labeled embeddings (unsupervised, fit on TP+FP together).
    # Orbitrap OOD: out-of-fold (never own neighbour). Blank/TTOF OOD: vs all Orbitrap.
    ood_o = np.empty(len(oi)); rng = np.arange(len(oi)) % 5
    for f in range(5):
        te = rng == f; ood_o[te] = _ood_distance(O[te], O[~te], k=K_NN)
    ood_b = _ood_distance(emb[bi], O, k=K_NN)

    # TEST 1 (reframed) — does OOD separate should-abstain BLANKS from genuine confident TP?
    # Match confidence ranges so the score can't trivially separate: TP restricted to the
    # blank confidence span [min,max]; all blanks are positives.
    lo, hi = float(np.min(conf_b)), float(np.max(conf_b))
    tp_band = is_tp & (conf_o >= lo) & (conf_o <= hi)
    tp_idx = np.where(tp_band)[0]
    print(f'\n=== TEST 1 (reframed) — OOD flags high-confidence BLANKS vs genuine TP ===')
    print(f'  confidence span matched to blanks: [{lo:.3f}, {hi:.3f}]')
    print(f'  pool: genuine confident-TP {len(tp_idx)}   should-abstain blanks {len(bi)}')

    pooled_ood  = np.concatenate([ood_o[tp_idx], ood_b])
    pooled_conf = np.concatenate([conf_o[tp_idx], conf_b])
    is_blank    = np.concatenate([np.zeros(len(tp_idx), bool), np.ones(len(bi), bool)])

    auc_ood  = _auc(pooled_ood, is_blank)
    auc_conf = _auc(-pooled_conf, is_blank)   # lower confidence ⇒ more likely blank (score baseline)
    print(f'  median OOD: blanks {np.median(ood_b):.3f}  vs confident-TP {np.median(ood_o[tp_idx]):.3f}')
    print(f'  AUC(flag blank): OOD distance {auc_ood:.3f}   |   (1-confidence) baseline {auc_conf:.3f}')

    # risk-coverage: retain most in-distribution (low OOD); how fast does retained blank-rate drop?
    print('  risk-coverage (blank-rate among retained, keep lowest-criterion first):')
    base = is_blank.mean()
    for nm, crit in [('OOD', pooled_ood), ('1-conf', -pooled_conf)]:
        order = np.argsort(crit)
        line = [f'cov{int(c*100)}={is_blank[order[:max(1,int(c*len(order)))]].mean()*100:4.1f}%'
                for c in (0.6, 0.7, 0.8, 0.9, 1.0)]
        print(f'      {nm:8s} {"  ".join(line)}')
    print(f'      {"base":8s} {"  ".join(f"cov{int(c*100)}={base*100:4.1f}%" for c in (0.6,0.7,0.8,0.9,1.0))}')

    # TEST 2 — TTOF (real platform shift) must be more OOD than Orbitrap
    if len(ti):
        ood_t = _ood_distance(emb[ti], O, k=K_NN)
        p95 = np.percentile(ood_o, 95)
        print('\n=== TEST 2 — distribution-shift control (TTOF vs Orbitrap) ===')
        print(f'  median OOD   orbi {np.median(ood_o):.3f}   ttof {np.median(ood_t):.3f}')
        print(f'  TTOF above orbi-95th pct ({p95:.3f}): {np.mean(ood_t > p95)*100:.1f}%')

    print('\n--- verdict gates ---')
    print(f'  OOD AUC(flag blank) > 0.60 ?        {auc_ood:.3f}  {"PASS" if auc_ood>0.60 else "FAIL"}')
    print(f'  OOD beats (1-conf) by > +0.05 ?     {auc_ood-auc_conf:+.3f}  {"PASS" if auc_ood-auc_conf>0.05 else "FAIL"}')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--prep', action='store_true')
    ap.add_argument('--embed', action='store_true')
    ap.add_argument('--analyze', action='store_true')
    a = ap.parse_args()
    if a.prep: prep()
    elif a.embed: embed()
    elif a.analyze: analyze()
    else: ap.print_help()
