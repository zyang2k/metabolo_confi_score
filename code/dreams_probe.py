"""dreams_probe.py — LABEL-FREE probe of DreaMS embeddings (NOT a GBM feature test).

Why label-free: a richer spectral representation as a feature on curator labels hits the
~0.92 Q1/Q2 ceiling (MS2DeepScore +0.0002, Set Transformer Δ≈0). So we test DreaMS only
where the curator-label GBM is structurally helpless:

  TEST B (primary) — CROSS-PLATFORM transfer. For compounds present in BOTH Orbitrap and
    TTOF, do DreaMS embeddings retrieve the same compound across platforms better than raw
    entropy_similarity? If yes → DreaMS is platform-invariant → candidate engine to close
    the Orbitrap→TTOF gap (~0.85 vs 0.91) our hand-features can't.
  TEST C (secondary) — REACH. Do uncurated (blank) bins get sensible curated neighbours in
    embedding space (not precursor-locked, unlike our 6.8%-reach propagation)?
  TEST A (deferred) — Q0 era separation; needs the glucose-1-P replicate pull (not cached).

DreaMS: transformer foundation model, self-supervised on ~700M GeMS spectra, 1024-d
embedding per MS². API: `from dreams.api import dreams_embeddings; embs = dreams_embeddings(mgf)`.
Requires Python 3.11 + torch + Zenodo weights — CANNOT run in the 3.9 sandbox; run --embed
in a DreaMS env.

STAGES (argparse):
  --prep     build the cross-platform probe set → MGF + meta + peaks sidecar  (rdkit/ms_entropy env)
  --embed    dreams_embeddings(MGF) → data/dreams_probe/emb.npy               (DreaMS py3.11 env)
  --analyze  cross-platform retrieval: DreaMS cosine vs entropy_similarity     (any env)

Install DreaMS (separate env):
  git clone https://github.com/pluskal-lab/DreaMS.git && cd DreaMS
  conda create -n dreams python==3.11 --yes && conda activate dreams && pip install -e .
  pip install rdkit ms_entropy        # for --prep/--analyze in the same env
"""
from __future__ import annotations
import os, json, csv, argparse
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
OUT  = ROOT/'data'/'dreams_probe'
FT   = ROOT/'data'/'feature_table_v2.csv'
QP_ORBI = ROOT/'data'/'query_peaks_cache_v2.json'
TTOF_NEG_SPEC = ROOT/'data'/'TTOF_HILIC_negESI_uncurated_041326.csv'
TTOF_POS_SPEC = ROOT/'data'/'TTOF_HILIC_posESI_uncurated_041326.csv'
TTOF_NEG_PK = ROOT/'data'/'ttof_neg_query_peaks_cache.json'
TTOF_POS_PK = ROOT/'data'/'ttof_pos_query_peaks_cache.json'

MGF   = OUT/'probe.mgf'
META  = OUT/'meta.csv'
PEAKS = OUT/'peaks.json'
EMB   = OUT/'emb.npy'
K_PER = 5          # max spectra per compound per platform
MS2_TOL = 0.02


def fnum(x, d=None):
    try: return float(x)
    except (TypeError, ValueError): return d

def ik14_from_smiles(smi):
    from rdkit import Chem, RDLogger
    from rdkit.Chem.inchi import InchiToInchiKey, MolToInchi
    RDLogger.DisableLog('rdApp.*')
    if not isinstance(smi, str) or not smi.strip(): return ''
    m = Chem.MolFromSmiles(smi)
    if m is None: return ''
    ik = InchiToInchiKey(MolToInchi(m))
    return ik[:14] if ik else ''


# ───────────────────────── PREP ─────────────────────────
def prep():
    OUT.mkdir(parents=True, exist_ok=True)
    qp_orbi = {k: v for k, v in json.load(open(QP_ORBI)).items() if v}

    # Orbitrap: labeled TP bins with a clean IK14 (one row per wiki_id)
    orbi = {}
    for r in csv.DictReader(open(FT)):
        if r['spectrum_label'] != 'TP': continue
        ik = (r.get('anno_ik14') or '').strip()
        w = r['wiki_id']
        if not ik or w in orbi or w not in qp_orbi: continue
        mz = fnum(r.get('precursor_mz'))
        if mz is None: continue
        orbi[w] = {'ik14': ik, 'mz': mz, 'pol': r.get('polarity'), 'peaks': qp_orbi[w]}
    print(f'Orbitrap TP bins w/ ik14+peaks: {len(orbi):,}')

    # TTOF: annotation-smiles -> ik14, peaks from ttof caches
    ttof = {}
    for spec_f, pk_f, pol in [(TTOF_NEG_SPEC, TTOF_NEG_PK, 'neg'), (TTOF_POS_SPEC, TTOF_POS_PK, 'pos')]:
        pk = {k: v for k, v in json.load(open(pk_f)).items() if v}
        for r in csv.DictReader(open(spec_f)):
            w = r.get('wiki_id'); smi = r.get('annotation-smiles')
            if not w or w in ttof or w not in pk: continue
            ik = ik14_from_smiles(smi)
            mz = fnum(r.get('precursor_mz'))
            if not ik or mz is None: continue
            ttof[w] = {'ik14': ik, 'mz': mz, 'pol': pol, 'peaks': pk[w]}
    print(f'TTOF bins w/ ik14+peaks: {len(ttof):,}')

    shared = sorted(set(d['ik14'] for d in orbi.values()) & set(d['ik14'] for d in ttof.values()))
    print(f'SHARED compounds (ik14 in both platforms): {len(shared):,}')

    # collect up to K spectra per compound per platform
    def take(src, ik):
        hits = [(w, d) for w, d in src.items() if d['ik14'] == ik][:K_PER]
        return hits
    rows = []; peaks_out = {}
    for ik in shared:
        for plat, src in [('orbi', orbi), ('ttof', ttof)]:
            for w, d in take(src, ik):
                sid = f'{plat}|{ik}|{w}'
                rows.append({'id': sid, 'platform': plat, 'ik14': ik,
                             'precursor_mz': d['mz'], 'polarity': d['pol'], 'wiki_id': w})
                peaks_out[sid] = d['peaks']
    print(f'probe spectra: {len(rows):,}  (orbi {sum(r["platform"]=="orbi" for r in rows)}, '
          f'ttof {sum(r["platform"]=="ttof" for r in rows)})')

    # write MGF (DreaMS reads PEPMASS + peaks)
    with open(MGF, 'w') as f:
        for r in rows:
            f.write('BEGIN IONS\n')
            f.write(f'TITLE={r["id"]}\n')
            f.write(f'PEPMASS={r["precursor_mz"]}\n')
            f.write('CHARGE=1-\n' if r['polarity'] in ('neg', '0') else 'CHARGE=1+\n')
            for mz, it in peaks_out[r['id']]:
                f.write(f'{mz} {it}\n')
            f.write('END IONS\n')
    with open(META, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['id', 'platform', 'ik14', 'precursor_mz', 'polarity', 'wiki_id'])
        w.writeheader(); w.writerows(rows)
    json.dump(peaks_out, open(PEAKS, 'w'))
    print(f'\nWrote {MGF}\n      {META}\n      {PEAKS}')
    print(f'\nNext: in a DreaMS env →  python code/dreams_probe.py --embed')


# ───────────────────────── EMBED (DreaMS env) ─────────────────────────
def embed():
    # uv editable-install .pth import-hook doesn't always register at startup; the
    # `dreams` package is a flat pkg at DreaMS/dreams, so add the repo root as fallback.
    import sys
    repo = ROOT/'DreaMS'
    if repo.is_dir() and str(repo) not in sys.path:
        sys.path.insert(0, str(repo))
    from dreams.api import dreams_embeddings   # requires py3.11 + torch + weights
    embs = dreams_embeddings(str(MGF))         # (N, 1024), MGF order
    np.save(EMB, np.asarray(embs))
    print(f'Wrote {EMB}  shape={np.asarray(embs).shape}')


# ───────────────────────── ANALYZE ─────────────────────────
def _cos_matrix(A, B):
    An = A / (np.linalg.norm(A, axis=1, keepdims=True) + 1e-9)
    Bn = B / (np.linalg.norm(B, axis=1, keepdims=True) + 1e-9)
    return An @ Bn.T

def analyze():
    import ms_entropy as me
    meta = list(csv.DictReader(open(META)))
    emb = np.load(EMB)
    peaks = json.load(open(PEAKS))
    assert len(meta) == len(emb), f'meta {len(meta)} vs emb {len(emb)}'
    plat = np.array([m['platform'] for m in meta])
    ik = np.array([m['ik14'] for m in meta])
    oi = np.where(plat == 'orbi')[0]; ti = np.where(plat == 'ttof')[0]
    print(f'orbi {len(oi)}  ttof {len(ti)}  shared-compounds {len(set(ik))}')

    # TEST B — cross-platform retrieval: each Orbitrap spectrum → nearest TTOF spectrum
    # (a) by DreaMS embedding cosine
    C = _cos_matrix(emb[oi], emb[ti])
    nn_dreams = ti[C.argmax(axis=1)]
    hit_dreams = np.mean(ik[oi] == ik[nn_dreams])
    # (b) by entropy_similarity on raw peaks (baseline)
    ids = [m['id'] for m in meta]
    def P(i): return np.asarray(peaks[ids[i]], float)
    hit_ent = 0
    for a in oi:
        best, bj = -1, -1
        qa = P(a)
        for b in ti:
            s = me.calculate_entropy_similarity(qa, P(b), ms2_tolerance_in_da=MS2_TOL, clean_spectra=True)
            if s > best: best, bj = s, b
        if ik[a] == ik[bj]: hit_ent += 1
    hit_ent /= max(len(oi), 1)

    print('\n=== TEST B — cross-platform same-compound retrieval (Orbitrap → TTOF) ===')
    print(f'  DreaMS embedding cosine : top-1 same-IK14 = {hit_dreams*100:.1f}%')
    print(f'  entropy_similarity      : top-1 same-IK14 = {hit_ent*100:.1f}%')
    print(f'  → DreaMS {"WINS" if hit_dreams>hit_ent else "does NOT beat"} entropy on platform transfer')

    # within-platform sanity: same-ik14 vs diff-ik14 mean cosine
    for name, idx in [('orbi', oi), ('ttof', ti)]:
        if len(idx) < 3: continue
        Cm = _cos_matrix(emb[idx], emb[idx]); iks = ik[idx]
        same = []; diff = []
        for a in range(len(idx)):
            for b in range(a+1, len(idx)):
                (same if iks[a] == iks[b] else diff).append(Cm[a, b])
        print(f'  [{name}] mean cosine  same-IK14 {np.mean(same):.3f}  vs diff-IK14 {np.mean(diff):.3f}')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--prep', action='store_true')
    ap.add_argument('--embed', action='store_true')
    ap.add_argument('--analyze', action='store_true')
    a = ap.parse_args()
    if a.prep: prep()
    elif a.embed: embed()
    elif a.analyze: analyze()
    else: ap.print_help()
