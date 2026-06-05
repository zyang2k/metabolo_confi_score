"""dreams_uncurated_probe.py — does the DreaMS embedding ORGANIZE the uncurated bins?

Tests the premise behind "borrow self-supervision to use the uncurated bins": before any
fine-tuning, do the (frozen) DreaMS embeddings already place uncurated blank bins sensibly
relative to the curated reference? For each uncurated bin we ask its nearest curated
neighbour in embedding space and bucket:

  near-known      nn_curated cos >= TAU_HIGH  -> looks like a known compound (annotation
                  candidate / carryover) — a curator-reviewable hit
  coherent-unknown  far from curated BUT close to other uncurated bins -> recurrent unknown
  isolated/noise  far from everything

Reference: curated bins' nearest *different-compound* curated cosine, to calibrate what
"near a known compound" looks like for genuine compounds. If uncurated bins land like real
compounds, the representation is already useful on the unlabeled population (→ SSL fine-tune
is the upgrade). If they're all isolated noise, fine-tuning probably won't rescue it.

STAGES: --prep (3.9 env) | --embed (DreaMS env) | --analyze (any)
"""
from __future__ import annotations
import os, sys, json, csv, argparse, random
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / 'data' / 'dreams_uncurated'
FT = ROOT / 'data' / 'feature_table_v2.csv'
QP_CUR = ROOT / 'data' / 'query_peaks_cache_v2.json'
BLANK_PK = {'pos': ROOT / 'data' / 'blanks_query_peaks_cache_pos.json',
            'neg': ROOT / 'data' / 'blanks_query_peaks_cache_neg.json'}
BLANK_HITS = {'pos': ROOT / 'data' / 'pos_blanks_hits.csv',
              'neg': ROOT / 'data' / 'neg_blanks_hits.csv'}

MGF = OUT / 'probe.mgf'
META = OUT / 'meta.csv'
PEAKS = OUT / 'peaks.json'
EMB = OUT / 'emb.npy'
N_CUR = 2500          # curated reference sample
N_BLANK = 2500        # uncurated query sample
TAU_HIGH = 0.85       # "near a known compound" cosine
SEED = 42


def fnum(x, d=None):
    try: return float(x)
    except (TypeError, ValueError): return d


def prep():
    OUT.mkdir(parents=True, exist_ok=True)
    rng = random.Random(SEED)

    # curated reference: TP bins with anno_ik14 + precursor + peaks
    cur_peaks = {k: v for k, v in json.load(open(QP_CUR)).items() if v}
    cur = {}
    for r in csv.DictReader(open(FT)):
        if r['spectrum_label'] != 'TP':
            continue
        w = r['wiki_id']
        if w in cur or w not in cur_peaks:
            continue
        ik = (r.get('anno_ik14') or '').strip()
        mz = fnum(r.get('precursor_mz'))
        if not ik or mz is None:
            continue
        cur[w] = {'ik14': ik, 'name': (r.get('name') or '').strip(),
                  'mz': mz, 'pol': r.get('polarity'), 'peaks': cur_peaks[w]}
    cur_ids = rng.sample(list(cur), min(N_CUR, len(cur)))
    print(f'curated TP pool {len(cur):,} → sampled {len(cur_ids):,}')

    # uncurated blanks: precursor estimated as median lib_precursor_mz across the bin's hits
    rows, peaks_out = [], {}
    for w in cur_ids:
        d = cur[w]
        sid = f'cur|{w}'
        rows.append({'id': sid, 'kind': 'curated', 'ik14': d['ik14'], 'name': d['name'],
                     'precursor_mz': d['mz'], 'polarity': d['pol'], 'wiki_id': w})
        peaks_out[sid] = d['peaks']

    blank_ids_all = []
    for pol in ('pos', 'neg'):
        bpk = {k: v for k, v in json.load(open(BLANK_PK[pol])).items() if v}
        # median observed precursor per blank from its library hits
        acc = {}
        for r in csv.DictReader(open(BLANK_HITS[pol])):
            w = r.get('wiki_id'); mz = fnum(r.get('lib_precursor_mz'))
            if w in bpk and mz is not None:
                acc.setdefault(w, []).append(mz)
        cand = [w for w in bpk if w in acc]
        take = rng.sample(cand, min(N_BLANK // 2, len(cand)))
        print(f'  {pol} blanks with peaks+precursor {len(cand):,} → sampled {len(take):,}')
        for w in take:
            sid = f'blk|{w}'
            rows.append({'id': sid, 'kind': 'uncurated', 'ik14': '', 'name': '',
                         'precursor_mz': float(np.median(acc[w])), 'polarity': pol, 'wiki_id': w})
            peaks_out[sid] = bpk[w]
            blank_ids_all.append(sid)
    print(f'probe spectra: {len(rows):,} '
          f'(curated {sum(r["kind"]=="curated" for r in rows)}, '
          f'uncurated {sum(r["kind"]=="uncurated" for r in rows)})')

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
        wr = csv.DictWriter(f, fieldnames=['id', 'kind', 'ik14', 'name', 'precursor_mz', 'polarity', 'wiki_id'])
        wr.writeheader(); wr.writerows(rows)
    json.dump(peaks_out, open(PEAKS, 'w'))
    print(f'\nWrote {MGF}\n      {META}\n      {PEAKS}')
    print('Next: in a DreaMS env →  python code/dreams_uncurated_probe.py --embed')


def embed():
    repo = ROOT / 'DreaMS'
    if repo.is_dir() and str(repo) not in sys.path:
        sys.path.insert(0, str(repo))
    from dreams.api import dreams_embeddings
    embs = np.asarray(dreams_embeddings(str(MGF)))
    np.save(EMB, embs)
    print(f'Wrote {EMB}  shape={embs.shape}')


def _cos(A, B):
    An = A / (np.linalg.norm(A, axis=1, keepdims=True) + 1e-9)
    Bn = B / (np.linalg.norm(B, axis=1, keepdims=True) + 1e-9)
    return An @ Bn.T


def _hist(x, edges):
    x = np.asarray(x)
    return '  '.join(f'{lo:.2f}-{hi:.2f}:{((x>=lo)&(x<hi)).mean()*100:4.1f}%'
                     for lo, hi in zip(edges[:-1], edges[1:]))


def analyze():
    meta = list(csv.DictReader(open(META)))
    emb = np.load(EMB)
    assert len(meta) == len(emb)
    kind = np.array([m['kind'] for m in meta])
    ik = np.array([m['ik14'] for m in meta])
    name = np.array([m['name'] for m in meta])
    ci = np.where(kind == 'curated')[0]
    bi = np.where(kind == 'uncurated')[0]
    print(f'curated {len(ci)}  uncurated {len(bi)}')

    # blanks → nearest curated
    Cb = _cos(emb[bi], emb[ci])
    nn_cur = Cb.max(axis=1)
    nn_cur_j = ci[Cb.argmax(axis=1)]
    # blanks → nearest other blank
    Cbb = _cos(emb[bi], emb[bi]); np.fill_diagonal(Cbb, -1)
    nn_blk = Cbb.max(axis=1)
    # reference: curated → nearest DIFFERENT-compound curated
    Cc = _cos(emb[ci], emb[ci]); np.fill_diagonal(Cc, -1)
    ikc = ik[ci]
    nn_diff = np.array([Cc[a, ikc != ikc[a]].max() if (ikc != ikc[a]).any() else -1
                        for a in range(len(ci))])

    edges = np.array([0, .3, .5, .7, .85, .95, 1.001])
    print('\nnearest-neighbour cosine distributions:')
    print(f'  curated→diff-compound curated : {_hist(nn_diff[nn_diff>=0], edges)}')
    print(f'  uncurated→nearest curated     : {_hist(nn_cur, edges)}')

    near = nn_cur >= TAU_HIGH
    coherent = (~near) & (nn_blk >= TAU_HIGH)
    isolated = (~near) & (~coherent)
    print(f'\nuncurated bin buckets (TAU={TAU_HIGH}):')
    print(f'  near-known (annotation candidate) : {near.sum():4d}  ({near.mean()*100:.1f}%)')
    print(f'  coherent-unknown (blank cluster)  : {coherent.sum():4d}  ({coherent.mean()*100:.1f}%)')
    print(f'  isolated / noise                  : {isolated.sum():4d}  ({isolated.mean()*100:.1f}%)')

    print('\nexample near-known matches (uncurated blank → nearest curated compound):')
    order = np.argsort(-nn_cur)
    shown = 0
    for k in order:
        if nn_cur[k] < TAU_HIGH:
            break
        nm = name[nn_cur_j[k]] or ik[nn_cur_j[k]]
        print(f'  {meta[bi[k]]["wiki_id"]:16s}  cos={nn_cur[k]:.3f}  →  {nm[:40]}')
        shown += 1
        if shown >= 12:
            break


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
