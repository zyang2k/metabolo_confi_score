"""dreams_ssl_corpus.py — assemble the domain-adaptation corpus for continued SSL.

Borrows DreaMS's self-supervised recipe to exploit our UNCURATED bins (the label-ceiling
escape). SSL needs no labels, so the corpus is EVERY spectrum we have peaks for — curated
(query_peaks_cache_v2) + uncurated blanks (blanks_query_peaks_cache_{pos,neg}) — to teach
the model our Orbitrap-HILIC spectral structure (the OOD gap the frozen probe exposed).

  --prep  build corpus.mgf (+ meta)              (3.9 env; precursor for blanks = median
                                                  lib_precursor_mz across the bin's hits)
  --hdf5  MSData.from_mgf(corpus.mgf) -> .hdf5    (DreaMS env; the training input format)

Then continued pretraining runs on GPU (Colab) from ssl_model.ckpt — see notebook.
"""
from __future__ import annotations
import os, sys, json, csv, argparse
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / 'data' / 'dreams_ssl'
FT = ROOT / 'data' / 'feature_table_v2.csv'
QP_CUR = ROOT / 'data' / 'query_peaks_cache_v2.json'
BLANK_PK = {'pos': ROOT / 'data' / 'blanks_query_peaks_cache_pos.json',
            'neg': ROOT / 'data' / 'blanks_query_peaks_cache_neg.json'}
BLANK_HITS = {'pos': ROOT / 'data' / 'pos_blanks_hits.csv',
              'neg': ROOT / 'data' / 'neg_blanks_hits.csv'}
MGF = OUT / 'corpus.mgf'
META = OUT / 'meta.csv'
HDF5 = OUT / 'corpus.hdf5'


def fnum(x, d=None):
    try: return float(x)
    except (TypeError, ValueError): return d


def prep():
    OUT.mkdir(parents=True, exist_ok=True)
    rows, n_written = [], 0
    with open(MGF, 'w') as f:
        # curated: precursor + polarity from feature_table (one row per wiki_id)
        cur_peaks = {k: v for k, v in json.load(open(QP_CUR)).items() if v}
        seen = set()
        for r in csv.DictReader(open(FT)):
            w = r['wiki_id']
            if w in seen or w not in cur_peaks:
                continue
            mz = fnum(r.get('precursor_mz'))
            if mz is None:
                continue
            seen.add(w)
            pol = r.get('polarity')
            f.write('BEGIN IONS\n'); f.write(f'TITLE=cur|{w}\n'); f.write(f'PEPMASS={mz}\n')
            f.write('CHARGE=1-\n' if pol in ('neg', '0') else 'CHARGE=1+\n')
            for m, it in cur_peaks[w]:
                f.write(f'{m} {it}\n')
            f.write('END IONS\n')
            rows.append({'id': f'cur|{w}', 'kind': 'curated', 'polarity': pol}); n_written += 1
        n_cur = n_written
        print(f'curated spectra: {n_cur:,}')

        # uncurated blanks: precursor = median lib_precursor_mz across the bin's hits
        for pol in ('pos', 'neg'):
            bpk = {k: v for k, v in json.load(open(BLANK_PK[pol])).items() if v}
            acc = {}
            for r in csv.DictReader(open(BLANK_HITS[pol])):
                w = r.get('wiki_id'); mz = fnum(r.get('lib_precursor_mz'))
                if w in bpk and mz is not None:
                    acc.setdefault(w, []).append(mz)
            for w in bpk:
                if w not in acc:
                    continue
                mz = float(np.median(acc[w]))
                f.write('BEGIN IONS\n'); f.write(f'TITLE=blk|{w}\n'); f.write(f'PEPMASS={mz}\n')
                f.write('CHARGE=1-\n' if pol == 'neg' else 'CHARGE=1+\n')
                for m, it in bpk[w]:
                    f.write(f'{m} {it}\n')
                f.write('END IONS\n')
                rows.append({'id': f'blk|{w}', 'kind': 'uncurated', 'polarity': pol}); n_written += 1
            print(f'  {pol} blanks added: {n_written - (n_cur if pol=="pos" else 0):,}' if False else
                  f'  {pol} blanks cumulative total: {n_written:,}')

    with open(META, 'w', newline='') as fh:
        wr = csv.DictWriter(fh, fieldnames=['id', 'kind', 'polarity']); wr.writeheader(); wr.writerows(rows)
    print(f'\nCORPUS: {n_written:,} spectra ({n_cur:,} curated + {n_written-n_cur:,} uncurated)')
    print(f'Wrote {MGF}\n      {META}')
    print('Next (DreaMS env): python code/dreams_ssl_corpus.py --hdf5')


def hdf5():
    repo = ROOT / 'DreaMS'
    if repo.is_dir() and str(repo) not in sys.path:
        sys.path.insert(0, str(repo))
    from dreams.utils.data import MSData
    print(f'Converting {MGF} → hdf5 ...')
    msd = MSData.from_mgf(str(MGF))            # auto-writes corpus.hdf5 alongside
    print(f'MSData: {msd.num_spectra if hasattr(msd, "num_spectra") else "?"} spectra')
    print(f'Expect {HDF5} (or .hdf5 next to the mgf). Listing data/dreams_ssl:')
    for p in sorted(OUT.glob('*')):
        print(f'  {p.name}  {p.stat().st_size/1e6:.1f} MB')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--prep', action='store_true')
    ap.add_argument('--hdf5', action='store_true')
    a = ap.parse_args()
    if a.prep: prep()
    elif a.hdf5: hdf5()
    else: ap.print_help()
