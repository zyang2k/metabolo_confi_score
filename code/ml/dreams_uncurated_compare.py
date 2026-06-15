"""dreams_uncurated_compare.py — Phase 3: did domain-adaptive SSL help?

Re-runs the uncurated-organization buckets on the SAME probe set, for the baseline backbone
vs the domain-adapted backbone (both produced on Colab via dreams_intermediates). The test:
do the rich-but-isolated uncurated blanks move into near-known / coherent-unknown after
adaptation? If isolated% drops meaningfully → self-supervision on the uncurated bins works.

Usage (after downloading the two npy from Colab into data/dreams_uncurated/):
  python code/dreams_uncurated_compare.py
"""
import csv
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
D = ROOT / 'data' / 'dreams_uncurated'
META = D / 'meta.csv'
EMB = {'baseline backbone': D / 'emb_base_backbone.npy',
       'adapted backbone':  D / 'emb_adapted_backbone.npy'}
TAU = 0.85


def cos(A, B):
    An = A / (np.linalg.norm(A, axis=1, keepdims=True) + 1e-9)
    Bn = B / (np.linalg.norm(B, axis=1, keepdims=True) + 1e-9)
    return An @ Bn.T


def buckets(emb, kind):
    ci = np.where(kind == 'curated')[0]
    bi = np.where(kind == 'uncurated')[0]
    nn_cur = cos(emb[bi], emb[ci]).max(axis=1)
    Cbb = cos(emb[bi], emb[bi]); np.fill_diagonal(Cbb, -1)
    nn_blk = Cbb.max(axis=1)
    near = nn_cur >= TAU
    coh = (~near) & (nn_blk >= TAU)
    iso = (~near) & (~coh)
    return dict(near=near.mean(), coh=coh.mean(), iso=iso.mean(),
                med_nncur=float(np.median(nn_cur)))


def main():
    meta = list(csv.DictReader(open(META)))
    kind = np.array([m['kind'] for m in meta])
    print(f'probe: {(kind=="curated").sum()} curated  {(kind=="uncurated").sum()} uncurated   τ={TAU}\n')
    print(f'{"representation":20s} {"near-known":>11s} {"coherent":>10s} {"isolated":>10s} {"med nn-cur":>11s}')
    res = {}
    for name, p in EMB.items():
        if not p.exists():
            print(f'{name:20s}   (missing {p.name} — run the Colab step first)')
            continue
        emb = np.load(p)
        assert len(emb) == len(meta), f'{name}: emb {len(emb)} vs meta {len(meta)}'
        b = buckets(emb, kind); res[name] = b
        print(f'{name:20s} {b["near"]*100:10.1f}% {b["coh"]*100:9.1f}% {b["iso"]*100:9.1f}% {b["med_nncur"]:11.3f}')
    if len(res) == 2:
        base, adpt = res['baseline backbone'], res['adapted backbone']
        d_iso = (adpt['iso'] - base['iso']) * 100
        d_org = ((adpt['near'] + adpt['coh']) - (base['near'] + base['coh'])) * 100
        print(f'\nΔ isolated   {d_iso:+.1f} pp   (negative = adaptation pulled blanks into structure)')
        print(f'Δ organized  {d_org:+.1f} pp   (near-known + coherent)')
        verdict = 'SSL adaptation HELPS — uncurated bins gained structure' if d_org > 5 else \
                  'no meaningful gain — adaptation did not organize the uncurated population'
        print(f'\nVERDICT: {verdict}')


if __name__ == '__main__':
    main()
