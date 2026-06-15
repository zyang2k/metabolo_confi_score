"""
build_tier_a_pdf.py — Same content as the PPTX review deck, rendered directly
to a multi-page PDF via matplotlib (no office suite required).

Each page: title/header block, mirror plot (query vs library), metadata footer.
Re-fetches query + reference peaks via MassWiki — set MASSWIKI_TOKEN env var.
"""

import os
import sys
import time
import json
import requests
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TIER_A_CSV    = os.path.join(ROOT, 'data', 'blank_high_confidence.csv')
NEG_CUR       = os.path.join(ROOT, 'data', 'Orbitrap_HILIC_negESI_curated_042126.csv')
POS_CUR       = os.path.join(ROOT, 'data', 'Orbitrap_HILIC_posESI_curated_042126.csv')
PDF_OUT       = os.path.join(ROOT, 'reports', 'tier_a_high_confidence_blanks_review.pdf')

BASE = 'https://masswiki.us-west-2.elasticbeanstalk.com'
GETDATA_URL = f'{BASE}/analysis/get_data'
LIB_URL     = f'{BASE}/reference_library/get_spectra_data'


def fetch_query_and_match(wiki_id: str, target_name: str, target_adduct: str, token: str):
    headers = {'Authorization': f'Bearer {token}'}
    for src, pub in [('binbase', 'false'), ('zyang2k', 'true')]:
        r = requests.get(GETDATA_URL,
                         params={'wiki_id': wiki_id, 'source': src, 'isPublic': pub},
                         headers=headers, timeout=20)
        if r.status_code == 200:
            break
    else:
        return {'error': f'HTTP {r.status_code}'}

    p = r.json()
    spec = p.get('spectrum') or {}
    isr = (p.get('analysis') or {}).get('reference_library', {}).get('identity_search') or []

    target_adduct_clean = str(target_adduct or '').strip().lower()
    target_name_clean = str(target_name or '').strip().lower()
    matched = None
    for e in isr:
        nm = str(e.get('name') or '').strip().lower()
        ad = str(e.get('adduct') or '').strip().lower()
        if nm == target_name_clean and ad == target_adduct_clean:
            matched = e; break
    if matched is None:
        for e in isr:
            if str(e.get('name') or '').strip().lower() == target_name_clean:
                matched = e; break
    if matched is None and isr:
        matched = isr[0]

    return {
        'query_peaks': spec.get('peaks') or [],
        'splash': spec.get('splash') or '',
        'precursor_mz': spec.get('precursor_mz'),
        'rt': spec.get('rt'),
        'matched_ref': matched or {},
    }


def fetch_lib_peaks_batch(library_wiki_ids):
    out = {}
    ids = [x for x in library_wiki_ids if isinstance(x, str) and x]
    for i in range(0, len(ids), 50):
        batch = ids[i:i+50]
        r = requests.post(LIB_URL, json={
            'id_list': batch, 'get_details': False, 'include_fields': ['peaks']
        }, timeout=30)
        if r.status_code == 200:
            for e in r.json():
                wid = e.get('wiki_id'); pk = e.get('peaks', [])
                if wid: out[wid] = pk
        time.sleep(0.1)
    return out


def draw_page(ax_mirror, q_peaks, r_peaks, title, tol=0.01):
    """Draw a mirror plot on `ax_mirror`. Query on top (blue), ref on bottom (orange=matched, gray=unmatched)."""
    def _norm(pks):
        if not pks: return np.array([]), np.array([])
        mz = np.asarray([p[0] for p in pks], dtype=float)
        it = np.asarray([p[1] for p in pks], dtype=float)
        mx = it.max() if len(it) else 1.0
        return mz, (it / mx * 100.0)

    q_mz, q_int = _norm(q_peaks)
    r_mz, r_int = _norm(r_peaks)
    q_mz_s = np.sort(q_mz) if len(q_mz) else np.array([])

    def matched(lmz):
        if len(q_mz_s) == 0: return False
        j = np.searchsorted(q_mz_s, lmz)
        d = min(
            abs(q_mz_s[j] - lmz) if j < len(q_mz_s) else np.inf,
            abs(q_mz_s[j-1] - lmz) if j > 0 else np.inf,
        )
        return d <= tol

    for m, i in zip(q_mz, q_int):
        ax_mirror.vlines(m, 0, i, color='steelblue', linewidth=1.2)
    for m, i in zip(r_mz, r_int):
        c = 'darkorange' if matched(m) else 'lightgray'
        ax_mirror.vlines(m, 0, -i, color=c, linewidth=1.2)

    ax_mirror.axhline(0, color='black', lw=0.6)
    all_mz = np.concatenate([q_mz, r_mz]) if (len(q_mz) or len(r_mz)) else np.array([0, 1])
    lo = max(0, all_mz.min() - 10); hi = all_mz.max() + 10
    ax_mirror.set_xlim(lo, hi); ax_mirror.set_ylim(-110, 110)
    ax_mirror.set_yticks([-100, -50, 0, 50, 100])
    ax_mirror.set_yticklabels(['100', '50', '0', '50', '100'])
    ax_mirror.set_xlabel('m/z')
    ax_mirror.set_ylabel('Rel. intensity\n(query ↑ / library ↓)')
    ax_mirror.set_title(title, fontsize=11, pad=4)

    ax_mirror.text(0.99, 0.96, f'query n={len(q_mz)}', transform=ax_mirror.transAxes,
                   ha='right', va='top', fontsize=8, color='steelblue')
    ax_mirror.text(0.99, 0.04, f'library n={len(r_mz)}', transform=ax_mirror.transAxes,
                   ha='right', va='bottom', fontsize=8, color='darkorange')


def make_page(pdf, row, info, lib_peaks_map, idx):
    """Create one PDF page for a bin."""
    ref = info.get('matched_ref') or {}
    lwid = ref.get('library_wiki_id', '')
    q_peaks = info.get('query_peaks') or []
    r_peaks = lib_peaks_map.get(lwid, [])
    splash = info.get('splash') or ''
    obs_mz = info.get('precursor_mz')
    obs_rt = info.get('rt')

    # Page: 13.33 x 7.5 inches (matches PPTX widescreen)
    fig = plt.figure(figsize=(13.33, 7.5))
    # Layout: title (top ~1.2"), mirror plot (middle ~4.5"), metadata (bottom ~1.5")
    # Use GridSpec to get proper vertical layout
    from matplotlib import gridspec
    gs = gridspec.GridSpec(3, 1, height_ratios=[1.2, 4.6, 1.4], hspace=0.35,
                           left=0.04, right=0.97, top=0.97, bottom=0.04)

    # Header
    ax_h = fig.add_subplot(gs[0]); ax_h.axis('off')
    ax_h.text(0.0, 0.85, f'#{idx+1:02d}   {row["pipeline_annotation"]}',
              fontsize=20, fontweight='bold', transform=ax_h.transAxes)
    ax_h.text(0.0, 0.50,
              f'Adduct: {row["adduct"]}    ·    Confidence: {row["confidence_pct"]:.1f}%    ·    '
              f'wiki_id: {row["wiki_id"]}',
              fontsize=13, color='#333333', transform=ax_h.transAxes)
    ax_h.text(0.0, 0.18,
              f'SPLASH (query): {splash or "(unknown)"}',
              fontsize=10, color='#666666', transform=ax_h.transAxes)

    # Mirror plot
    ax_m = fig.add_subplot(gs[1])
    draw_page(ax_m, q_peaks, r_peaks,
              title='Query vs. library reference spectrum')

    # Metadata footer
    ax_f = fig.add_subplot(gs[2]); ax_f.axis('off')
    obs_mz_s = f'{obs_mz:.4f}' if isinstance(obs_mz, (int, float)) else '-'
    obs_rt_s = f'{obs_rt:.1f}s' if isinstance(obs_rt, (int, float)) else '-'
    lib_mz = ref.get('precursor_mz', '-')
    rt_pred_hilic = ref.get('predicted_rt_hilic', '-')
    d_rt_pred = ref.get('delta_predicted_rt', '-')
    ref_splash = ref.get('splash', '') or '-'

    # Safer delta_mda formatting
    try:
        dmda = f'{float(row.get("delta_mda")):.2f}'
    except (TypeError, ValueError):
        dmda = str(row.get('delta_mda', '?'))

    lines = [
        f'Library: {row["db"]}    ·    IK14: {row["hit_ik14"]}    ·    library_wiki_id: {lwid or "(unknown)"}',
        f'Observed: m/z {obs_mz_s}, RT {obs_rt_s}    |    Library: m/z {lib_mz}, '
        f'RT_pred(HILIC) {rt_pred_hilic}, ΔRT_pred {d_rt_pred}s',
        f'entropy_sim {row["entropy_sim"]:.3f}    ·    sim_gap {row["sim_gap"]:.3f}    ·    '
        f'delta_mda {dmda}    ·    n_close_alternatives {int(row.get("n_close_alternatives", 0))}',
        f'Library SPLASH: {ref_splash}',
    ]
    y = 0.88
    for line in lines:
        ax_f.text(0.0, y, line, fontsize=10, color='#333333', transform=ax_f.transAxes)
        y -= 0.24

    pdf.savefig(fig, dpi=150)
    plt.close(fig)


def title_page(pdf):
    fig = plt.figure(figsize=(13.33, 7.5))
    fig.text(0.5, 0.66, 'Tier-A High-Confidence Blank Bins',
             ha='center', fontsize=36, fontweight='bold')
    fig.text(0.5, 0.55,
             '26 pos-mode bins scoring ≥ 90% confidence without a curator annotation',
             ha='center', fontsize=18, color='#555555')
    fig.text(0.5, 0.49,
             'Potential missed annotations for Oliver to review',
             ha='center', fontsize=18, color='#555555')
    fig.text(0.5, 0.28, f'Ziyue Yang · {time.strftime("%Y-%m-%d")}',
             ha='center', fontsize=14, color='#888888')
    pdf.savefig(fig, dpi=150)
    plt.close(fig)


def main():
    token = os.environ.get('MASSWIKI_TOKEN')
    if not token:
        sys.exit('Set MASSWIKI_TOKEN env var')

    print('Loading Tier A review list...')
    bh = pd.read_csv(TIER_A_CSV)
    tier_a = bh[bh['review_tier'] == 'A_confident_≥90%'].copy()
    tier_a = tier_a.sort_values('confidence_pct', ascending=False).reset_index(drop=True)
    print(f'  {len(tier_a)} bins')

    # SPLASH map from curated CSVs (raw_splash column)
    neg = pd.read_csv(NEG_CUR, low_memory=False)[['wiki_id', 'raw_splash']]
    pos = pd.read_csv(POS_CUR, low_memory=False)[['wiki_id', 'raw_splash']]
    splash_map = pd.concat([neg, pos], ignore_index=True).set_index('wiki_id')['raw_splash'].to_dict()

    print(f'Fetching query spectra + matched ref metadata for {len(tier_a)} bins...')
    bin_info = {}
    with ThreadPoolExecutor(max_workers=6) as ex:
        futs = {
            ex.submit(fetch_query_and_match, r['wiki_id'], r['pipeline_annotation'], r['adduct'], token): r['wiki_id']
            for _, r in tier_a.iterrows()
        }
        for f in as_completed(futs):
            wid = futs[f]
            try:
                bin_info[wid] = f.result()
            except Exception as e:
                bin_info[wid] = {'error': str(e)}
    n_err = sum('error' in v for v in bin_info.values())
    print(f'  done ({n_err} errors)')

    lwids = [ (bin_info.get(r['wiki_id'], {}).get('matched_ref') or {}).get('library_wiki_id', '')
              for _, r in tier_a.iterrows() ]
    lwids = [l for l in lwids if isinstance(l, str) and l]
    print(f'Fetching {len(lwids)} library peak lists (batched)...')
    lib_peaks = fetch_lib_peaks_batch(lwids)
    print(f'  retrieved {len(lib_peaks)} / {len(lwids)}')

    print(f'Writing {PDF_OUT}...')
    with PdfPages(PDF_OUT) as pdf:
        title_page(pdf)
        for idx, row in tier_a.iterrows():
            info = bin_info.get(row['wiki_id'], {})
            if 'error' in info:
                print(f'  SKIP {row["wiki_id"]}: {info["error"]}'); continue
            # prefer query splash from API; fall back to curated raw_splash
            if not info.get('splash'):
                info['splash'] = splash_map.get(row['wiki_id'], '')
            make_page(pdf, row, info, lib_peaks, idx)
            print(f'  page #{idx+1}: {row["pipeline_annotation"][:40]} ({row["confidence_pct"]:.0f}%)')

        # PDF metadata
        d = pdf.infodict()
        d['Title'] = 'Tier-A High-Confidence Blank Bins Review'
        d['Author'] = 'Ziyue Yang'
        d['Subject'] = 'metabolo_confi_score — 26 bins ≥ 90% confidence, no curator annotation'
        d['Keywords'] = 'LC-MS, MassWiki, confidence scoring, NoTA'
        d['CreationDate'] = None  # let matplotlib fill it

    print(f'\nSaved {PDF_OUT}')


if __name__ == '__main__':
    main()
