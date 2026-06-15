"""
build_tier_a_pptx.py — Generate a PPTX review deck for the 26 Tier-A
high-confidence blank bins (confidence ≥ 90%, no curator annotation yet).

Each slide has:
  • Header: pipeline annotation + adduct + confidence + wiki_id + SPLASH
  • Mirror plot: query spectrum on top (blue, up), reference spectrum on bottom (red, down)
  • Metadata footer: DB, IK14, observed precursor/RT, entropy_sim, sim_gap, n_close_alternatives

Data pulled fresh from MassWiki:
  • Query peaks + SPLASH via /analysis/get_data (per wiki_id, authenticated)
  • Reference peaks via /reference_library/get_spectra_data (batched, unauth)

Token via MASSWIKI_TOKEN env var. Never written to disk.
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

from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_SHAPE

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TIER_A_CSV    = os.path.join(ROOT, 'data', 'blank_high_confidence.csv')
NEG_CUR       = os.path.join(ROOT, 'data', 'Orbitrap_HILIC_negESI_curated_042126.csv')
POS_CUR       = os.path.join(ROOT, 'data', 'Orbitrap_HILIC_posESI_curated_042126.csv')
OUT_DIR       = os.path.join(ROOT, 'data', 'tier_a_pptx_assets')
PPTX_OUT      = os.path.join(ROOT, 'reports', 'tier_a_high_confidence_blanks_review.pptx')

BASE = 'https://masswiki.us-west-2.elasticbeanstalk.com'
GETDATA_URL = f'{BASE}/analysis/get_data'
LIB_URL     = f'{BASE}/reference_library/get_spectra_data'

os.makedirs(OUT_DIR, exist_ok=True)


def fetch_query_and_match(wiki_id: str, target_name: str, target_adduct: str, token: str):
    """Fetch one bin's query spectrum + identify the matched reference entry.

    Returns dict with query_peaks, splash, precursor_mz, rt, matched_ref (dict).
    """
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

    # Find matched ref entry by name + adduct (flexible match)
    target_adduct_clean = str(target_adduct or '').strip().lower()
    target_name_clean = str(target_name or '').strip().lower()
    matched = None
    for e in isr:
        nm = str(e.get('name') or '').strip().lower()
        ad = str(e.get('adduct') or '').strip().lower()
        if nm == target_name_clean and ad == target_adduct_clean:
            matched = e
            break
    if matched is None:
        # fall back to name-only match
        for e in isr:
            if str(e.get('name') or '').strip().lower() == target_name_clean:
                matched = e
                break
    if matched is None and isr:
        matched = isr[0]  # last resort

    return {
        'query_peaks': spec.get('peaks') or [],
        'splash': spec.get('splash') or '',
        'precursor_mz': spec.get('precursor_mz'),
        'rt': spec.get('rt'),
        'matched_ref': matched or {},
    }


def fetch_lib_peaks_batch(library_wiki_ids):
    """POST to /reference_library/get_spectra_data to fetch peaks for a batch of lwids."""
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


def plot_mirror(query_peaks, ref_peaks, title, subtitle, out_png, tol=0.01):
    """Two-pane mirror plot: query up, reference down. Both L1-normalized for visual fairness."""
    fig, ax = plt.subplots(figsize=(9.5, 4.5))

    def _norm(pks):
        if not pks: return np.array([]), np.array([])
        mz = np.asarray([p[0] for p in pks], dtype=float)
        it = np.asarray([p[1] for p in pks], dtype=float)
        mx = it.max() if len(it) else 1.0
        return mz, (it / mx * 100.0)

    q_mz, q_int = _norm(query_peaks)
    r_mz, r_int = _norm(ref_peaks)

    # Matched pairs: library peak with a query peak within tol
    q_mz_s = np.sort(q_mz) if len(q_mz) else np.array([])

    def matched(lmz):
        if len(q_mz_s) == 0: return False
        j = np.searchsorted(q_mz_s, lmz)
        d = min(
            abs(q_mz_s[j] - lmz) if j < len(q_mz_s) else np.inf,
            abs(q_mz_s[j-1] - lmz) if j > 0 else np.inf,
        )
        return d <= tol

    # Top (query, positive)
    if len(q_mz):
        for m, i in zip(q_mz, q_int):
            ax.vlines(m, 0, i, color='steelblue', linewidth=1.2)
    # Bottom (reference, negative)
    if len(r_mz):
        for m, i in zip(r_mz, r_int):
            c = 'darkorange' if matched(m) else 'lightgray'
            ax.vlines(m, 0, -i, color=c, linewidth=1.2)

    ax.axhline(0, color='black', lw=0.6)
    all_mz = np.concatenate([q_mz, r_mz]) if (len(q_mz) or len(r_mz)) else np.array([0, 1])
    lo = max(0, all_mz.min() - 10); hi = all_mz.max() + 10
    ax.set_xlim(lo, hi); ax.set_ylim(-110, 110)
    ax.set_yticks([-100, -50, 0, 50, 100])
    ax.set_yticklabels(['100', '50', '0', '50', '100'])
    ax.set_xlabel('m/z')
    ax.set_ylabel('Rel. intensity\n(query ↑ / library ↓)')
    ax.set_title(title, fontsize=11, pad=4)
    if subtitle:
        ax.text(0.5, 1.01, subtitle, transform=ax.transAxes, ha='center', va='bottom',
                fontsize=8, color='gray')

    # Annotate # peaks
    ax.text(0.99, 0.96, f'query n={len(q_mz)}', transform=ax.transAxes,
            ha='right', va='top', fontsize=8, color='steelblue')
    ax.text(0.99, 0.04, f'library n={len(r_mz)}', transform=ax.transAxes,
            ha='right', va='bottom', fontsize=8, color='darkorange')

    plt.tight_layout()
    fig.savefig(out_png, dpi=150, bbox_inches='tight')
    plt.close(fig)


def main():
    token = os.environ.get('MASSWIKI_TOKEN')
    if not token:
        sys.exit('Set MASSWIKI_TOKEN env var')

    # Load tier A bins + per-bin metadata
    print('Loading Tier A review list...')
    bh = pd.read_csv(TIER_A_CSV)
    tier_a = bh[bh['review_tier'] == 'A_confident_≥90%'].copy()
    tier_a = tier_a.sort_values('confidence_pct', ascending=False).reset_index(drop=True)
    print(f'  {len(tier_a)} bins')

    # Get SPLASH from curated CSVs
    neg = pd.read_csv(NEG_CUR, low_memory=False)[['wiki_id','raw_splash']]
    pos = pd.read_csv(POS_CUR, low_memory=False)[['wiki_id','raw_splash']]
    splash_map = pd.concat([neg, pos], ignore_index=True).set_index('wiki_id')['raw_splash'].to_dict()

    # Fetch query spectra + find matched ref entries
    print(f'Fetching query spectra + ref metadata for {len(tier_a)} bins...')
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
    print(f'  done ({sum("error" in v for v in bin_info.values())} errors)')

    # Collect library_wiki_ids to fetch peaks
    lwids = []
    for wid, info in bin_info.items():
        lwid = (info.get('matched_ref') or {}).get('library_wiki_id')
        if isinstance(lwid, str) and lwid:
            lwids.append(lwid)
    print(f'Fetching {len(lwids)} library peak lists (batched)...')
    lib_peaks = fetch_lib_peaks_batch(lwids)
    print(f'  retrieved {len(lib_peaks)} / {len(lwids)}')

    # Build PPTX
    print('Building PPTX...')
    prs = Presentation()
    prs.slide_width = Inches(13.333)
    prs.slide_height = Inches(7.5)
    blank = prs.slide_layouts[6]  # blank layout

    # Title slide
    s = prs.slides.add_slide(blank)
    tx = s.shapes.add_textbox(Inches(0.5), Inches(2.5), Inches(12), Inches(2))
    tf = tx.text_frame
    tf.word_wrap = True
    p = tf.paragraphs[0]; p.alignment = 1
    r = p.add_run()
    r.text = 'Tier-A High-Confidence Blank Bins'
    r.font.size = Pt(36); r.font.bold = True
    p2 = tf.add_paragraph(); p2.alignment = 1
    r2 = p2.add_run()
    r2.text = (f'26 pos-mode bins scoring ≥ 90% confidence without a curator annotation\n'
               f'Potential missed annotations for Oliver to review')
    r2.font.size = Pt(18); r2.font.color.rgb = RGBColor(0x55, 0x55, 0x55)
    p3 = tf.add_paragraph(); p3.alignment = 1
    r3 = p3.add_run()
    r3.text = f'Ziyue Yang · {time.strftime("%Y-%m-%d")}'
    r3.font.size = Pt(14); r3.font.color.rgb = RGBColor(0x88, 0x88, 0x88)

    for idx, row in tier_a.iterrows():
        wid = row['wiki_id']
        info = bin_info.get(wid, {})
        if 'error' in info:
            print(f'  SKIP {wid}: {info["error"]}'); continue

        ref = info.get('matched_ref') or {}
        lwid = ref.get('library_wiki_id', '')
        ref_peaks = lib_peaks.get(lwid, [])
        q_peaks = info.get('query_peaks') or []
        splash = info.get('splash') or splash_map.get(wid, '')

        # Plot mirror
        title = f"{row['pipeline_annotation']}   [{row['adduct']}]"
        subtitle = (f"Query spectrum (top, blue) vs. library match (bottom, orange=matched, gray=unmatched)")
        png = os.path.join(OUT_DIR, f'{wid.replace("/","_")}.png')
        try:
            plot_mirror(q_peaks, ref_peaks, title, subtitle, png)
        except Exception as e:
            print(f'  plot error {wid}: {e}'); continue

        # Build slide
        s = prs.slides.add_slide(blank)

        # Header box (top)
        hx = s.shapes.add_textbox(Inches(0.35), Inches(0.2), Inches(12.6), Inches(0.95))
        hf = hx.text_frame; hf.word_wrap = True
        p = hf.paragraphs[0]
        r = p.add_run()
        r.text = f'#{idx+1:02d}   {row["pipeline_annotation"]}   ·   [{row["adduct"]}]'
        r.font.size = Pt(22); r.font.bold = True
        p2 = hf.add_paragraph()
        r2 = p2.add_run()
        r2.text = f'Confidence: {row["confidence_pct"]:.1f}%   ·   wiki_id: {wid}'
        r2.font.size = Pt(14); r2.font.color.rgb = RGBColor(0x33, 0x33, 0x33)
        p3 = hf.add_paragraph()
        r3 = p3.add_run()
        r3.text = f'SPLASH (query): {splash}'
        r3.font.size = Pt(11); r3.font.color.rgb = RGBColor(0x88, 0x88, 0x88)

        # Mirror plot (middle)
        s.shapes.add_picture(png, Inches(0.35), Inches(1.4), width=Inches(12.6))

        # Metadata footer (bottom)
        fx = s.shapes.add_textbox(Inches(0.35), Inches(6.1), Inches(12.6), Inches(1.3))
        ff = fx.text_frame; ff.word_wrap = True
        entry = info.get('matched_ref') or {}
        ref_splash = entry.get('splash', '')
        lines = [
            f'Library source: {row["db"]}   ·   IK14: {row["hit_ik14"]}   ·   Library entry: {lwid or "(unknown)"}',
            f'Observed  m/z {info.get("precursor_mz", "-"):.4f}' +
            (f'   RT {info.get("rt", 0):.1f}s' if info.get("rt") is not None else '   RT -') +
            f'   ·   Library m/z {entry.get("precursor_mz", "-")}' +
            f'   RT_pred(HILIC) {entry.get("predicted_rt_hilic", "-")}' +
            f'   ΔRT_pred {entry.get("delta_predicted_rt", "-")}s',
            f'entropy_sim: {row["entropy_sim"]:.3f}   ·   sim_gap: {row["sim_gap"]:.3f}   ·   '
            f'delta_mda: {row.get("delta_mda", "?")}   ·   n_close_alts: {row.get("n_close_alternatives", "?")}',
            f'Library SPLASH: {ref_splash}',
        ]
        for i, line in enumerate(lines):
            if i > 0:
                pf = ff.add_paragraph()
            else:
                pf = ff.paragraphs[0]
            rf = pf.add_run()
            rf.text = str(line)
            rf.font.size = Pt(10)
            rf.font.color.rgb = RGBColor(0x33, 0x33, 0x33)

        print(f'  slide #{idx+1}: {row["pipeline_annotation"][:40]} ({row["confidence_pct"]:.0f}%)')

    prs.save(PPTX_OUT)
    print(f'\nSaved {PPTX_OUT}')
    print(f'Mirror-plot PNGs in {OUT_DIR}')


if __name__ == '__main__':
    main()
