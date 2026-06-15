"""build_sibling_conflict_pptx.py — Curator review deck for sibling-conflict triage queue.

Minimal version. For each Tier-1 collision group:
  • Title: compound name + adduct + N bins + RT range
  • Stack of all bins' MS² spectra, labeled by RT only
  • Verdict box at the bottom

No anchors, no per-spectrum confidence labels, no cosine matrix in the slide
(internal cosine still computed for the figure title only as a single number).
"""

import os
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_SHAPE
from pptx.enum.text import PP_ALIGN

ROOT = Path(__file__).resolve().parent.parent
GROUPS = ROOT / 'data' / 'oliver_review_groups_summary.csv'
FLAGS  = ROOT / 'data' / 'bench_joint_inference_flags.csv'
PEAKS  = ROOT / 'data' / 'query_peaks_cache_v2.json'
# Curated CSVs hold the curator-given name with yy_ prefix preserved on raw_name
CUR_NEG = ROOT / 'data' / 'Orbitrap_HILIC_negESI_curated_042126.csv'
CUR_POS = ROOT / 'data' / 'Orbitrap_HILIC_posESI_curated_042126.csv'
OUT    = ROOT / 'reports' / 'sibling_conflict_review.pptx'
PLOT_TMP = ROOT / 'data' / '_sibling_conflict_plots'

PLOT_TMP.mkdir(parents=True, exist_ok=True)
OUT.parent.mkdir(parents=True, exist_ok=True)

SLIDE_W = Inches(13.333)
SLIDE_H = Inches(7.5)
TITLE_H = Inches(0.6)

C_HEADER_BG = RGBColor(0x2C, 0x3E, 0x50)
C_HEADER_FG = RGBColor(0xFF, 0xFF, 0xFF)
C_TIER1_BG  = RGBColor(0xC0, 0x39, 0x2B)
C_TIER2_BG  = RGBColor(0xE6, 0x7E, 0x22)


def normalize_peaks(peaks):
    if peaks is None: return np.array([]).reshape(0, 2)
    if hasattr(peaks, '__len__') and len(peaks) == 0:
        return np.array([]).reshape(0, 2)
    p = np.asarray(peaks, dtype=float)
    if p.size == 0 or p.ndim == 1: return np.array([]).reshape(0, 2)
    p = p[p[:, 0].argsort()]
    mx = p[:, 1].max()
    if mx > 0:
        p[:, 1] /= mx
    return p


def peak_cosine(p_a, p_b, tol=0.01):
    a = normalize_peaks(p_a)
    b = normalize_peaks(p_b)
    if len(a) == 0 or len(b) == 0: return 0.0
    dot = 0.0
    used = np.zeros(len(b), dtype=bool)
    for mz_a, int_a in a:
        diffs = np.abs(b[:, 0] - mz_a)
        diffs[used] = np.inf
        j = int(diffs.argmin())
        if diffs[j] <= tol:
            dot += int_a * b[j, 1]
            used[j] = True
    norm_a = np.sqrt((a[:, 1] ** 2).sum())
    norm_b = np.sqrt((b[:, 1] ** 2).sum())
    if norm_a == 0 or norm_b == 0: return 0.0
    return float(dot / (norm_a * norm_b))


def avg_pairwise_cosine(rows, peaks_map):
    """Mean pairwise cosine across all bins in the group."""
    peaks = [peaks_map.get(r['wiki_id'], []) for _, r in rows.iterrows()]
    n = len(peaks)
    if n < 2: return float('nan')
    vals = []
    for i in range(n):
        for j in range(i + 1, n):
            vals.append(peak_cosine(peaks[i], peaks[j]))
    return float(np.mean(vals)) if vals else float('nan')


# tab10 palette as hex for matching python-pptx text color to plot color
TAB10_HEX = [
    '1F77B4', 'FF7F0E', '2CA02C', 'D62728', '9467BD',
    '8C564B', 'E377C2', '7F7F7F', 'BCBD22', '17BECF',
]


def make_spectra_panel(rows, peaks_map, out_path):
    """Stack N spectra vertically, label each only with its RT + confidence.

    Returns (avg_pairwise_cosine, rows_sorted_by_rt) — the sort order is
    needed to align colors between the figure and the slide's SPLASH list.
    """
    rows = rows.sort_values('rt').reset_index(drop=True)
    n = len(rows)
    # Compact figure: ~0.6 in per spectrum, capped 2.5..4.5 in total
    h = min(4.5, max(2.5, 0.6 * n))
    fig, axes = plt.subplots(n, 1, figsize=(11, h),
                              sharex=True, squeeze=False, dpi=120)
    axes = axes.flatten()

    all_mzs = []
    for _, r in rows.iterrows():
        p = normalize_peaks(peaks_map.get(r['wiki_id'], []))
        if len(p): all_mzs.extend(p[:, 0])
    xmax = (max(all_mzs) * 1.05) if all_mzs else 200

    colors = plt.get_cmap('tab10')
    for i, (ax, (_, r)) in enumerate(zip(axes, rows.iterrows())):
        p = normalize_peaks(peaks_map.get(r['wiki_id'], []))
        c = colors(i % 10)
        if len(p):
            ax.vlines(p[:, 0], 0, p[:, 1], color=c, lw=1.5)
        ax.set_xlim(0, xmax)
        ax.set_ylim(0, 1.10)
        ax.set_yticks([0, 1])
        ax.set_yticklabels(['0', '1'], fontsize=8)
        # No per-subplot ylabel — would overlap when subplots are compact.
        # RT + confidence label inside figure. SPLASH lives in the slide text box below.
        conf = r['gbm_cal'] if pd.notna(r['gbm_cal']) else float('nan')
        ax.text(0.99, 0.85,
                f'RT = {r["rt"]:.1f}s   conf {conf:.3f}',
                transform=ax.transAxes, ha='right', va='top',
                fontsize=12, fontweight='bold', color=c)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    axes[-1].set_xlabel('m/z', fontsize=10)
    # Single shared y-axis label on the middle-ish subplot
    axes[n // 2].set_ylabel('relative intensity', fontsize=9)
    avg_cos = avg_pairwise_cosine(rows, peaks_map)
    fig.suptitle(f'MS² query spectra for the {n} bins  ·  avg pairwise cosine = {avg_cos:.3f}',
                 fontsize=11, y=0.99)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out_path, dpi=140, bbox_inches='tight')
    plt.close(fig)
    return avg_cos, rows


def add_splash_list(slide, rows_rt_sorted, top, height):
    """Add a copy-pastable textbox listing each bin's SPLASH, color-matched to its spectrum.

    `rows_rt_sorted` must be in the same order as the spectra panel (sorted by RT
    ascending), so row i in the textbox corresponds to color TAB10_HEX[i].
    """
    n = len(rows_rt_sorted)
    tb = slide.shapes.add_textbox(Inches(0.4), top, Inches(12.5), height)
    tf = tb.text_frame
    tf.word_wrap = False
    tf.margin_left = Pt(8); tf.margin_right = Pt(8)
    tf.margin_top = Pt(4); tf.margin_bottom = Pt(4)
    for i, (_, r) in enumerate(rows_rt_sorted.iterrows()):
        splash = str(r['splash']) if pd.notna(r['splash']) else 'splash unknown'
        rt = r['rt']
        conf = r['gbm_cal'] if pd.notna(r['gbm_cal']) else float('nan')
        # Build one paragraph per bin, colored to match its plot color
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        # Prefix: "RT 11.4s  conf 0.618    "
        p.text = ''
        run_meta = p.add_run()
        run_meta.text = f'RT {rt:>6.1f}s   conf {conf:.3f}    '
        run_meta.font.size = Pt(10)
        run_meta.font.bold = True
        run_meta.font.name = 'Consolas'
        run_meta.font.color.rgb = RGBColor.from_string(TAB10_HEX[i % 10])
        # SPLASH: monospace, copy-pastable
        run_splash = p.add_run()
        run_splash.text = splash
        run_splash.font.size = Pt(10)
        run_splash.font.name = 'Consolas'
        run_splash.font.color.rgb = RGBColor(0x20, 0x20, 0x20)


def add_blank_slide(prs):
    return prs.slides.add_slide(prs.slide_layouts[6])


def add_title_bar(slide, title_text, subtitle_text=None, tier_color=C_HEADER_BG):
    bar = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, 0, 0, SLIDE_W, TITLE_H)
    bar.fill.solid()
    bar.fill.fore_color.rgb = tier_color
    bar.line.fill.background()
    tf = bar.text_frame
    tf.margin_left = Pt(14); tf.margin_right = Pt(14)
    tf.margin_top = Pt(5); tf.margin_bottom = Pt(2)
    p = tf.paragraphs[0]; p.text = title_text
    r = p.runs[0]
    r.font.size = Pt(20); r.font.bold = True
    r.font.color.rgb = C_HEADER_FG
    if subtitle_text:
        p2 = tf.add_paragraph(); p2.text = subtitle_text
        r2 = p2.runs[0]
        r2.font.size = Pt(11); r2.font.color.rgb = C_HEADER_FG


# Front matter

def slide_cover(prs, n_primary, n_appendix):
    s = add_blank_slide(prs)
    bg = s.shapes.add_shape(MSO_SHAPE.RECTANGLE, 0, 0, SLIDE_W, SLIDE_H)
    bg.fill.solid(); bg.fill.fore_color.rgb = C_HEADER_BG
    bg.line.fill.background()
    tb = s.shapes.add_textbox(Inches(0.7), Inches(1.7), Inches(12), Inches(1.5))
    tf = tb.text_frame; tf.word_wrap = True
    p = tf.paragraphs[0]; p.text = 'Sibling-conflict review'
    p.runs[0].font.size = Pt(44); p.runs[0].font.bold = True
    p.runs[0].font.color.rgb = C_HEADER_FG
    p2 = tf.add_paragraph()
    p2.text = 'Bins labeled as the same compound but at chemically-incompatible retention times'
    p2.runs[0].font.size = Pt(18)
    p2.runs[0].font.color.rgb = RGBColor(0xBD, 0xC3, 0xC7)
    tb2 = s.shapes.add_textbox(Inches(0.7), Inches(4.0), Inches(12), Inches(3))
    tf2 = tb2.text_frame; tf2.word_wrap = True
    lines = [
        f'• {n_primary} groups in PRIMARY review queue (no prior yy_ rejection)',
        f'• {n_appendix} groups in APPENDIX — already include a yy_-rejected bin, likely confirms your prior calls',
        '• Each primary slide shows the bins\' MS² spectra and asks one question',
        '• Date: 2026-05-29',
    ]
    for ln in lines:
        pp = tf2.add_paragraph() if tf2.paragraphs[0].text else tf2.paragraphs[0]
        pp.text = ln
        pp.runs[0].font.size = Pt(16); pp.runs[0].font.color.rgb = C_HEADER_FG


def slide_method(prs):
    s = add_blank_slide(prs)
    add_title_bar(s, 'How we found these')
    tb = s.shapes.add_textbox(Inches(0.7), Inches(1.0), Inches(12), Inches(5))
    tf = tb.text_frame; tf.word_wrap = True

    def add(text, bold=False, size=12, color=None):
        p = tf.paragraphs[0] if not tf.paragraphs[0].text else tf.add_paragraph()
        p.text = text
        if p.runs:
            p.runs[0].font.size = Pt(size); p.runs[0].font.bold = bold
            if color: p.runs[0].font.color.rgb = color
    add('The setup', bold=True, size=15, color=C_HEADER_BG)
    add('Under fixed HILIC column and method, a compound elutes at one retention time. '
        'If multiple bins are curator-labeled as the same compound + adduct but at very '
        'different RTs, something is inconsistent — either the labels or the chromatography.')
    add('')
    add('The filter', bold=True, size=15, color=C_HEADER_BG)
    add('• Group bins by (compound, adduct) where the same identity is the top-1 call across ≥2 bins')
    add('• Restrict to RT range > 30s within a group (chemistry-incompatible under fixed method)')
    add('• Restrict to groups where at least one bin is curator-annotated AND model-confident (≥ 0.5)')
    add('')
    add('Deck structure', bold=True, size=15, color=C_HEADER_BG)
    add('• Primary detail slides: groups with NO prior yy_-rejected bins — potentially new findings to decide on')
    add('• Appendix table: groups that already include a yy_-rejected bin — model is re-discovering peaks you\'ve already flagged; skim and let us know if any look new')
    add('')
    add('What each detail slide shows', bold=True, size=15, color=C_HEADER_BG)
    add('• Compound name + adduct + number of bins + RT range')
    add('• Each bin\'s MS² query spectrum, stacked, labeled by RT + confidence')
    add('• Copy-pastable SPLASH for each bin below the figure')
    add('• Average pairwise MS² cosine across the bins (one number)')


def slide_action(prs):
    s = add_blank_slide(prs)
    add_title_bar(s, 'Action requested')
    tb = s.shapes.add_textbox(Inches(0.7), Inches(1.0), Inches(12), Inches(5.5))
    tf = tb.text_frame; tf.word_wrap = True

    def add(text, bold=False, size=12, color=None):
        p = tf.paragraphs[0] if not tf.paragraphs[0].text else tf.add_paragraph()
        p.text = text
        if p.runs:
            p.runs[0].font.size = Pt(size); p.runs[0].font.bold = bold
            if color: p.runs[0].font.color.rgb = color
    add('What we need from you', bold=True, size=16, color=C_HEADER_BG)
    add('• Walk through the Tier-1 slides and let us know which groups are real conflicts')
    add('• If a Tier-2 group looks interesting, flag it and we\'ll add a detail slide')
    add('')
    add('How to read the spectra', bold=True, size=16, color=C_HEADER_BG)
    add('• If all spectra look essentially the same → they really are the same compound; the RT discrepancy is the puzzle (curator label issue or method drift)')
    add('• If spectra differ → at least one of the labels is wrong; the model is over-extending one good call to a different peak')
    add('• The average pairwise cosine (single number) is a fast quantitative check: ≥ 0.7 means spectra match; < 0.3 means they don\'t')


def slide_appendix_summary(prs, appendix_groups, n_yy_bins_by_group):
    """Single slide(s) listing the yy_-containing groups. May split across slides
    if there are too many to fit comfortably (max ~30 rows/slide)."""
    if not len(appendix_groups):
        return
    ag = appendix_groups.sort_values(['severity', 'n_bins'],
                                      ascending=[False, False]).reset_index(drop=True)
    ag['n_yy'] = ag['collision_group_id'].map(n_yy_bins_by_group).fillna(0).astype(int)
    rows_per_slide = 28
    n_slides = (len(ag) + rows_per_slide - 1) // rows_per_slide
    for k in range(n_slides):
        s = add_blank_slide(prs)
        title = (f'Appendix — yy_-containing groups   ({k+1} / {n_slides})  ·  '
                 f'{len(ag)} groups total')
        subtitle = 'These already include a curator-rejected (yy_) bin. Most likely confirm your prior calls — skim and flag anything that surprises you.'
        add_title_bar(s, title, subtitle)
        chunk = ag.iloc[k*rows_per_slide:(k+1)*rows_per_slide]
        cols = ['#', 'Compound (top-1 name)', 'Adduct', 'n_bins', 'n_yy_bins', 'RT range (s)', 'max_conf', 'Polarity']
        n_rows = len(chunk) + 1
        tbl = s.shapes.add_table(n_rows, len(cols),
                                  Inches(0.3), Inches(0.95),
                                  Inches(12.7), Inches(6.3)).table
        cw = [Inches(0.4), Inches(5.7), Inches(1.3), Inches(0.7),
              Inches(0.9), Inches(1.4), Inches(0.9), Inches(1.0)]
        for i, w in enumerate(cw): tbl.columns[i].width = w
        for i, c in enumerate(cols):
            cell = tbl.cell(0, i)
            cell.fill.solid(); cell.fill.fore_color.rgb = C_HEADER_BG
            p = cell.text_frame.paragraphs[0]; p.text = c
            p.runs[0].font.size = Pt(9); p.runs[0].font.bold = True
            p.runs[0].font.color.rgb = C_HEADER_FG
        for i, (_, r) in enumerate(chunk.iterrows(), start=1):
            cells = [str(k*rows_per_slide + i),
                     str(r['name'])[:60], str(r['adduct']),
                     str(int(r['n_bins'])), str(int(r['n_yy'])),
                     f'{r["rt_range"]:.1f}', f'{r["max_conf"]:.3f}',
                     str(r['polarities'])]
            for j, val in enumerate(cells):
                p = tbl.cell(i, j).text_frame.paragraphs[0]; p.text = val
                p.runs[0].font.size = Pt(8.5)


def main():
    print(f'Loading {GROUPS}')
    groups = pd.read_csv(GROUPS)
    flags = pd.read_csv(FLAGS)
    qualifying = set(zip(groups['hit_ik14'], groups['adduct']))
    flags['_grp_key'] = list(zip(flags['hit_ik14'], flags['adduct']))
    rows = flags[flags['_grp_key'].isin(qualifying)].copy()
    rows['collision_group_id'] = (rows['hit_ik14'].fillna('') + '|' + rows['adduct'].fillna(''))
    print(f'  {len(groups)} groups, {len(rows)} bins after expanding to all members')

    # Join curator-given names (raw_name preserves yy_ prefix) so we can label
    # each group yy_-containing or yy-free.
    cur_neg = pd.read_csv(CUR_NEG, low_memory=False,
                          usecols=['wiki_id', 'name', 'raw_name'])
    cur_pos = pd.read_csv(CUR_POS, low_memory=False,
                          usecols=['wiki_id', 'name', 'raw_name'])
    cur = pd.concat([cur_neg, cur_pos], ignore_index=True).drop_duplicates(subset='wiki_id')
    cur = cur.rename(columns={'name': 'curator_name', 'raw_name': 'curator_raw_name'})
    rows = rows.merge(cur, on='wiki_id', how='left')
    rows['is_yy'] = (
        rows['curator_name'].fillna('').str.match(r'^\s*yy_', case=False) |
        rows['curator_raw_name'].fillna('').str.match(r'^\s*yy_', case=False)
    )
    yy_per_group = rows.groupby('collision_group_id')['is_yy'].sum()

    print(f'Loading peaks from {PEAKS}')
    peaks_map = json.load(open(PEAKS))
    print(f'  {len(peaks_map):,} bins cached')

    # Partition: primary = no yy_ in any bin, appendix = ≥1 yy_ bin
    groups['n_yy_bins'] = groups['collision_group_id'].map(yy_per_group).fillna(0).astype(int)
    primary = groups[groups['n_yy_bins'] == 0].sort_values(
        ['severity', 'n_bins'], ascending=[False, False]).reset_index(drop=True)
    appendix = groups[groups['n_yy_bins'] >= 1].sort_values(
        ['severity', 'n_bins'], ascending=[False, False]).reset_index(drop=True)
    print(f'  Primary (no yy_):   {len(primary)} groups → detail slides')
    print(f'  Appendix (has yy_): {len(appendix)} groups → summary table')

    prs = Presentation()
    prs.slide_width = SLIDE_W
    prs.slide_height = SLIDE_H

    slide_cover(prs, n_primary=len(primary), n_appendix=len(appendix))
    slide_method(prs)

    print(f'\nBuilding {len(primary)} primary detail slides...')
    for i, (_, gr) in enumerate(primary.iterrows(), 1):
        gid = gr['collision_group_id']
        members = rows[rows['collision_group_id'] == gid].copy()
        if len(members) < 2:
            print(f'  [{i}/{len(primary)}] SKIP {gid[:30]}…  only {len(members)} member')
            continue
        print(f'  [{i}/{len(primary)}] {gr["name"][:50]} | {gr["adduct"]} | n={len(members)}')

        plot_path = PLOT_TMP / f'primary_{i:03d}.png'
        avg_cos, rt_sorted = make_spectra_panel(members, peaks_map, plot_path)

        # Recompute count + RT range from the bins actually plotted (broader
        # collision, including any unannotated / low-confidence members shown
        # for context). The group-summary values reflect only the curator-
        # annotated subset and would mismatch the figure otherwise.
        plotted_n = len(members)
        plotted_rt_range = members['rt'].max() - members['rt'].min()

        s = add_blank_slide(prs)
        tier_color = C_TIER1_BG if gr['severity'] >= 330 else C_TIER2_BG
        title = f'{gr["name"]}  ·  {gr["adduct"]}'
        subtitle = (f'{plotted_n} bins where the model picked this compound  ·  '
                    f'RT range {plotted_rt_range:.1f} s  ·  polarity {gr["polarities"]}')
        add_title_bar(s, title, subtitle, tier_color=tier_color)
        pic = s.shapes.add_picture(str(plot_path),
                                    Inches(0.4), Inches(0.85),
                                    width=Inches(12.5))
        # Place SPLASH textbox right below the picture; pic.height auto-computed
        splash_top = Inches(0.85) + pic.height + Inches(0.1)
        # Cap the bottom at slide height
        remaining = SLIDE_H - splash_top - Inches(0.1)
        if remaining > Inches(0.4):
            add_splash_list(s, rt_sorted, splash_top, remaining)

    slide_appendix_summary(prs, appendix, yy_per_group.to_dict())
    slide_action(prs)

    prs.save(OUT)
    print(f'\n✓ Saved {OUT}')
    print(f'  Total slides: {len(prs.slides)}')


if __name__ == '__main__':
    main()
