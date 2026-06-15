"""build_reverse_score_deck.py — Methodology deck: how reverse score was probed (multiple ways)
and how the 'everything-not-in-library = noise' assumption was tested.

Numbers are the validated results from the 2026-06-09..12 reverse-score probes
(see reports/reverse_score_revisited.md and project_reverse_search_isf_probe_20260609 memory).
Two figures are generated from those numbers; the rest are text/table slides.
"""
import os
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / 'reports' / 'reverse_score_probes.pptx'
TMP = ROOT / 'data' / '_revscore_deck_plots'
TMP.mkdir(parents=True, exist_ok=True)

SLIDE_W, SLIDE_H = Inches(13.333), Inches(7.5)
NAVY   = RGBColor(0x1F, 0x2D, 0x3D)
BLUE   = RGBColor(0x2C, 0x3E, 0x50)
TEAL   = RGBColor(0x16, 0xA0, 0x85)
RED    = RGBColor(0xC0, 0x39, 0x2B)
ORANGE = RGBColor(0xE6, 0x7E, 0x22)
GREY   = RGBColor(0x5D, 0x6D, 0x7E)
LGREY  = RGBColor(0xEC, 0xF0, 0xF1)
WHITE  = RGBColor(0xFF, 0xFF, 0xFF)
DARK   = RGBColor(0x2B, 0x2B, 0x2B)

# ---------- figures ----------
def fig_injection():
    frac = np.array([0, 25, 50, 75])
    pred = np.array([0.769, 0.577, 0.384, 0.192])
    obs  = np.array([0.769, 0.567, 0.382, 0.193])
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    ax.plot(frac, pred, '--', color=GREY_hex(), lw=2, label='predicted  baseline×(1−frac)')
    ax.plot(frac, obs, 'o-', color='#16A085', lw=2.5, ms=10, label='observed purity (median)')
    for x, y in zip(frac, obs):
        ax.annotate(f'{y:.3f}', (x, y), textcoords='offset points', xytext=(8, 8), fontsize=11)
    ax.set_xlabel('injected contaminant intensity (%)', fontsize=12)
    ax.set_ylabel('query-side purity', fontsize=12)
    ax.set_title('Injection test: purity is a calibrated contamination gauge', fontsize=12.5, weight='bold')
    ax.set_ylim(0, 0.85); ax.legend(fontsize=11, frameon=False); ax.grid(alpha=0.25)
    fig.tight_layout(); p = TMP / 'injection.png'; fig.savefig(p, dpi=200); plt.close(fig)
    return str(p)

def fig_assumption():
    labels = ['single best\nentry', 'UNION same-\ncompound (IK14)', 'UNION random\n(control)']
    vals = [0.769, 0.854, 0.023]
    cols = ['#5D6D7E', '#16A085', '#C0392B']
    fig, ax = plt.subplots(figsize=(6.6, 4.4))
    bars = ax.bar(labels, vals, color=cols, width=0.6)
    for b, v in zip(bars, vals):
        ax.text(b.get_x()+b.get_width()/2, v+0.015, f'{v:.3f}', ha='center', fontsize=12, weight='bold')
    ax.annotate('', xy=(1, 0.86), xytext=(0, 0.775),
                arrowprops=dict(arrowstyle='->', color='#16A085', lw=2))
    ax.text(0.30, 0.96, '+8.5 pp real fragments', color='#16A085', fontsize=11.5, ha='center', weight='bold')
    ax.set_ylabel('median query-intensity explained', fontsize=12)
    ax.set_title('Is "not in the library" really contamination?', fontsize=12.5, weight='bold')
    ax.set_ylim(0, 1.08); ax.grid(axis='y', alpha=0.25)
    fig.tight_layout(); p = TMP / 'assumption.png'; fig.savefig(p, dpi=200); plt.close(fig)
    return str(p)

def fig_isf_controls():
    labels = ['real\n(neutral losses)', 'RT-shuffled\n(break co-elution)', 'random Δm/z\n(not a loss)']
    vals = [14.9, 7.6, 0.0]
    cols = ['#16A085', '#E67E22', '#C0392B']
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    bars = ax.bar(labels, vals, color=cols, width=0.6)
    for b, v in zip(bars, vals):
        ax.text(b.get_x()+b.get_width()/2, v+0.4, f'{v:.1f}%', ha='center', fontsize=13, weight='bold')
    ax.set_ylabel('% of co-eluting Δm/z pairs with containment ≥ 0.5', fontsize=11)
    ax.set_title('Relational ISF: controls (random Δm/z → 0%)', fontsize=12.5, weight='bold')
    ax.set_ylim(0, 17); ax.grid(axis='y', alpha=0.25)
    fig.tight_layout(); p = TMP / 'isf_controls.png'; fig.savefig(p, dpi=200); plt.close(fig)
    return str(p)

def GREY_hex(): return '#5D6D7E'

# ---------- slide helpers ----------
def _set_bg(slide, color):
    slide.background.fill.solid(); slide.background.fill.fore_color.rgb = color

def _box(slide, l, t, w, h, fill=None, line=None):
    sp = slide.shapes.add_shape(1, l, t, w, h)
    sp.fill.solid(); sp.fill.fore_color.rgb = fill if fill else WHITE
    if line is None: sp.line.fill.background()
    else: sp.line.color.rgb = line
    sp.shadow.inherit = False
    return sp

def _txt(slide, l, t, w, h, text, size=18, color=DARK, bold=False, align=PP_ALIGN.LEFT, font='Calibri', anchor=MSO_ANCHOR.TOP):
    tb = slide.shapes.add_textbox(l, t, w, h); tf = tb.text_frame
    tf.word_wrap = True; tf.vertical_anchor = anchor
    lines = text.split('\n')
    for i, ln in enumerate(lines):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.alignment = align
        r = p.add_run(); r.text = ln
        r.font.size = Pt(size); r.font.bold = bold; r.font.color.rgb = color; r.font.name = font
    return tb

def header(slide, kicker, title):
    _box(slide, 0, 0, SLIDE_W, Inches(1.15), fill=BLUE)
    if kicker:
        _txt(slide, Inches(0.6), Inches(0.12), Inches(12), Inches(0.35), kicker.upper(),
             size=13, color=RGBColor(0x9F,0xB6,0xC9), bold=True)
    _txt(slide, Inches(0.6), Inches(0.42), Inches(12.1), Inches(0.62), title, size=26, color=WHITE, bold=True)

def bullets(slide, items, top=Inches(1.5), left=Inches(0.7), width=Inches(12), size=18, gap=0.06):
    tb = slide.shapes.add_textbox(left, top, width, Inches(5.4)); tf = tb.text_frame; tf.word_wrap = True
    for i, it in enumerate(items):
        if isinstance(it, tuple): text, lvl, color = it
        else: text, lvl, color = it, 0, DARK
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.level = lvl; p.space_after = Pt(gap*72)
        r = p.add_run(); r.text = ('•  ' if lvl == 0 else '–  ') + text
        r.font.size = Pt(size - 2*lvl); r.font.color.rgb = color; r.font.name = 'Calibri'
    return tb

def table(slide, data, left, top, width, col_w=None, header_fill=BLUE, fsize=14, row_h=Inches(0.42)):
    rows, cols = len(data), len(data[0])
    gt = slide.shapes.add_table(rows, cols, left, top, width, row_h*rows).table
    if col_w:
        for j, w in enumerate(col_w): gt.columns[j].width = w
    for i in range(rows):
        for j in range(cols):
            c = gt.cell(i, j); c.text = str(data[i][j])
            para = c.text_frame.paragraphs[0]; para.alignment = PP_ALIGN.LEFT if j == 0 else PP_ALIGN.CENTER
            run = para.runs[0]; run.font.size = Pt(fsize); run.font.name = 'Calibri'
            c.margin_top = Pt(2); c.margin_bottom = Pt(2)
            if i == 0:
                c.fill.solid(); c.fill.fore_color.rgb = header_fill
                run.font.color.rgb = WHITE; run.font.bold = True
            else:
                c.fill.solid(); c.fill.fore_color.rgb = WHITE if i % 2 else LGREY
                run.font.color.rgb = DARK
    return gt

def verdict_chip(slide, text, color, left, top, w=Inches(2.0)):
    b = _box(slide, left, top, w, Inches(0.5), fill=color)
    _txt(slide, left, top+Inches(0.04), w, Inches(0.42), text, size=14, color=WHITE, bold=True, align=PP_ALIGN.CENTER, anchor=MSO_ANCHOR.MIDDLE)

# ---------- build ----------
prs = Presentation(); prs.slide_width = SLIDE_W; prs.slide_height = SLIDE_H
blank = prs.slide_layouts[6]
def newslide(bg=WHITE):
    s = prs.slides.add_slide(blank); _set_bg(s, bg); return s

inj = fig_injection(); asm = fig_assumption(); isfc = fig_isf_controls()

# 1 — title
s = newslide(NAVY)
_txt(s, Inches(0.9), Inches(2.2), Inches(11.5), Inches(1.2), 'Reverse Score, Probed Multiple Ways',
     size=44, color=WHITE, bold=True)
_txt(s, Inches(0.95), Inches(3.5), Inches(11.5), Inches(0.7),
     'Testing reverse similarity as identity vs. noise — and the "everything else is contamination" assumption',
     size=20, color=RGBColor(0x9F,0xB6,0xC9))
_txt(s, Inches(0.95), Inches(6.3), Inches(11), Inches(0.5),
     'Ziyue Yang · West Coast Metabolomics Center · 2026-06-12', size=15, color=RGBColor(0x7F,0x8C,0x9D))

# 2 — the setup
s = newslide(); header(s, 'Why revisit this', 'Two prompts, one old idea, one strong assumption')
bullets(s, [
    ('Oliver (deployed-score review): "we really need reverse-similarity scoring that ignores noise and only looks at which library fragments actually match."', 0, DARK),
    ('GNPS poster (Xing et al.): enhanced reverse cosine recovers +37% TPs at 90% purity, +252% at 50% — by NOT penalizing query peaks the reference lacks.', 0, DARK),
    ('NIST, verbatim: "in a reverse search, nonmatching unknown peaks are assumed to be contaminants."', 0, GREY),
    ('That last line is a STRONG assumption. This deck: I probed reverse score 5 ways, then tested the assumption directly.', 0, RED),
], top=Inches(1.6), size=18)

# 3 — what reverse score is
s = newslide(); header(s, 'Framing', 'Reverse score = asymmetric / containment similarity')
_txt(s, Inches(0.7), Inches(1.4), Inches(12), Inches(0.6),
     'It asks "is the reference explained by the query?" — not "are the two identical?"  In precision/recall terms:', size=17)
table(s, [
    ['metric', 'measures', 'side'],
    ['forward_cosine', 'query AND reference peaks both penalize', 'both'],
    ['reverse_cosine', 'anchored on reference; query-only peaks dropped', 'reference (≈ recall)'],
    ['cov_count / cov_int', 'matched-library fraction', 'reference'],
    ['query-side PURITY', 'matched-query / total-query intensity', 'QUERY (≈ precision)'],
], Inches(0.7), Inches(2.2), Inches(11.9),
   col_w=[Inches(3.2), Inches(5.7), Inches(3.0)], fsize=14)
_txt(s, Inches(0.7), Inches(5.3), Inches(12), Inches(1.1),
     'Every metric we ship is reference-side ("is the library entry present?"). The one quantity none of them captures is the\nquery side — how much of the OBSERVED spectrum the reference explains. That discarded quantity IS the contamination.',
     size=16, color=BLUE, bold=True)

# 4 — probe map
s = newslide(); header(s, 'Roadmap', 'I probed reverse score six ways')
table(s, [
    ['#', 'probe', 'objective', 'verdict'],
    ['1', 'reverse_cosine as a GBM feature', 'identity', 'redundant'],
    ['2', 'forward−reverse GAP as ISF detector', 'identity', 'null (chance)'],
    ['3', 'reverse-only ranking after an RT gate', 'identity', 'precision collapses'],
    ['4', 'query-side PURITY meter (+ injection)', 'noise', 'valid but ≈ the gap'],
    ['5', 'the assumption: single vs union vs random', 'assumption', 'wrong, bounded'],
    ['6', 'RELATIONAL ISF: co-elute + Δm/z + containment', 'detection', 'WORKS (experiment)'],
], Inches(0.7), Inches(1.5), Inches(11.9),
   col_w=[Inches(0.6), Inches(5.2), Inches(2.6), Inches(3.5)], fsize=14, row_h=Inches(0.5))
_txt(s, Inches(0.7), Inches(5.9), Inches(12), Inches(0.9),
     'Probes 1–3: can reverse identify? (no).   4: can it measure noise? (yes, but not new).   5: is the core assumption true? (no).\n6: point containment at the RIGHT pair — fragment ⊂ parent — and it works.',
     size=14, color=GREY)

# 5 — probe 1
s = newslide(); header(s, 'Probe 1 — identity', 'Reverse cosine as a GBM feature')
bullets(s, [
    ('Added reverse_cosine (and forward_cosine, cov_count, cov_int) to the feature table.', 0, DARK),
    ('Marginal AUC over the shipped feature set: ΔAUC ≈ 0.', 0, DARK),
    ('entropy_similarity already dominates (AUC 0.90+); reverse is strictly more lenient, so it adds no separating power.', 0, GREY),
], top=Inches(1.6), size=19)
verdict_chip(s, 'REDUNDANT', GREY, Inches(0.7), Inches(5.6))

# 6 — probe 2
s = newslide(); header(s, 'Probe 2 — identity', 'Forward−reverse GAP as an ISF / chimera detector')
bullets(s, [
    ('Hypothesis: reverse ≫ forward ⇒ query carries extra peaks ⇒ in-source fragment / contamination.', 0, DARK),
    ('Detects adduct-notation ISF weakly (AUC 0.62) …', 0, DARK),
    ('… but separating TP from FP WITHIN the ISF band: AUC 0.52–0.54 = coin flip.', 0, RED),
    ('entropy_similarity still holds 0.84 in that same band.', 0, GREY),
], top=Inches(1.6), size=19)
verdict_chip(s, 'NULL (CHANCE)', RED, Inches(0.7), Inches(5.7))

# 7 — probe 3
s = newslide(); header(s, 'Probe 3 — identity', 'Reverse-only ranking after an RT gate')
_txt(s, Inches(0.7), Inches(1.35), Inches(12), Inches(0.5),
     'Gate on RT-plausibility, then rank by reverse score. It gets WORSE as the gate tightens:', size=17)
table(s, [
    ['RT gate', 'reverse_cosine AUC'],
    ['all rows with peaks', '0.804'],
    ['|ΔRT| ≤ 30 s', '0.736'],
    ['|ΔRT| ≤ 15 s', '0.716'],
    ['|ΔRT| ≤ 5 s', '0.677'],
], Inches(0.7), Inches(2.0), Inches(6.2), col_w=[Inches(3.2), Inches(3.0)], fsize=15)
_txt(s, Inches(7.4), Inches(2.0), Inches(5.4), Inches(3),
     'The "rescue zone"\n(RT≤15 s, esim<0.7, reverse≥0.7):\n\n• recovers 318 / 390 missed TPs\n• but drags in 5,526 FPs\n• precision 5.4%  (base rate 16%)',
     size=17, color=DARK)
verdict_chip(s, 'PRECISION COLLAPSES', RED, Inches(0.7), Inches(6.0), w=Inches(3.2))

# 8 — why identity fails
s = newslide(BLUE)
_txt(s, Inches(0.8), Inches(0.7), Inches(11.7), Inches(1.0),
     'Why all three identity probes fail — one structural reason', size=30, color=WHITE, bold=True)
bullets(s, [
    ('The poster\'s reverse helps when the extra peaks are random co-isolated UNRELATED compounds.', 0, WHITE),
    ('In OUR confusable band the extra peaks are the true compound\'s own ISF fragments, or a co-eluting isomer.', 0, WHITE),
    ('So the reference is present in BOTH the TP and the FP — reverse fires identically for both.', 0, RGBColor(0xF5,0xB7,0xB1)),
    ('The real question is not "is the reference present" but "which compound owns these peaks" —', 0, WHITE),
    ('which only RT + the curator\'s cross-bin logic answers, i.e. the labeling rule itself → circular.', 1, RGBColor(0xF5,0xB7,0xB1)),
], top=Inches(2.0), size=19)

# 9 — reframe
s = newslide(); header(s, 'The pivot', 'Reverse score for NOISE, not identity')
bullets(s, [
    ('The property that kills reverse as an identity signal — it ignores query-only peaks — is exactly what a noise meter wants.', 0, DARK),
    ('Because the query-only peaks ARE the noise.', 1, TEAL),
    ('So flip the objective: measure query-side PURITY = fraction of observed intensity the reference explains.', 0, DARK),
    ('Decomposition that falls out:  forward_cosine = identity × cleanliness.', 0, BLUE),
    ('Forward conflates "right compound?" with "clean spectrum?"; purity isolates the cleanliness factor.', 1, GREY),
], top=Inches(1.55), size=18)

# 10 — noise test 1 (injection)
s = newslide(); header(s, 'Probe 4 — noise · construct validity', 'Inject known contamination → does purity track it?')
s.shapes.add_picture(inj, Inches(0.6), Inches(1.4), height=Inches(4.7))
_txt(s, Inches(7.4), Inches(1.7), Inches(5.4), Inches(4),
     'Add X% contaminant intensity to clean\nspectra; purity should fall as (1−X).\n\nObserved tracks predicted almost\nexactly (0.567 / 0.382 / 0.193).\n\n→ purity is a CALIBRATED\n   contamination gauge.',
     size=18, color=DARK)
verdict_chip(s, 'CONSTRUCT-VALID', TEAL, Inches(7.4), Inches(5.7), w=Inches(2.7))

# 11 — noise test 1b orthogonality + honesty
s = newslide(); header(s, 'Probe 4 — noise · is it new?', 'Orthogonal to identity — but not a new number')
table(s, [
    ['Spearman ρ of purity with …', 'value', 'reading'],
    ['entropy_similarity', '0.43', 'orthogonal to identity ✓'],
    ['reverse_cosine', '0.11', 'reverse alone is BLIND to noise'],
    ['forward_cosine', '0.78', 'forward partly sees it (confounded)'],
    ['L2 forward/reverse GAP (probe 2)', '0.95', 'same construct, just L1 vs L2'],
], Inches(0.7), Inches(1.6), Inches(11.9),
   col_w=[Inches(5.0), Inches(1.6), Inches(5.3)], fsize=15, row_h=Inches(0.5))
_txt(s, Inches(0.7), Inches(5.4), Inches(12), Inches(1.3),
     'Honest verdict: the scalar contamination axis is REAL and calibrated, but ρ=0.95 with the gap we already tried.\nThe number is not new — only the OBJECTIVE is (spectrum-level QC, not pair-level ID). 13% of well-matched\nspectra are "identifiable-but-contaminated" (esim≥0.7, purity<0.4) — a real population, but the scalar is exhausted.',
     size=15, color=BLUE, bold=True)

# 12 — the assumption test design
s = newslide(); header(s, 'Probe 5 — the assumption', '"Everything not in the library = contamination" — is it?')
bullets(s, [
    ('Single best entry explains 76.9% of median TP query intensity → assumption claims the other 23.1% is ALL noise.', 0, DARK),
    ('Ground truth for "real fragment": a peak that appears in ANOTHER library spectrum of the SAME compound (IK14)', 0, DARK),
    ('cannot be contamination from something else.', 1, GREY),
    ('Test: re-explain the query with (a) single entry, (b) UNION of same-compound entries, (c) UNION of random entries (inflation control).', 0, DARK),
], top=Inches(1.55), size=18)
s.shapes.add_picture(asm, Inches(6.9), Inches(3.0), height=Inches(4.0))

# 13 — the assumption result
s = newslide(); header(s, 'Probe 5 — result', 'Wrong — by a measurable, BOUNDED amount')
bullets(s, [
    ('Same-compound siblings recover +8.5 pp the single reference called "noise" — provably real fragments.', 0, DARK),
    ('Random size-matched references recover only 2.3% → the +8.5 pp is NOT pooling inflation; it is compound-specific.', 0, DARK),
], top=Inches(1.5), size=18)
table(s, [
    ['quantity', 'bound'],
    ['real signal mislabeled as contamination', '[ 8.5 , 23.1 ] pp'],
    ['genuine contamination', '[ 0 , 14.6 ] pp'],
], Inches(0.7), Inches(3.0), Inches(7.4), col_w=[Inches(5.0), Inches(2.4)], fsize=16, row_h=Inches(0.55))
_txt(s, Inches(0.7), Inches(4.9), Inches(12), Inches(1.6),
     'Bounded, not a point: 8.5 pp is a LOWER bound (real fragments deposited in none of our finite references stay\nhidden); 14.6% residual is an UPPER bound on true contamination. Only MS1 co-isolation ground truth (BinBase)\ncollapses the interval to a point.',
     size=15, color=GREY)
verdict_chip(s, 'ASSUMPTION FALSE (BOUNDED)', RED, Inches(8.6), Inches(3.0), w=Inches(4.0))

# 14 — cross-domain
s = newslide(); header(s, 'Sanity check', 'We re-derived a universal pattern')
_txt(s, Inches(0.7), Inches(1.3), Inches(12), Inches(0.6),
     'One-sided / containment similarity is the right tool for "signal + clutter" data — and every field pairs it with a counter-penalty:', size=15)
table(s, [
    ['field', 'their reverse score', 'their counter-penalty'],
    ['Metagenomics', 'containment index (sourmash / Mash Screen)', 'size-aware significance / ANI'],
    ['NLP / IR', 'ROUGE (recall) vs BLEU (precision)', 'brevity penalty; MULTIPLE references'],
    ['Cheminformatics', 'Tversky similarity (RDKit; substructure)', 'β weight on distinctive features'],
    ['Audio (Shazam)', 'landmark-peak fingerprinting', 'peak pairing for distinctiveness'],
], Inches(0.7), Inches(2.0), Inches(11.9),
   col_w=[Inches(2.6), Inches(5.3), Inches(4.0)], fsize=14, row_h=Inches(0.55))
_txt(s, Inches(0.7), Inches(5.6), Inches(12), Inches(1.0),
     'We independently re-derived both halves: reverse-alone floods FPs (over-reward), and the fix is the purity axis +\nreference-sparsity control (counter-penalty). BLEU\'s "multiple references" = our compensation in §next.',
     size=15, color=BLUE, bold=True)

# 15 — correction: what hit_is_isf really is
s = newslide(); header(s, 'Probe 6 — setup', 'First, a correction: our ISF flag is just the adduct STRING')
bullets(s, [
    ('hit_is_isf = (hit_adduct_cat == "isf") — from a 684-entry taxonomy CSV + an "M+H−…" regex.', 0, DARK),
    ('It reads what the ANNOTATION claims the ion is. No RT, no second bin, no m/z arithmetic, no spectrum.', 0, RED),
    ('So Probe 2 was "detecting" the label, semi-circularly. The physical event was never tested.', 0, GREY),
    ('The real, de-novo question (what CAMERA / RAMClust do, and Oliver\'s suggestion #1):', 0, DARK),
    ('fragment bin F is in-source fragment of parent bin P  ⇔  co-elute (|ΔRT|≤4s)  +  Δm/z = a neutral loss  +  F\'s MS² CONTAINED in P\'s', 1, BLUE),
], top=Inches(1.55), size=18)
_txt(s, Inches(0.7), Inches(5.9), Inches(12), Inches(0.8),
     'That last clause is the reverse / containment score — finally pointed at the right pair: fragment ⊂ parent, within-sample (not query ↔ library).',
     size=15, color=TEAL, bold=True)

# 16 — relational ISF result + controls
s = newslide(); header(s, 'Probe 6 — result  (experiment, not production)', 'It works — and survives hard controls')
s.shapes.add_picture(isfc, Inches(0.6), Inches(1.4), height=Inches(4.6))
_txt(s, Inches(7.3), Inches(1.5), Inches(5.6), Inches(4.5),
     'Anchor (2-aminoadipic):\n  containment 0.976 via H₂O ✓\n\n865 flagged → 558 STRONG set\n(+ fragment precursor present\nas a peak in parent: 64.5%;\nmedian 27 peaks → not sparse)\n\nControls (left):\n• random Δm/z → 0% contained\n• co-elution doubles the hits',
     size=17, color=DARK)
verdict_chip(s, 'DETECTOR VALIDATED', TEAL, Inches(7.3), Inches(6.0), w=Inches(3.0))

# 17 — relational ISF validation + payoff
s = newslide(); header(s, 'Probe 6 — why I believe it', 'Validation no string-flag could fake')
bullets(s, [
    ('The detected loss matches the annotated loss, INDEPENDENTLY — levetiracetam [M+H−NH3]+ via NH₃,', 0, DARK),
    ('ephedrine / GABA [M+H−H2O]+ via H₂O — from observed masses + RT + spectra only.', 1, GREY),
    ('It re-discovers curator free-text BLIND — bin named "yy_in source to ethylbenzyl cation" flagged ISF without reading the name.', 0, DARK),
    ('Payoff: 30% of the strong set (167 bins) are yy_ curator-REJECTED → the detector gives the MECHANISM for the rejection.', 0, TEAL),
    ('Also surfaces candidate mis-annotations (proline [M+H]+, 9 peaks, fully inside a co-eluting +H₂O parent).', 0, DARK),
    ('Label-free / reference-free → escapes the circularity that killed every identity route.', 0, BLUE),
], top=Inches(1.5), size=17)

# 17b — this is ISFrag (prior art)
s = newslide(); header(s, 'Probe 6 — prior art', 'This is ISFrag — we re-derived a published method')
bullets(s, [
    ('ISFrag (Guo, Shen, XING, Yu, Huan — Anal. Chem. 2021; Shipei Xing = the reverse-search-poster author).', 0, DARK),
    ('Three tiers: L3 co-elution (MS1 PEAK-SHAPE correlation) → L2 fragment-m/z in precursor MS² → L1 fragmentation match.', 0, DARK),
    ('ISFrag\'s Level-1 uses a REVERSE DOT PRODUCT (reverse_dp) + containment ratio (>0.7) — exactly our containment(F⊂P) at 0.5.', 0, BLUE),
    ('So our detector ≈ ISFrag\'s strongest tier. Strong validation — and a directive: DON\'T rebuild ISFrag.', 0, TEAL),
    ('What ISFrag has that we don\'t: MS1 peak-shape co-elution (needs raw EIC → Gert). Our novelty = the APPLICATION:', 0, DARK),
    ('ISFrag cleans a raw feature table upstream; we audit ALREADY-ANNOTATED bins for label errors + confidence triage.', 1, GREY),
], top=Inches(1.5), size=17)

# 17c — what we ruled out: continuous noise-budget
s = newslide(); header(s, 'Probe 6 — and what we ruled out', 'The continuous "noise budget" feature — tested, negative')
bullets(s, [
    ('Tempting idea: decompose query intensity {reference / siblings / co-eluting neighbors / UNEXPLAINED}, use unexplained-residual as a calibration/abstention signal.', 0, DARK),
    ('Descriptively true: median query = 0.73 reference + 0.06 neighbors + only ~0.07 truly unexplained → most noise is STRUCTURED.', 0, DARK),
    ('But as a signal it FAILS: unexplained-residual TP/FP AUC 0.555 (vs purity 0.689); within low-purity band 0.539 ≈ chance.', 0, RED),
    ('The "low-purity TP = benign explained noise" hypothesis is mildly INVERTED (FP explained-away 0.560 > TP 0.542).', 0, RED),
    ('Expected, not a bug: it is a NOISE axis, and noise ⊥ identity — testing it against TP/FP is the wrong yardstick.', 0, GREY),
    ('Only legit role = calibration/abstention — UNVALIDATABLE without a noise-quality ground truth (BinBase MS1 purity / curator junk-label). PARKED.', 0, BLUE),
], top=Inches(1.5), size=16)

# 18 — verdict / next
s = newslide(NAVY)
_txt(s, Inches(0.8), Inches(0.6), Inches(11.7), Inches(0.9), 'Verdict & next step', size=32, color=WHITE, bold=True)
bullets(s, [
    ('Reverse as identity, standalone scalar, OR continuous noise-budget → exhausted/ruled out. Stop investing.', 0, RGBColor(0xF5,0xB7,0xB1)),
    ('What survives: the DISCRETE relational-ISF call (= ISFrag Level-1, applied to annotated bins). Don\'t rebuild ISFrag.', 0, WHITE),
    ('NEXT (still experiment): curator-review deck of the 558 strong-set ISF + the proline-type mis-annotations.', 0, RGBColor(0xA3,0xE4,0xD7)),
    ('THEN: wire as detection-side mirror of P1 compute_confirmation(); expand peak coverage; (later) ISFrag peak-shape co-elution via Gert\'s EIC.', 0, WHITE),
    ('PARKED behind a noise-quality ground truth (BinBase MS1 purity): the continuous calibration/abstention signal.', 0, RGBColor(0xA3,0xE4,0xD7)),
    ('Nothing here is production yet — these are probes.', 0, RGBColor(0xF5,0xB7,0xB1)),
], top=Inches(1.7), size=18)

prs.save(str(OUT))
print('wrote', OUT)
print('slides:', len(prs.slides._sldIdLst))
