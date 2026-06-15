"""
Weekly report — slide-style PDF.
Generates reports/weekly_report_20260522.pdf with one slide per page.
"""

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

# Style
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Helvetica', 'Arial', 'DejaVu Sans']
plt.rcParams['pdf.fonttype'] = 42

INK = '#1f2937'      # slate-800
MUTED = '#6b7280'    # gray-500
ACCENT = '#0e7490'   # cyan-700
GOOD = '#15803d'     # green-700
BAD = '#b91c1c'      # red-700
BG_BOX = '#f3f4f6'   # gray-100

OUT_PATH = 'reports/weekly_report_20260522.pdf'


def new_slide():
    fig, ax = plt.subplots(figsize=(13.33, 7.5))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')
    return fig, ax


def slide_header(ax, title, sub=None, num=None, total=None):
    ax.text(0.05, 0.92, title, fontsize=24, fontweight='bold', color=INK, va='top')
    if sub:
        ax.text(0.05, 0.86, sub, fontsize=13, color=MUTED, va='top')
    ax.plot([0.05, 0.95], [0.83, 0.83], color=ACCENT, linewidth=2)
    if num is not None:
        ax.text(0.95, 0.04, f'{num} / {total}', fontsize=9, color=MUTED, ha='right')
    ax.text(0.05, 0.04, 'Ziyue · 2026-05-22', fontsize=9, color=MUTED)


def bullet(ax, x, y, text, size=12, color=INK, weight='normal', bullet_char='•'):
    ax.text(x, y, bullet_char, fontsize=size, color=ACCENT, va='top', weight='bold')
    ax.text(x + 0.018, y, text, fontsize=size, color=color, va='top',
            weight=weight, wrap=True)


def callout(ax, x, y, w, h, text, color=BG_BOX, edge=ACCENT, fontsize=11,
            text_color=INK, weight='normal'):
    from matplotlib.patches import FancyBboxPatch
    box = FancyBboxPatch((x, y - h), w, h,
                         boxstyle="round,pad=0.01,rounding_size=0.012",
                         linewidth=1.2, edgecolor=edge, facecolor=color)
    ax.add_patch(box)
    ax.text(x + 0.015, y - 0.02, text, fontsize=fontsize, color=text_color,
            va='top', weight=weight, wrap=True)


with PdfPages(OUT_PATH) as pdf:

    # ---------------- Slide 1 — Title ----------------
    fig, ax = new_slide()
    ax.text(0.5, 0.72, 'Weekly update', fontsize=40, fontweight='bold',
            color=INK, ha='center')
    ax.plot([0.3, 0.7], [0.66, 0.66], color=ACCENT, linewidth=3)
    ax.text(0.5, 0.58, 'Distilling Oliver · ensemble_sd shipped · Q1/Q2 framing',
            fontsize=16, color=MUTED, ha='center', style='italic')
    ax.text(0.5, 0.52, 'EIC for zz_ · Bayesian wind-down',
            fontsize=16, color=MUTED, ha='center', style='italic')
    ax.text(0.5, 0.20, 'Ziyue Yang', fontsize=14, color=INK, ha='center')
    ax.text(0.5, 0.16, '2026-05-22', fontsize=12, color=MUTED, ha='center')
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)

    # ---------------- Slide 2 — Distilling Oliver ----------------
    fig, ax = new_slide()
    slide_header(ax, '1. "Distilling" Oliver — curation reviewer',
                 sub='Automated reviewer that returns Oliver\'s likely verdict + 2-3 chemistry reasons. No training needed.',
                 num=2, total=7)
    bullet(ax, 0.05, 0.76, 'Rules from ~90 doubly-reviewed Orbi HILIC bins + handcrafted heuristics', size=12)
    bullet(ax, 0.05, 0.71, 'Multi-adduct anchor · ISF-no-partner · Δm thresholds · reverse cosine · noise per Leon', size=11, color=MUTED)
    bullet(ax, 0.05, 0.64, 'Validated on 535 doubly-reviewed bins (171 names where Oliver actively chose between RT siblings)', size=12)

    # Headline numbers
    callout(ax, 0.05, 0.55, 0.42, 0.10,
            'TP bins — would Oliver keep?\n88% agreement',
            color='#ecfdf5', edge=GOOD, fontsize=13, text_color=GOOD, weight='bold')
    callout(ax, 0.53, 0.55, 0.42, 0.10,
            'yy_ bins — would the reviewer catch?\n50% agreement — structural',
            color='#fef2f2', edge=BAD, fontsize=13, text_color=BAD, weight='bold')

    bullet(ax, 0.05, 0.39, 'Why the gap? My hypothesis: the bottleneck is upstream.', size=13, weight='bold')
    ax.text(0.07, 0.34,
            'Same reason GBM plateaus at AUC ~0.92 regardless of added features (Morgan, learned embeddings),\n'
            'and every neural approach converges to the same ceiling.',
            fontsize=11, color=INK, va='top')

    callout(ax, 0.05, 0.21, 0.90, 0.16,
            'Features answer Q1 — "is the library match plausible?"\n'
            'Labels (TP / yy_) encode Q2 — "should this bin be annotated this way?"\n'
            'Q2 needs cross-bin context (sibling bins, per-compound RT prior, co-elution). Models can\'t see it.',
            color=BG_BOX, edge=ACCENT, fontsize=12)
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)

    # ---------------- Slide 3 — Why the model isn't improving ----------------
    fig, ax = new_slide()
    slide_header(ax, '2. Why the model isn\'t improving',
                 sub='Every architecture lands at the same ceiling. That\'s the diagnostic — the limit is data, not models.',
                 num=3, total=7)

    # Left: Q1/Q2 table
    ax.text(0.05, 0.76, 'The two questions', fontsize=14, fontweight='bold', color=INK)
    ax.text(0.05, 0.72,
            '             Q1 (row-level)             Q2 (bin + cross-bin)',
            fontsize=10, color=MUTED, family='monospace')
    ax.plot([0.05, 0.48], [0.71, 0.71], color=MUTED, linewidth=0.5)
    q1q2 = [
        ('Asks', '"library match plausible?"', '"should bin be annotated?"'),
        ('Inputs', 'spectrum, Δm, ΔRT, adduct', '+ sibling bins, RT prior, co-elution'),
        ('Who sees it', 'features, GBM, skill', 'Oliver (only)'),
        ('In our pipeline', 'covered', 'MISSING'),
    ]
    y = 0.66
    for label, q1, q2 in q1q2:
        ax.text(0.05, y, label, fontsize=10, color=INK, weight='bold')
        ax.text(0.17, y, q1, fontsize=10, color=INK)
        ax.text(0.36, y, q2, fontsize=10, color=INK)
        y -= 0.035

    # Right: convergence table
    ax.text(0.55, 0.76, 'Every architecture to ~0.92', fontsize=14, fontweight='bold', color=INK)
    ax.plot([0.55, 0.95], [0.71, 0.71], color=MUTED, linewidth=0.5)
    conv = [
        ('GBM (production)',         '0.913'),
        ('+ Morgan fingerprints',    '~0.913'),
        ('+ learned MS2 embeddings', '0.9133'),
        ('Random Forest',            '0.907'),
        ('MLP',                      '0.894'),
        ('Contrastive learning',     '~0.87'),
        ('Set Transformer',          '0.867 – 0.913'),
    ]
    y = 0.66
    for name, auc in conv:
        ax.text(0.55, y, name, fontsize=11, color=INK)
        ax.text(0.92, y, auc, fontsize=11, color=INK, weight='bold', ha='right')
        y -= 0.035

    # Putrescine example
    callout(ax, 0.05, 0.32, 0.90, 0.12,
            'Concrete: a yy_ putrescine bin had spectral match 1.00, tight Δm, multi-adduct anchor.\n'
            'Every row-level signal said keep. Oliver rejected it because the REAL putrescine lives in a\n'
            'different bin at the correct RT. Q1 = perfect. Q2 = rejected. No row-level feature recovers that.',
            color=BG_BOX, edge=ACCENT, fontsize=11)

    callout(ax, 0.05, 0.16, 0.90, 0.10,
            'Hypothesis (one line): our features and our labels are answering different questions.\n'
            'Until we close that gap upstream, no model is going to do meaningfully better than 0.92.',
            color='#fef3c7', edge='#a16207', fontsize=12, text_color=INK, weight='bold')
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)

    # ---------------- Slide 4 — ensemble_sd ----------------
    fig, ax = new_slide()
    slide_header(ax, '3. Curator triage signal — shipped ensemble_sd',
                 sub='K=10 IK14-group-block bootstrap of the production GBM. New column for curator-side triage.',
                 num=4, total=7)
    bullet(ax, 0.05, 0.76, 'Splits "intrinsically ambiguous bin" from "model lacks training neighbors"', size=12)
    bullet(ax, 0.05, 0.71, 'high sd + high confidence to candidate FP worth a second look', size=11, color=MUTED)
    bullet(ax, 0.05, 0.67, 'low sd + high confidence to trust the call', size=11, color=MUTED)
    bullet(ax, 0.05, 0.61, 'Two falsifiable tests pre-committed (ship if ≥1.5×, drop if <1.2×)', size=12)
    ax.text(0.075, 0.56,
            '· [0.6, 0.9] band-error ratio:  1.65× to 1.20× (after SMILES filter)\n'
            '· FP @ ≥0.7 concentration:     1.70× to 2.88×  (passes)',
            fontsize=11, color=INK, family='monospace', va='top')
    bullet(ax, 0.05, 0.45, 'Also: trustworthy isotonic refit (curator-verified slice) + 93-row Oliver override CSV', size=12)

    # Headline
    callout(ax, 0.05, 0.38, 0.90, 0.13,
            'HEADLINE — deck discrimination on Oliver\'s 23 cases\n'
            'Spread between OK and yy_ candidates: −3.9 pp  to  +40.9 pp\n'
            'Blank ≥0.7 suppression: 1.1% to 0.4%   ·   Blank ≥0.9: 0.3% to 0% (zero)   ·   Min Level-1: 108/108 unchanged',
            color='#ecfdf5', edge=GOOD, fontsize=12, text_color=GOOD, weight='bold')

    callout(ax, 0.05, 0.18, 0.90, 0.10,
            'Connection to Q1/Q2: ensemble_sd is the model saying "I\'m uncertain on this row." It correlates with hard\n'
            'Q2 cases but doesn\'t structurally fix the Q1/Q2 gap — a useful triage axis, not a ceiling-breaker.',
            color=BG_BOX, edge=ACCENT, fontsize=11)
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)

    # ---------------- Slide 5 — EIC / zz_ ----------------
    fig, ax = new_slide()
    slide_header(ax, '4. EIC access for catching zz_',
                 sub='Talked to Gert and Fanzhou — chromatographic data, not spectral, is what distinguishes zz_.',
                 num=5, total=7)
    bullet(ax, 0.05, 0.72, 'Tried training a model to flag zz_ (no-peak / column-level junk) from MS2 alone', size=12)
    bullet(ax, 0.05, 0.66, 'Doesn\'t work — the problem isn\'t spectral, it\'s chromatographic', size=12)
    bullet(ax, 0.075, 0.61, 'EIC peak shape distinguishes a real bin from zz_', size=11, color=MUTED, bullet_char='–')
    bullet(ax, 0.075, 0.56, 'We need access to the raw chromatogram', size=11, color=MUTED, bullet_char='–')

    callout(ax, 0.05, 0.45, 0.90, 0.18,
            'OPEN PROBLEM\n'
            'mzML is strictly needed, but processing mzML in batch locally has no clean path yet.\n\n'
            'Any ideas welcome.',
            color='#fef3c7', edge='#a16207', fontsize=13, weight='bold', text_color=INK)
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)

    # ---------------- Slide 6 — Statistical side / Bayesian wind-down ----------------
    fig, ax = new_slide()
    slide_header(ax, '5. Statistical side — closing out the (b) track',
                 sub='Huang τ² penalty lifted Bayesian +1 pt. Proposal: wind down further Bayesian investment.',
                 num=6, total=7)
    bullet(ax, 0.05, 0.76, 'τ² (between-channel disagreement variance) — signal is real:', size=12, weight='bold')
    ax.text(0.075, 0.72,
            '· FP candidates: ~2× the inverse-variance-weighted disagreement of TP\n'
            '· mean 1.31 (FP) vs 0.70 (TP)   ·   τ² alone AUC = 0.668',
            fontsize=11, color=INK, family='monospace', va='top')

    callout(ax, 0.05, 0.62, 0.42, 0.06,
            'Bayesian scorer + τ² penalty (α=0.5)\nAUC 0.873 to 0.882  — SHIPPED',
            color='#ecfdf5', edge=GOOD, fontsize=11, text_color=GOOD, weight='bold')
    callout(ax, 0.53, 0.62, 0.42, 0.06,
            'GBM + τ² as feature\nΔAUC −0.0003 — wash',
            color='#fef2f2', edge=BAD, fontsize=11, text_color=BAD, weight='bold')

    bullet(ax, 0.05, 0.49, 'Why? Trees already encode equivalent interactions (entropy_sim × sim_gap × signed_delta_rt).', size=11)
    bullet(ax, 0.05, 0.45, 'Bayesian sum-of-logLRs structurally cannot represent τ²-style interactions; trees can.', size=11, weight='bold')

    callout(ax, 0.05, 0.36, 0.90, 0.13,
            'PROPOSAL — wind down further Bayesian investment\n'
            'The (b) track is shipped at 0.882. Remaining gap to GBM (~3 pt) is the structural ceiling.\n'
            'Only unique advantage: per-channel logLR explainability (GBM SHAP covers most of it).\n'
            'Redirect effort to upstream feature engineering per Q1/Q2.',
            color='#fef3c7', edge='#a16207', fontsize=11, weight='bold', text_color=INK)

    bullet(ax, 0.05, 0.18, 'Keep (model-agnostic):  per-channel DDI as deck eval metric  ·  NoTA machinery preserved', size=11, color=MUTED)
    bullet(ax, 0.05, 0.13, 'Park, not kill:  Inforpower for principled family selection (no AUC lift, no Q1/Q2 fix)', size=11, color=MUTED)
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)

    # ---------------- Slide 7 — Where this leaves us ----------------
    fig, ax = new_slide()
    slide_header(ax, '6. Where this leaves us — next round of work',
                 sub='Modeling side is saturated. The unblocked lift is on the data side.',
                 num=7, total=7)

    ax.text(0.05, 0.77, 'Tier 1 — high-impact, can start now', fontsize=13, fontweight='bold', color=ACCENT)
    bullet(ax, 0.07, 0.72, 'Sibling-aware signal — per-bin, within-name comparison. Directly hits Test B failure mode.', size=11)
    bullet(ax, 0.07, 0.68, 'Per-compound RT prior — empirical RT distribution from confident prior TPs.', size=11)
    bullet(ax, 0.07, 0.64, 'Snorm noise gate for un-annotated bins — feature exists; data plumbing pending.', size=11)
    bullet(ax, 0.07, 0.60, 'Re-label split of yy_ into Q2 vs Q3 (curator-rejected wrong-anno vs no-anno).', size=11)

    ax.text(0.05, 0.52, 'Tier 2 — high-impact, but harder', fontsize=13, fontweight='bold', color=ACCENT)
    bullet(ax, 0.07, 0.47, 'EIC / chromatographic features for zz_ — blocked on mzML batch processing.', size=11)
    bullet(ax, 0.07, 0.43, 'Cross-replicate consistency — same compound, same RT across runs.', size=11)
    bullet(ax, 0.07, 0.39, 'Molecular networking / co-elution (GNPS-style).', size=11)

    ax.text(0.05, 0.31, 'Adjacent — supporting features (won\'t close Q1/Q2 alone)', fontsize=13, fontweight='bold', color=MUTED)
    bullet(ax, 0.07, 0.26, 'MolRex structure embedding · within-spectrum features · OOD detection (paused)', size=11, color=MUTED)

    callout(ax, 0.05, 0.18, 0.90, 0.10,
            'My read: start with the sibling-aware signal. Pure within-dataset groupby, no blockers,\n'
            'directly addresses the structural Test B gap. Highest ROI of any item on this list.',
            color='#ecfdf5', edge=GOOD, fontsize=12, weight='bold', text_color=GOOD)
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)

print(f'Wrote {OUT_PATH}')
