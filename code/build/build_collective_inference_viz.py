"""
Visualization: NLP-style collective inference for bin-level annotation.

Two-panel comparison using the putrescine FP example:
  TOP    — current pipeline: each bin scored in isolation, model picks the FP.
  BOTTOM — collective inference: chemistry-aware edges between bins flip the call.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from matplotlib.lines import Line2D

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Helvetica', 'Arial', 'DejaVu Sans']
plt.rcParams['pdf.fonttype'] = 42

INK = '#1f2937'
MUTED = '#6b7280'
ACCENT = '#0e7490'
GOOD = '#15803d'
BAD = '#b91c1c'
ORANGE = '#c2410c'
PURPLE = '#7c3aed'
BG = '#f9fafb'

# ---- Define the bins (RT in seconds, m/z) and their Q1 (row-level) scores ----
bins = {
    # The FP putrescine — Oliver yy_'d. Row-level model loves it.
    'A': dict(rt=120, mz=89.1, name='putrescine', q1=0.92,
              label='Bin A — putrescine?\nQ1: 0.92  (model: KEEP)',
              truth='yy_'),
    # The TP putrescine — moderate spectral match, but right RT.
    'B': dict(rt=145, mz=89.1, name='putrescine', q1=0.65,
              label='Bin B — putrescine?\nQ1: 0.65  (model: maybe)',
              truth='TP'),
    # Larger compound co-eluting with Bin A — possible ISF parent.
    'C': dict(rt=120, mz=204.1, name='unknown parent', q1=0.55,
              label='Bin C — large compound\nco-elutes with A',
              truth='?'),
    # Putrescine [M+Na]+ — would co-elute with B at the +Na mass offset.
    'D': dict(rt=144, mz=111.1, name='putrescine [M+Na]+', q1=0.70,
              label='Bin D — putrescine [M+Na]+\nat B\'s RT',
              truth='TP'),
}

# ---- Edges (relationship type, from, to, label) ----
edges = [
    ('sibling',   'A', 'B', 'sibling repulsion\n(same compound, different RT)'),
    ('multi_adduct', 'B', 'D', 'multi-adduct attraction\n(B is M+H, D is M+Na, same RT)'),
    ('isf_parent', 'A', 'C', 'ISF parent-child\n(A could be ISF of C)'),
]


def setup_axis(ax, title, subtitle=None):
    ax.set_xlim(80, 175)
    ax.set_ylim(50, 230)
    ax.set_xlabel('Retention time (s)', fontsize=11, color=INK)
    ax.set_ylabel('Precursor m/z', fontsize=11, color=INK)
    ax.grid(True, linestyle=':', alpha=0.35, color=MUTED)
    ax.set_facecolor(BG)
    ax.tick_params(colors=MUTED, labelsize=9)
    for spine in ax.spines.values():
        spine.set_edgecolor(MUTED)
        spine.set_linewidth(0.8)
    ax.set_title(title, fontsize=14, fontweight='bold', color=INK, loc='left', pad=8)
    if subtitle:
        ax.text(0.0, 1.05, subtitle, transform=ax.transAxes,
                fontsize=10, color=MUTED, va='bottom')


def draw_bin(ax, key, b, decision=None, highlight=None):
    # Color by Q1 confidence (blue scale)
    intensity = b['q1']
    face = plt.cm.Blues(0.3 + 0.5 * intensity)
    edge_color = INK
    edge_width = 1.0

    # Override edge based on decision (for the bottom panel)
    if decision == 'keep':
        edge_color = GOOD
        edge_width = 2.5
    elif decision == 'reject':
        edge_color = BAD
        edge_width = 2.5
    elif decision == 'demote':
        edge_color = ORANGE
        edge_width = 2.0

    circ = patches.Circle((b['rt'], b['mz']), radius=4.8, facecolor=face,
                          edgecolor=edge_color, linewidth=edge_width, zorder=4)
    ax.add_patch(circ)
    ax.text(b['rt'], b['mz'], key, ha='center', va='center',
            fontsize=12, fontweight='bold', color=INK, zorder=5)

    # Label box
    label_x = b['rt'] + 7
    label_y = b['mz']
    ha = 'left'
    if key == 'D':
        label_x = b['rt'] + 6
        label_y = b['mz'] - 14
    if key == 'C':
        label_x = b['rt'] - 7
        label_y = b['mz'] + 8
        ha = 'right'
    if key == 'B':
        label_x = b['rt'] + 7
        label_y = b['mz'] + 10
    if key == 'A':
        label_x = b['rt'] + 7
        label_y = b['mz'] + 10

    ax.text(label_x, label_y, b['label'], fontsize=8.5, color=INK,
            ha=ha, va='center', zorder=5,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                      edgecolor=MUTED, linewidth=0.5))


EDGE_STYLE = {
    'sibling':      dict(color=BAD, linestyle='--', label='Sibling repulsion (same name)'),
    'multi_adduct': dict(color=GOOD, linestyle='-',  label='Multi-adduct attraction (cross-m/z, same RT)'),
    'isf_parent':   dict(color=ORANGE, linestyle=':', label='ISF parent-child relationship'),
}


def draw_edge(ax, b1, b2, kind):
    style = EDGE_STYLE[kind]
    ax.annotate('', xy=(b2['rt'], b2['mz']), xytext=(b1['rt'], b1['mz']),
                arrowprops=dict(arrowstyle='<->', color=style['color'],
                                linestyle=style['linestyle'],
                                linewidth=2.0, alpha=0.85,
                                connectionstyle='arc3,rad=0.0'),
                zorder=3)


# ---- Build the figure ----
fig, (ax_top, ax_bot) = plt.subplots(2, 1, figsize=(13, 11),
                                      gridspec_kw=dict(hspace=0.35))

# === TOP: current pipeline ===
setup_axis(ax_top,
           'CURRENT — Each bin scored in isolation (row-level Q1 only)',
           'Model picks the bin with the highest row-level confidence. No cross-bin awareness.')
for key, b in bins.items():
    draw_bin(ax_top, key, b)

# A box highlighting what the model picks
ax_top.annotate('Model\'s pick: Bin A\n(Q1 confidence = 0.92)\nBut Oliver labels this yy_',
                xy=(120, 89.1), xytext=(85, 175),
                fontsize=10, color=BAD, fontweight='bold',
                ha='left', va='center',
                bbox=dict(boxstyle='round,pad=0.5', facecolor='#fef2f2',
                          edgecolor=BAD, linewidth=1.5),
                arrowprops=dict(arrowstyle='->', color=BAD, linewidth=1.5))

# === BOTTOM: collective inference ===
setup_axis(ax_bot,
           'PROPOSED — Collective inference: chemistry-aware edges between bins',
           'Each edge type encodes a Q2-level relationship. The joint assignment flips the call.')

# Draw edges first (under nodes)
for kind, b1key, b2key, _ in edges:
    draw_edge(ax_bot, bins[b1key], bins[b2key], kind)

# Decisions under collective inference:
#   A: demoted/rejected — same-compound sibling with B; also potential ISF of C
#   B: promoted — has multi-adduct partner D at compatible RT/m/z
#   C: kept as parent context
#   D: corroborates B
decisions = {'A': 'reject', 'B': 'keep', 'C': 'demote', 'D': 'keep'}
for key, b in bins.items():
    draw_bin(ax_bot, key, b, decision=decisions[key])

# Decision callouts
ax_bot.annotate('Collective decision:\nBin B is the real putrescine\n(corroborated by D multi-adduct)',
                xy=(145, 89.1), xytext=(150, 200),
                fontsize=10, color=GOOD, fontweight='bold',
                ha='left', va='center',
                bbox=dict(boxstyle='round,pad=0.5', facecolor='#ecfdf5',
                          edgecolor=GOOD, linewidth=1.5),
                arrowprops=dict(arrowstyle='->', color=GOOD, linewidth=1.5))

ax_bot.annotate('Bin A demoted:\nsibling of B at wrong RT,\npossible ISF of C',
                xy=(120, 89.1), xytext=(82, 60),
                fontsize=9.5, color=BAD, fontweight='bold',
                ha='left', va='center',
                bbox=dict(boxstyle='round,pad=0.5', facecolor='#fef2f2',
                          edgecolor=BAD, linewidth=1.5),
                arrowprops=dict(arrowstyle='->', color=BAD, linewidth=1.5))

# Legend for edge types
legend_lines = [
    Line2D([0], [0], color=EDGE_STYLE['sibling']['color'],
           linestyle='--', linewidth=2.5, label='Sibling repulsion (same compound name)'),
    Line2D([0], [0], color=EDGE_STYLE['multi_adduct']['color'],
           linestyle='-',  linewidth=2.5, label='Multi-adduct attraction (mass/RT-coherent adduct partners)'),
    Line2D([0], [0], color=EDGE_STYLE['isf_parent']['color'],
           linestyle=':',  linewidth=2.5, label='ISF parent-child (co-elution + mass relationship)'),
]
ax_bot.legend(handles=legend_lines, loc='lower right', fontsize=9,
              framealpha=0.95, edgecolor=MUTED)

# Overall title
fig.suptitle('Collective inference for LC-MS annotation — NLP-inspired Q2 layer',
             fontsize=16, fontweight='bold', color=INK, y=0.995)

# Footer / framing
fig.text(0.5, 0.005,
         'Analog of NLP entity linking: don\'t decide each mention in isolation; jointly optimize over all mentions in a document. '
         'Here: don\'t decide each bin\'s identity from row features alone — jointly assign identities such that sibling, '
         'multi-adduct, and ISF relationships are mutually consistent.',
         ha='center', fontsize=9, color=MUTED, style='italic', wrap=True)

out = 'reports/collective_inference_viz.pdf'
fig.savefig(out, bbox_inches='tight', dpi=200)
print(f'Wrote {out}')

# Also save as PNG for slides
out_png = 'reports/collective_inference_viz.png'
fig.savefig(out_png, bbox_inches='tight', dpi=200)
print(f'Wrote {out_png}')

plt.close(fig)
