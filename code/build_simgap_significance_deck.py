"""build_simgap_significance_deck.py — proposal deck: "When is a similarity gap significant?"
A spectrum-adaptive (ion-counting bootstrap) significance test for sim_gap.
Numbers: methylcytidine result from project_ms2_uq_literature_20260609 (0.006 gap, 94% direction, 1.6sigma).
"""
import os
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / 'reports' / 'simgap_significance_proposal.pptx'
TMP = ROOT / 'data' / '_simgap_deck_plots'; TMP.mkdir(parents=True, exist_ok=True)

SLIDE_W, SLIDE_H = Inches(13.333), Inches(7.5)
NAVY=RGBColor(0x1F,0x2D,0x3D); BLUE=RGBColor(0x2C,0x3E,0x50); TEAL=RGBColor(0x16,0xA0,0x85)
RED=RGBColor(0xC0,0x39,0x2B); ORANGE=RGBColor(0xE6,0x7E,0x22); GREY=RGBColor(0x5D,0x6D,0x7E)
LGREY=RGBColor(0xEC,0xF0,0xF1); WHITE=RGBColor(0xFF,0xFF,0xFF); DARK=RGBColor(0x2B,0x2B,0x2B)

def fig_significance():
    x=np.linspace(-0.010,0.020,400)
    def g(mu,s): return np.exp(-0.5*((x-mu)/s)**2)/(s*np.sqrt(2*np.pi))
    mu=0.006; s_sparse=mu/1.6; s_dense=mu/5.0          # methylcytidine is 1.6 sigma
    fig,ax=plt.subplots(figsize=(8.6,4.9))
    ax.fill_between(x, g(mu,s_sparse), color='#16A085', alpha=0.18)
    ax.plot(x, g(mu,s_sparse), color='#16A085', lw=2.4, label='sparse 4-ion spectrum (methylcytidine): σ large → gap = 1.6σ')
    ax.plot(x, g(mu,s_dense), color='#2C3E50', lw=2.4, label='dense spectrum (illustrative): σ small → same gap = 5σ')
    ax.axvline(0, color=RED.__str__() if False else '#C0392B', ls='-', lw=1.6)
    ax.text(0.0002,ax.get_ylim()[1]*0.04,'gap = 0\n(tie boundary)',color='#C0392B',fontsize=9,va='bottom')
    ax.axvline(mu, color='#555', ls='--', lw=1.2); ax.text(mu+0.0004, ax.get_ylim()[1]*0.92,'observed gap = 0.006',fontsize=10,color='#333')
    # shade sparse mass below 0 (where the runner-up would win)
    xs=x[x<=0]; ax.fill_between(xs, g(mu,s_sparse)[x<=0], color='#C0392B', alpha=0.25)
    ax.set_xlabel('similarity gap  (top-1 − next-best entropy similarity)', fontsize=12)
    ax.set_ylabel('bootstrap density', fontsize=12); ax.set_yticks([])
    ax.set_title('Same 0.006 gap, opposite verdicts —\nsignificance depends on the spectrum\'s ion-counting noise', fontsize=12.5, weight='bold')
    ax.legend(fontsize=9.5, loc='upper right', framealpha=0.92)
    fig.tight_layout(); p=TMP/'significance.png'; fig.savefig(p,dpi=200); plt.close(fig); return str(p)

# ---- slide helpers ----
def _bg(s,c): s.background.fill.solid(); s.background.fill.fore_color.rgb=c
def _box(s,l,t,w,h,fill=None):
    sp=s.shapes.add_shape(1,l,t,w,h); sp.fill.solid(); sp.fill.fore_color.rgb=fill or WHITE
    sp.line.fill.background(); sp.shadow.inherit=False; return sp
def _txt(s,l,t,w,h,text,size=18,color=DARK,bold=False,align=PP_ALIGN.LEFT,anchor=MSO_ANCHOR.TOP):
    tb=s.shapes.add_textbox(l,t,w,h); tf=tb.text_frame; tf.word_wrap=True; tf.vertical_anchor=anchor
    for i,ln in enumerate(text.split('\n')):
        p=tf.paragraphs[0] if i==0 else tf.add_paragraph(); p.alignment=align
        r=p.add_run(); r.text=ln; r.font.size=Pt(size); r.font.bold=bold; r.font.color.rgb=color; r.font.name='Calibri'
    return tb
def header(s,kick,title):
    _box(s,0,0,SLIDE_W,Inches(1.15),fill=BLUE)
    if kick: _txt(s,Inches(0.6),Inches(0.12),Inches(12),Inches(0.35),kick.upper(),size=13,color=RGBColor(0x9F,0xB6,0xC9),bold=True)
    _txt(s,Inches(0.6),Inches(0.42),Inches(12.1),Inches(0.62),title,size=26,color=WHITE,bold=True)
def bullets(s,items,top=Inches(1.55),left=Inches(0.7),width=Inches(12),size=18):
    tb=s.shapes.add_textbox(left,top,width,Inches(5.4)); tf=tb.text_frame; tf.word_wrap=True
    for i,it in enumerate(items):
        text,lvl,color=it if isinstance(it,tuple) else (it,0,DARK)
        p=tf.paragraphs[0] if i==0 else tf.add_paragraph(); p.level=lvl; p.space_after=Pt(5)
        r=p.add_run(); r.text=('•  ' if lvl==0 else '–  ')+text
        r.font.size=Pt(size-2*lvl); r.font.color.rgb=color; r.font.name='Calibri'
def chip(s,text,color,left,top,w=Inches(2.4)):
    _box(s,left,top,w,Inches(0.5),fill=color)
    _txt(s,left,top+Inches(0.04),w,Inches(0.42),text,size=14,color=WHITE,bold=True,align=PP_ALIGN.CENTER,anchor=MSO_ANCHOR.MIDDLE)

prs=Presentation(); prs.slide_width=SLIDE_W; prs.slide_height=SLIDE_H; blank=prs.slide_layouts[6]
def newslide(bg=WHITE): s=prs.slides.add_slide(blank); _bg(s,bg); return s
sig=fig_significance()

# 1 title
s=newslide(NAVY)
_txt(s,Inches(0.9),Inches(2.3),Inches(11.6),Inches(1.2),'When is a similarity gap significant?',size=42,color=WHITE,bold=True)
_txt(s,Inches(0.95),Inches(3.6),Inches(11.5),Inches(0.8),'A spectrum-adaptive significance test for sim_gap — replacing the flat floor with the spectrum\'s own ion-counting noise',size=19,color=RGBColor(0x9F,0xB6,0xC9))
_txt(s,Inches(0.95),Inches(6.4),Inches(11),Inches(0.5),'Ziyue Yang · West Coast Metabolomics Center · proposal',size=15,color=RGBColor(0x7F,0x8C,0x9D))

# 2 problem
s=newslide(); header(s,'The problem','sim_gap decides isomer calls — but the threshold ignores the spectrum')
bullets(s,[
    ('When a bin has competing candidates, the model leans on sim_gap = (sim to top-1) − (sim to next-best).',0,DARK),
    ('Today sim_gap is a RAW entropy-similarity difference, and "big enough to trust" is a flat floor (~0.05).',0,DARK),
    ('Oliver: "a 0.006 gap means something very different on a dense, well-counted spectrum than on a single-ion MS². The amount of the gap has to go into the scoring."',0,BLUE),
    ('A fixed floor cannot do that — it ignores how noisy THIS spectrum is.',0,RED),
])

# 3 motivating failure
s=newslide(); header(s,'Why it matters','Methylcytidine: a 0.006 gap drives a 30-point swing')
bullets(s,[
    ('2′-C- vs 2′-O-methylcytidine: entropy similarity 0.836 vs 0.831 — gap 0.006, on a 4-peak spectrum.',0,DARK),
    ('Deployed score: top-1 (2′-C) 22% vs 2′-O 6.5% — a ~30-point confidence swing from a 0.006 coin-flip.',0,RED),
    ('RT actually favors 2′-O (2′-C has no RT prediction) — but the spectral near-tie overrode it.',0,DARK),
    ('This is the bug: the model treats 0.006 as decisive when the spectrum cannot resolve these two.',0,BLUE),
])

# 4 current state
s=newslide(); header(s,'Current state','How sim_gap is calibrated today — it isn\'t')
bullets(s,[
    ('sim_gap = raw difference of two entropy-similarity numbers (build_features_v2.py:749; same in deploy).',0,DARK),
    ('No significance test, no spectrum-quality scaling, no soft-threshold.',1,RED),
    ('A 0.006 gap on a 4-ion spectrum and a 200-ion spectrum produce the IDENTICAL feature value.',1,GREY),
    ('The only calibration is downstream isotonic on the FINAL GBM probability — not on sim_gap.',0,DARK),
    ('The GBM could learn an implicit quality interaction, but methylcytidine proves it doesn\'t.',1,GREY),
])

# 5 insight
s=newslide(); header(s,'The insight','MS² peak intensities are ion counts → they carry error bars')
bullets(s,[
    ('Each peak\'s intensity is, at bottom, a count of ions hitting the detector.',0,DARK),
    ('Counts follow shot statistics: a peak from N ions has uncertainty ≈ √N.',0,DARK),
    ('So a spectrum is a NOISY measurement — and it tells you its own noise:',0,BLUE),
    ('many ions (dense) → measured precisely;  few ions (sparse) → measured sloppily.',1,GREY),
])

# 5b introduce the idea — how resampling works
s=newslide(); header(s,'The idea','Bootstrap: re-draw the spectrum within its ion-counting noise')
_img=str(ROOT/'reports'/'simgap_resampling_explainer.png')
if os.path.exists(_img):
    s.shapes.add_picture(_img, Inches(1.9), Inches(1.35), height=Inches(5.9))

# 6 proposal
s=newslide(); header(s,'The proposal','Bootstrap the gap within the spectrum\'s counting noise')
bullets(s,[
    ('Repeat B times: re-draw the spectrum within its ion-counting noise (Poisson / multinomial resample of peak intensities).',0,DARK),
    ('Recompute the similarity to top-1 AND next-best; record the gapᵦ.',0,DARK),
    ('The spread of {gapᵦ} is how much the gap wobbles under re-measurement. Two readings:',0,BLUE),
    ('direction consistency — does top-1 stay top-1?   AND   magnitude — gap / σ_gap (a z-score).',1,DARK),
    ('SIGNIFICANT if gap / σ_gap ≥ ~2 (≈97.5% direction consistency).',0,TEAL),
])

# 7 figure
s=newslide(); header(s,'Why it works','Spectrum-adaptive: the same gap, opposite verdicts')
s.shapes.add_picture(sig, Inches(0.7), Inches(1.35), height=Inches(4.9))
_txt(s,Inches(9.5),Inches(1.7),Inches(3.4),Inches(4.5),'Dense spectrum:\nσ small → 0.006 clears 2σ\n→ trust top-1.\n\nSparse spectrum:\nσ large → 0.006 buried in noise\n(red mass < 0 = runner-up\ncould win)\n→ TIE.',size=15,color=DARK)

# 8 methylcytidine under the gate
s=newslide(); header(s,'The fix in action','Methylcytidine under the significance gate')
bullets(s,[
    ('Bootstrap the 0.006 gap on the 4-peak spectrum:',0,DARK),
    ('direction consistency 94% — 2′-C usually stays ahead (gap is real in sign),',1,DARK),
    ('but only ~1.6σ in magnitude — BELOW the 2σ bar.',1,RED),
    ('Verdict: 0.836 ≈ 0.831 is a statistical TIE. MS² cannot separate these two here.',0,BLUE),
    ('→ call it a tie and defer to RT (which favors 2′-O) — instead of a 0.006 coin-flip driving 30 points.',0,TEAL),
])
chip(s,'TIE → DEFER TO RT',TEAL,Inches(0.7),Inches(5.9),w=Inches(3.4))

# 9 decision rule
s=newslide(NAVY)
_txt(s,Inches(0.8),Inches(0.7),Inches(11.7),Inches(0.9),'The decision rule',size=30,color=WHITE,bold=True)
bullets(s,[
    ('gap / σ_gap ≥ 2  →  the spectrum CAN separate the candidates  →  trust top-1\'s confidence.',0,RGBColor(0xA3,0xE4,0xD7)),
    ('gap / σ_gap < 2   →  MS² is "mute" between them  →  call a tie, defer to RT, raise the identity-ambiguity flag.',0,WHITE),
    ('This is exactly Oliver\'s two-number idea: a per-bin confidence AND a separate "can we differentiate" (Δ) confidence.',0,RGBColor(0x9F,0xB6,0xC9)),
],top=Inches(2.0),size=19)

# 10 caveats
s=newslide(); header(s,'Honest caveats','The bootstrap σ is a LOWER bound')
bullets(s,[
    ('It models only SHOT noise. It ignores chemical noise (co-isolation), the electronic floor, and peak correlations.',0,DARK),
    ('Biggest: it ignores collision-energy / run-to-run variability. Our glucose-1-P data: within-batch cosine 0.95+, cross-batch ~0 — real spectra vary FAR more than shot noise.',0,RED),
    ('So σ_gap is optimistic → the test is CONSERVATIVE: a FAIL (tie) is robust; a PASS is necessary, not sufficient.',0,BLUE),
    ('Methylcytidine FAILS even this optimistic test (1.6σ) → the tie verdict is rock-solid. Trust ties; don\'t over-trust passes.',0,TEAL),
    ('Also needs intensity → ion-count calibration (Orbitrap injection time / AGC). Gold-standard upgrade: empirical replicate resampling.',0,GREY),
])

# 11 deploy
s=newslide(); header(s,'Deployment','Cheap at scale + train/deploy parity')
bullets(s,[
    ('Bootstrap (500× per spectrum) is expensive in production.',0,DARK),
    ('Closed-form shortcut: propagate per-peak Poisson variance through the similarity (delta method) → σ_gap analytically, no loop.',0,TEAL),
    ('Parity is mandatory: the transform must be identical in build_features_v2.py (training) AND deploy/masswiki_gbm (inference) — or train/infer skew (the soft-threshold lesson).',0,RED),
    ('Slots into compute_sim_gaps: feed gap / σ_gap (or a tie flag) instead of the raw gap.',0,DARK),
])

# 12 plan
s=newslide(NAVY)
_txt(s,Inches(0.8),Inches(0.7),Inches(11.7),Inches(0.9),'Proposed plan',size=30,color=WHITE,bold=True)
bullets(s,[
    ('1. Prototype the gate in compute_sim_gaps (bootstrap first; delta-method for deploy).',0,WHITE),
    ('2. Validate on methylcytidine (tie → RT picks 2′-O) + a significance-vs-flat-floor comparison on the confusable band.',0,RGBColor(0xA3,0xE4,0xD7)),
    ('3. GUARD: must beat base + n_candidate_adducts + rel_nn_sim, not just base (don\'t re-litigate rank-vs-confidence).',0,WHITE),
    ('4. Ship with train/deploy parity + re-freeze; surface as the "Δ confidence" of Oliver\'s two-number scheme.',0,RGBColor(0xA3,0xE4,0xD7)),
    ('Novelty: ion-counting bootstrap on ENTROPY similarity for the annotation gap (SpecReBoot did cosine for networking).',0,RGBColor(0x9F,0xB6,0xC9)),
],top=Inches(1.9),size=18)

prs.save(str(OUT)); print('wrote',OUT,'| slides',len(prs.slides._sldIdLst))
