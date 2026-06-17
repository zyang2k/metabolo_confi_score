"""build_unconfirmed_deorphan_deck.py — short deck on the UNCONFIRMED candidate-bin de-orphan.

Numbers from the 2026-06-16 run (project_reverse_search_isf_probe_20260609 memory);
examples pulled live from data/lcb_unconfirmed_denoise.csv + data/lcb_hilicneg_bins.csv.
"""
import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR

ROOT = Path(__file__).resolve().parent.parent.parent
OUT = ROOT / "reports" / "unconfirmed_deorphan_slides.pptx"
TMP = ROOT / "data" / "_deorphan_deck"; TMP.mkdir(parents=True, exist_ok=True)

SW, SH = Inches(13.333), Inches(7.5)
NAVY=RGBColor(0x1F,0x2D,0x3D); TEAL=RGBColor(0x16,0xA0,0x85); RED=RGBColor(0xC0,0x39,0x2B)
ORANGE=RGBColor(0xE6,0x7E,0x22); GREY=RGBColor(0x5D,0x6D,0x7E); LGREY=RGBColor(0xEC,0xF0,0xF1)
WHITE=RGBColor(0xFF,0xFF,0xFF); DARK=RGBColor(0x2B,0x2B,0x2B)
TEALh, GREYh, ORANGEh, REDh = "#16A085", "#5D6D7E", "#E67E22", "#C0392B"

prs = Presentation(); prs.slide_width=SW; prs.slide_height=SH
BLANK = prs.slide_layouts[6]


def _box(s, l, t, w, h):
    tb = s.shapes.add_textbox(l, t, w, h); tb.text_frame.word_wrap = True; return tb.text_frame

def _para(tf, text, size, color=DARK, bold=False, bullet=False, space=6, first=False):
    p = tf.paragraphs[0] if first else tf.add_paragraph()
    p.space_after = Pt(space)
    r = p.add_run(); r.text = ("•  " if bullet else "") + text
    r.font.size = Pt(size); r.font.bold = bold; r.font.color.rgb = color; r.font.name = "Calibri"
    return p

def bar(s, t, l, top):
    bg = s.shapes.add_shape(1, 0, 0, SW, Inches(1.15)); bg.fill.solid()
    bg.fill.fore_color.rgb = NAVY; bg.line.fill.background()
    tf = bg.text_frame; tf.margin_left=Inches(0.5); tf.vertical_anchor=MSO_ANCHOR.MIDDLE
    p=tf.paragraphs[0]; r=p.add_run(); r.text=t; r.font.size=Pt(26); r.font.bold=True; r.font.color.rgb=WHITE

def slide_title():
    s = prs.slides.add_slide(BLANK)
    bg = s.shapes.add_shape(1,0,0,SW,SH); bg.fill.solid(); bg.fill.fore_color.rgb=NAVY; bg.line.fill.background()
    tf=_box(s, Inches(0.8), Inches(2.4), Inches(11.7), Inches(2.8))
    _para(tf,"De-orphaning UNCONFIRMED candidate bins",34,WHITE,True,first=True)
    _para(tf,"A reverse / containment match: which candidates are artifacts of compounds we already have?",18,LGREY,space=18)
    _para(tf,"LCBinBase · HILIC-negative · reverse-score pipeline · 2026-06-16",14,TEAL)

def slide_bullets(title, items, sub=None):
    s = prs.slides.add_slide(BLANK); bar(s,title,0,0)
    tf=_box(s, Inches(0.7), Inches(1.5), Inches(11.9), Inches(5.6))
    if sub: _para(tf, sub, 16, GREY, True, space=12, first=True); first=False
    else: first=True
    for it in items:
        lvl=it[0] if isinstance(it,tuple) else 0; txt=it[1] if isinstance(it,tuple) else it
        _para(tf, txt, 18 if lvl==0 else 15, DARK if lvl==0 else GREY, bullet=True, space=8, first=first)
        first=False
    return s

def slide_fig(title, fig_path, caption=None):
    s = prs.slides.add_slide(BLANK); bar(s,title,0,0)
    s.shapes.add_picture(fig_path, Inches(1.4), Inches(1.4), height=Inches(5.2))
    if caption:
        tf=_box(s, Inches(0.7), Inches(6.7), Inches(11.9), Inches(0.7))
        _para(tf, caption, 14, GREY, first=True)
    return s

# ---------- figures ----------
def fig_breakdown():
    fig, ax = plt.subplots(figsize=(7.2,4.6))
    seg = [("in-source fragments",2171,TEALh),("adducts",354,"#1F8F76"),
           ("¹³C-isotope peaks",105,"#0E6B57"),("no confirmed parent\n(candidate novels)",23897,GREYh)]
    labels=[x[0] for x in seg]; vals=[x[1] for x in seg]; cols=[x[2] for x in seg]
    bars=ax.barh(range(len(seg)), vals, color=cols)
    for i,(v) in enumerate(vals):
        ax.text(v+250, i, f"{v:,} ({100*v/26527:.1f}%)", va="center", fontsize=11, weight="bold")
    ax.set_yticks(range(len(seg))); ax.set_yticklabels(labels, fontsize=11)
    ax.invert_yaxis(); ax.set_xlim(0,27000); ax.set_xlabel("UNCONFIRMED candidate bins (26,527 total)",fontsize=11)
    ax.set_title("9.9% are a relational ion of a co-eluting CONFIRMED compound",fontsize=12.5,weight="bold")
    ax.spines[["top","right"]].set_visible(False)
    fig.tight_layout(); p=TMP/"breakdown.png"; fig.savefig(p,dpi=200); plt.close(fig); return str(p)

def fig_controls():
    labels=["real","random Δm/z\n(break chemistry)","RI-shuffle\n(break co-elution)","decoy losses\n(non-physical)"]
    vals=[9.9,1.5,7.6,0.0]; cols=[TEALh,REDh,ORANGEh,REDh]
    fig,ax=plt.subplots(figsize=(7.2,4.6))
    b=ax.bar(labels,vals,color=cols,width=0.62)
    for r,v in zip(b,vals): ax.text(r.get_x()+r.get_width()/2, v+0.25, f"{v:.1f}%", ha="center", fontsize=12, weight="bold")
    ax.set_ylabel("explained rate (%)",fontsize=11); ax.set_ylim(0,15)
    ax.set_title("Null controls: the chemistry carries the signal",fontsize=12.5,weight="bold")
    ax.spines[["top","right"]].set_visible(False); ax.grid(axis="y",alpha=0.25)
    fig.tight_layout(); p=TMP/"controls.png"; fig.savefig(p,dpi=200); plt.close(fig); return str(p)

def fig_crossmethod():
    groups = ["HILIC\n(polar metab.)", "Lipids C18\nnegative", "Lipids C18\npositive"]
    flagged = [9.9, 18.2, 10.2]; rish = [7.6, 3.8, 2.2]; rmz = [1.5, 0.9, 1.1]
    dom = ["ISF", "ISF+iso", "adducts"]
    x = np.arange(3); w = 0.25
    fig, ax = plt.subplots(figsize=(8.2, 4.6))
    ax.bar(x - w, flagged, w, label="flagged (real)", color=TEALh)
    ax.bar(x,     rish,    w, label="RI-shuffle (null)", color=ORANGEh)
    ax.bar(x + w, rmz,     w, label="random Δm/z (null)", color=REDh)
    for i, v in enumerate(flagged):
        ax.text(x[i] - w, v + 0.3, f"{v}%", ha="center", fontsize=11.5, weight="bold")
        ax.text(x[i] - w, -1.4, dom[i], ha="center", fontsize=9, color=GREYh, style="italic")
    ax.set_xticks(x); ax.set_xticklabels(groups, fontsize=10.5)
    ax.set_ylabel("% of UNCONFIRMED candidate bins", fontsize=11); ax.set_ylim(0, 21)
    ax.set_title("Lipids: tighter null controls; relation profile tracks ionization chemistry",
                 fontsize=12, weight="bold")
    ax.legend(fontsize=10, frameon=False); ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout(); p = TMP / "crossmethod.png"; fig.savefig(p, dpi=200); plt.close(fig); return str(p)

# ---------- examples table ----------
def example_rows():
    o = pd.read_csv(ROOT/"data"/"lcb_unconfirmed_denoise.csv")
    b = pd.read_csv(ROOT/"data"/"lcb_hilicneg_bins.csv")[["wiki_id","name","precursor_mz"]]
    b = b.rename(columns={"wiki_id":"parent","name":"pn","precursor_mz":"pmz"})
    e = o[o.explained].merge(b,on="parent",how="left")
    e = e[~e["pn"].fillna("").str.startswith(("unknown_","yy_"))]
    e = e[e["precursor_mz"] >= 100]   # cleaner illustrations: skip bare low-mass marker ions
    e["rel"]=e["relation"].fillna(""); picks=[]
    def take(mask,n):
        sub=e[mask].sort_values("containment",ascending=False)
        for _,r in sub.head(n).iterrows(): picks.append(r)
    take(e.rel.str.contains("CO2"),1); take(e.rel.str.contains("NH3"),1)
    take(e.rel.str.contains("hexose|pentose"),1); take(e.rel.str.contains("H2O")& ~e.rel.str.contains("CO2"),1)
    take(e.rel.str.startswith("isotope"),2); take(e.rel.str.startswith("adduct"),2)
    rows=[]
    seen=set()
    for r in picks:
        key=(round(r.precursor_mz,3),r.rel)
        if key in seen: continue
        seen.add(key)
        rows.append((f"{r.precursor_mz:.4f}", r.rel.replace("isotope:","iso ").replace("ISF:","−").replace("adduct:","adduct "),
                     f"{r.containment:.2f}", str(r["pn"])[:34], f"{r.pmz:.4f}"))
    return rows[:8]

def slide_examples():
    s = prs.slides.add_slide(BLANK); bar(s,"Examples — candidate → confirmed parent it's a related ion of",0,0)
    rows = example_rows()
    cols = ["candidate m/z","relation","contain.","confirmed parent","parent m/z"]
    tbl = s.shapes.add_table(len(rows)+1, len(cols), Inches(0.6), Inches(1.5), Inches(12.1), Inches(0.5)).table
    widths=[Inches(2.2),Inches(2.4),Inches(1.4),Inches(4.3),Inches(1.8)]
    for j,w in enumerate(widths): tbl.columns[j].width=w
    for j,c in enumerate(cols):
        cell=tbl.cell(0,j); cell.text=c; cell.fill.solid(); cell.fill.fore_color.rgb=NAVY
        pr=cell.text_frame.paragraphs[0]; pr.runs[0].font.size=Pt(13); pr.runs[0].font.bold=True; pr.runs[0].font.color.rgb=WHITE
    for i,row in enumerate(rows,1):
        for j,val in enumerate(row):
            cell=tbl.cell(i,j); cell.text=val
            cell.fill.solid(); cell.fill.fore_color.rgb = LGREY if i%2 else WHITE
            r0=cell.text_frame.paragraphs[0].runs[0]; r0.font.size=Pt(12.5); r0.font.color.rgb=DARK
            if j==1: r0.font.bold=True; r0.font.color.rgb=TEAL
    tf=_box(s, Inches(0.6), Inches(6.7), Inches(12), Inches(0.7))
    _para(tf,"Textbook in-source losses (carboxylic acids shedding CO₂, etc.); ¹³C peaks; same fragment often spawns several candidates.",13,GREY,first=True)

# ---------- build ----------
slide_title()
slide_bullets("The problem", [
    "BinBase generates candidate bins; the UNCONFIRMED ones passed every QC gate (ion count, S/N, scans, ISTD, clean-spectra) and await acceptance as real compounds.",
    "But some candidates aren't new — they're in-source fragments, adducts, or isotope peaks of compounds we've ALREADY confirmed.",
    "A recurring in-source fragment recurs as often as its parent, so it clears the threshold and becomes its own candidate.",
    "Question: which candidates should NOT be promoted (artifact of an existing compound) vs are genuinely novel?",
], sub="UNCONFIRMED = the promote / don't-promote pile (26,527 candidates, HILIC-negative)")
slide_bullets("The method — reverse / containment", [
    "Flag a candidate O as a related ion of a confirmed bin C when ALL hold:",
    (1,"they co-elute (Δ retention-index ≤ 4)"),
    (1,"prec(C) − prec(O) is a neutral-loss / adduct / ¹³C-isotope spacing"),
    (1,"O's ions are CONTAINED in C's spectrum (reverse score)"),
    "Why containment, not forward cosine: a fragment is a SUBSPECTRUM of its parent — forward cosine is dragged down by the parent's extra ions and misses the link; one-sided containment stays high.",
    "Isotopes restricted to ¹³C (³⁴S/³⁷Cl need a parent formula the table doesn't carry).",
])
slide_fig("Result — HILIC (polar metabolites)", fig_breakdown(),
          "2,630 / 26,527 (9.9%) flagged as relational artifacts → don't-promote list; 1,576 link to a NAMED confirmed parent. Dominated by in-source fragments.")
slide_fig("Generalizes across methods — and tracks the chemistry", fig_crossmethod(),
          "C18-neg 18.2% (ISF-dominated: water/CO₂/acetate); C18-pos 10.2% (ADDUCT-dominated: TG [M+NH₄]⁺, PC [M+Na/K]⁺ — exactly what positive-mode lipids form). RP co-elution is strongly discriminating both polarities (RI-shuffle 2–4% vs HILIC's 8%). The relation profile shifting neg→pos is itself a validation.")
slide_examples()
slide_fig("How far to trust it", fig_controls(),
          "Scramble masses → 1.5% (chemistry is real); non-physical losses → 0%. Not sparse spectra (median 9 peaks).")
slide_bullets("Honest caveats + a side finding", [
    "Method-level co-elution is WEAK on its own — the confirmed-bin retention axis is crowded, so most candidates have some neighbour within a few RI units (RI-shuffle 7.6%).",
    (1,"→ firmest calls are the named-parent in-source fragments and adducts."),
    "Isotopes are a small slice (105) — they MUST be the heavier ion of a confirmed bin; an earlier two-sided match inflated them ~10×, now fixed.",
    "The 'novel' 90% residual is candidate novels — still needs the spectrum-quality (S_norm) check before being called real.",
    "Side finding: the same test among CONFIRMED bins flags ~13% (858/6,398) as related-ions of another confirmed bin — a separate library-cleanup list (re-finds cases already yy_-flagged).",
])
slide_bullets("Suggested use + next steps", [
    "Per pending pile: auto-label each candidate 'in-source fragment / ¹³C-isotope / adduct of [compound] — don't promote' vs 'no confirmed parent — review'.",
    "Stops redundant artifacts competing for promotion; surfaces the genuinely-new candidates.",
    "Next: (a) spot-check ~20 flagged candidates; (b) ship the confirmed-bin cleanup list; (c) run on HILIC-pos / C18 / lipidomics to confirm it generalizes.",
])

prs.save(str(OUT))
print("saved", OUT, "—", len(prs.slides.__iter__.__self__._sldIdLst), "slides")
