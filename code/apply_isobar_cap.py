"""apply_isobar_cap.py — post-hoc isobar-ambiguity cap on the GBM confidence.

WHY post-hoc, not a feature: sim_gap is already a GBM feature (6.8% gain) yet the model
stays ~96% confident on near-tied isobars, because Oliver's labels are RT-driven and never
penalised isomer ambiguity (project_oliver_curation_method). A signal absent from the labels
can't be learned — so we IMPOSE it as a transparent rule on the output instead.

Rule (user-chosen "cap, don't crush", 2026-06-04):
  For each bin, in its MassWiki candidate list find the best DIFFERENT compound (different
  hit_ik14 AND not a stereo-variant of the top name). If that competitor is within ISOBAR_GAP
  of the top AND the top is a real match (>= ISOBAR_MIN_SIM), cap confidence at CAP.
  Otherwise unchanged. Downward only.

This SHIFTS the meaning of confidence from "valid annotation at the right RT" (label-calibrated)
to "valid annotation AND the specific structure is MS2-distinguishable". The capped value is a
principled heuristic, NOT label-calibrated (no isomer ground truth). Raw score kept as
`confidence_raw_cal`. Authentic-standard-verified bins (Min) are EXEMPT — a standard confirms
identity regardless of MS2 ambiguity.

Output: data/deliverable_scores_capped.csv
"""
import re, csv
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DELIV = ROOT/'data'/'deliverable_scores_v2.csv'
CAND  = ROOT/'data'/'candidate_scores_v2.csv'
MIN   = ROOT/'benchmark'/'_internal'/'min_verified_entries.csv'
OUT   = ROOT/'data'/'deliverable_scores_capped.csv'

ISOBAR_GAP = 0.05
ISOBAR_MIN_SIM = 0.50
CAP = 0.70             # MS² ambiguous, RT inconclusive
CAP_CONFLICT = 0.50    # MS² ambiguous AND RT fits the competitor better (strongest flag)
RT_MARGIN = 10.0       # s — |ΔRT| difference for RT to break the tie

_META = ('collisionenergy', 'massbank', ';', ' lc-', '@', ' - ', ' ev', 'energy')
def norm_name(n):
    """Compound-identity key: cut library metadata + collapse stereo (KEEP locants/sequence
    so real isomers stay distinct). 'GLUTARIC ACID CollisionEnergy0' -> 'glutaric acid'."""
    s = (n or '').lower().strip()
    for m in _META:
        i = s.find(m)
        if i > 2: s = s[:i]      # keep at least a stub; don't blank short names
    s = re.sub(r'^(d|l|dl|rac|\(r\)|\(s\)|\(\+\)|\(-\)|\(\+/-\)|\(±\))[- ]', '', s)
    s = re.sub(r'[()±]', '', s)
    return re.sub(r'\s+', ' ', s).strip()

def fnum(x, d=None):
    try: return float(x)
    except (TypeError, ValueError): return d


def main():
    cand_by_wiki = {}
    for r in csv.DictReader(open(CAND)):
        if r.get('hit_ik14'):
            cand_by_wiki.setdefault(r['wiki_id'], []).append(r)

    min_wiki = set()
    if MIN.exists():
        min_wiki = {r['wiki_id'] for r in csv.DictReader(open(MIN)) if r.get('wiki_id')}

    rows = list(csv.DictReader(open(DELIV)))
    out_cols = rows[0].keys()
    extra = ['confidence_raw_cal','isobar_gap','isobar_competitor','isobar_rt_top','isobar_rt_comp',
             'isobar_rt_verdict','isobar_capped','confidence_capped','confidence_capped_pct']
    n_cap=0; drops=[]; vcount={'rt_supports':0,'rt_conflict':0,'rt_inconclusive':0,'rt_missing':0}
    for r in rows:
        w=r['wiki_id']; conf=fnum(r.get('confidence'),0)
        r['confidence_raw_cal']=r.get('confidence')
        r['isobar_gap']=''; r['isobar_competitor']=''; r['isobar_rt_top']=''; r['isobar_rt_comp']=''
        r['isobar_rt_verdict']=''; r['isobar_capped']='False'
        r['confidence_capped']=conf; r['confidence_capped_pct']=r.get('confidence_pct')
        cl=sorted(cand_by_wiki.get(w,[]), key=lambda x:-fnum(x['entropy_similarity'],0))
        if len(cl)>=2 and w not in min_wiki:
            top=cl[0]; topsim=fnum(top['entropy_similarity'],0); tn=norm_name(top.get('name'))
            comp=next((x for x in cl[1:] if x['hit_ik14']!=top['hit_ik14'] and norm_name(x.get('name'))!=tn), None)
            if comp and topsim>=ISOBAR_MIN_SIM:
                gap=topsim-fnum(comp['entropy_similarity'],0)
                r['isobar_gap']=round(gap,3)
                r['isobar_competitor']=f"{comp.get('name')} ({fnum(comp['entropy_similarity'],0):.3f})"
                if gap<=ISOBAR_GAP:
                    # RT tie-breaker: |ΔRT| of annotated top vs competitor (per-candidate MolRex)
                    rt_t=fnum(top.get('signed_delta_rt')); rt_c=fnum(comp.get('signed_delta_rt'))
                    if rt_t is not None and rt_c is not None:
                        at,ac=abs(rt_t),abs(rt_c)
                        r['isobar_rt_top']=round(at,1); r['isobar_rt_comp']=round(ac,1)
                        if ac-at>RT_MARGIN:   verdict='rt_supports'      # observed RT fits the annotation
                        elif at-ac>RT_MARGIN: verdict='rt_conflict'      # observed RT fits the competitor
                        else:                 verdict='rt_inconclusive'
                    else:
                        verdict='rt_missing'
                    r['isobar_rt_verdict']=verdict; vcount[verdict]+=1
                    if verdict=='rt_supports':
                        capped=conf                       # RT breaks the tie → keep confidence
                        r['isobar_capped']='rt_resolved'
                    elif verdict=='rt_conflict':
                        capped=min(conf,CAP_CONFLICT); r['isobar_capped']='rt_conflict'
                    else:
                        capped=min(conf,CAP); r['isobar_capped']='True'
                    r['confidence_capped']=round(capped,4)
                    r['confidence_capped_pct']=round(capped*100,1)
                    if capped<conf: n_cap+=1; drops.append((r.get('annotation'),conf,capped,gap,comp.get('name'),verdict))

    with open(OUT,'w',newline='') as f:
        wr=csv.DictWriter(f, fieldnames=list(out_cols)+extra); wr.writeheader(); wr.writerows(rows)

    print(f'bins: {len(rows):,} | confidence reduced: {n_cap} | Min-exempt: {len(min_wiki)}')
    print(f'RT verdicts among isobar-flagged: {vcount}')
    print(f'  rt_supports → kept (RT breaks tie) | rt_conflict → cap {CAP_CONFLICT} | inconclusive/missing → cap {CAP}')
    conflicts=[d for d in drops if d[5]=='rt_conflict']
    print(f'\nRT-CONFLICT cases (MS² ambiguous AND RT fits competitor — strongest flags): {len(conflicts)}')
    for a,c0,c1,g,cn,v in sorted(conflicts,key=lambda x:x[3])[:8]:
        print(f"  {str(a)[:28]:28s} {c0*100:5.1f}% → {c1*100:4.0f}%  vs {str(cn)[:24]}")
    print(f'\nWrote {OUT}')


if __name__ == '__main__':
    main()
