"""build_review_app.py — Path-1 interactive curator-review app (self-contained HTML).

Two queues for Oliver (project_two_path_roadmap_20260604):
  SUSPICIOUS  — confident annotations the MassWiki library search can BARELY distinguish
                from a different compound at the SAME precursor (small sim_gap to a
                different-IK14 candidate). The GBM confidence is driven by adduct + RT,
                NOT spectral isomer discrimination — and Oliver's labels never penalised
                isomer ambiguity, so the model stays confident on genuinely ambiguous
                isobars. This queue surfaces those. +/- : + = real problem, - = fine.
  PROMISING   — missed-annotation candidates (spectral propagation worklist).
                +/- : + = annotate it, - = no.

Confusability is LIBRARY-NATIVE (from the MassWiki candidate list, apples-to-apples with
the 0.849 the curator sees) — NOT the local bin-to-bin metric (different engine; that
stays a GBM feature, not a displayed reason).

Card: confidence + WHY-suspicious caveat, the competing candidate table, query MS2,
neighbour context, SPLASH, polarity. Mirror library trace for SUSPICIOUS is PENDING the
MassWiki reference re-fetch (local library_peaks_cache is unreliable — 28% of entries
≤5 peaks, doesn't reproduce the API sim). PROMISING mirror = query vs the confirmed
reference bin (both real observed spectra, trustworthy). species/organ/instrument/EIC
pending the scoped pull / FlashEIC service. localStorage autosave + Export.

Output: reports/curator_review_app.html
"""
import json, csv, re
from collections import Counter, defaultdict
from pathlib import Path

_META = ('collisionenergy','massbank',';',' lc-','@',' - ',' ev','energy')
def norm_name(n):
    """Compound-identity key — matches apply_isobar_cap.norm_name (cut library metadata + stereo,
    keep locants/sequence) so the queue and the cap agree on what counts as a competitor."""
    s=(n or '').lower().strip()
    for m in _META:
        i=s.find(m)
        if i>2: s=s[:i]
    s=re.sub(r'^(d|l|dl|rac|\(r\)|\(s\)|\(\+\)|\(-\)|\(\+/-\)|\(±\))[- ]','',s)
    s=re.sub(r'[()±]','',s)
    return re.sub(r'\s+',' ',s).strip()

ROOT = Path(__file__).resolve().parent.parent
DELIV = ROOT/'data'/'deliverable_scores_capped.csv'   # isobar-capped confidence (apply_isobar_cap.py)
CAND  = ROOT/'data'/'candidate_scores_v2.csv'
LIBPK = ROOT/'data'/'library_peaks_cache.json'
WORK  = ROOT/'data'/'missed_annotations_worklist.csv'
CUR_N = ROOT/'data'/'Orbitrap_HILIC_negESI_curated_042126.csv'
CUR_P = ROOT/'data'/'Orbitrap_HILIC_posESI_curated_042126.csv'
QP_CUR= ROOT/'data'/'query_peaks_cache_v2.json'
QP_POS= ROOT/'data'/'blanks_query_peaks_cache_pos.json'
QP_NEG= ROOT/'data'/'blanks_query_peaks_cache_neg.json'
META  = ROOT/'data'/'review_metadata.csv'      # per-sample species/organ/instrument (Phase-2 pull)
OUT   = ROOT/'reports'/'curator_review_app.html'

N_SUSPICIOUS = 60
SUSP_GAP = 0.05       # top-1 vs best DIFFERENT-compound candidate (library-native confusability)
SUSP_MIN_SIM = 0.50   # require a real top-1 match
TOPK_PEAKS = 40


def rd(p): return list(csv.DictReader(open(p)))
def fnum(x, d=None):
    try: return float(x)
    except (TypeError, ValueError): return d
def topk(arr):
    if not arr: return []
    a = sorted(arr, key=lambda p: -p[1])[:TOPK_PEAKS]
    m = max((p[1] for p in a), default=1) or 1
    return [[round(float(mz),4), round(float(it)/m,4)] for mz,it in a]

def load_metadata():
    """splash -> {matrix(dominant species·organ), n_matrix, instrument, n_inst, rep_sample(max-intensity), n_samples}.
    Reads review_metadata.csv + any review_metadata_*.csv supplements (later pulls for queue changes)."""
    files=sorted((ROOT/'data').glob('review_metadata*.csv'))
    if not files: return {}
    agg=defaultdict(lambda:{'mx':Counter(),'inst':Counter(),'rep':(None,-1.0),'n':0})
    rows=(r for f in files for r in csv.DictReader(open(f)))
    for r in rows:
        d=agg[r['splash']]; d['n']+=1
        sps=(r.get('species') or '').strip().lower(); org=(r.get('organ') or '').strip().lower()
        if sps or org: d['mx'][f"{sps or '?'} · {org or '?'}"]+=1
        inst=(r.get('instrument') or '').strip()
        if inst: d['inst'][re.sub(r'\s*#?\d+$','',inst)]+=1   # strip slot numbers
        it=fnum(r.get('intensity'),-1.0)
        if it>d['rep'][1]: d['rep']=(r.get('sample'),it)
    out={}
    for sp,d in agg.items():
        mx=d['mx'].most_common(); inst=d['inst'].most_common()
        out[sp]={'matrix':mx[0][0] if mx else '?','n_matrix':len(mx),
                 'instrument':inst[0][0] if inst else '?','n_inst':len(inst),
                 'rep_sample':d['rep'][0],'n_samples':d['n']}
    return out


def main():
    deliv = {r['wiki_id']: r for r in rd(DELIV)}
    meta = {}
    for f, pol in [(CUR_N,'neg'), (CUR_P,'pos')]:
        for r in rd(f):
            w=r['wiki_id']
            if w not in meta:
                meta[w]={'splash': r.get('raw_splash',''), 'polarity': pol}
    cand_by_wiki = {}
    for r in rd(CAND):
        cand_by_wiki.setdefault(r['wiki_id'], []).append(r)
    qp = {}
    for f in (QP_CUR, QP_POS, QP_NEG):
        for k,v in json.load(open(f)).items():
            if k not in qp and v: qp[k]=v
    libpk = json.load(open(LIBPK))   # MassWiki reference cache (76% faithful; 30% sparse → flagged)
    bio = load_metadata()            # splash -> dominant matrix/instrument/rep-sample/spread

    # ---------- SUSPICIOUS: library-native isobaric confusability ----------
    susp=[]
    for w,d in deliv.items():
        rows=[r for r in cand_by_wiki.get(w,[]) if r.get('hit_ik14')]
        if len(rows)<2: continue
        rows=sorted(rows,key=lambda r:-fnum(r['entropy_similarity'],0))
        top=rows[0]; topsim=fnum(top['entropy_similarity'],0)
        if topsim<SUSP_MIN_SIM: continue
        tn=norm_name(top.get('name'))
        comp=next((r for r in rows[1:] if r['hit_ik14']!=top['hit_ik14'] and norm_name(r.get('name'))!=tn), None)
        if not comp: continue
        gap=topsim-fnum(comp['entropy_similarity'],0)
        if gap>SUSP_GAP: continue
        verdict=d.get('isobar_rt_verdict','')
        if verdict=='rt_supports': continue   # RT breaks the tie toward the annotation → not suspicious
        conf=fnum(d.get('confidence'),0)
        # competing candidates = distinct compounds (dedup by ik14), top few
        seen=set(); comps=[]
        for r in rows:
            ik=r['hit_ik14']
            if ik in seen: continue
            seen.add(ik)
            comps.append({'name':r.get('name'),'esim':round(fnum(r['entropy_similarity'],0),3),
                          'adduct':r.get('adduct')})
            if len(comps)>=4: break
        boost=2.0 if verdict=='rt_conflict' else 1.0
        susp.append((conf*(1-gap)*boost, w, d, conf, topsim, gap, comps, top.get('library_wiki_id'), verdict))
    susp.sort(reverse=True, key=lambda x:x[0])
    susp=susp[:N_SUSPICIOUS]
    print(f'suspicious queue: {len(susp)}  (rt_conflict {sum(1 for x in susp if x[8]=="rt_conflict")})')

    suspicious=[]
    for prio,w,d,conf,topsim,gap,comps,top_lib,verdict in susp:
        top_name=comps[0]['name']; comp_name=comps[1]['name']; comp_sim=comps[1]['esim']
        capped=str(d.get('isobar_capped')) in ('True','rt_conflict')
        rt_t=fnum(d.get('isobar_rt_top')); rt_c=fnum(d.get('isobar_rt_comp'))
        cap_pct=round(fnum(d.get('confidence_capped'),conf)*100,1)
        raw_pct=round(conf*100,1)
        lib_raw=libpk.get(top_lib) if top_lib else None
        lib_sparse = (not lib_raw) or len(lib_raw)<=5
        suspicious.append({
            'wiki_id':w,'annotation':d.get('annotation'),'adduct':d.get('adduct'),
            'confidence_pct':cap_pct,'raw_pct':raw_pct,'capped':capped,
            'ensemble_sd_pct':round((fnum(d.get('ensemble_sd'),0) or 0)*100,1),
            'top_sim':round(topsim,3),'gap':round(gap,3),
            'adduct_cat':d.get('hit_adduct_cat'),
            'splash':meta.get(w,{}).get('splash',''),'polarity':meta.get(w,{}).get('polarity',''),
            'bio':bio.get(meta.get(w,{}).get('splash',''),{}),
            'competitors':comps,'rt_verdict':verdict,
            'rt_top':rt_t,'rt_comp':rt_c,
            'reasons':[
              f'MS² barely separates two DIFFERENT compounds at this precursor: '
              f'<b>{top_name}</b> ({comps[0]["esim"]}) vs <b>{comp_name}</b> ({comp_sim}) — Δ only {gap:.3f}.',
              (f'Model raw confidence {raw_pct}% (adduct + RT driven, NOT spectral isomer discrimination) '
               f'— <b>capped to {cap_pct}%</b>.') if capped else f'Confidence {cap_pct}%.',
              ] + ([
              f'⚠ <b>Retention time fits the COMPETITOR better</b> (observed |ΔRT| {rt_t}s for {top_name} '
              f'vs {rt_c}s for {comp_name}) — the annotation may be misassigned to the wrong isomer.'
              ] if verdict=='rt_conflict' else
              [f'Retention time cannot break the tie (|ΔRT| {rt_t}s vs {rt_c}s) — likely true co-eluting isomers.']
              if verdict in ('rt_inconclusive','rt_missing') and rt_t is not None else []),
            'query_peaks':topk(qp.get(w)),
            'lib_peaks':[] if lib_sparse else topk(lib_raw),
            'lib_name':top_name,'lib_sparse':lib_sparse,
        })

    # ---------- PROMISING: missed-annotation worklist ----------
    promising=[]
    for r in rd(WORK):
        w=r['wiki_id']; refw=r.get('ref_wiki_id')
        promising.append({
            'wiki_id':w,'tier':r.get('tier'),'proposed_identity':r.get('proposed_identity'),
            'spectral_sim':fnum(r.get('spectral_sim')),'delta_rt_to_ref':fnum(r.get('delta_rt_to_ref')),
            'lipid_class_only':r.get('lipid_class_only'),'generator_name':r.get('generator_name'),
            'splash':r.get('splash',''),'polarity':r.get('polarity',''),
            'bio':bio.get(r.get('splash',''),{}),
            'reasons':[f"Spectral match sim {fnum(r.get('spectral_sim'),0):.2f} to a CONFIRMED "
                       f"{r.get('proposed_identity')} (Δrt {fnum(r.get('delta_rt_to_ref'),0):+.0f}s), "
                       f"but never annotated."
                       + (' — LIPID class-only (acyl isomer uncertain).' if str(r.get('lipid_class_only')).lower()=='true' else '')],
            'query_peaks':topk(qp.get(w)),
            'lib_peaks':topk(qp.get(refw)) if refw else [],   # confirmed reference bin (real, trustworthy)
            'lib_name':f"confirmed {r.get('proposed_identity')}",
        })

    data={'suspicious':suspicious,'promising':promising,'generated':'2026-06-04',
          'n_suspicious':len(suspicious),'n_promising':len(promising)}
    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(HTML.replace('/*DATA*/', json.dumps(data)))
    print(f'Wrote {OUT}  (suspicious {len(suspicious)}, promising {len(promising)})')


HTML = r"""<!DOCTYPE html><html><head><meta charset="utf-8"><title>Curator Review</title>
<style>
 body{font-family:-apple-system,Segoe UI,Arial;margin:0;background:#f4f5f7;color:#1a1a1a}
 header{background:#222;color:#fff;padding:10px 18px;display:flex;gap:14px;align-items:center}
 .tab{padding:6px 14px;border-radius:6px;cursor:pointer;background:#444}
 .tab.active{background:#2d7}.tab.active.susp{background:#e85}
 .wrap{max-width:780px;margin:18px auto;padding:0 12px}
 .card{background:#fff;border-radius:10px;box-shadow:0 1px 4px #0002;padding:16px;margin-bottom:14px}
 .hd{display:flex;justify-content:space-between;align-items:baseline}
 .name{font-size:19px;font-weight:600}.conf{font-size:22px;font-weight:700}
 .sub{color:#666;font-size:13px;margin:2px 0 8px}
 .reasons{background:#fff6e5;border-left:4px solid #e85;padding:8px 10px;border-radius:4px;font-size:14px;margin:8px 0;line-height:1.5}
 .promising .reasons{background:#eaf7ee;border-left-color:#2d7}
 table.cand{width:100%;border-collapse:collapse;font-size:13px;margin:8px 0}
 table.cand td,table.cand th{text-align:left;padding:3px 6px;border-bottom:1px solid #eee}
 table.cand tr.top{font-weight:600;background:#fbfbe8}
 .meta{display:flex;flex-wrap:wrap;gap:6px 16px;font-size:12px;color:#555;margin:8px 0}
 .meta b{color:#222}.pending{color:#aaa;font-style:italic}
 svg{background:#fafafa;border:1px solid #eee;border-radius:6px;width:100%}
 .pendbox{background:#fafafa;border:1px dashed #ccc;border-radius:6px;padding:14px;text-align:center;color:#999;font-size:13px}
 .btns{display:flex;gap:12px;margin-top:12px}
 button.v{flex:1;padding:12px;border:0;border-radius:8px;font-size:16px;font-weight:600;cursor:pointer}
 .yes{background:#2d7;color:#fff}.no{background:#bbb;color:#fff}
 .done{opacity:.5}.chosen{outline:3px solid #1a1a1a}
 .bar{position:sticky;top:0;background:#fff;padding:8px 12px;border-bottom:1px solid #ddd;font-size:13px;display:flex;justify-content:space-between;z-index:9}
 .exp{background:#2d7;color:#fff;border:0;padding:6px 12px;border-radius:6px;cursor:pointer}
</style></head><body>
<header><b>Curator Review</b>
 <span class="tab susp active" onclick="show('suspicious')" id="t_s"></span>
 <span class="tab" onclick="show('promising')" id="t_p"></span>
 <span style="flex:1"></span><button class="exp" onclick="exportCSV()">Export decisions</button>
</header>
<div class="bar"><span id="prog"></span><span>+ / − keys · autosaves</span></div>
<div class="wrap" id="wrap"></div>
<script>
const DATA=/*DATA*/;
const KEY='curator_decisions_v2';
let dec=JSON.parse(localStorage.getItem(KEY)||'{}');
let cur='suspicious';
const SEM={suspicious:{yes:'Real problem (reject)',no:'Fine (keep)'},promising:{yes:'Annotate it',no:'No'}};
function mirror(q,l,pendMsg){
 if((!l||!l.length)&&pendMsg) return `<div class="pendbox">query MS² shown below · reference: ${pendMsg}</div>`+onlyTop(q);
 const W=740,H=150,mid=H/2; let mx=1; q.concat(l).forEach(p=>mx=Math.max(mx,p[0]));
 const X=mz=>40+(mz/mx)*(W-50); let s=`<svg viewBox="0 0 ${W} ${H}">`;
 s+=`<line x1="40" y1="${mid}" x2="${W}" y2="${mid}" stroke="#999"/>`;
 q.forEach(p=>{s+=`<line x1="${X(p[0]).toFixed(1)}" y1="${mid}" x2="${X(p[0]).toFixed(1)}" y2="${(mid-p[1]*(mid-8)).toFixed(1)}" stroke="#2563eb"/>`;});
 l.forEach(p=>{s+=`<line x1="${X(p[0]).toFixed(1)}" y1="${mid}" x2="${X(p[0]).toFixed(1)}" y2="${(mid+p[1]*(mid-8)).toFixed(1)}" stroke="#e85d2a"/>`;});
 s+=`<text x="44" y="12" font-size="10" fill="#2563eb">query (top)</text>`;
 s+=`<text x="44" y="${H-4}" font-size="10" fill="#e85d2a">reference (bottom)</text></svg>`; return s;
}
function onlyTop(q){
 const W=740,H=90; let mx=1; q.forEach(p=>mx=Math.max(mx,p[0])); const X=mz=>40+(mz/mx)*(W-50);
 let s=`<svg viewBox="0 0 ${W} ${H}"><line x1="40" y1="${H-12}" x2="${W}" y2="${H-12}" stroke="#999"/>`;
 q.forEach(p=>{s+=`<line x1="${X(p[0]).toFixed(1)}" y1="${H-12}" x2="${X(p[0]).toFixed(1)}" y2="${(H-12-p[1]*(H-24)).toFixed(1)}" stroke="#2563eb"/>`;});
 return s+`<text x="44" y="12" font-size="10" fill="#2563eb">query MS²</text></svg>`;
}
function candTable(c){
 let s='<table class="cand"><tr><th>candidate (MassWiki search)</th><th>adduct</th><th>esim</th></tr>';
 c.forEach((r,i)=>{s+=`<tr class="${i==0?'top':''}"><td>${r.name||'?'}</td><td>${r.adduct||''}</td><td>${r.esim}</td></tr>`;});
 return s+'</table>';
}
function card(it,kind){
 const id=it.wiki_id, v=dec[id]&&dec[id].verdict;
 const b=it.bio||{};
 const matrixS=b.matrix?`<span>matrix <b>${b.matrix}</b>${b.n_matrix>1?` (+${b.n_matrix-1} more)`:''}</span>`:'';
 const instS=b.instrument?`<span>instrument <b>${b.instrument}</b>${b.n_inst>1?` (+${b.n_inst-1})`:''}</span>`:'';
 const nS=b.n_samples?`<span>${b.n_samples} samples</span>`:'';
 const meta=`<div class="meta"><span>SPLASH <b>${(it.splash||'?').slice(0,32)}</b></span>`+
   `<span>mode <b>${it.polarity||'?'}</b></span>`+(it.adduct?`<span>adduct <b>${it.adduct}</b></span>`:'')+
   matrixS+instS+nS+`<span class="pending">EIC: pending FlashEIC</span></div>`;
 if(kind==='suspicious'){
   const confHtml=`${it.confidence_pct}%`+(it.capped?` <span style="font-size:12px;color:#e85;font-weight:600">(capped from ${it.raw_pct}%)</span>`:'');
   return `<div class="card suspicious ${v?'done':''}" id="c_${id}">
     <div class="hd"><span class="name">${it.annotation||'(no name)'}</span><span class="conf">${confHtml}</span></div>
     <div class="sub">top MS² sim ${it.top_sim} · gap to next compound ${it.gap} · sd ±${it.ensemble_sd_pct}%</div>
     <div class="reasons">${it.reasons.join('<br><br>')}</div>
     ${candTable(it.competitors)}
     ${mirror(it.query_peaks,it.lib_peaks, it.lib_sparse?'reference sparse in local cache — confirm in MassWiki':'MassWiki reference')}
     ${meta}
     <div class="btns">
       <button class="v yes ${v==='yes'?'chosen':''}" onclick="vote('${id}','yes')">+ ${SEM[kind].yes}</button>
       <button class="v no ${v==='no'?'chosen':''}" onclick="vote('${id}','no')">− ${SEM[kind].no}</button></div></div>`;
 }
 return `<div class="card promising ${v?'done':''}" id="c_${id}">
   <div class="hd"><span class="name">${it.proposed_identity||'(no name)'}</span><span class="sub">tier ${it.tier}</span></div>
   <div class="sub">spectral sim ${it.spectral_sim?.toFixed?.(2)} · ΔRT ${it.delta_rt_to_ref?.toFixed?.(0)}s · generator said "${it.generator_name||'?'}"</div>
   <div class="reasons">${it.reasons.join('<br>')}</div>
   ${mirror(it.query_peaks,it.lib_peaks,null)}
   ${meta}
   <div class="btns">
     <button class="v yes ${v==='yes'?'chosen':''}" onclick="vote('${id}','yes')">+ ${SEM[kind].yes}</button>
     <button class="v no ${v==='no'?'chosen':''}" onclick="vote('${id}','no')">− ${SEM[kind].no}</button></div></div>`;
}
function render(){
 const items=DATA[cur];
 document.getElementById('wrap').innerHTML=items.map(it=>card(it,cur)).join('');
 const done=items.filter(it=>dec[it.wiki_id]).length;
 document.getElementById('prog').textContent=`${cur}: ${done}/${items.length} reviewed`;
}
function vote(id,v){dec[id]={verdict:v,queue:cur,ts:Date.now()};localStorage.setItem(KEY,JSON.stringify(dec));render();}
function show(q){cur=q;document.getElementById('t_s').className='tab susp'+(q==='suspicious'?' active':'');
 document.getElementById('t_p').className='tab'+(q==='promising'?' active':'');render();}
function exportCSV(){let rows=[['wiki_id','queue','verdict','timestamp']];
 for(const[id,o]of Object.entries(dec))rows.push([id,o.queue,o.verdict,new Date(o.ts).toISOString()]);
 const blob=new Blob([rows.map(r=>r.join(',')).join('\n')],{type:'text/csv'});
 const a=document.createElement('a');a.href=URL.createObjectURL(blob);a.download='curator_decisions.csv';a.click();}
document.addEventListener('keydown',e=>{if(!'+-=_'.includes(e.key))return;
 const next=DATA[cur].find(it=>!dec[it.wiki_id]);
 if(next)vote(next.wiki_id,(e.key==='+'||e.key==='=')?'yes':'no');});
document.getElementById('t_s').textContent='Suspicious ('+DATA.n_suspicious+')';
document.getElementById('t_p').textContent='Promising ('+DATA.n_promising+')';
render();
</script></body></html>"""

if __name__ == '__main__':
    main()
