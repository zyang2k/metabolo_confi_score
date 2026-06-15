"""build_missed_annotations.py — Consolidated missed-annotation worklist for Oliver.

Objective (chosen 2026-06-03): find uncurated bins that are real, unannotated instances
of compounds the curator already knows. NOT novel discovery (those are outside the
curated manifold — project_knowledge_graph_pivot 93% ceiling, which is the CORRECT scope
boundary for this objective).

Engine = layer-3 spectral propagation (bench_spectral_propagate.py): identity derived
from the SPECTRUM (trustworthy), not the candidate generator's name (97% junk on
uncurated bins). Layer-2 RT-anchor propagation (bench_name_rt_propagate.py) is folded in
only as corroboration. Layer 1 is not used here (it was scaffolding / a diagnostic).

A bin is a missed-annotation candidate when its spectrum matches a curated TRUE-POSITIVE
bin (entropy sim ≥ 0.70) at a consistent RT (|Δrt| ≤ 30 s) — i.e. it elutes where that
compound is known to elute AND fragments like it. Precision tiers by spectral sim + RT:

  A  sim ≥ 0.90 and |Δrt| ≤ 15 s        (near-certain replicate of a known compound)
  B  sim ≥ 0.80 and |Δrt| ≤ 30 s
  C  sim ≥ 0.70 and |Δrt| ≤ 30 s        (CONFIRMED floor)

Enrichment flags:
  name_corroborated  — layer-2 also flagged this bin (claimed name + RT agree)
  pipeline_override  — proposed IK14 differs from the bin's own GBM top-1 IK14
                       (the generator named it something else → a real catch, or an
                       isobaric confusion; ~20% of sim≥0.7 matches are wrong by IK14,
                       so tier A first)

Output: data/missed_annotations_worklist.csv
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
SPEC = ROOT / 'data' / 'spectral_propagated.csv'
L2 = ROOT / 'data' / 'missed_annotation_candidates.csv'
GRAPH = ROOT / 'data' / 'bench_name_rt_graph.csv'
OUT = ROOT / 'data' / 'missed_annotations_worklist.csv'


def tier(sim, drt):
    a = abs(drt)
    if sim >= 0.90 and a <= 15:
        return 'A'
    if sim >= 0.80 and a <= 30:
        return 'B'
    return 'C'


def main():
    spec = pd.read_csv(SPEC)
    core = spec[spec['verdict'] == 'CONFIRMED'].copy()   # spectral edge + RT consistent
    print(f'spectral CONFIRMED (engine): {len(core):,}')

    # bin's own GBM top-1 IK14 + stable splash id, to detect pipeline overrides
    g = pd.read_csv(GRAPH, usecols=['wiki_id', 'splash', 'hit_ik14'])
    g = g.rename(columns={'hit_ik14': 'gen_ik14'})
    core = core.merge(g, on='wiki_id', how='left')

    # layer-2 corroboration (name-anchor agreed this bin is a missed annotation)
    l2 = set(pd.read_csv(L2)['wiki_id'])
    core['name_corroborated'] = core['wiki_id'].isin(l2)

    both_ik = (core['prop_ik14'].notna() & (core['prop_ik14'] != '')
               & core['gen_ik14'].notna() & (core['gen_ik14'] != ''))
    core['pipeline_override'] = both_ik & (core['prop_ik14'] != core['gen_ik14'])

    core['tier'] = [tier(s, d) for s, d in zip(core['spectral_sim'], core['delta_rt_to_ref'])]
    core['priority'] = core['spectral_sim'] * np.exp(-(core['delta_rt_to_ref'] ** 2) / (2 * 20.0 ** 2))

    # Lipid-class caveat: lipids share headgroup fragments, so a high entropy-similarity
    # match confirms the CLASS but not the acyl-chain isomer (project_lipid_hilic_check).
    # These are real (bin IS a lipid of that class) but the exact identity needs the
    # curator's eye — so we surface, not rank, them at the top.
    lipid = r'(?:\b(?:lpc|lpe|pc|pe|pg|pi|ps|car|acylcarnitine|tg|dg|mg|sm|cer|fa)\b|\d+:\d+)'
    core['lipid_class_only'] = core['prop_identity'].fillna('').str.contains(lipid, case=False, regex=True)
    # Sort: compound-level first within each tier (more trustworthy identity), then priority
    core = core.sort_values(['tier', 'lipid_class_only', 'priority'],
                            ascending=[True, True, False])

    out = core[[
        'tier', 'priority', 'wiki_id', 'splash', 'polarity', 'precursor_mz', 'rt',
        'prop_identity', 'prop_ik14', 'lipid_class_only', 'ref_wiki_id', 'ref_rt', 'delta_rt_to_ref',
        'spectral_sim', 'name_corroborated', 'pipeline_override',
        'name', 'gen_ik14', 'gbm_cal',
    ]].rename(columns={
        'name': 'generator_name', 'rt': 'obs_rt', 'gbm_cal': 'generator_conf',
        'prop_identity': 'proposed_identity',
    })
    out.to_csv(OUT, index=False)

    print('\n=== Missed-annotation worklist ===')
    cl = out['lipid_class_only']
    print(f'  COMPOUND-LEVEL (trustworthy identity): {int((~cl).sum()):,}')
    print(f'  LIPID-CLASS only (class ok, isomer uncertain): {int(cl.sum()):,}')
    for t in ['A', 'B', 'C']:
        sub = out[out['tier'] == t]
        print(f'  Tier {t}: {len(sub):,}  (compound-level {int((~sub["lipid_class_only"]).sum()):,}) '
              f'| name-corrob {int(sub["name_corroborated"].sum()):,} '
              f'| override {int(sub["pipeline_override"].sum()):,}')
    print(f'  TOTAL: {len(out):,}')
    print(f'  both engines agree (name_corroborated):  {int(out["name_corroborated"].sum()):,}')
    print(f'  overrides the pipeline top-1 call:        {int(out["pipeline_override"].sum()):,}')
    print(f'\nWrote {OUT}')

    print('\n=== Tier A preview (near-certain missed annotations) ===')
    for _, b in out[out['tier'] == 'A'].head(12).iterrows():
        ov = 'OVERRIDE' if b['pipeline_override'] else ''
        co = 'name✓' if b['name_corroborated'] else ''
        print(f'   {str(b["wiki_id"]):<15} sim={b["spectral_sim"]:.2f} Δrt={b["delta_rt_to_ref"]:+5.1f}s '
              f'gen_conf={b["generator_conf"]:.2f} pol={b["polarity"]} '
              f'"{str(b["generator_name"])[:18]}" -> "{str(b["proposed_identity"])[:24]}" {ov} {co}')


if __name__ == '__main__':
    main()
