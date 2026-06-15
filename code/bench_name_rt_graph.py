"""bench_name_rt_graph.py — Same-name / RT-disagreement graph (KG prototype layer 1).

First layer of the knowledge-graph pivot (project_knowledge_graph_pivot). Generalizes
the collision-group triage in bench_joint_inference_flags.py in three ways:

  1. Re-keys from (hit_ik14, adduct) of the top-1 call to the *claimed identity name*
     (top-1 candidate name, normalized). Oliver curates by NAME — RT-disagreement-
     within-name is the actual labeling rule (project_oliver_curation_method).
  2. Extends to UNCURATED bins (drops the is_manual_annotated filter). This is where
     NOVEL findings live — curated-only can only re-confirm what Oliver already saw.
  3. Builds an explicit graph (node = bin, edge = same-name + co-eluting) so later
     edge types (adduct / ISF / cross-polarity) and propagation plug in. networkx is
     not installed; we use a portable edge-list + scipy-free single-linkage on the RT
     axis. Swap to networkx/igraph later by loading data/name_rt_graph_edges.csv.

Graph / metric definition
-------------------------
Within each name-clique (all bins sharing a normalized claimed identity), bins are
clustered along the RT axis by SINGLE LINKAGE at RT_TIGHT (30 s): sort by RT, start a
new RT cluster whenever the gap to the previous same-name bin exceeds 30 s. (Single-
linkage threshold, not GMM — feedback_gmm_failed: clustering gives a binary split on
pre-filtered data. A hard 30 s gap is the lab's "same chromatographic peak" tolerance.)

  - coherent compound  -> one RT cluster spanning the clique
  - RT disagreement    -> ≥2 RT clusters; bins outside the dominant cluster are outliers

Per bin: rt_cluster_size, n_rt_clusters_in_clique, in_dominant_rt_cluster,
rt_coherence = rt_cluster_size / name_clique_size, rt_outlier (in a disagreement clique
and not in the dominant cluster).

CIRCULARITY CAVEAT
------------------
This REPRODUCES the labeling rule, so we do NOT report "AUC vs yy_" as an independent
test (project_oliver_curation_method: RT dominance is partly a labeling artifact).
Success is two descriptive numbers:
  SANITY — fraction of known yy_/FP bins this re-derives as RT outliers (should be high).
  VALUE  — count of NOVEL RT-disagreement flags among bins the curator has NOT rejected,
           split curated-not-yet-yy vs uncurated. Per project_curator_triage_v1, expect
           most flags to already be known yy_; the novel slice is the deliverable.

Outputs
  data/bench_name_rt_graph.csv   — per-bin metrics on all bins
  data/name_rt_review.csv        — RT-outlier shortlist (conf ≥ 0.5), novelty-tagged
  data/name_rt_graph_edges.csv   — same-name edge list (portable graph object)
"""
import re
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
CAND = ROOT / 'data' / 'candidate_scores_v2.csv'
SPEC_NEG = ROOT / 'data' / 'Orbitrap_HILIC_negESI_curated_042126.csv'
SPEC_POS = ROOT / 'data' / 'Orbitrap_HILIC_posESI_curated_042126.csv'
OUT_BINS = ROOT / 'data' / 'bench_name_rt_graph.csv'
OUT_REVIEW = ROOT / 'data' / 'name_rt_review.csv'
OUT_EDGES = ROOT / 'data' / 'name_rt_graph_edges.csv'

RT_TIGHT = 30.0   # seconds — "same chromatographic peak" tolerance (shared w/ bench_joint_inference_flags)
RT_STRONG = 60.0  # seconds — strong incompatibility
MIN_CONF_FOR_REVIEW = 0.5


def normalize_name(name) -> str:
    """Claimed-identity key. Lowercase, collapse whitespace, strip a stray yy_ prefix.

    Grouping is by the CLAIMED identity (candidate top-1 name), independent of the
    curator's verdict — so a legit 'glucose' bin and a 'yy_glucose' rejection both
    land in the same clique and the rejection surfaces as the RT outlier.
    """
    if not isinstance(name, str):
        return ''
    s = name.strip().lower()
    s = re.sub(r'^yy_', '', s)
    s = re.sub(r'\s+', ' ', s)
    return s


def main():
    print(f'Loading {CAND.name}')
    cand = pd.read_csv(CAND)
    cand = cand[cand['gbm_cal'].notna()].copy()

    # GBM top-1 per bin = the claimed identity for that bin.
    cand_sorted = cand.sort_values(['wiki_id', 'gbm_cal'], ascending=[True, False])
    top1 = cand_sorted.groupby('wiki_id', as_index=False).first()
    print(f'  bins (top-1 rows): {len(top1):,}')

    # Observed RT (seconds), SPLASH, curator name + annotation status from curated tables.
    # curated `name` carries the curator verdict (yy_ prefix on rejections); candidate
    # `name` is the claimed identity. SPLASH is the stable id (reference_splash_vs_wiki_id).
    print(f'Loading curated metadata from {SPEC_NEG.name} + {SPEC_POS.name}')
    cur_cols = ['wiki_id', 'rt', 'raw_splash', 'is_manual_annotated', 'name']
    neg = pd.read_csv(SPEC_NEG, usecols=cur_cols, low_memory=False); neg['polarity'] = 'neg'
    pos = pd.read_csv(SPEC_POS, usecols=cur_cols, low_memory=False); pos['polarity'] = 'pos'
    spec = pd.concat([neg, pos], ignore_index=True).drop_duplicates(subset='wiki_id', keep='first')
    spec = spec.rename(columns={'raw_splash': 'splash', 'name': 'curator_name'})

    df = top1.merge(
        spec[['wiki_id', 'rt', 'polarity', 'splash', 'is_manual_annotated', 'curator_name']],
        on='wiki_id', how='left',
    )

    # Identity key = what the bin is CALLED. For curated bins use the curator's
    # assignment (curator_name) — that's the identity Oliver actually groups by; the
    # candidate top-1 can disagree (e.g. claims "pyruvic acid" while curator said
    # "guanidine"), and grouping on the guess pollutes the clique. Uncurated bins have
    # no curator name, so fall back to the candidate top-1 name to let them participate.
    cand_id = df['name'].apply(normalize_name)
    cur_id = df['curator_name'].apply(normalize_name)
    use_curator = df['is_manual_annotated'].fillna(False) & (cur_id != '')
    df['norm_name'] = np.where(use_curator, cur_id, cand_id)
    df['identity_source'] = np.where(use_curator, 'curator', 'candidate')

    # A bin can only join the name graph if it has an observed RT and a claimed identity.
    n_before = len(df)
    df = df[df['rt'].notna() & (df['norm_name'] != '')].copy()
    print(f'  bins with RT + claimed identity: {len(df):,}  (dropped {n_before - len(df):,})')

    # ---- Single-linkage RT clustering within each name-clique -------------------
    df = df.sort_values(['norm_name', 'rt']).reset_index(drop=True)
    name_change = df['norm_name'] != df['norm_name'].shift()
    rt_gap = df['rt'] - df['rt'].shift()
    new_cluster = name_change | (rt_gap > RT_TIGHT)
    df['rt_cluster_id'] = new_cluster.cumsum().astype(int)

    cl_size = df.groupby('rt_cluster_id')['wiki_id'].transform('size')
    df['rt_cluster_size'] = cl_size
    df['name_clique_size'] = df.groupby('norm_name')['wiki_id'].transform('size')
    df['n_rt_clusters_in_clique'] = df.groupby('norm_name')['rt_cluster_id'].transform('nunique')

    # Dominant RT cluster per clique = highest-CONFIDENCE cluster, tie-break by size.
    # NOT by size first: a real high-conf call is often outnumbered by low-conf junk
    # bins sharing the name, and size-first would flag the real call as the outlier.
    # The compound's true RT is where its confident evidence sits.
    cl = (df.groupby(['norm_name', 'rt_cluster_id'])
            .agg(sz=('wiki_id', 'size'), maxconf=('gbm_cal', 'max'))
            .reset_index()
            .sort_values(['norm_name', 'maxconf', 'sz'], ascending=[True, False, False]))
    dom = (cl.drop_duplicates('norm_name')[['norm_name', 'rt_cluster_id']]
             .rename(columns={'rt_cluster_id': 'dominant_cluster_id'}))
    df = df.merge(dom, on='norm_name', how='left')
    df['in_dominant_rt_cluster'] = df['rt_cluster_id'] == df['dominant_cluster_id']

    df['rt_coherence'] = df['rt_cluster_size'] / df['name_clique_size']
    df['rt_outlier'] = (df['n_rt_clusters_in_clique'] >= 2) & (~df['in_dominant_rt_cluster'])

    # ---- Verdict / novelty tagging ----------------------------------------------
    df['is_curated'] = df['is_manual_annotated'].fillna(False).astype(bool)
    already_yy = df['curator_name'].fillna('').str.lower().str.startswith('yy_')
    is_fp = df['hit_label'] == 0
    df['known_problem'] = already_yy | is_fp   # curator already rejected this call
    df['novel_flag'] = (
        df['rt_outlier']
        & (df['gbm_cal'] >= MIN_CONF_FOR_REVIEW)
        & (~df['known_problem'])
    )

    # ---- Graph edge list (portable: node=bin, edge=same-name consecutive in RT) --
    # Path approximation of the same-name clique along the RT axis: rt_compatible
    # edges (gap ≤ 30 s) are exactly the single-linkage cluster bonds; the broken
    # edges (gap > 30 s) are the disagreement boundaries. Sufficient scaffold for
    # RT clustering and for hanging adduct/ISF/cross-polarity edges off later.
    same = ~name_change.values
    src = df['wiki_id'].shift().values
    edges = pd.DataFrame({
        'source': src, 'target': df['wiki_id'].values,
        'norm_name': df['norm_name'].values,
        'delta_rt': rt_gap.values, 'same_name': same,
    })
    edges = edges[edges['same_name']].copy()
    edges['edge_type'] = 'same_name'
    edges['rt_compatible'] = edges['delta_rt'].abs() <= RT_TIGHT
    edges = edges[['source', 'target', 'edge_type', 'delta_rt', 'rt_compatible', 'norm_name']]

    # ---- Console: graph stats ---------------------------------------------------
    n_cliques = df['norm_name'].nunique()
    multi = df.drop_duplicates('norm_name')
    n_disagree = int((multi['n_rt_clusters_in_clique'] >= 2).sum())
    print(f'\n=== Graph ===')
    print(f'  bins (nodes):                 {len(df):,}')
    print(f'  name cliques:                 {n_cliques:,}')
    print(f'  disagreement cliques (≥2 RT clusters): {n_disagree:,}')
    print(f'  same-name edges:              {len(edges):,}  '
          f'({int((~edges["rt_compatible"]).sum()):,} RT-incompatible)')
    print(f'  RT-outlier bins:              {int(df["rt_outlier"].sum()):,}')

    # ---- Success number 1: SANITY (re-derive known yy_/FP) ----------------------
    print('\n=== SANITY (reproduces the labeling rule — NOT an independent test) ===')
    for label, mask in [('curator yy_ bins', already_yy), ('all FP (hit_label==0)', is_fp)]:
        m = mask & df['rt'].notna()
        n = int(m.sum())
        if n:
            caught = int((df.loc[m, 'rt_outlier']).sum())
            in_disagree = int((df.loc[m, 'n_rt_clusters_in_clique'] >= 2).sum())
            print(f'  {label}: {n:,} | RT-outlier {caught:,} ({100*caught/n:.0f}%) | '
                  f'in a disagreement clique {in_disagree:,} ({100*in_disagree/n:.0f}%)')
        else:
            print(f'  {label}: 0 present')

    # ---- Success number 2: VALUE (novel flags) ----------------------------------
    nov = df[df['novel_flag']]
    print('\n=== VALUE (novel RT-disagreement flags, conf ≥ 0.5, not curator-rejected) ===')
    print(f'  total novel:                  {len(nov):,}')
    print(f'    curated, not yet yy_:       {int(nov["is_curated"].sum()):,}')
    print(f'    uncurated (never reviewed): {int((~nov["is_curated"]).sum()):,}')

    # ---- Write artifacts --------------------------------------------------------
    bin_cols = [
        'wiki_id', 'splash', 'polarity', 'rt', 'is_curated',
        'name', 'norm_name', 'curator_name', 'hit_ik14', 'adduct', 'hit_label',
        'gbm_cal', 'sim_gap', 'entropy_similarity',
        'name_clique_size', 'rt_cluster_id', 'rt_cluster_size',
        'n_rt_clusters_in_clique', 'in_dominant_rt_cluster',
        'rt_coherence', 'rt_outlier', 'known_problem', 'novel_flag',
    ]
    bin_cols = [c for c in bin_cols if c in df.columns]
    df[bin_cols].to_csv(OUT_BINS, index=False)
    print(f'\nWrote {OUT_BINS} ({len(df):,} rows)')

    review = df[df['rt_outlier'] & (df['gbm_cal'] >= MIN_CONF_FOR_REVIEW)].copy()
    review = review.sort_values(
        ['novel_flag', 'is_curated', 'gbm_cal'], ascending=[False, True, False])
    review[bin_cols].to_csv(OUT_REVIEW, index=False)
    print(f'Wrote {OUT_REVIEW} ({len(review):,} rows; {int(review["novel_flag"].sum()):,} novel)')

    edges.to_csv(OUT_EDGES, index=False)
    print(f'Wrote {OUT_EDGES} ({len(edges):,} edges)')

    # ---- Console: a few novel disagreement cliques ------------------------------
    print('\n=== Top novel disagreement cliques (preview) ===')
    novel_names = (nov.sort_values('gbm_cal', ascending=False)
                      .drop_duplicates('norm_name')['norm_name'].head(5))
    for nm in novel_names:
        clq = df[df['norm_name'] == nm].sort_values('rt')
        print(f'\nname="{nm}"  clique={len(clq)}  RT clusters={clq["rt_cluster_id"].nunique()}')
        for _, b in clq.head(8).iterrows():
            tag = 'OUTLIER' if b['rt_outlier'] else ('dom' if b['in_dominant_rt_cluster'] else '')
            nv = 'NOVEL' if b['novel_flag'] else ('yy_/FP' if b['known_problem'] else '')
            print(f'   wiki_id={str(b["wiki_id"]):<16} rt={b["rt"]:7.1f}s conf={b["gbm_cal"]:.3f} '
                  f'pol={b["polarity"]} cur="{str(b["curator_name"])[:24]}" {tag:7} {nv}')


if __name__ == '__main__':
    main()
