"""bench_relational.py — Do GNPS-style relational spectral-network features help the GBM?

Edge 3 of the knowledge-graph pivot (project_knowledge_graph_pivot) as scoring input,
not as a worklist. Bare NN cosine is redundant with entropy_similarity
(project_learned_similarity_mvp); these are the *relational* columns — neighbourhood
purity / density / claim-corroboration — computed leakage-safe (neighbours blocked from
sharing anno_ik14). Expect the win in triage/disagreement, not AUC.
"""
from bench_harness import BenchSpec, run_bench

SPEC = BenchSpec(
    name='relational',
    # rel_claim_corrob DROPPED — structural leakage: the same-anno_ik14 exclusion makes
    # it fire only on FPs (a TP's claimed IK14 == its own excluded anno_ik14), so
    # corrob==1 => 100% FP (145/145). It inflated the delta to +0.0133 (label readout).
    # rel_nbr_purity DROPPED — auditor drop-one: -0.0009 (84% NaN, ~0 own-corr, noise).
    # Honest arm = the 2 clean confusability features (auditor: +0.0070, dominates 3-feat).
    extra_numeric=['rel_nn_sim', 'rel_n_confirmed_nbr'],
    extra_table='data/relational_features.csv',
    merge_on=('wiki_id',),
)

if __name__ == '__main__':
    run_bench(SPEC)
