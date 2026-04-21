
# mvp_run_example.py
# Example glue (edit paths and run in your environment)
import pandas as pd
from mvp_ingest import load_spectra, load_hits, join_hits_to_spectra, build_confusable_sets, collapse_duplicates_by_structure
from mvp_score import fit_channel_params, per_hit_logLR, aggregate_to_structure, local_rank_probability, compute_global_posterior
from mvp_export import make_assertion_table, make_top_calls

# --- Edit these paths ---
SPECTRA_PATH = "spectrum.parquet"
HITS_PATH = "combined_df.parquet"

spectra = load_spectra(SPECTRA_PATH)
hits = load_hits(HITS_PATH)

joint = join_hits_to_spectra(hits, spectra)
cset = build_confusable_sets(joint, ms2_min=0.75, zrt_k=3.0, rt_sigma_sec=2.0)

# Build a quick labeled subset for parameter fitting if available
# Expect columns: is_manual_annotated (True for TP spectra)
# For MVP, create labels at hit level by inheriting spectrum's label (crude but workable)
labeled = cset.dropna(subset=["entropy_similarity"]).copy()
labeled["label"] = labeled["is_manual_annotated"].astype(int)

params = fit_channel_params(labeled)

cset["logLR_hit"] = per_hit_logLR(cset, params)

struct_scores = aggregate_to_structure(cset, quality_cols=("entropy_similarity",), topK=None, lambda_var=0.0)
struct_scores = local_rank_probability(struct_scores)

post = compute_global_posterior(struct_scores, prior_policy="uniform", lr_nota=0.5, prior_nota=0.5)

assertions = make_assertion_table(post, features_df=cset)
top_calls = make_top_calls(assertions, nota_threshold=0.6)

assertions.to_parquet("assertions.parquet")
top_calls.to_csv("top_calls.csv", index=False)
print("Wrote assertions.parquet and top_calls.csv")
