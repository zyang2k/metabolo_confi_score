# Archived Notebooks

Past attempts that informed — or were superseded by — the current pipeline. Kept for reference, not part of the active workflow. If you're looking for current production code, see `code/bayesian_3channel_explorer.ipynb` (current 3-channel Bayesian, validated 2026-04-17, Orbi AUC 0.841), the `*_gbm*.ipynb` notebooks (GBM v2, AUC 0.882), and the `fdr_benchmark_*.ipynb` notebooks.

Each entry below covers: **intent → outcome → why archived / lesson learned**.

---

## Superseded by the current Bayesian + GBM pipeline

### `confidence_score_gmm.ipynb`
- **Intent:** Fit an unsupervised 2-component GMM on 7 features (identity, entropy, sim_gap, delta_rt_pred, n_adducts, delta_ppm, id_fuzzy_ratio) and use `P(confident | features)` as the score. No labels needed.
- **Outcome:** ~50/50 binary split (Very Low vs High); only ~0.6% landed in the Low/Medium middle. No smooth gradient.
- **Why it failed — two compounding reasons:**
  1. *No latent bimodality exists in the data.* Annotated spectra are already pre-filtered (they passed a threshold to be named). There's no natural "correct vs incorrect" split for GMM to find, so it fell back to splitting along the axis of maximum variance (`delta_ppm`) — an artifact, not confidence.
  2. *Wrong tool for a smooth score.* Clustering produces discrete partitions; GMM posteriors saturate to 0/1 away from the decision boundary.
- **Lesson:** For pre-filtered homogeneous data, use rank-/percentile-based or likelihood-ratio methods that produce gradation by construction. Drove the move to the LR-based Bayesian pipeline.

### `confidence_score.ipynb`
- **Intent:** Combine 4 component scores (MS2 match, uniqueness gap, RT agreement, spectrum quality) via geometric mean → confidence.
- **Outcome:** Correlates with Oliver's probability ratings but flagged ~12% as "absurd" (score < 0.3).
- **Why archived:** Ad-hoc combination has no probabilistic meaning; a geometric mean treats scores on different scales as equivalent. Replaced by the principled Bayesian LR formulation.

### `proposal_mvp.ipynb` (758 KB)
- **Intent:** March 2026 4-channel Bayesian proposal — null-adjusted MS2 LR (entropy_similarity via KDE TP/FP), mass accuracy, RT agreement, spectrum quality prior. Included detailed parameter fitting and library-gap handling.
- **Outcome:** Methodology validated on Oliver data, but channel count was unstable — RT often unavailable, and MS1 quality prior was weakly informative.
- **Why archived:** Superseded on 2026-04-17 by the cleaner 3-channel model (entropy_sim, sim_gap, delta_rt) in `bayesian_3channel_explorer.ipynb`.

### `new_entropy_sim.ipynb`
- **Intent:** Adjust entropy similarity by p-value weighting: `S' = -log10(p) × S`, where `p` comes from per-entropy-bin bootstrap null pairs at ±10 ppm.
- **Outcome:** Multiple partial implementations; never finalized or validated end-to-end.
- **Why archived:** The core idea (stratify nulls by entropy) was kept, but implemented cleanly in the current pipeline. This notebook was a working sketch.

### `entropy_vs_similarity.ipynb` (Sep 2024)
- **Intent:** Test whether simple (low-peak) spectra artificially inflate similarity. Compared three null-construction methods: (A) stratify by peak count, (B) stratify by entropy, (C) stratify by peak-matching fraction.
- **Outcome:** Entropy-stratified nulls had the best behavior (tightest FP tail for high-entropy queries).
- **Why archived:** Finding adopted into current null calibration. The exploratory notebook itself isn't needed.

### `entropy_null_exploration.ipynb` (23 MB)
- **Intent:** Foundational large-scale exploration of entropy-stratified null distributions, same-formula isomer pairs, and isomer-specific nulls on NIST23.
- **Outcome:** Produced the evidence base for stratifying by entropy and handling isomers specially.
- **Why archived:** Methods adopted into the production null calibration. File is bloated with repeated plotting cells.

### `peak_null.ipynb`
- **Intent:** Build null by ±10 ppm mass-matched random pairs stratified by entropy and peak count.
- **Outcome:** Same theme as `entropy_null_exploration.ipynb`, cleaner but narrower scope.
- **Why archived:** Subsumed into current null-distribution module.

### `peak-by-peak.ipynb`
- **Intent:** Propose a peak-by-peak weighted similarity (`weight = m/z^α × intensity^β`, tolerance 0.01 Da, missing-peak penalty 0.1) as an alternative to entropy similarity.
- **Outcome:** ROC-AUC 0.6854 (peak-by-peak) vs 0.6923 (entropy); correlation 0.7575. Methods are complementary, entropy marginally better.
- **Why archived:** Not adopted. Entropy similarity won; peak-by-peak may still inform future ensemble work but isn't in the active pipeline.

---

## Failed / abandoned experiments

### `transfer_null.ipynb`
- **Intent:** Transfer TP/FP distributions learned on NIST23 public data to in-house Orbitrap via λ-weighted mixing of public and in-house KDEs. Included EM and a λ sweep.
- **Outcome:** Explicitly marked **"不work" (doesn't work)** in the notebook. NIST and Orbitrap distributions differ too much (different instrument physics, fragmentation energies, population of compounds).
- **Lesson:** **NIST → Orbitrap transfer gap is real.** Don't assume public TP/FP distributions transfer to in-house data without re-calibration. This motivates the current approach of fitting TP/FP distributions *on the target platform*.

### `em_soft_label.ipynb`
- **Intent:** Soft-label EM refining KDE-based scores via iterative E/M on annotated + unannotated spectra. Tracked KDE convergence, per-iteration calibration, and `mix_alpha` sensitivity.
- **Outcome:** EM converged, but high-similarity FP mass remained stubbornly high after iteration (KDE couldn't pull those FPs down).
- **Why archived:** Superseded by `tn_initialized_em.ipynb`, which itself didn't beat the 3-channel Bayesian.

### `tn_initialized_em.ipynb`
- **Intent:** Fix the naive-EM FP bias by initializing with 342 curator-confirmed hard FPs (Oliver's "yy_" tier, scoring 0.59–1.0). Compared TN-only vs TN + in-house FP init; re-scored TN against NIST23.
- **Outcome:** TN-init kept the FP tail more realistic (mean 0.212 vs 0.177 naive). LR at sim=1.0 stayed reasonable (17.6 vs a wildly overconfident 45.6 from naive EM).
- **Why archived:** Genuine improvement over naive EM but still lost to the 3-channel Bayesian on validation AUC. Keep the TN-augmentation idea in mind if EM is revisited.

### `reversed_ml_classifier.ipynb`
- **Intent:** Reverse-classify spectra into their BinBase bins using MSMS data — treat "which bin does this spectrum belong to?" as the primary supervised target.
- **Outcome:** Author's note in the notebook: *"not super doable because inner-bin spectra lack MSMS."* Many bins have too few MSMS-bearing spectra for a classifier.
- **Why archived:** Data infeasibility. Bin-level ML needs either more MSMS coverage or a different problem formulation.

---

## Early exploration, subsumed

### `benchmark_mar2024.ipynb`
- **Intent:** Early (March 2024) benchmark of confidence-scoring approaches on the Orbitrap dataset before methodology stabilized.
- **Why archived:** Superseded by `benchmark_oct2024.ipynb` (larger scope) and the current `fdr_benchmark_composite.ipynb` / `fdr_benchmark_scores.ipynb` suite.

### `example_for_benchmark.ipynb`
- **Intent:** Template / skeleton notebook showing how to generate benchmark pairs.
- **Why archived:** No substantive analysis; workflow templates now live in the active benchmarking notebooks.

### `bond_similarity_1_2.ipynb`
- **Intent:** Analyze bond-level structural similarity vs entropy similarity for isomer pairs. Test whether bond-based metrics improve confidence for near-isomers.
- **Outcome:** Exploratory; two bond-similarity metrics computed and compared, but no integration path surfaced.
- **Why archived:** Not adopted. Structural similarity is handled via InChIKey-14 dedup and MCS filtering upstream instead.

### `structural_similarity_mcs.ipynb`
- **Intent:** Use maximum common substructure (MCS, RDKit) to filter "one-bond-different isomers" from null pairs — build cleaner isomer-aware nulls.
- **Outcome:** MCS filtering worked; produced cleaner nulls in isolation.
- **Why archived:** Not integrated into the production pipeline. The simpler InChIKey-14 dedup (see `feedback_null_ik14_dedup.md` in memory — AUC 0.832 → 0.844) gave most of the win with less complexity.

---

## If you ever revisit

If any of the above ideas come back up, key files to re-read:
- GMM / clustering → start with `confidence_score_gmm.ipynb` for the failure mode.
- EM refinement → `tn_initialized_em.ipynb` has the best TN-augmented formulation.
- Null-distribution design → `entropy_null_exploration.ipynb` for the full exploration, `peak_null.ipynb` for the cleaner version.
- Peak-by-peak similarity → `peak-by-peak.ipynb` has parameter sweeps.
- NIST transfer → `transfer_null.ipynb` documents why it fails (useful if someone proposes it again).
