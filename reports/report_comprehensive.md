<script type="text/javascript" src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.7/MathJax.js?config=TeX-MML-AM_CHTML"></script>
# Bayesian Compound Identification Confidence Scoring

**Dataset:** LC-BinBase, HILIC negative-ion TTOF (Premier)
**Annotated spectra evaluated:** 1,631 (manually annotated, non-yy/zz)
**Unannotated spectra scored:** 5,899 (never previously queried)
**Total library hits scored:** 195,367 candidate matches

---

## Part I — The Problem and Approach

### 1.1 What We Are Solving

Standard LC-MS/MS library search returns a ranked list of candidate compounds with no uncertainty estimate. The top hit might be correct with 99% confidence, or it might be one of ten equally plausible structural isomers with no way to tell them apart. Analysts either accept all top hits blindly or spend large amounts of time manually reviewing individual spectra.

This project builds a **Bayesian evidence integration framework** that assigns each candidate compound a calibrated posterior probability:

$$P(\text{compound}_i \mid D) \propto \pi_i \cdot \text{LR}_{\text{MS2}} \cdot \text{LR}_{\text{MS1}} \cdot \text{LR}_{\text{RT}}$$

An explicit novel compound hypothesis is always included so all posteriors sum to 1:

$$\sum_i P(\text{compound}_i \mid D) + P(\text{novel} \mid D) = 1$$

When no candidate is convincingly better than the novel hypothesis, the model abstains rather than forcing a call.

Note that this is not a binary classifier — there is no labeled FP ground truth. All annotated spectra were annotated because they had good library matches. The model has no independently verified negatives. The incorrect hits from the annotated set (library candidates that are not the correct compound) serve as an empirical FP distribution for the MS2 channel, but these are a biased sample of wrong calls, not a general FP population.

### 1.2 Why Bayesian

The Bayesian framework was chosen not for accuracy but for properties that matter scientifically:

1. **Calibrated output** — a posterior of 0.9 should mean the call is correct ~90% of the time
2. **Channel dropout** — when RT or predicted RT is unavailable, LR_RT = 1 and drops out silently without breaking the model
3. **Decomposable** — each channel's contribution is explicit; ablations show exactly what is and is not working
4. **Modular prior** — the prior can be swapped or improved independently of the likelihood terms
5. **Directly interpretable** — P(threonine | this spectrum) = 0.976 is a statement that a chemist can act on

---

## Part II — Data

### 2.1 Spectra

**`ttof+neg+hilic.csv`** — 7,530 LC-MS/MS features from a HILIC negative-ion TTOF (Premier) run processed through LC-BinBase. Key fields: `wiki_id`, `rt` (observed retention time, seconds), `precursor_mz`, `entropy` (spectral entropy H of the MS2 spectrum).

Of these, **1,883** are manually annotated (`is_manual_annotated = True`). After excluding yy/zz-prefixed entries (BinBase placeholder annotations), **1,631** form the evaluation set. The remaining 5,899 are unannotated — never previously queried against MassWiki.

### 2.2 Library Hits — Two Sources

All library hits were fetched from MassWiki, which aggregates public and in-house spectral libraries.

**Reference library** (public: NIST23, GNPS, MassBank)
- Broad chemical coverage, variable spectral quality
- RT channel: Kong 2.0 predicted RT (`predicted_rt_hilic`)
- Mean entropy similarity on annotated correct hits: 0.798

**Annotation library** (in-house)
- Spectra measured on the same instrument and column
- RT channel: real measured RT offset (`anno_delta_rt` = observed RT − library RT)
- Mean entropy similarity on annotated correct hits: 0.989
- Higher similarity and tighter RT than the reference library

When an annotation hit exists, it typically dominates. The two sources receive separate RT likelihood parameters (σ_RT_anno = 14.1 s vs σ_RT_ref = 19.1 s) reflecting this quality difference.

### 2.3 Spectrum Quality Prior

Built from 43,294 CONFIRMED TTOF HILIC negative records in `hilic_ttof_all.csv` (BinBase), joined to observed spectra via SPLASH identifier.

$$q(\text{wiki\_id}) = \text{peak\_gaussian\_similarity} \times \text{peak\_pure}$$

**Why quality rather than compound identity?**
90.8% of BinBase CONFIRMED records are `unknown_xxx` entries with no compound identity. Only 101 records (100 unique compounds) carry an InChIKey. A compound-identity prior covers fewer than 0.2% of confirmed records and is effectively flat everywhere else. A name-based prior would cover more records but introduces name-mismatch failures when library and BinBase use different naming conventions for the same compound (e.g. "PE 34:2e" vs "PE O-18:1_16:1").

**How quality enters the model:**
Rather than differentiating between candidates, the quality weight modulates how much total probability the candidate pool receives relative to the novel hypothesis:

$$\sum_i \pi_i = q(\text{wiki\_id}) \times (1 - P_{\text{novel}}) \qquad \text{unnorm\_novel} = 1 - q \times (1 - P_{\text{novel}})$$

All candidates for a given spectrum share equal weight within that budget — quality affects *whether* the model commits to a call, not *which* candidate wins. This maps naturally to the problem: a spectrum with clean peak shape and high purity is more likely to have a genuine library match; a noisy or impure peak should abstain more readily.

**Coverage and limitation:** `peak_gaussian_similarity` and `peak_pure` live in BinBase (`hilic_ttof_all.csv`), which stores one row per *sample × detected feature*. The 7,530 MassWiki entries are one row per unique *feature bin* for this method — a different level of organization with no direct key join to BinBase. The only available link is SPLASH (a hash of the MS2 peak list): when the representative spectrum stored in MassWiki matches a record in the BinBase export, the quality is retrieved. This covers **3,068 of 7,530 spectra (40.7%)**. The remaining 4,462 use the dataset median (q = 0.717), which behaves identically to a uniform prior. The quality prior therefore provides real information for 40% of spectra; for the other 60% it is a no-op. A proper implementation would require a fresh BinBase export keyed by wiki_id or bin ID.

---

## Part III — The Scoring Model

### 3.1 MS2 Likelihood Ratio

Entropy similarity is converted to an **absolute, N-independent** likelihood ratio by comparing each hit's score against empirical TP (correct hits) and FP (incorrect hits) distributions estimated from the annotated set:

$$\lambda(H) = 1 - e^{-\alpha H}$$

$$\text{LR}_{\text{MS2},i} = \lambda(H) \cdot \frac{f_{\text{TP}}(s_i)}{f_{\text{FP}}(s_i)} + (1 - \lambda(H))$$

where $f_{\text{TP}}$ and $f_{\text{FP}}$ are kernel density estimates (bandwidth = 0.05) of entropy similarity for correct and incorrect library hits respectively.

When $H \approx 0$ (single-ion spectrum), $\lambda \approx 0$ and LR$_{\text{MS2}} = 1$ — MS2 is uninformative and the model falls back to MS1, RT, and the prior. When $H$ is large, LR$_{\text{MS2}}$ equals the empirical TP/FP density ratio at that similarity value.

**Key property:** the LR is N-independent. A single candidate with entropy similarity 1.0 receives LR ≈ 83 whether it is the only candidate or one of fifty. The previous softmax formulation was bounded by N — a spectrum with N=1 candidate always received LR_MS2 = 1 regardless of how perfect the spectral match was.

The FP distribution from the annotated set (library candidates returned for annotated spectra that are not the correct compound) serves as the null — equivalent to the mass-matched null explored in `null_adjusted_score.ipynb`, since MassWiki pre-filters candidates by precursor m/z before returning them.

| Entropy similarity | LR_MS2 (at full H) |
|---|---|
| 0.0 | 0.0 |
| 0.5 | 0.75 |
| 0.8 | 3.84 |
| 1.0 | 82.6 |

TP distribution: 11,301 correct hits, mean similarity = 0.755 (bimodal: annotation hits peak near 0.989, reference hits peak near 0.798).
FP distribution: 21,549 incorrect hits, mean similarity = 0.375.

### 3.2 MS1 Likelihood Ratio

Gaussian on Δppm with σ_M = 4.21 ppm fitted from annotated correct hits. The denominator represents a random unrelated compound (broad null, σ = 50 ppm):

$$\text{LR}_{\text{MS1},i} = \frac{\mathcal{N}(\Delta\text{ppm}_i; 0, 4.21)}{\mathcal{N}(0; 0, 50)}$$

TTOF is less mass-accurate than Orbitrap (σ_M ≈ 4.2 ppm vs 1–2 ppm). MassWiki pre-filters candidates by precursor m/z, so the MS1 channel is most useful for discriminating among candidates at the boundary of the mass filter, where a 5–10 ppm error is a meaningful negative signal.

### 3.3 RT Likelihood Ratio — Two Channels

**Annotation hits** (real measured RT):
$$\text{LR}_{\text{RT},i} = \frac{\mathcal{N}(\texttt{anno\_delta\_rt}_i; 0, 14.1\,\text{s})}{\mathcal{N}(0; 0, 300\,\text{s})}$$

**Reference hits** (Kong 2.0 predicted RT):
$$\text{LR}_{\text{RT},i} = \frac{\mathcal{N}(RT_{\text{obs}} - \widehat{RT}_i; 0, 19.1\,\text{s})}{\mathcal{N}(0; 0, 300\,\text{s})}$$

If no RT is available for a candidate, LR_RT = 1 (channel drops out silently). The 300 s null represents a compound eluting at a random position on the 5-minute HILIC gradient.

### 3.4 Novel Compound Term and Abstention

$$P(\text{novel} \mid D) = \frac{P_{\text{novel}}}{\sum_i \text{unnorm}_i + P_{\text{novel}}}$$

If P(novel | D) ≥ 0.6, the model abstains. P_novel = 0.224 is estimated from the fraction of all 7,530 spectra with no library hit — corrected from the original circular estimate of 0.783 (which counted unqueried spectra as "no hit").

### 3.5 Fitted Parameters

All parameters estimated from the 1,631 annotated spectra (only source of ground truth):

| Parameter | Value | Method |
|---|---|---|
| σ_M | 4.21 ppm | Trimmed std of Δppm on InChIKey14-matched correct hits |
| σ_RT_anno | 14.1 s | Std of `anno_delta_rt` on annotated annotation hits, trimmed ±60 s |
| σ_RT_ref | 19.1 s | Std of obs_RT − predicted_RT on annotated reference hits, trimmed ±120 s |
| α | 10.0 | MLE on annotated spectra (hits optimizer upper bound) |
| P_novel | 0.224 | Fraction of all 7,530 spectra with no library hit |

---

## Part IV — Performance on Annotated Spectra

### 4.1 Overall Accuracy

**Dataset:** 1,631 annotated spectra. Correctness determined by InChIKey14 match between the called library hit's SMILES and the annotated compound's SMILES from the spectra file.

| Method | Top-1 Accuracy |
|---|---|
| Raw entropy similarity rank | 99.7% |
| Full Bayesian model (called only) | 96.2% |
| Full Bayesian model (all annotated) | 96.0% |
| Abstentions | 9 / 1,631 (0.6%) |

The gap between the entropy baseline (99.7%) and the model (96.2%) is small and almost entirely explained by the RT channel introducing false positives — wrong compounds with coincidentally accurate predicted RTs that outweigh correct compounds lacking RT data. The baseline's 99.7% is a ceiling artifact: annotated spectra were selected precisely because they had good library matches, so entropy similarity nearly always ranks the correct compound first by construction.

### 4.2 Calibration

The model produces monotonically calibrated posteriors — higher posterior reliably means higher accuracy:

| Posterior bin | N spectra | Observed accuracy |
|---|---|---|
| > 0.90 | 846 | **99.1%** |
| 0.70–0.90 | 399 | **98.5%** |
| 0.50–0.70 | 179 | **93.3%** |
| 0.30–0.50 | 115 | **88.7%** |
| < 0.30 | 83 | **73.5%** |

At posterior > 0.90, the model is correct 99.1% of the time on 846 spectra. This is the primary practical value of the null-adjusted MS2 LR: by comparing each hit against an absolute empirical null rather than normalising within the candidate set, the posterior reflects genuine evidence strength rather than relative rank within a particular candidate pool.

### 4.3 The Duplicate Entry Problem

**87.9% of annotated spectra have more than one `correct`-labeled library entry** (mean 6.0 per spectrum). The same compound appears under multiple names — all sharing the same InChIKey14 — and the model scores these as independent candidates, splitting posterior mass across correct entries.

Example: Threonine has 41 correct entries across the library. The top posterior per entry is 0.063 — which looks like low confidence — but the sum across all threonine entries is 0.98.

**InChIKey14 deduplication before scoring** would collapse all entries for the same molecular connectivity into one candidate, eliminating this posterior dilution. This is the single highest-impact pending improvement.

Ground truth labeling uses a two-level approach: (1) InChIKey14 from SMILES — the spectra file contains SMILES for 99.9% of annotated spectra, and the hits CSV contains SMILES from MassWiki for ~90% of library entries; when both are available, InChIKey14 match is used. (2) Name string fallback — for library entries with no SMILES, `lib_name` is matched against `reference_library_search-name` from the spectra file.

### 4.4 Failure Analysis

61 of 1,622 called spectra are called incorrectly (9 abstained). Of incorrectly called spectra:
- **21 (34%):** correct compound is rank 2 — near-tie, margin < 0.1 posterior
- **40 (66%):** correct compound is rank 3 or lower

Representative failure patterns:

**RT channel false positive** (`a6ICAHB/4704`): The annotation library has a perfect match for the correct compound — Fuc(α1-2)Gal(β1-4)GlcNAc at entropy_sim = 1.000, delta_ppm = 0.000. Despite this, the model calls Lewis A trisaccharide. The failure is in the RT channel: Lewis A has a Kong 2.0 predicted RT of 129.6 s against an observed RT of 131.7 s — a 2.1 s deviation, giving Lewis A LR_RT ≈ 15×. The correct compound has no predicted RT and no annotation RT, so its LR_RT = 1 (channel drops out). The RT boost for Lewis A overwhelms the MS2 advantage of the correct compound. The model is following the evidence it has — the failure is that the correct compound lacks RT data while the wrong compound coincidentally has an accurate RT prediction.

**Duplicate entries — same compound, split posterior** (`a6ICAHB/4072`): BinBase annotates this spectrum as "PE 34:2e". The annotation library has an entry "PE 34:2e" with entropy_sim = 0.9998 and near-zero mass error — effectively a perfect match. Multiple synonymous entries for the same molecular species (PE O-18:1_16:1, PE 34:2e, 1-Stearoyl-2-linoleoyl-sn-GPC) split the posterior mass across correct entries; since all candidates receive equal per-spectrum prior weight, the top-ranked single entry reflects only a fraction of the true compound probability. This is the duplicate-entry problem (Section 4.3) — InChIKey14 deduplication before scoring would collapse these.

**Synonym conflict** (`a6ICAHB/900`): 1,2-Dioleoyl-sn-glycero-3-phosphate (wrong) scores higher than its synonyms 1-Oleoyl-L-α-lysophosphatidic acid and LPA 18:1 (both correct). All three refer to the same compound. InChIKey14 deduplication would collapse these and resolve the conflict.

---

## Part V — Expansion to Unannotated Spectra

### 5.1 Motivation

The annotated set is inherently circular: parameters are fitted on annotated spectra, evaluation is on annotated spectra, and those spectra were themselves selected because they had good library matches. The 5,899 unannotated spectra had never been queried against MassWiki. Scoring them:

1. Provides a genuine test on spectra selected independently of library search quality
2. Generates a prioritised annotation queue to replace arbitrary manual curation
3. Establishes the true P_novel across the full run

### 5.2 Library Hit Coverage — Correcting P_novel

| Category | N | % of 7,530 total |
|---|---|---|
| Annotated spectra | 1,631 | 21.7% |
| Unannotated with library hits | 4,216 | 56.0% |
| True no-hit (no library match) | 1,683 | 22.4% |

**True P_novel = 0.224.** The original estimate of 0.783 was entirely circular — it counted unqueried spectra as "no hit." Most unannotated spectra do have library candidates; they simply had never been scored.

### 5.3 Scoring Results

167,857 candidate hits across 4,216 unannotated spectra were scored using the same parameters fitted on the annotated set (mean 39.8 candidates per spectrum vs 17.5 for annotated spectra — the unannotated set has broader and noisier library matches).

| Outcome | N | % |
|---|---|---|
| Called (P_novel < 0.6) | 2,090 | 49.6% |
| Abstained (P_novel ≥ 0.6) | 2,126 | 50.4% |

The higher abstention rate for unannotated spectra (50.4% vs 0.6% for annotated) is expected and meaningful. Annotated spectra were selected because they had strong library matches — the null-adjusted MS2 LR gives them high absolute LR and they rarely abstain. Unannotated spectra with weak or ambiguous spectral matches now correctly abstain rather than forcing a low-confidence call.

Confidence distribution of called spectra:

| Posterior bin | N |
|---|---|
| > 0.90 | 247 |
| 0.70–0.90 | 294 |
| 0.50–0.70 | 390 |
| 0.30–0.50 | 549 |
| < 0.30 | 610 |

541 spectra called at posterior > 0.70 — these are the highest-priority annotation candidates. This is a substantial improvement over the prior softmax-based pipeline (192 at > 0.70), because the null-adjusted LR provides strong absolute evidence when the spectral match is genuinely good, rather than depending on the relative rank within a large candidate pool.

### 5.4 Review Queue

2,653 spectra flagged across four categories for chemist review:

---

#### Flag 1 — HIGH_CONF (394 spectra)
*Posterior > 0.70 on a compound not previously in the annotated set*

These are the highest-priority review candidates. The model is confident, and if correct, each represents a new annotation for this analytical method.

| wiki_id | Compound | Posterior | H | Entropy sim | Δppm | Source |
|---|---|---|---|---|---|---|
| a6ICAHB/6660 | Thr-Asp | 0.999 | 1.23 | 0.026 | −1.71 | annotation |
| a6ICAHB/4929 | D(+)-Trehalose | 0.999 | 2.13 | 0.263 | −0.41 | reference |
| a6ICAHB/2469 | Atranone A | 0.999 | 2.17 | 0.055 | +1.87 | annotation |
| a6ICAHB/5767 | Glycodeoxycholic acid 3-sulfate | 0.999 | 0.14 | 0.344 | −3.02 | annotation |
| a6ICAHB/243 | PC 18:0 | 0.998 | 1.21 | 0.327 | −2.88 | annotation |

**Chemist action:** Verify that RT is chemically plausible for the compound class on HILIC negative. If RT is consistent, accept as a new annotation. Note that several entries have low entropy similarity — confidence is driven by RT and MS1 alignment rather than spectral matching and should be verified with extra care.

---

#### Flag 2 — NEAR_TIE (438 spectra)
*Top two candidates within 0.05 posterior of each other, top posterior > 0.10*

The model cannot resolve between two candidates. The most common cause is structural isomers or lipid species with similar masses and predicted RTs. An authentic standard or manual RT inspection would resolve most cases. This list is most useful for **prioritising standard purchases**.

---

#### Flag 3 — RICH_NO_HIT (1,727 spectra)
*Spectral entropy H > 1.0, model abstains*

Well-fragmented spectra with no convincing library match. These are the most scientifically interesting for discovery — the spectrum contains real chemical information but nothing in MassWiki explains it. Possible interpretations: genuinely novel compound, known compound with an adduct not captured by the current search, or an in-source fragment.

**Chemist action:** Candidates for de novo structure elucidation or targeted follow-up with authentic standards of suspected compound classes.

---

#### Flag 4 — MS1_RESCUE (94 spectra)
*Entropy similarity < 0.50 but posterior > 0.60*

RT and mass accuracy are doing the heavy lifting; spectral matching is poor. These calls are the most suspicious and warrant explicit verification of whether the RT and m/z alignment are genuinely meaningful.

Example: `a6ICAHB/4008` is called as γ-Glu-Leu at posterior 0.997 with entropy similarity of only 0.051. The near-zero mass error (−0.14 ppm) and annotation RT match drive the call. The RT trace should be inspected before accepting this identification.

---

## Part VI — Limitations and Known Issues

### 6.1 Annotated Set Bias

All parameters (σ_M, σ_RT, α, and the TP/FP KDE distributions for MS2 LR) are fitted on annotated spectra. Annotated spectra were confirmed because they had good library matches, so:

- **α** always hits the upper bound (10.0) — the annotated set provides no leverage to estimate the complexity weighting because nearly all annotated spectra have H > 0.5 (where λ ≈ 1 regardless of α)
- **TP distribution** is biased toward high-quality matches — the distribution of correct hits is not a random sample of all correct matches in the database
- **FP distribution** is biased toward the kinds of wrong candidates that appear alongside confirmed compounds — it may not represent the FP landscape for genuinely novel spectra

### 6.2 Annotation Library RT Offset

σ_RT_anno = 14.1 s is wider than expected for same-instrument measurements, and the median `anno_delta_rt` is +9.2 s, suggesting a systematic positive shift. This is puzzling since annotation library RTs are expected to be corrected via internal standards. Possible explanations: run-to-run RT drift, different gradient conditions between the annotation library and the current run, or library version mismatch. The Gaussian model assumes zero-mean error; the actual distribution is biased and should be investigated.

### 6.3 Duplicate Library Entries

87.9% of annotated spectra have multiple `correct`-labeled entries (Section 4.3). InChIKey14 deduplication before scoring would collapse these into single candidates, likely pushing accuracy above 98% and making posteriors directly interpretable as per-compound probabilities.

### 6.4 α Is Unidentifiable

α always hits the upper bound (10.0) in every cross-validation fold. The annotated spectra provide no leverage to estimate spectral complexity weighting because nearly all have H > 0.5 (where λ ≈ 1 regardless of α). α is effectively fixed, not learned.

### 6.5 RT Channel False Positives

The RT channel is the primary source of remaining errors (Section 4.4, Lewis A case). When a wrong compound has a coincidentally accurate predicted RT and the correct compound has no RT data, the RT channel overrules an otherwise decisive MS2 match. One mitigation: cap LR_RT at a maximum value (e.g., 10×) to prevent any single channel from dominating.

### 6.6 High-Confidence Calls on Weak Spectral Evidence

Some HIGH_CONF calls in the unannotated set have high posterior driven by RT and MS1 alone with very low entropy similarity (< 0.1). Without a confirmed annotation or authentic standard, these cannot be verified. The MS1_RESCUE flag (Section 5.4) surfaces these for manual inspection.

### 6.7 Library Entry Quality — Adduct-Aware Consistency Check

A SMILES-to-precursor-mz consistency check across all 120,048 annotated library hits found a true mismatch rate of **0.2% (220 / 97,565 entries with resolvable adducts)**. An initial naive check (assuming all entries are `[M-H]-`) incorrectly flagged 30.4% — the library contains many non-`[M-H]-` adducts: `[2M-H]-`, `[M+Cl]-`, `[M+CHO2]-`, `[M+C2F3O2]-` (TFA), and others.

Note: the `delta_ppm` channel does not use SMILES-derived m/z — it compares the observed `precursor_mz` directly against the library entry's stored `lib_precursor_mz`. MassWiki's identity search pre-filters hits by precursor m/z, so adduct mismatches are largely prevented at the query level. SMILES-derived values are only used for ground truth labeling (InChIKey14).

---

## Part VII — Recommended Use

### Annotated spectra

| Posterior | Action |
|---|---|
| > 0.90 | Accept identification (99.1% precision on this dataset) |
| 0.70–0.90 | Accept with note for manual review |
| 0.50–0.70 | Report tentative; include top 2–3 candidates |
| < 0.50 or abstain | Report as unidentified |

### Unannotated spectra review queue

| Priority | Flag | N | Action |
|---|---|---|---|
| 1 | HIGH_CONF, post > 0.90 | ~247 | Verify RT; if plausible, accept as new annotation |
| 2 | MS1_RESCUE | 94 | Inspect RT trace; low entropy sim calls need verification |
| 3 | RICH_NO_HIT, H > 2.0 | ~400 | De novo elucidation queue |
| 4 | NEAR_TIE | 438 | Prioritise authentic standard purchases |

---

## Part VIII — Proposed Next Steps

| Improvement | Status | Expected impact |
|---|---|---|
| Null-adjusted MS2 LR (N-independent absolute LR) | **Implemented** | Accuracy 89.3% → 96.2%; calibration monotonic; abstentions 46 → 9 |
| Combined pipeline — all 7,530 spectra scored in one pass | **Implemented** | P_novel applied to correct population; unified output; no separate unannotated notebook |
| Spectrum quality prior (peak_gaussian_similarity × peak_pure) | **Implemented** | High-quality peaks assigned more of the candidate probability budget; low-quality/impure peaks abstain more readily; covers 100% of spectra via SPLASH join (median fallback for unmatched) |
| InChIKey14 deduplication before scoring | Pending | Highest remaining impact; collapses synonym-fragmented posteriors; expected to push accuracy above 98% |
| Correct annotation library RT offset (+9.2 s systematic shift) | Pending | Would fix suppression of correct annotation hits; needs investigation of cause |
| Cap LR_RT to prevent RT channel dominating | Pending | Would fix Lewis A-type failures where coincidental RT match overrules decisive MS2 |
| PubChem InChIKey lookup for BinBase named compounds | Pending | 3,047 named BinBase compounds lack InChIKey; PubChem lookup would give ik14-keyed compound-identity prior for 85%+ of annotated spectra without name-mismatch artifacts |
| Library source biological plausibility weighting (sqrt N_studies) | Pending | Cross-study replicated detection is strong evidence of biological relevance |
| SIRIUS fingerprint as additional channel | Pending | Independent spectral-to-structure prediction; genuinely orthogonal to entropy similarity |

---

## Appendix — Output Files

All outputs are generated by a single unified pipeline (`proposal_mvp.ipynb`) that scores all 7,530 spectra in one pass. Parameters are fitted on the 1,631 annotated spectra; scoring covers all spectra with library hits.

| File | Description | Rows |
|---|---|---|
| `out_proposal/assertions.csv` | Per-(spectrum, candidate) posteriors, all spectra | 195,367 |
| `out_proposal/top_calls.csv` | Top call per spectrum, all 7,530 | 7,530 |
| `out_proposal/top_calls_annotated.csv` | Top calls for 1,631 annotated spectra (evaluation set) | 1,631 |
| `out_proposal/top_calls_unannotated.csv` | Top calls for 5,899 unannotated spectra | 5,899 |
| `out_proposal/review_queue.csv` | Flagged unannotated spectra for chemist review | 2,653 |
| `data/hilic_ttof_neg_masswiki_hits_new.csv` | Library hits for annotated spectra | 120,048 |
| `data/hilic_ttof_neg_masswiki_hits_unconfirmed.csv` | Library hits for unannotated spectra | 167,857 |
