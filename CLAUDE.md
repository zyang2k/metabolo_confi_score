# CLAUDE.md

## Project Overview

**metabolo_confi_score** is a metabolomics compound identification confidence scoring project. The core problem: given an unknown LC-MS/MS spectrum, how confidently can we identify the compound?

The project builds a multi-component probabilistic framework combining:
- **MS2 spectral library matching** via entropy similarity
- **Statistical null distribution adjustment** for FDR control
- **Spectrum quality assessment** via ML classifier
- **Measurement uncertainty quantification** via bootstrap
- **Literature-based chemical priors** (PubMed/patent frequency)
- **Bayesian posterior scoring** with a "No Target Available" (NoTA) abstention mechanism

Primary language is **Python 3** with 25 Jupyter notebooks containing the core research. R scripts in the repo root are early exploratory analysis only.

## Technology Stack

- **Primary Language:** Python 3
- **Notebooks:** Jupyter (25 notebooks in `code/`)
- **Secondary:** R (early EDA scripts in repo root, not part of main pipeline)

**Key Python Dependencies** (no `requirements.txt` — install manually):
```
pandas, numpy, scipy, scikit-learn
ms_entropy       # spectral entropy similarity scoring
rdkit            # SMILES/MCS chemistry (isomer detection)
requests, urllib3  # MassWiki API
matplotlib       # visualization
```

**R packages** (exploratory only): `dplyr`, `ggplot2`, `readr`, `httr`, `jsonlite`, `purrr`, `car`, `fitdistrplus`

## Project Structure

```
metabolo_confi_score/
├── CLAUDE.md
├── README.md
├── metabolo_confi_score.Rproj        # RStudio project (R EDA only)
├── nist_entropy.py                   # NIST MSP parsing & entropy null distributions
├── code/                             # PRIMARY: all research notebooks & scripts
│   ├── 0918mvp/                      # Core MVP Bayesian scoring pipeline
│   │   ├── mvp_ingest.py             # Data loading & normalization
│   │   ├── mvp_score.py              # Bayesian scoring engine (3-channel logLR)
│   │   ├── mvp_export.py             # Output formatting
│   │   └── mvp_run_example.py        # End-to-end pipeline example
│   │
│   ├── — Research Notebooks —
│   ├── null_adjusted_score.ipynb     # Core pipeline: null dist + FDR-adjusted scoring
│   ├── reliability_metric.ipynb      # MS1/MS2 quality metrics framework
│   ├── spectrum_quality.ipynb        # Random Forest spectrum quality classifier
│   ├── peak_null.ipynb               # Null distribution from NIST23 (p-value calc)
│   ├── nist_null.ipynb               # Null distributions from NIST23 reference spectra
│   ├── peak-by-peak.ipynb            # Alternative MS/MS peak matching algorithm
│   ├── entropy_vs_similarity.ipynb   # Entropy metric exploration
│   ├── new_entropy_sim.ipynb         # Refined entropy similarity metrics
│   ├── bootstrap.ipynb               # Bootstrap CIs for mass/RT measurement uncertainty
│   ├── pubmed_prior.ipynb            # Literature-based chemical priors (PubChem API)
│   ├── reliability_metric.ipynb      # Comprehensive spectrum quality framework
│   ├── bond_similarity_1_2.ipynb     # Bond-level structural similarity
│   ├── istd_bin_distribution.ipynb   # Internal standard (ISTD) distribution analysis
│   ├── FDR_confidence.ipynb          # False Discovery Rate confidence calculations
│   ├── reversed_ml_classifier.ipynb  # Reverse ML classification experiments
│   ├── get_wiki_matches.ipynb        # MassWiki match retrieval client
│   ├── query_wiki.ipynb              # MassWiki API query utilities
│   ├── example_for_benchmark.ipynb   # Benchmarking workflow example
│   ├── benchmarking_0304.ipynb       # Benchmark run (Mar 2024)
│   ├── benchmarking_1030.ipynb       # Benchmark run (Oct 2024)
│   ├── goodspectrum.ipynb            # Good spectrum definition/criteria
│   │
│   ├── — Analysis Scripts —
│   ├── wikimatch_analysis0918.py     # Label fitting & logLR computation
│   ├── wikimatch_analysis0919.py     # Entropy similarity > 0.99 analysis
│   ├── masswiki_pipeline_with_token.py  # Multi-threaded MassWiki API fetcher
│   ├── masswiki_loop_all.py          # Sequential MassWiki fetcher (resumable)
│   ├── masswiki_min_tester.py        # Single wiki_id API test
│   ├── masswiki_min_from_csv.py      # CSV-based batch wiki_id fetcher
│   ├── masswiki_binbase_private.py   # Private BinBase connector
│   ├── diagnostics_channels.py       # QQ plots, PIT, AIC/BIC model diagnostics
│   ├── gof_channels.py               # Goodness-of-fit with ground-truth labels
│   └── hitscore_analysis.py          # High-confidence hit analysis
├── data/
│   ├── ttof+neg+hilic.csv            # Negative ion HILIC LC-MS spectra
│   ├── ttof+pos+rp.csv               # Positive ion RP LC-MS spectra
│   ├── masswiki_result.csv           # wiki_id + SPLASH identifiers
│   └── hilic_masswiki_reference_hits.csv  # MassWiki library hits
├── out_0918/                         # Pipeline outputs (Sep 18 run)
├── out_1017/                         # Pipeline outputs (Oct 17 run)
├── out_gof/                          # Goodness-of-fit diagnostic plots
├── viz/                              # Visualizations
├── MegaMoNA/                         # MassBank of North America reference data
├── cleaned_code/                     # Archived/cleaned code versions
└── 102424_bnb/                       # Meeting/presentation materials
```

## Research Summary (by Notebook)

### Null Distribution & FDR Control (`null_adjusted_score.ipynb`)
The core statistical engine. Builds null distributions of entropy similarity scores from random cross-compound spectral comparisons, then converts raw library match scores to p-values.

- **Data:** MoNA QTOF LC-MS/MS database (47,726 spectra), NIST23
- **Null construction:** 50,000+ random pairs from different compounds (different InChIKey), optionally mass-matched (±10 mDa)
- **P-value:** `p = mean(null_distribution >= observed_similarity)`
- **Key finding:** 95th percentile cutoff ~0.110 entropy similarity; raw score of 0.42 → p=0.0023 (highly significant)
- **FDR control at 5%** using adjusted thresholds
- **Mass-matched null** outperforms random null — accounts for isotope/fragment patterns

### Spectrum Quality Classifier (`spectrum_quality.ipynb`)
Random Forest classifier to filter low-quality spectra before library matching.

- **Training data:** 5,135 manually labeled spectra (91 negative, 5,044 positive)
- **Features (11):** corrected intensity, peak shapeness, standards found %, peak purity, raw intensity, SNR, entropy, gaussian similarity, dynamic range, etc.
- **Top features:** corrected intensity (0.214), peak shapeness (0.174), standards found % (0.115)
- **Performance:** 5-fold CV ROC-AUC = 0.9526; threshold=0.8 passes ~96% of spectra
- **External validation:** HILIC Orbitrap ROC-AUC = 0.8917; TTOF ROC-AUC = 0.4581 (platform generalization gap)

### Reliability Metrics Framework (`reliability_metric.ipynb`)
Defines composite quality scores for MS1 and MS2 measurements.

**MS1 quality metrics:**
- Ion Statistics Reliability (ISR): `1 - 1/sqrt(raw_intensity)`
- SNR reliability: `1 - exp(-SNR/threshold)`
- Peak quality: `cubic_root(gaussian_similarity × peak_purity × peak_shapeness)`

**MS2 quality metrics:**
- Entropy quality: `sqrt(normalized_entropy × entropy)`
- Run quality: standards found % (threshold = 80%)

**Final reliability score:** `base_reliability × (0.6 + 0.2×quality_modifiers + 0.2×run_quality)`

Database-level results: mass accuracy quality mean=0.741, standards quality mean=0.933, overall quality mean=0.710.

### Peak-by-Peak Matching (`peak-by-peak.ipynb`)
Alternative MS/MS matching algorithm evaluated against entropy similarity.

- Normalizes spectra (sum=1), finds best m/z match within tolerance
- Intensity similarity: `min(i1,i2) / max(i1,i2)`
- Weighted: `weight = mz^alpha × intensity^beta`, missing peak penalty = 0.1
- Best params: tolerance=0.01, alpha=0.5, beta=0.5
- **ROC-AUC:** peak-by-peak = 0.6854 vs entropy = 0.6923; correlation = 0.7575
- Conclusion: entropy similarity marginally better; methods are complementary

### NIST Null Distributions (`peak_null.ipynb`, `nist_null.ipynb`)
Builds statistical background from NIST23 negative mode library (50,000+ spectra).

- Stratified by number of peaks and spectral entropy
- Cross-compound comparisons only (different InChIKey)
- Mass tolerance ±10 mDa, MS2 tolerance 0.01 Da
- Outputs percentile-based p-value lookup tables

### Bootstrap Measurement Uncertainty (`bootstrap.ipynb`)
Quantifies instrument-level uncertainty for mass and retention time.

- 1,000 bootstrap iterations, 95% CI
- 1,322 compounds with ≥50 samples each (sample sizes 31–62)
- RT match score: `exp(-(Δrt²) / (2 × total_tolerance²))` where `total_tolerance = max(rt_uncertainty, method_tolerance)`
- Tighter bootstrap CI → higher confidence in identification

### Literature Chemical Priors (`pubmed_prior.ipynb`)
Incorporates chemical knowledge from literature into Bayesian priors.

- PubChem API: InChIKey → CID → PubMed reference count + patent count
- Prior score: `0.5 × min(1, log1p(pubmed_count)/6) + 0.5 × min(1, log1p(patent_count)/5)`
- Examples: caffeine prior=0.771, aspirin prior=0.861
- Well-studied metabolites get higher baseline confidence

### MassWiki Data Retrieval (`get_wiki_matches.ipynb`, `query_wiki.ipynb`)
Robust clients for the MassWiki spectral database API.

- Session pooling with retry (HTTP 429, 500–504)
- Throttled at 6 RPS, 8 concurrent workers (ThreadPoolExecutor)
- Outputs tidy DataFrame: `db, id, name, adduct, precursor_mz, entropy_similarity, rt, ri`
- Filters: manual annotations only, excludes `yy`/`zz` prefixed entries

### Bayesian MVP Pipeline (`code/0918mvp/`)
Three-channel log-likelihood ratio scoring with Bayesian posterior and NoTA abstention.

See [Scoring Architecture](#scoring-architecture) below.

## Scoring Architecture

### Full Confidence Score

The overall confidence is a product of independent components:

```
confidence = reliability × MS2_match_score × statistical_adjustment × prior
```

### Component 1: Spectrum Reliability
- **MS1:** ion statistics, SNR, peak shape quality
- **MS2:** spectral entropy, run quality (ISTD pass rate)
- Multiplicative modifier (0–1) applied to all downstream scores

### Component 2: MS2 Library Matching (Entropy Similarity)
- Primary metric: entropy similarity (ms_entropy library)
- Range [0, 1]; logit-transformed for Bayesian channel: `logit(entropy_similarity)`
- Compared against null distribution to yield p-value

### Component 3: Statistical Null Adjustment
- Null built from NIST23 or MoNA cross-compound random pairs
- P-value: fraction of null ≥ observed similarity
- Mass-matched null preferred (accounts for fragment overlap)
- FDR controlled at 5%

### Component 4: Three-Channel Bayesian LogLR (MVP pipeline)

Three independent channels combined via log-sum:

| Channel | TP Model | FP Model |
|---|---|---|
| **delta_ppm** (mass accuracy) | Student-t (df=3, scale ~3–5 ppm) | Laplace (scale ~8 ppm) |
| **entropy_similarity** (MS2) | Normal on logit(sim), μ≈1.2, σ≈0.6 | Normal, μ≈−0.2, σ≈1.0 |
| **zRT** (retention time) | N(0, 1) | N(0, σ_fp≈2) — optional |

**Per-hit log LR:** sum of per-channel log likelihood ratios
**Structure-level aggregation:** log-sum-exp over multiple hits for same compound
**Local posterior:** softmax over candidates per spectrum
**Global posterior:** includes NoTA (No Target Available) hypothesis — abstain if P_NoTA > 0.6

### Component 5: Literature Prior (optional)
```python
prior = 0.5 * min(1, log1p(pubmed_count)/6) + 0.5 * min(1, log1p(patent_count)/5)
```

## Running the Pipeline

### MVP Bayesian pipeline
```python
# See code/0918mvp/mvp_run_example.py
from mvp_ingest import load_spectra, load_hits, join_hits_to_spectra, build_confusable_sets
from mvp_score import fit_channel_params, per_hit_logLR, aggregate_to_structure, compute_global_posterior
from mvp_export import make_assertion_table, make_top_calls

spectra = load_spectra("data/ttof+neg+hilic.csv")
hits    = load_hits("data/hilic_masswiki_reference_hits.csv")
joint   = join_hits_to_spectra(hits, spectra)         # computes delta_ppm
cset    = build_confusable_sets(joint)                 # filter by entropy_sim ≥ 0.75, |zRT| ≤ 3
params  = fit_channel_params(label_df)                 # fit TP/FP distributions
scores  = per_hit_logLR(cset, params)
structs = aggregate_to_structure(scores)
post    = compute_global_posterior(structs)
```

### Fetch MassWiki reference hits
```python
# code/masswiki_pipeline_with_token.py
# Multi-threaded (8 workers), 6 RPS throttle, auto-fallback binbase → zyang2k
# Requires Bearer token (set externally)
```

### Build null distributions
```python
# code/peak_null.ipynb or nist_null.ipynb
# Requires NIST23 MSP file or MoNA JSON
# nist_entropy.py: parse_msp(), compute_random_similarity()
```

## Core Functions (`code/0918mvp/`)

### `mvp_ingest.py`
| Function | Description |
|---|---|
| `load_spectra(path)` | Load spectrum file; required: `wiki_id`, `rt`, `precursor_mz` |
| `load_hits(path)` | Load library hits; required: `wiki_id`, `library_id`, `name`, `adduct`, `lib_precursor_mz`, `entropy_similarity` |
| `join_hits_to_spectra(hits, spectra)` | Merge and compute `delta_ppm` |
| `build_confusable_sets(joint, ms2_min=0.75, zrt_k=3.0)` | Filter by entropy threshold and RT deviation |
| `collapse_duplicates_by_structure(cset, key_cols)` | Aggregate hits per structure via median/max |

### `mvp_score.py`
| Function | Description |
|---|---|
| `fit_channel_params(label_df)` | Fit TP/FP distributions per channel from labeled data |
| `per_hit_logLR(df, params)` | Per-hit log likelihood ratio across all channels |
| `aggregate_to_structure(df, quality_cols, topK, lambda_var)` | Log-sum-exp to structure level |
| `local_rank_probability(struct_scores)` | Softmax over candidates per spectrum |
| `compute_global_posterior(struct_scores, prior_policy, lr_nota, prior_nota)` | Bayesian posterior with NoTA |

### `mvp_export.py`
| Function | Description |
|---|---|
| `make_assertion_table(post)` | Full per-(spectrum, candidate) posteriors + features |
| `make_top_calls(post)` | Best call per spectrum, or abstain if P_NoTA > threshold |

## Data Sources

**MassWiki API:**
- `https://masswiki.us-west-2.elasticbeanstalk.com/analysis/get_data?wiki_id=<id>`
- Auth: Bearer token (not stored in repo — configure externally)
- Sources: BinBase (primary), zyang2k (fallback)

**PostgreSQL (AWS RDS):**
- Host: `lcb-standalone-cluster.cluster-czbqhgrlaqbf.us-west-2.rds.amazonaws.com`
- Database: `carrot-prod`
- Credentials: not stored in repo

**Reference spectral databases:**
- NIST23 (MSP format) — for null distribution construction
- MassBank of North America (MoNA) — 47,726 QTOF spectra used for benchmarking
- MegaMoNA/ folder in repo

## Output Files

| File | Description |
|---|---|
| `assertions.csv` | Per-(spectrum, candidate) posteriors + all features |
| `top_calls.csv` | Best identification per spectrum (or abstain if P_NoTA > threshold) |
| `hit_scores.csv` | Hit-level breakdowns with per-hit posteriors |
| `out_gof/` | Goodness-of-fit diagnostic plots (QQ, PIT, PDF overlays) |

## Key Findings & Benchmarks

- **Entropy similarity vs peak-by-peak:** ROC-AUC 0.6923 vs 0.6854; entropy preferred, methods are complementary (correlation 0.7575)
- **Null distribution p-values:** 95th percentile cutoff ~0.110 entropy similarity; raw score 0.42 → p=0.0023
- **Mass-matched null** outperforms random null for FDR control
- **Spectrum quality classifier:** ROC-AUC 0.9526 in-domain; drops to 0.4581 on TTOF (platform generalization gap needs attention)
- **Literature priors:** log1p-scaled PubMed + patent counts provide meaningful compound-level prior

## Known Issues / Open Questions

- No `requirements.txt` — dependencies must be installed manually
- Spectrum quality classifier does not generalize well across instrument platforms (TTOF AUC=0.458)
- Library search results often missing retention time data (zRT channel frequently unavailable)
- De novo identification workflow not implemented
- Entries without passing MSMS match thresholds need defined handling
- Many spectra lack compound names (large number of NAs in annotation data)

## Configuration

**`.gitignore`** excludes: `.Rproj.user`, `.Rhistory`, `.RData`, `.Ruserdata`

API tokens and database credentials are not stored in the repo and must be configured externally.
