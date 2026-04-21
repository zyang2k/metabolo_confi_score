# Annotation Confidence Score v2

**Date:** 2026-04-08
**Author:** Ziyue Yang
**Dataset:** Orbitrap HILIC negESI, 1,513 spectra (1,298 TP + 215 FP)
**Notebook:** `code/annotation_confidence_v2.ipynb`

## Goal

For each of Oliver's 1,298 manually reviewed annotations, produce a confidence score (0–1) that flags which annotations are most likely to be wrong. The score combines 9 evidence channels spanning retention time, spectral matching, peak coverage, and candidate uniqueness.

## Training Data

| Label | Count | Source |
|---|---|---|
| TP | 1,298 | Rows 0–1297 in Oliver's spreadsheet (reviewed twice) |
| FP | 215 | `yy_` prefix entries (confirmed wrong by Oliver) |

Within TP, Oliver color-coded the `annotation-smiles` column:
- **Green** (178): highest confidence — annotation confirmed golden
- **Light blue** (1,008): regular reviewed TP
- **Yellow/orange** (112): solid TPs

The model trains on all 1,298 TPs vs 215 FPs. Colors are used for analysis, not for training labels.

## The 9 Evidence Channels

All spectral channels (3–9) compare the query spectrum against the **annotation compound's** library entry specifically, matched by InChIKey14. This ensures every channel evaluates the same compound pair.

### Channel 1: Delta RT (52.5% importance)
- **What:** `|observed RT - Kong predicted RT|`
- **Source:** Spreadsheet column `DRTpred`, predicted from annotation SMILES
- **Coverage:** 97.7%
- **Signal:** Large delta RT = compound shouldn't elute here. FP mean = 34.4s vs TP mean = 10.2s.

### Channel 2: Sim gap (13.5%)
- **What:** Annotation compound's entropy similarity minus the best competitor's entropy similarity
- **Source:** Computed from deduplicated reference library hits (IK14-matched)
- **Coverage:** 74.5% (needs annotation compound in external library)
- **Signal:** Positive = annotation is clearly the best match. Negative = a different compound matches better. FP annotations have smaller gaps.

### Channel 3: Annotation entropy similarity (9.6%)
- **What:** Entropy similarity between query spectrum and annotation compound's best external library entry
- **Source:** Computed from IK14-matched reference hits
- **Coverage:** 74.5%
- **Signal:** How well the query matches the annotation compound. Unlike the spreadsheet's `identity_score_reference_library_` (which is the top hit regardless of compound), this always measures the correct compound.

### Channel 4: Reverse score (8.7%)
- **What:** Fraction of the query spectrum's entropy-weighted signal explained by the annotation compound's library entry
- **Source:** Peak-level computation using `ms_entropy.clean_spectrum()` + `apply_weight_to_intensity()`
- **Coverage:** 74.5%
- **Signal:** Low reverse = the query has peaks the annotation compound can't explain (wrong compound, co-elution, ISF). FP mean = 0.42 vs TP mean = 0.60.

### Channel 5: Max deviation (5.0%)
- **What:** `max(|log2(query_intensity_i / library_intensity_i)|)` across matched peaks
- **Source:** Peak-level computation
- **Coverage:** 62.8% (needs ≥2 matched peaks)
- **Signal:** The worst single-peak intensity mismatch. One badly wrong peak is more diagnostic than uniformly mediocre matching.

### Channel 6: Delta PPM (4.3%)
- **What:** `|observed precursor m/z - library precursor m/z|` in ppm
- **Source:** Annotation hit's `lib_precursor_mz`
- **Coverage:** 74.5%
- **Signal:** Mass accuracy against the annotation compound specifically.

### Channel 7: Spectral entropy (3.8%)
- **What:** Shannon entropy of the query spectrum's intensity distribution
- **Source:** Spreadsheet column `entropy` (computed by BinBase)
- **Coverage:** 100%
- **Signal:** How informative the MS2 is. S < 1 = poor fragmentation, unreliable. S = 1–3 = sweet spot (per Yuanyue Li, Nature Methods 2021).

### Channel 8: Forward score (1.9%)
- **What:** Fraction of the annotation compound's library entry explained by the query spectrum (entropy-weighted)
- **Source:** Peak-level computation
- **Coverage:** 74.5%
- **Signal:** Low forward = the query is missing peaks the library predicts.

### Channel 9: Annotation rank (0.7%)
- **What:** Rank of annotation compound among all reference library hits by entropy similarity (1 = best match)
- **Source:** Computed from deduplicated reference hits
- **Coverage:** 74.5%
- **Signal:** Low importance because sim gap already captures this more precisely.

## Model

- **Algorithm:** GradientBoostingClassifier (100 trees, max_depth=3, learning_rate=0.1)
- **Labels:** TP=1 (1,298 spectra), FP=0 (215 spectra)
- **NaN handling:** Fill with training set median per feature. 74.5% of spectra have all 9 channels; 25.5% are scored on delta_rt + spectral_entropy + median-filled spectral channels.
- **Output:** P(annotation is correct) between 0 and 1

### 5-fold CV AUC: 0.855 ± 0.023

## Flagged Annotations

**22 out of 1,298 TPs (1.7%)** score below the FP 75th percentile (confidence < 0.603). These annotations have weak evidence across multiple channels and are worth re-examining.

| wiki_id | Tier | Score | Delta RT | Entropy sim | Reverse | Sim gap | Name |
|---|---|---|---|---|---|---|---|
| aPUDE1U/2848 | blue | 0.374 | 32.6s | 0.67 | 0.46 | 0.08 | nonenedioic acid |
| aPUDE1U/5278 | blue | 0.387 | 60.5s | N/A | N/A | N/A | carglumic acid |
| aPUDE1U/3482 | blue | 0.415 | 38.5s | 0.75 | 0.28 | 0.12 | chloramphenicol |
| aPUDE1U/809 | golden | 0.464 | 45.3s | 0.72 | 0.34 | 0.00 | N-acetylleucine |
| aPUDE1U/4598 | blue | 0.468 | 25.1s | 0.67 | 0.39 | 0.07 | octyl sulfate |
| aPUDE1U/3235 | golden | 0.494 | 46.6s | 0.99 | 0.69 | 0.16 | 4-hydroxyphenyl sulfate |
| aPUDE1U/748 | golden | 0.496 | 30.7s | N/A | N/A | N/A | imidazolepropionic acid |
| aPUDE1U/5057 | golden | 0.513 | 58.2s | 0.60 | 0.00 | 0.02 | dihydrouracil |
| aPUDE1U/676 | blue | 0.516 | 52.0s | 0.83 | 0.00 | 0.05 | 5-aminonicotinic acid |
| aPUDE1U/3249 | blue | 0.532 | 60.1s | 0.71 | 0.47 | 0.02 | N-acetylcysteine |
| aPUDE1U/5179 | golden | 0.534 | 41.7s | 0.85 | 0.57 | 0.00 | N-acetylleucine_minor |
| aPUDE1U/4256 | golden | 0.540 | 41.7s | 0.90 | 0.50 | 0.07 | N-acetylvaline |
| aPUDE1U/3653 | blue | 0.551 | 22.1s | 0.71 | 0.00 | 0.18 | ascorbic acid |
| aPUDE1U/5526 | blue | 0.563 | 37.1s | 0.76 | 0.21 | 0.01 | erythrose |
| aPUDE1U/3857 | yellow | 0.568 | 103.1s | N/A | N/A | N/A | tetraconazole |
| aPUDE1U/3269 | blue | 0.568 | 33.8s | 0.80 | 0.38 | 0.01 | carmofur |
| aPUDE1U/1059 | blue | 0.571 | 22.4s | 0.70 | 0.24 | -0.05 | 5-(methoxymethyl)pyridine-2,3-dicar |
| aPUDE1U/5808 | blue | 0.590 | 21.5s | 0.71 | 0.53 | 0.08 | 3-ureidopropionic acid_to check |
| aPUDE1U/4069 | blue | 0.598 | 37.5s | 0.87 | 0.65 | 0.07 | isonicotinylglycine |
| aPUDE1U/422 | blue | 0.599 | 36.7s | 0.87 | 0.16 | 0.07 | purine |

### Common patterns in flagged spectra

- **Large delta RT** (20–103s) — all 22 flagged spectra have delta RT > 20s
- **reverse = 0.000** — dihydrouracil, 5-aminonicotinic acid, ascorbic acid have zero query peaks matching the annotation library entry
- **No library hit** — carglumic acid, imidazolepropionic acid, tetraconazole: annotation compound not found in external libraries
- **sim_gap ≤ 0** — N-acetylleucine, N-acetylleucine_minor: a competitor compound matches equally well or better
- **"_to check" in name** — 3-ureidopropionic acid_to check: Oliver himself flagged uncertainty in the name

## Coverage Note

25.5% of spectra (mostly those whose annotation compound is not in NIST23/GNPS/MassBank/MB-EU) lack IK14-matched spectral features. These are scored on delta_rt + spectral_entropy only, with median fill for the missing channels. If a fresh MassWiki API token is obtained, some of these gaps can be closed.

## Output Files

| File | Rows | Description |
|---|---|---|
| `results/annotation_confidence_v2/annotation_confidence_scores.csv` | 1,513 | All spectra: confidence score + all 9 channels + tier + Oliver's score |
| `results/annotation_confidence_v2/flagged_annotations.csv` | 22 | Suspicious TPs sorted by confidence |
| `data/green_golden_tp.csv` | 277 | Green-highlighted wiki_ids (golden annotations) |
| `data/spreadsheet_colors.csv` | 5,971 | All spreadsheet rows with extracted annotation-smiles column color |
| `results/reverse_similarity/identity_score_mismatches.csv` | 201 | Spectra where identity_score measures a different compound |
