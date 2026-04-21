# Annotation Confidence Scoring — Forward/Reverse Similarity & Multi-Channel Evidence

**Date:** 2026-04-07
**Dataset:** Orbitrap HILIC negESI, 1,988 named entries (362 golden TP + 936 regular TP + 475 first-pass + 215 FP)
**Notebook:** `code/annotation_confidence_v2.ipynb`

## Goal

For each of Oliver's 1,988 named annotations, score how trustworthy the annotation is given all available spectral evidence — using **correct IK14-matched compound pairs** for every channel.

## The Spreadsheet

`data/masswiki_Orbitrap HILIC negESI_2026-03-19.xlsx` — 5,971 total entries, 1,988 with non-blank names.

| Tier | Count | Definition |
|---|---|---|
| **Golden TP** | 362 | Solid TPs (112, yellow/orange) ∪ commented (82) ∪ multi-adduct confirmed (227), all within rows 0–1297 |
| **Regular TP** | 936 | Rows 0–1297, reviewed twice by Oliver, not golden |
| **First-pass** | 475 | After row 1297, named but never re-examined |
| **FP** | 215 | `yy_` prefix, confirmed wrong |

## Finding: identity_score Measures the Wrong Compound 13.3% of the Time

`identity_score_reference_library_` in the MassWiki output is the entropy similarity of the **top reference library hit** — whichever compound scores highest, regardless of whether it matches the annotation.

| | Count | % |
|---|---|---|
| Top hit IS the annotation compound (IK14 match) | 1,309 | 86.7% |
| Top hit is a DIFFERENT compound | 201 | 13.3% |
| → annotation still in results, just not #1 | 173 | |
| → annotation not in results at all | 28 | |

**Impact:** When Oliver evaluates an annotation, `delta_rt` is computed against the annotation compound (correct), but `identity_score` may be computed against a different compound. For 13.3% of spectra, these channels evaluate different compounds without Oliver knowing.

**Example:** Spectrum aPUDE1U/1006 ("maleic acid"):
- `delta_rt`: evaluated against **maleic acid** ✓
- `identity_score = 0.857`: evaluated against **DL-Malic acid** (different compound!) ✗
- Maleic acid's actual sim = 0.814 (rank #2)

See `results/reverse_similarity/identity_score_mismatches.csv` (201 rows) for all cases.

## Evidence Channels

### Always available from spreadsheet

| Channel | Coverage | Source | What it measures |
|---|---|---|---|
| **Spectral entropy** | 1,988 (100%) | BinBase | Is the spectrum informative? S<1 unreliable, S=1-3 sweet spot |
| **Delta RT** | 1,943 (97.7%) | Observed RT − Kong-predicted RT from annotation SMILES | Does the compound elute when expected? |

### IK14-matched annotation hit channels

These require finding the annotation compound in external reference libraries (NIST23, GNPS, MassBank, MB-EU) by InChIKey14 matching. **All channels use the correct compound pair.**

| Channel | Coverage | What it measures |
|---|---|---|
| **anno_entropy_sim** | 1,482 (74.5%) | Entropy similarity of annotation compound's library entry |
| **anno_forward** | 1,482 (74.5%) | Query signal explained by annotation compound (entropy-weighted) |
| **anno_reverse** | 1,482 (74.5%) | Library signal found in query (entropy-weighted) |
| **anno_delta_ppm_abs** | 1,482 (74.5%) | Mass error against annotation compound's library entry |
| **sim_gap** | 1,482 (74.5%) | Annotation sim minus best competitor's sim |
| **anno_rank** | 1,482 (74.5%) | Rank of annotation compound among all reference hits |
| **max_deviation** | 1,249 (62.8%) | Worst single-peak intensity mismatch |

**Coverage gap:** 506 spectra lack IK14-matched channels. 475 are first-pass annotations whose hits haven't been fetched (need fresh API token). 31 are spectra where the annotation compound isn't in any external library.

### Why not use identity_score?

`identity_score_reference_library_` is available for 99.6% of spectra but measures the **top reference hit**, not the annotation compound. For 13.3% of spectra, this is a different compound. Our `anno_entropy_sim` always measures the correct pair.

## Model

**Training:** GradientBoostingClassifier on 362 golden TPs vs 215 FPs.
**Features:** 9 channels, all using correct IK14-matched pairs (except delta_rt and spectral_entropy which don't need matching).
**NaN handling:** Fill with training set median.

### Cross-validation: AUC = 0.833

### Feature importance

| Feature | Importance | Per-feature AUC |
|---|---|---|
| delta_rt_abs | 37.8% | 0.770 |
| sim_gap | 15.9% | 0.602 |
| anno_forward | 12.8% | 0.697 |
| anno_entropy_sim | 12.1% | 0.685 |
| anno_delta_ppm_abs | 7.8% | 0.556 |
| max_deviation | 6.5% | 0.655 |
| spectral_entropy | 4.3% | 0.554 |
| anno_reverse | 2.5% | 0.592 |
| anno_rank | 0.3% | 0.550 |

Delta RT dominates (37.8%), but spectral channels collectively contribute 62.2% — with reverse score (12.8%), sim gap (15.9%), and annotation entropy sim (12.1%) as the main spectral drivers.

## Results

### Score distribution by tier

| Tier | n | Mean | Median | P10 | P90 |
|---|---|---|---|---|---|
| Golden TP | 362 | 0.868 | 0.918 | 0.694 | 0.974 |
| Regular TP | 936 | 0.794 | 0.882 | 0.453 | 0.974 |
| First-pass | 475 | 0.354 | 0.265 | 0.120 | 0.753 |
| FP | 215 | 0.222 | 0.168 | 0.043 | 0.491 |

**First-pass entries score low (median 0.265)** because their spectral channels are NaN (hits not yet fetched). They're scored on delta_rt and spectral_entropy alone. Their true scores will be higher once hits are fetched.

### Flagged annotations

**49 regular TPs (5.2%)** score below the FP 75th percentile (0.319). These were reviewed twice by Oliver but the evidence is weak:

| Spectrum | Score | Delta RT | Entropy sim | Reverse | Name |
|---|---|---|---|---|---|
| aPUDE1U/676 | 0.115 | 52.0s | 0.831 | 0.000 | 5-aminonicotinic acid |
| aPUDE1U/5509 | 0.154 | 53.4s | 0.807 | 0.556 | 4-acetamidobutanoic acid |
| aPUDE1U/5808 | 0.101 | 21.5s | 0.707 | 0.525 | 3-ureidopropionic acid_to check |
| aPUDE1U/591 | 0.119 | 1.2s | 0.483 | 0.000 | DL-2-Hydroxyvaleric acid |
| aPUDE1U/3209 | 0.119 | 1.3s | 0.500 | 0.381 | olivetolic acid |

Common patterns: large delta RT (compound elutes at wrong time), reverse=0.000 (no query peaks match the library entry), or low entropy sim with the annotation compound.

## Output Files

| File | Rows | Description |
|---|---|---|
| `annotation_confidence_scores.csv` | 1,988 | All named entries with confidence score + all evidence channels |
| `identity_score_mismatches.csv` | 201 | Spectra where identity_score ≠ annotation compound |
| `identity_score_audit.csv` | 1,510 | All fetched spectra with ik14_match flag |
| `deduped_hits_validation.csv` | 13,505 | All deduped hits with hit_correct/is_anno_hit |
| `annotations_missing_library_match.csv` | 30 | Spectra where annotation compound not in external libraries |

## Next Steps

1. **Fetch hits for 475 first-pass entries** — need fresh MassWiki API token. This will raise their coverage from delta_rt-only to full multi-channel scoring.
2. **Send identity_score mismatch finding to Oliver** — the 201-row CSV shows where MassWiki's identity_score measures a different compound than the annotation.
3. **Validate on holdout** — 112 solid TPs are held out. Check that they score consistently high.
4. **Calibration** — convert raw model scores to calibrated probabilities (reliability diagram).
