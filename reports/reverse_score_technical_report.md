# Reverse (containment) spectral matching for bin de-orphaning — technical report

**Author:** Ziyue Yang · **Date:** 2026-06-17 · **Status:** internal, drafts-on-disk

---

## Abstract

We evaluated the NIST-style **reverse (containment)** spectral score for two jobs in the LCBinBase
confidence pipeline. As an **identity** discriminator (rescuing co-isolation-contaminated library
matches) it failed — it recovers contaminated true positives but rates retention-time-wrong false
positives equally, so its rescue-zone precision equals the base rate. As a **de-orphaning** tool —
deciding whether an unaccepted candidate bin is merely an in-source fragment / adduct / isotope of a
compound already in the library — it works, because that task is pure *containment* (a fragment is a
subspectrum of its parent). On the production `compound` table we flag **9.9%** of HILIC candidate
bins, **18.2%** of C18-negative (lipid) candidates, and **10.2%** of C18-positive — each validated by
permutation null controls (signal/null ≈ 6–20×), with the relation profile tracking ionization mode
(losses in negative, adducts in positive). The deliverable is a per-method "don't-promote / link-to-
parent" list.

---

## 1. Background & objective

The reverse score normalizes the matched dot product by the *query's* matched intensity only, ignoring
unmatched query peaks. The literature (Xing et al., reverse spectral search) reports it rescues
contaminated/chimeric annotations. Objective: determine where, if anywhere, it adds value in our
pipeline. Two candidate jobs were tested — **identity rescue** (§5.1) and **candidate-bin
de-orphaning** (§5.2–5.4).

---

## 2. Data

**Primary source.** carrot-prod (BinBase) PostgreSQL, `compound` table, read-only, exported per method
via DBeaver (the RDS is VPN-gated; not reachable headlessly). Relevant columns:

| column | role |
|---|---|
| `id` | row key |
| `accurate_mass` | precursor m/z |
| `retention_index` (RI) | cross-sample-aligned retention axis (co-elution; preferred over raw `retention_time`) |
| `msms` | peak list, `"mz:intensity …"`, absolute intensity |
| `target_type` | population label (below) |
| `fragment_of` | CARROT's own in-source-fragment link (found empty — see §7) |

**Populations** (`target_type`):
- `CONFIRMED` — accepted library bins → **reference** set
- `UNCONFIRMED` — generated candidate bins, QC-passed, awaiting acceptance; **method-level, one row per splash** → **query** set (the de-orphan target)
- `INVALID_TARGET` — **per-sample** detections rejected by bin-generation rules, dominated by run-level ISTD-coverage failures (i.e. real compounds from QC-failed injections, *not* spectral noise) → initially used, then demoted (§5.5)

**Per-method counts:**

| method | UNCONFIRMED | CONFIRMED | INVALID_TARGET |
|---|---|---|---|
| HILIC amide · Orbitrap · neg | 26,527 | 6,398 | 286,011 |
| C18 "splash one" · Orbitrap · neg (lipid) | 360,231 | 23,059 | 336,133 |
| C18 "splash one" · Orbitrap · pos (lipid) | 385,346 | 37,263 | 749,951 |

(C18-pos UNCONFIRMED was sampled to 50,000 via `ORDER BY random() LIMIT 50000`; HILIC and C18-neg are
full.) The identity-rescue probe (§5.1) used a separate bin-level export — `lipidomics_{neg,pos}.csv`
(10,854 + 18,310 spectra) — with query/library peaks fetched from MoNA REST, GNPS, MassBank-EU, and a
local NIST23 MSP (MassWiki's by-library-id endpoint is dead). Labels there: TP = compound (IK14) with
≥2 adducts within a 10 s RT window; FP = `yy_`-prefixed; reverse/forward cosine via `compute_ms2_scores`.

---

## 3. Methods

**Unit of analysis (de-orphan):** candidate `O` (UNCONFIRMED) vs reference bins `C` (CONFIRMED),
matched method-level on RI.

**Decision rule** — flag `O` as a relational artifact of `C` if all hold:
1. **Co-elution:** `|RI(C) − RI(O)| ≤ RT_WIN` (sorted-RI window).
2. **Mass relation** on `Δ = precursor(C) − precursor(O)`:
   - neutral loss (fragment, C heavier): `|Δ − L| ≤ MTOL`
   - adduct (bidirectional): `|‖Δ‖ − L| ≤ MTOL`
   - isotope (**directional**, candidate is the heavier ¹³C peak): `|Δ + L| ≤ MTOL`
3. **Containment (reverse score):** `Σ(O intensities whose m/z is present in C within TOL) / Σ(all O intensity)`. One-sided — extra peaks in C are not penalized.
4. **Acceptance:** ISF requires `containment ≥ CONT_MIN` **and** O's precursor m/z present as a peak in C (strong signature); adduct/isotope require `containment ≥ CONT_MIN_RELATED`.

Output metric: **flagged rate** = fraction of candidates with ≥1 accepted parent. Each flagged
candidate carries `{parent_id, relation, containment}`.

### Thresholds

| param | value | gates |
|---|---|---|
| `RT_WIN` | 4.0 RI units | co-elution (= BinBase annotation RI window) |
| `MTOL` | 0.006 Da | precursor Δm/z ↔ loss/adduct/isotope |
| `TOL` | 0.01 Da | MS² peak match in containment |
| `CONT_MIN` | 0.50 | containment, ISF (+ precursor-in-parent) |
| `CONT_MIN_RELATED` | 0.70 | containment, adduct/isotope |

**Dictionaries.** ~20 generic neutral losses (H₂O, 2H₂O, 3H₂O, NH₃, CO, CO₂, HCOOH, CH₂O, CH₃OH,
CH₃COOH, C₂H₄, C₃H₆, C₂H₂O, H₂O+CO₂, SO₃, H₃PO₄, pentose, hexose, glucuronide, HCl); adducts (Na–H,
K–H, NH₄–H); isotopes (¹³C, 2×¹³C). ³⁴S/³⁷Cl excluded (valid only with parent formula, which the table
lacks). For `lipid`-tagged runs, +18 lipid losses merge in (6 headgroups: phosphocholine,
phosphoethanolamine, serine, glycerophosphate, inositol, TMA; 12 fatty-acyl as RCOOH and ketene for
16:0/18:0/18:1/18:2/20:4/22:6). Cutoffs are principled defaults (mirror ISFrag's reverse-dot 0.5/0.7
and BinBase deconvolution tolerances); no ROC sweep performed yet.

---

## 4. Validation design

**Permutation null controls** (specificity): break one necessary condition, hold the rest, re-run the
full denoise. Spectra stay bound to each candidate by id; only one column is relabeled. Fixed seed.
- **random-Δm/z:** permute `accurate_mass` across candidates → breaks the mass relation; keeps
  co-elution + spectra. Estimates the coincidental-match floor (precision proxy).
- **RI-shuffle:** permute `retention_index` → breaks co-elution; keeps mass + spectra. Measures
  co-elution's contribution.
- **decoy-loss:** swap the loss/adduct/isotope dictionaries for non-physical masses (41.30, 55.73,
  73.11, 88.62) on real data → tests whether the dictionary fabricates matches (expect ~0).

**Other checks:** chemical coherence of named-parent hits; per-injection stability (per-sample
variant); same-injection-presence (per-sample variant, via `sample_annotation_data`); cross-method
generalization; the confirmed-vs-confirmed library audit (§5.4).

---

## 5. Results

### 5.1 Identity rescue — negative

On lipidomics, criterion `entropy_sim < 0.7 ∧ reverse_cos > 0.7`, adduct-matched (n=4,908): **281 TP
rescued** (their containment ≈ 0.92 — the correct compound's peaks are present, entropy tanked by
co-isolation), but rescue-zone **precision = 0.646 ≈ base rate 0.670**. Rescued TP and rescued FP are
statistically identical because the FP are `yy_` (retention-time-wrong matches with good MS²), which a
spectral metric correctly rates high. An early 0.80 precision was a coverage artifact (the FP class was
under-covered). **Verdict:** reverse is a triage/explanation flag, not an identity discriminator.

### 5.2 De-orphaning UNCONFIRMED — HILIC

**9.9%** (2,630/26,527) flagged: ISF 2,171 / adduct 354 / ¹³C-isotope 105; 1,576 link to a named
parent (e.g. glutamate→102 −CO₂, β-alanine→71 −NH₃). Null controls: real 9.9% vs **random-Δm/z 1.5%**,
RI-shuffle 7.6%, decoy 0%; median 9 peaks. The high RI-shuffle (7.6%) shows method-level co-elution is
weak discrimination on the crowded HILIC RI axis → firmest calls are isotopes/adducts/named-parent
fragments; random-Δm/z (1.5%) is the tighter precision floor (~85% non-coincidental).

### 5.3 Cross-method — lipidomics (C18)

| method | flagged | dominant relation | real / random-Δm/z / RI-shuffle / decoy |
|---|---|---|---|
| HILIC-neg | 9.9% (2,630/26,527) | in-source fragments | 9.9 / 1.5 / 7.6 / 0 |
| C18-neg (lipid) | 18.2% (65,573/360,231) | ISF + ¹³C-isotope | 18.2 / 0.9 / 3.8 / 0 |
| C18-pos (lipid) | 10.2% (5,106/50,000) | **adducts (NH₄/Na/K)** | 10.2 / 1.1 / 2.2 / 0 |

Two structural findings: **(a)** on reverse-phase the RI axis is wide and uncrowded, so RI-shuffle
collapses to 2–4% — co-elution, the weak link in HILIC, is strongly discriminating in lipids, making
the lipid numbers the most trustworthy. **(b)** the relation profile **tracks ionization mode** —
negative is led by losses (water/CO₂/acetate, the `[M+OAc]⁻`/`[M−H]⁻` family), positive by adducts
(TG [M+NH₄]⁺, PC [M+Na/K]⁺) — an independent confirmation the rule reads chemistry. Adding
lipid-specific headgroup/acyl losses changed C18-neg by only +0.6 pp (17.6→18.2%), so the artifacts are
mostly generic adducts/isotopes/small losses, not classic acyl fragmentation.

### 5.4 Confirmed-library audit (side use)

Running the same test confirmed-vs-confirmed flags **13.4%** of HILIC confirmed bins (858/6,398) as a
relational ion of another confirmed bin — candidate mis-promoted fragments/isotopes. It re-discovers
curator `yy_`/"in source" flags blind. A separate library-QC deliverable, triage not verdict.

### 5.5 Methodological corrections

- **Isotope directionality bug.** The isotope matcher initially used `|‖Δ‖ − L|` (two-sided), flagging
  candidates *lighter* than the parent (e.g. "2×¹³C of mannitol", impossible). Requiring the candidate
  to be the heavier peak dropped HILIC ¹³C flags 1,013→105 and the headline **13.2%→9.9%**. Caught while
  building the slide deck.
- **INVALID_TARGET demoted.** Reading CARROT's `AddToLibrary` rules showed `INVALID_TARGET` is mostly
  run-level ISTD-coverage rejection (real compounds from bad runs), not noise; the per-sample run on it
  (HILIC 20.1% method-wide, 12.0% same-injection) was set aside in favor of `UNCONFIRMED`.
- **Coverage-parity lesson (§5.1).** Always check label-class coverage before quoting precision.

---

## 6. Caveats & limitations

- **No internal ground truth:** `fragment_of` is empty across all methods, so validation rests on the
  null controls + chemistry + (per-sample) same-injection presence; curator spot-check pending.
- **Residual ≠ novel:** the 82–90% unflagged is "no qualifying parent found" — true novels + dictionary
  misses + background; still needs the within-spectrum quality (S_norm) gate.
- **Method-level co-elution is weak in HILIC** (RI-shuffle 7.6%); strong in RP lipids (2–4%).
- **C18-pos is a 50k sample** (±~0.3% on the rate); no ROC threshold sweep yet; ³⁴S/³⁷Cl disabled
  pending a parent-formula source.

---

## 7. Recommendations / next steps

1. Ship the per-method **don't-promote / link-to-parent** list for curator review; spot-check ~20
   flagged candidates (named-parent fragments + low-mass marker ions).
2. Run the full C18-pos (385k) for an exact figure; add a **threshold ROC** once a labeled set exists.
3. Treat the **confirmed-library audit** as its own QC pass (the ¹³C-isotope and high-containment ISF
   cases are highest-priority).
4. Keep reverse score **out** of identity scoring — it is a completeness/containment tool only.

---

## Appendix — artifacts

- Pipeline: `code/analysis/binbase_orphan_denoise_run.py` (`--unconfirmed`, `--tag`, `--validate-bins`, `--same-sample`; null controls built in); kernel `code/analysis/isf_orphan_denoise.py`.
- Extraction SQL: `code/sql/binbase_extract_denoise.sql` (HILIC), `binbase_lipid_unconfirmed.sql` (C18-neg), `binbase_lipid_unconfirmed_pos.sql` (C18-pos).
- Identity-rescue harness: `code/bench_reverse_rescue_lipids.py`.
- Outputs: `data/{tag}_unconfirmed_denoise.csv`; deck `reports/unconfirmed_deorphan_slides.pptx`; memo `reports/oliver_isf_denoise_memo_20260616.md`.
- Branch `kg-frozen-graph-gate2`; reverse-score work in commits d1c41e3 … 3c1c8a1.
