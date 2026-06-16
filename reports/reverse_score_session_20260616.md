# Reverse Score, Two Ways: Identity-Rescue vs Relational-Denoise

**Session report — 2026-06-15/16**
Topic: where (if anywhere) does the NIST-style *reverse* spectral match earn its keep in our
confidence pipeline? We tested it on two distinct jobs and got two clean, opposite answers.

---

## TL;DR

| Job | Question | Verdict |
|---|---|---|
| **A. Identity rescue (lipidomics)** | Does reverse rescue contaminated-but-correct *library matches* that entropy/forward wrongly rejects? | **No, not as a discriminator.** It recovers real contaminated TPs, but rescues wrong-RT FPs at the same rate — precision ≈ base rate. Use only as a triage/explanation signal. |
| **B. Relational denoise (real LCB)** | Can reverse/containment cut the impact of noise ions among bin-unmatched MS² in regular studies? | **Yes.** ~12% (parent-presence-confirmed) to ~20% (method-wide) of floating MS² are explained as ISF/adduct/isotope of a co-eluting confirmed bin — triple-validated. |

One throughline: **reverse score's home is containment/completeness, not identity.** Confirmed on
both fronts — it fails at "which compound is this?" and succeeds at "is this a fragment-subset of
that co-eluting compound?"

---

## A. Identity rescue on lipidomics

**Setup.** Local `data/lipidomics_{neg,pos}.csv` (bin-level annotations; no reverse score, no peaks).
Labels per user: TP = compound (IK14) with ≥2 distinct adducts within 10 s RT; FP = `yy_`; ignore
`zz_`; rest unlabelled → 6,048 annotated working bins (TP 1,214 / FP 610 / unlabelled 4,224).

**Data plumbing built.** Query peaks via MassWiki `get_data` (`spectrum.peaks`); library peaks per
source by accession — **MoNA REST**, **GNPS** (`gnps2/gnpsspectrum`), **MassBank EU** (GitHub raw
records), **NIST23** from the local 2.7 GB MSP. (`get_data`-by-library-id is dead — HTTP 400.)
Reverse/forward cosine via `compute_ms2_scores`. Coverage reached **88%** (5,323/6,048).

**Result (criterion: `entropy_sim < 0.7` AND `reverse_cos > 0.7`, adduct-matched n=4,908):**
- **281 TP rescued** (of 334 entropy-rejected TP), mechanism real: rescued TP have entropy median
  0.57 but reverse 0.96 and library-coverage `cov_int` 0.92 — the correct compound's peaks are all
  present, entropy was just dragged down by co-isolated extras.
- **But precision = 0.646 ≈ the base TP rate (0.670).** Reverse does *not* enrich correct over
  incorrect. Rescued TP and rescued FP are statistically identical (reverse 0.96 vs 0.97, cov_int
  0.92 vs 0.89) because the FP here are `yy_` = RT-disagreement rejects with genuinely *good* MS²,
  which a spectral metric correctly rates high. Measuring "precision vs FP" conflates Q1 (spectral
  rescue) with Q2 (curator's RT decision).

**Key correction.** A first pass on MoNA-only coverage (58%) reported precision 0.80 — an artifact:
MoNA coverage silently dropped ~80% of the FP class. Restoring coverage put precision at its true
~0.65. **Lesson: always check label-class coverage parity before quoting precision.**

**Verdict.** Far better than the prior HILIC rescue probe (5.4% precision), so lipid co-isolation is
genuinely where reverse *could* matter — but it remains a triage/explanation signal ("low entropy
here is contamination, not mismatch — don't auto-reject, verify RT"), **not an accept gate.**

---

## B. Relational denoise of bin-unmatched MS² (real LCB)

This is reverse's actual strength: an in-source-fragment/adduct/isotope orphan **O** is a *subset* of
its richer co-eluting parent **C**. Forward cosine fails (C's extra peaks tank it); one-sided
containment(O⊂C) stays high. = the CAMERA/ISFrag paradigm pointed at the orphan population.

**Data.** carrot-prod `compound` table (read-only via DBeaver — the RDS is VPN/private-VPC, not
reachable from home or by the agent). Method `5m hilic premier | orbitrap | beh amide | negative`.
- bins (reference) = `target_type='CONFIRMED'` (6,398, all with MS²)
- orphans (floating) = `INVALID_TARGET` (286,011) + `UNCONFIRMED` (26,527), named `unknown_*`
- co-elution axis = `retention_index` (cross-sample aligned), window 4 (= BinBase's RI window)
- decision rule: co-elution + Δm/z ∈ {neutral loss / adduct / isotope} + containment(O⊂C) ≥ 0.5 +
  (for ISF) the orphan precursor present as a peak in the parent.

**Headline (50 injections, 24,423 floating MS²):**
- **20.1% explained** as ISF/adduct/isotope of a co-eluting confirmed bin → collapsible noise.
- **Stable** across injections: median 20.4%, IQR 18.8–21.5%, range 15.2–25.1%.
- Relation mix: ISF 3,738 / isotope 784 / adduct 376. Top losses H₂O, CO₂, ¹³C, **SO₃ (sulfate)**,
  acetate, formate, **³⁴S** — chemically coherent for a bile-acid (`negBA`) matrix.
- containment median 0.92; **63% link to a *named* confirmed bin** (citric/salicylic→CO₂,
  phenylacetic→CO, dehydroascorbic→CO₂ — textbook in-source fragments).

**Validation.**
1. **Null controls (all pass):** real 20.1% vs random-Δm/z **3.2%** (chemistry is real → ~84% of
   calls non-coincidental), RI-shuffle 8.3%, decoy-loss (non-physical) **0.0%** (the dictionary
   doesn't manufacture matches).
2. **Same-injection parent** (parent must be detected in the orphan's own run, via
   `sample_annotation_data`): on 10 injections (5,344 orphans), method-wide 19.4% → **same-injection
   12.0%**; the RI-shuffle background **halves, 8.3% → 4.1%**. So ~7 pp of the method-wide rate were
   links to compounds not actually present (over-count, conservatively removed). **12% is the
   defensible floor** (the detection list counts only confirmed, non-gap-filled parents).
3. **CARROT `fragment_of` concordance — unavailable.** Zero `fragment_of` labels on bins or orphans
   in this method: CARROT does not annotate in-source fragments here, which is precisely the gap this
   tool fills. (A cross-method count is queued to see if any method populates it.)

**Verdict.** Reverse/containment removes **~12–20%** of bin-unmatched MS² as explainable relational
noise — validated by three null controls + same-injection presence + chemistry + cross-injection
stability. For the curator, each removed orphan is labeled "in-source fragment of [named compound],"
and the residual becomes a cleaner novel-candidate worklist (which still needs the within-spectrum
quality gate before anything is called novel — the residual is *not* "all novels").

**Explicit answer to the driving question** — *Can reverse score help in LCB to reduce the impact of
noise ions for the hundreds of MS/MS that don't match bins in regular studies?* **Yes:** ~20%
method-wide / ~12% same-injection of the unconfirmed MS² per run are relational artifacts of
co-eluting confirmed bins, auto-labelable and collapsible, stable across injections, and validated by
controls + chemistry.

### Side finding — auditing the CONFIRMED library (confirmed-vs-confirmed)
The main run uses confirmed bins as a clean reference. Running the same test *among* the confirmed
bins flags **1,108 / 6,398 (17.3%)** as an ISF/adduct/isotope of another co-eluting confirmed bin
(657 ISF / 290 isotope / 161 adduct). A recurring in-source fragment recurs as often as its parent →
can clear the recurrence threshold and be promoted to its own bin. Validation that it's real: it
re-discovers **159 already-curator-flagged** cases (`yy_` / "in source to…") blind, and the chemistry
is textbook (genistein ← genistein-4′-glucuronide via glucuronide loss; xanthine ← xanthosine via
pentose; 4-methylcatechol ← guaiacol sulfate via SO₃). **949 are not curator-flagged** (518 ISF / 281
isotope / 150 adduct) → a review list; the **281 isotope-of-another-bin** cases are the highest
suspicion (isotopes shouldn't be separate compounds). Caveat: triage, not verdict — mixes genuine
mis-bins, legitimate separate adducts, and coincidental isobar co-elutions; needs the same
same-injection + control + curator rigor before any specific bin is called an error. Arguably the
higher-stakes use, since a mis-promoted bin propagates into downstream biology/quantification.

---

## Infrastructure learnings (for next time)
- **Library peaks:** `get_data?wiki_id=<library_id>` is dead. Fetch by accession — MoNA REST
  (`MoNA*`), GNPS gnps2 (`CCMSLIB*`), MassBank EU GitHub raw (`MSBNK-*`); NIST23 only from local MSP.
- **Cache-by-index is broken:** `library_peaks_cache.json` numeric `index` is per-pull position, not
  a stable id (a cached GNPS index returned forward=0 vs server 0.58). Splash is the only stable key.
- **carrot-prod:** VPN/private-VPC (DNS fails from home; agent connection auto-blocked). Use DBeaver.
  `compound` holds bins+floating split by `target_type`; `sample_annotation_data` is a ~600-partition
  view that **must** be bounded by a literal `acquired` date to prune (study acquired 2024-03-05→08).
- **Security:** `code/r_eda/get_binbase_data.R` has the DB password in plaintext — rotate + env-var.

## Artifacts
- `code/analysis/binbase_orphan_denoise_run.py` — end-to-end runner (`--from-csv`, `--same-sample`,
  `--validate-bins`; per-injection stability + 3 null controls + concordance built in)
- `code/analysis/isf_orphan_denoise.py` — containment kernel + expanded loss/adduct dictionary
- `code/bench_reverse_rescue_lipids.py` — lipidomics identity-rescue harness
- `code/sql/binbase_extract_denoise.sql` — pure-SQL extraction (bins, orphans, detections)
- data: `reverse_rescue_lipids.csv`, `lcb_orphan_denoise.csv`, `lcb_hilicneg_{bins,orphans}.csv`,
  `lcb_sample_detections.csv`; peak caches `mona_lib_peaks_cache.json`, `lipid_web_libpeaks_cache.json`

## Open / next (all optional)
- Curator spot-check of ~20 explained calls (esp. low-containment) + residual → human precision.
- Cross-method run (HILIC-pos / C18 / lipid) → generalization + the bad-spectrum-vs-coverage reason
  split (≈absent here: 22,977 coverage rejects vs 2 bad-spectrum) + possibly revive `fragment_of` #1.
- Package as a repeatable "denoise a study" deliverable (method + acquired window → per-orphan
  {explained, parent, relation} + residual novel worklist).
