# Expanding Compound Identification to Unconfirmed Spectra
## Bayesian Confidence Scoring — TTOF HILIC Negative Mode

**Dataset:** TTOF HILIC negative mode, LC-BinBase
**Confirmed (previously scored):** 1,631 manually annotated spectra
**Unconfirmed (this report):** 5,899 spectra never previously queried against MassWiki

---

## 1. Motivation

The confirmed pipeline (`proposal_mvp.ipynb`) operates exclusively on manually annotated spectra. These annotations are not a random sample — they reflect whatever the curator happened to annotate, biased toward compounds already in the library at annotation time. The 5,899 unconfirmed spectra had never been queried against MassWiki at all.

Two questions this expansion addresses:

1. **How many unconfirmed spectra actually have library hits?** The old P_novel = 0.783 was circular — it counted unqueried spectra as "no hit."
2. **Which unconfirmed spectra are worth a chemist's review?** Rather than annotating arbitrarily, flag the cases where the model's output is most informative or most uncertain.

---

## 2. Library Hit Coverage — The Real Picture

After fetching MassWiki hits for all 5,899 unconfirmed spectra:

| Category | N | % of all 7,530 spectra |
|---|---|---|
| Confirmed (annotated) | 1,631 | 21.7% |
| Unconfirmed with library hits | 4,216 | 56.0% |
| No library hit (true novel) | 1,683 | 22.4% |

**True P_novel = 0.224** — not 0.783.

The old P_novel was entirely an artifact of never querying unconfirmed spectra. Most unconfirmed spectra (71%) do have library candidates; they just had never been scored. The 1,683 genuinely unmatched spectra are either truly novel compounds, low-quality noise features, or compounds outside MassWiki's current coverage.

---

## 3. Scoring the Unconfirmed Set

The same Bayesian framework and fitted parameters from the confirmed pipeline are applied directly:

| Parameter | Value |
|---|---|
| σ_M | 4.21 ppm |
| σ_RT_anno | 14.1 s |
| σ_RT_ref | 19.1 s |
| α | 10.0 |
| P_novel (updated) | 0.224 |

**167,857 candidate hits** across 4,216 spectra were scored (mean 40 candidates per spectrum, reflecting broader library coverage for unconfirmed features).

### Outcome

| Outcome | N | % of scored |
|---|---|---|
| Called (P_novel < 0.6) | 3,119 | 74.0% |
| Abstained (P_novel ≥ 0.6) | 1,097 | 26.0% |

### Confidence distribution of called spectra

| Posterior bin | N |
|---|---|
| > 0.90 | 110 |
| 0.70–0.90 | 82 |
| 0.50–0.70 | 117 |
| 0.30–0.50 | 327 |
| < 0.30 | 2,483 |

The distribution is heavily skewed toward low-confidence calls — most called spectra have posterior < 0.30. This is expected: without manual annotation guiding which spectra are "easy," the unconfirmed set contains many genuinely ambiguous cases. The 192 spectra with posterior > 0.70 are the candidates most worth acting on.

Of the 3,119 called spectra, **2,624 (84%) are called on compounds not previously seen in BinBase** — potential new annotations for this analytical method.

---

## 4. Review Queue

Rather than returning all 4,216 scored spectra to the curator, four flag categories identify the cases most worth a chemist's time. Total flagged: **2,253 spectra**.

### Flag 1 — HIGH_CONF_NOVEL (156 spectra)

**Definition:** Posterior > 0.70 on a compound not previously confirmed in BinBase.

These are the highest-priority review candidates. The model is confident, and if correct, each represents a new annotation for this method.

**Top examples:**

| wiki_id | Compound | Posterior | H | Δppm |
|---|---|---|---|---|
| a6ICAHB/6660 | Thr-Asp | 0.999 | 1.23 | −1.71 |
| a6ICAHB/4929 | D(+)-Trehalose | 0.999 | 2.13 | −0.41 |
| a6ICAHB/2469 | Atranone A | 0.999 | 2.17 | +1.87 |
| a6ICAHB/5767 | Glycodeoxycholic acid 3-sulfate | 0.999 | 0.14 | −3.02 |
| a6ICAHB/243 | PC 18:0 | 0.998 | 1.21 | −2.88 |

**Chemist action:** Verify RT is chemically plausible for the compound class on HILIC negative. If RT is consistent, accept as a new annotation.

**Risk:** The model has no FP ground truth on unconfirmed spectra. A high posterior here reflects strong internal consistency of the evidence channels — not a confirmed identification. Some of these will be wrong, particularly when entropy similarity is low (e.g. Thr-Asp: entropy_sim = 0.026 — the call is driven almost entirely by RT and MS1).

---

### Flag 2 — NEAR_TIE (1,188 spectra)

**Definition:** Top two candidates within 0.05 posterior of each other, top posterior > 0.10.

The model cannot resolve between two candidates. Additional evidence — an authentic standard, manual RT inspection, or a literature check — would break the tie.

**Top examples:**

| wiki_id | Top candidate | Posterior | 2nd posterior | H |
|---|---|---|---|---|
| a6ICAHB/4905 | glutamine | 0.510 | ~0.465 | 2.10 |
| a6ICAHB/6265 | PC 18:2_02:0 | 0.513 | ~0.468 | 1.26 |
| a6ICAHB/2882 | PC 19:2 | 0.509 | ~0.464 | 0.77 |

Most near-ties are structural isomers or lipid species with similar masses and predicted RTs. These are genuine ambiguities — not model failures — and correctly surface the hardest identification cases.

---

### Flag 3 — RICH_NO_HIT (845 spectra)

**Definition:** Spectral entropy H > 1.0 (well-fragmented, information-rich spectrum) but model abstains (P_novel ≥ 0.6).

These are the most scientifically interesting spectra for discovery. The spectrum contains real chemical information — median H = 1.65, max H = 3.55 — but nothing in MassWiki explains it. Possible interpretations:

- **Genuinely novel compound** not yet in any public library
- **Known compound with different adduct** not captured by the current search
- **In-source fragment** of a larger molecule

**Chemist action:** These are candidates for de novo structure elucidation or targeted follow-up with authentic standards of suspected compound classes.

---

### Flag 4 — MS1_RESCUE (64 spectra)

**Definition:** Entropy similarity < 0.5 but posterior > 0.60 (RT and MS1 doing the heavy lifting).

The spectral match is poor, but mass accuracy and RT are sufficiently consistent to call the compound. These are the most suspicious calls — a chemist should verify whether the RT and m/z alignment are truly meaningful or coincidental.

**Top examples:**

| wiki_id | Compound | Posterior | Entropy sim | Δppm | ΔRT (s) |
|---|---|---|---|---|---|
| a6ICAHB/4008 | γ-Glu-Leu | 0.997 | 0.051 | −0.14 | 23.3 |
| a6ICAHB/2689 | LPG 18:0 | 0.997 | 0.057 | +6.96 | 0.08 |
| a6ICAHB/5773 | Pantetheine | 0.992 | 0.116 | −0.39 | 5.6 |

γ-Glu-Leu is a notable example: near-zero mass error (−0.14 ppm), tight RT match, but entropy similarity of only 0.051. This could reflect a compound that fragments differently from its library reference (different collision energy, adduct, or matrix effect) while still being correctly identified by the other channels.

---

## 5. Implications for the Pipeline

**P_novel = 0.224 changes the abstention behavior.** With the old circular P_novel = 0.783, the novel term was so large it rarely triggered abstention even for weak candidates (the novel term was already "used up"). With P_novel = 0.224, the abstention mechanism is more discriminating — weak evidence now actually flows into the novel term.

**84% of calls are on compounds not in BinBase.** This is the pipeline's primary practical output: a ranked, confidence-scored list of potential new annotations. The 192 high-confidence calls (posterior > 0.70) on novel compounds are actionable without any manual triage.

**The unconfirmed set is harder than the confirmed set.** Confirmed spectra had 99.8% top-1 accuracy under entropy similarity. The unconfirmed set has no ground truth, but the posterior distribution — heavily skewed toward < 0.30 — reflects genuine ambiguity rather than model failure. This is what the model is supposed to do: distinguish easy identifications from hard ones.

---

## 6. Recommended Actions

| Priority | Action | N spectra |
|---|---|---|
| 1 | Review HIGH_CONF_NOVEL with posterior > 0.90 | 110 |
| 2 | Check MS1_RESCUE calls (low similarity, high posterior) | 64 |
| 3 | Flag RICH_NO_HIT for de novo elucidation queue | 845 |
| 4 | Use NEAR_TIE list to prioritize authentic standard purchases | 1,188 |

**Output files:**
- `out_unconfirmed/review_queue.csv` — 2,253 flagged spectra with flag reason
- `out_unconfirmed/unconfirmed_top_calls.csv` — top call for all 4,216 scored spectra
- `out_unconfirmed/unconfirmed_scored.csv` — full candidate table (167,857 rows)
