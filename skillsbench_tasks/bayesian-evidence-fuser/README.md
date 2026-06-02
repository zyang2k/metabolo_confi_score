# bayesian-evidence-fuser

**Domain**: Untargeted negative-mode LC-MS/MS metabolomics (human blood plasma)
**Shape**: single-feature evidence fusion — combine 3 independent channels into a normalized posterior, then rank
**Scored on**: 5 per-candidate posteriors + normalization + ranking + 2 format gates (9 tests, reward = passed/total)

## Task Summary

Given 5 putative metabolite identities for one detected feature — each with a
`Prior_Probability`, `Observed_mz`, `Mass_Error_ppm`, and `MSMS_Similarity_Score`
— the agent must:

1. Convert each evidence channel to a likelihood **per the lab QC SOP**.
2. Fuse them via discrete Bayes' rule: `unnormalized = prior * mass_lik * msms_lik`.
3. Normalize across candidates so posteriors sum to 1.0.
4. Rank descending and emit `/root/ranked_posteriors.json`.

## Why This Is Hard Without Skills

The Bayesian fusion math is textbook and fully derivable — that is **not** the
crux. The skill delta lives entirely in two **non-derivable lab QC conventions**:

- **Channel B — m/z-zoned mass-accuracy sigma.** The Gaussian sigma is not a
  single constant; it is zoned by observed m/z (`<250 -> 4.0`, `250..450 -> 2.0`,
  `>450 -> 3.0`), reflecting measured calibration drift at the mass extremes. A
  no-skill agent picks one global sigma.
- **Channel C — MS/MS reliability floor.** Cosines below `0.50` are replaced by a
  fixed `0.01` floor (background-dominated, no compound-specific information). A
  no-skill agent uses the raw cosine.

Neither can be inferred from chemistry or the data alone. A textbook fusion gets
every posterior wrong by far more than the `1e-4` tolerance, and **mis-ranks the
bottom two candidates** — the `LysoPC_16:0 <-> D-Glucose` flip is driven entirely
by the two cruxes (Glucose's `0.48` cosine is floored to `0.01`; LysoPC's `4 ppm`
error is judged against `sigma=3.0`, not `2.0`).

## Benchmark Results (probed 2026-06-02, claude-opus-4-8)

| Condition | Reward | posterior values | normalization | ranking |
|---|---|---|---|---|
| Without skills (live agent) | **0.33** (3/9) | 0/5 | pass | wrong |
| Without skills (assumption sweep) | **0.33–0.44** (3–4/9) | 0/5 | pass | varies |
| With skills (live agent) | **1.00** (9/9) | 5/5 | pass | correct |
| Oracle (`solution/solve.sh`, in-container) | **1.00** (9/9) | 5/5 | pass | correct |

The live no-skill agent performed textbook fusion (chose a single global
`sigma=5 ppm`, used the raw MS/MS cosine) and **explicitly noted that no SOP was
available to it**, getting every posterior wrong and mis-ranking the candidates.

Because a single agent run is one stochastic draw, the no-skill ceiling is also
**bounded analytically** over the assumptions a skill-less agent could plausibly
adopt (`dev/build.py` reproduces this):

- any constant `sigma` in 1–10 ppm + raw cosine → **3/9** (0/5 values, wrong rank)
- m/z-zoned sigma but **no** MS/MS floor (half-skill) → **3/9**
- MS/MS floor but constant sigma (half-skill) → **4/9** (0/5 values; ranking
  happens to land correct but no posterior is within tolerance)

Only applying **both** cruxes reaches 9/9 — the per-candidate values are 0/5 in
every partial case because normalization couples all five candidates, so a single
posterior cannot be correct to 1e-4 unless both conventions are applied
throughout. **Skill delta ≈ +0.6**, no-skill capped well below 50% — a genuine C3
pass, not a knife-edge.

> **If the no-skill arm ever creeps up** (a model guessing "calibration degrades
> at mass extremes"), sharpen the sigma zones to non-monotonic / counter-intuitive
> values so they cannot be reverse-engineered from chemistry. Re-probe with
> `dev/build.py` + an isolated no-skill agent run.

## Files

- `instruction.md` — runtime prompt (math-free; points at the lab SOP).
- `environment/skills/bayesian-evidence-fuser/SKILL.md` — the SOP (the crux).
- `environment/data/input.csv` — the 5-candidate input.
- `environment/Dockerfile` — Ubuntu + python3 + pandas/scipy + pytest.
- `solution/solve.sh` — oracle reference implementation.
- `tests/test_outputs.py`, `tests/test.sh` — pytest suite + reward wrapper.
- `eval.py` — standalone scorer (mirrors the pytest checks).
- `dev/build.py` — regenerates the input CSV and prints gold + textbook answers.

## Gold Values

```
L-Phenylalanine 0.3548799017   (rank 1)
Citric_acid     0.3274002778   (rank 2)
Adenosine       0.1988559215   (rank 3)
LysoPC_16:0     0.1143642591   (rank 4)
D-Glucose       0.0044996399   (rank 5)
```
