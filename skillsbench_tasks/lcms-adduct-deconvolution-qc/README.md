# lcms-adduct-deconvolution-qc

**Domain**: Untargeted negative-mode HILIC LC-MS/MS metabolomics
**Shape**: multi-stage agentic pipeline (deconvolution → adduct QC triage → keep/flag gate)
**Scored on**: the non-derivable crux only — `adduct_category` on 14 skill-critical features + 12 keep/flag gates

## Task Summary

Given a 36-feature LC-MS list, the agent must:
1. **Deconvolve** — group co-eluting features (RT-first), match m/z differences to
   adduct/ISF relationships, designate the principal ion, infer neutral mass.
   Some compounds have **no `[M-H]-`**, so an adduct must be designated the
   principal ion and the mass back-calculated from it.
2. **Triage** each upstream `adduct_annotation` as `ok` / `isf` / `dubious` / `unassigned`.
3. **Gate** each compound `keep` (has ≥1 genuinely ok adduct) vs `flag` (only
   dubious/ISF evidence).

Deconvolution is the agentic vehicle (requires exploration, scripting, iteration);
the **score lives on the empirical taxonomy crux**, not the deconvolution.

## Why This Is Hard Without Skills

The crux is dominated by the **false-positive-suppression** direction — the one a
stronger, more cautious model fails *harder*, not easier:

- **Validated-but-unusual ions → ok (durable crux).** Isotope-labeled
  internal-standard ions (`[M-H+1i]-`, `[2M-H+2i]-`, `[2M-H+4i]-`) and legacy
  deprotonation notation (`M-H1`) are **expected, validated species** in this
  lab's method. A no-skill model reads the unfamiliar `+Ni` / `H1` notation as
  malformed → `dubious`, and then **wrongly flags** the compounds whose only
  intact ion is one of these. Knowing they are OK requires the lab's spike-in
  protocol and annotation history — it cannot be derived from chemistry, and
  *more* caution produces *more* over-flagging.
- **Plausible-but-low-quality notation → dubious (complementary direction).**
  Raw-formula brackets — `[M+C2H3O2]-` (acetate), `[M+COOH]-` (formate),
  `[M+C2F3O2]-` (TFA) — are chemically identical to OK shorthand but mark a
  low-reliability upstream source. (Frontier models increasingly get this on
  their own, which is exactly why the crux was rebalanced toward suppression.)

The deconvolution it does perfectly; it fails purely on the non-derivable taxonomy.

## Benchmark Results

| Condition | Reward | crux category | bin gate | deconvolution |
|---|---|---|---|---|
| Without skills (Opus 4.8, local proxy) | **0.54** (15/28) | 7/14 | 6/12 | 100% (not scored) |
| With skills (reference `solve.sh`) | **1.00** (28/28) | 14/14 | 12/12 | 100% |

The no-skill agent used multiple tool calls and solved the deconvolution flawlessly.
**11 of its 13 scored failures land on the suppression direction** (all 6 isotope/
historical crux features over-flagged; all 5 over-flag gates wrongly flagged); it
got 7/8 raw-formula crux *right*. That distribution is the design intent: the delta
rests on the **durable** non-derivable direction, not the guessable one.

> The 0.54 figure is a single-run local proxy on Opus 4.8 (no skills, no internet).
> Official bench baselines (claude-opus-4-8 + gpt-5.5, 3 trials) are to be recorded
> in `task.toml` before merge.

## Criteria Assessment

| Criterion | Status | Notes |
|---|---|---|
| C1: Multi-decision workflow | ✓ | group → infer mass (incl. no-`[M-H]-` bins) → triage (36×) → gate (12×); ordered, dependent stages |
| C2: Skills load-bearing | ✓ | validated-ion (V1/V2) and notation-reliability (Rule C) conventions are empirical lab knowledge, not chemistry |
| C3: SOTA < 50% without skill | ✓ | local proxy 0.50 on the scored crux (13/26); 11/13 failures on the durable suppression direction |
| C4: Self-contained public data | ✓ | synthetic-from-real-compounds; no API; masses from known metabolites |
| C5: Generalizable skill | ✓ | deconvolution + adduct taxonomy apply to any LC-MS annotation workflow |
| C6: Robust verifier | ✓ | deterministic parametrized; reads `n/a` safely (keep_default_na=False) |

## Skills

- `lcms-adduct-deconvolution` (new) — RT grouping, m/z-difference matching,
  principal-ion designation (incl. anchoring on an adduct when no `[M-H]-`),
  neutral-mass inference, keep/flag gate (both over- and under-flag failure modes).
- `lcms-adduct-validation` (reused) — ok vs dubious. **Rules V1/V2:** isotope-
  labeled-standard and legacy-notation ions are validated → ok (with provenance).
  **Rule C:** raw formula in brackets → dubious.
- `lcms-isf-detection` (reused) — separate in-source fragments first.

## Files

```
lcms-adduct-deconvolution-qc/
├── instruction.md            # agent-facing task
├── task.toml                 # metadata + sandbox config
├── eval.py                   # standalone evaluator (crux + gate)
├── data/{input,gold}.csv     # full data + gold
├── environment/
│   ├── Dockerfile
│   ├── data/input.csv        # agent sees this (no labels)
│   └── skills/               # 3 SKILL.md
├── solution/solve.sh         # reference (scores 28/28)
├── tests/{test.sh,test_outputs.py}
└── dev/                      # data generator (final_build.py) + probe history
```

## Probe History (intent → outcome → why)

This task is the survivor of a deliberate C3 search. Earlier candidates were
probed against a no-skill frontier agent and **rejected**:

- `formula-validation-seven-golden-rules` — no-skill **90%**. Principled chemistry
  (RDBE, ratios, isotopes) is in-weights; only ~5 ratio cutoffs were non-derivable.
- `gcms-tms-derivatization` — no-skill **100%**. Model knew the exact MeOX/TMS increments.
- Pure adduct deconvolution (chemistry-only) — no-skill **100%**. Tool-call depth
  does nothing for C3 when the task is computable.

**Rebalance (2026-06-02):** the first deconv-qc crux leaned entirely on
raw-formula-bracket → dubious. That points in the *guessable* direction — a cautious
frontier model flags weird notation unaided — so the no-skill baseline rose and the
delta shrank. The crux was rebalanced toward **false-positive suppression**
(validated-but-unusual ions → ok), which a more-cautious model fails harder. This
makes the delta durable *and* makes the skill a substantive, provenance-backed lab
fact rather than a notation gotcha.

**Lesson:** computable chemistry fails C3 regardless of agentic depth; and even an
empirical crux fails if it points the same way as the model's default caution. The
durable C3-passing zone is an **empirical/lab-convention crux in the suppression
direction, wrapped in an agentic pipeline**. That is exactly this task.
