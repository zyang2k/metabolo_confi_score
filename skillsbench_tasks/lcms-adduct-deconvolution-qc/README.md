# lcms-adduct-deconvolution-qc

**Domain**: Untargeted negative-mode HILIC LC-MS/MS metabolomics
**Shape**: multi-stage agentic pipeline (deconvolution → adduct QC triage → keep/flag gate)
**Scored on**: the non-derivable crux only — `adduct_category` on 14 skill-critical features + 12 keep/flag gates

## Task Summary

Given a 39-feature LC-MS list, the agent must:
1. **Deconvolve** — group co-eluting features (RT-first), match m/z differences to
   adduct/ISF relationships, designate the principal ion, infer neutral mass.
2. **Triage** each upstream `adduct_annotation` as `ok` / `isf` / `dubious` / `unassigned`.
3. **Gate** each compound `keep` (has ≥1 genuinely ok adduct) vs `flag` (only
   dubious/ISF evidence).

Deconvolution is the agentic vehicle (requires exploration, scripting, iteration);
the **score lives on the empirical taxonomy crux**, not the deconvolution.

## Why This Is Hard Without Skills

The crux is the lab's empirical **notation rule**, not chemistry: a raw molecular
formula written in brackets — `[M+C2H3O2]-` (acetate), `[M+COOH]-` (formate),
`[M+C2F3O2]-` (TFA) — is chemically identical to an OK shorthand adduct but is
flagged **dubious** (it signals a low-reliability upstream source). A frontier
model reads these as valid adducts → `ok`, and therefore **keeps** the flag-bins
that are built only from them. The deconvolution it does perfectly.

## Benchmark Results (tested)

| Condition | Reward | crux category | bin gate | deconvolution |
|---|---|---|---|---|
| Without skills (frontier model) | **0.36** (10/28) | 0/14 | 8/12 | 100% (not scored) |
| With skills (reference) | **1.00** (28/28) | 14/14 | 12/12 | 100% |

The no-skill agent used multiple tool calls and solved the deconvolution flawlessly —
it fails purely on the non-derivable taxonomy crux. This is the design intent:
**agentic depth (answers the "too few tool calls" critique) + a load-bearing skill
(C3 holds at well under 50%).**

## Criteria Assessment

| Criterion | Status | Notes |
|---|---|---|
| C1: Multi-decision workflow | ✓ | group → infer mass → triage (39×) → gate (12×); ordered, dependent stages |
| C2: Skills load-bearing | ✓ | raw-formula-vs-shorthand notation rule is empirical lab convention, not chemistry |
| C3: SOTA < 50% without skill | ✓ | Tested 0.36 (0/14 crux; all 4 flag-bins flipped to keep) |
| C4: Self-contained public data | ✓ | synthetic-from-real-compounds; no API; masses from known metabolites |
| C5: Generalizable skill | ✓ | deconvolution + adduct taxonomy apply to any LC-MS annotation workflow |
| C6: Robust verifier | ✓ | deterministic parametrized; reads `n/a` safely (keep_default_na=False) |

## Skills

- `lcms-adduct-deconvolution` (new) — RT grouping, m/z-difference matching,
  principal-ion designation, neutral-mass inference, keep/flag gate.
- `lcms-adduct-validation` (reused) — ok vs dubious via notation rules (Rule C: raw
  formula in brackets → dubious).
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
└── dev/                      # data generator + probe history (provenance)
```

## Probe History (intent → outcome → why)

This task is the survivor of a deliberate C3 search. Earlier candidates were
probed against a no-skill frontier agent and **rejected**:

- `formula-validation-seven-golden-rules` — no-skill **90%**. Principled chemistry
  (RDBE, ratios, isotopes) is in-weights; only ~5 ratio cutoffs were non-derivable
  (a transcription skill, not a real one).
- `gcms-tms-derivatization` — no-skill **100%**. Model knew the exact MeOX/TMS
  increments and group-counting.
- Pure adduct deconvolution (chemistry-only) — no-skill **100%**. Adding tool-call
  depth does nothing for C3 when the task is computable.

**Lesson:** computable chemistry fails C3 regardless of agentic depth; the
C3-passing zone is an **empirical/lab-convention crux wrapped in an agentic
pipeline**. That is exactly this task.
