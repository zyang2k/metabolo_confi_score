# Add task: lcms-adduct-deconvolution-qc

## Summary

A multi-stage **agentic** LC-MS metabolomics task: the agent deconvolves a raw
negative-mode HILIC feature list into compounds, triages each adduct annotation,
and applies a per-compound keep/flag quality gate. It is the agentic successor to
`lcms-adduct-taxonomy-lookup` (which reviewers noted was a single-shot lookup) —
same domain skill, now embedded in a pipeline that requires real exploration,
scripting, and iteration.

The task is scored **only on the non-derivable crux** (adduct QC triage + gate);
the deconvolution is the tool-call vehicle and is intentionally not scored.

## What the agent does

1. **Deconvolve** — group co-eluting features (retention-time-first), match m/z
   differences to adduct/in-source-fragment relationships, designate the principal
   ion, and infer each compound's neutral monoisotopic mass.
2. **Triage** — classify each upstream `adduct_annotation` as `ok` / `isf` /
   `dubious` / `unassigned`.
3. **Gate** — `keep` a compound if it has ≥1 genuinely `ok` adduct, else `flag`.

## Why a skill is required (C3)

The crux is an **empirical lab notation convention, not chemistry**: a raw
molecular formula written in brackets — `[M+C2H3O2]-` (acetate), `[M+COOH]-`
(formate), `[M+C2F3O2]-` (TFA) — is chemically identical to an OK shorthand adduct
but is flagged **dubious** (it marks a low-reliability upstream source). A frontier
model reads these as valid adducts → `ok`, and therefore **keeps** the flag-bins
built only from them. It performs the deconvolution perfectly; it fails purely on
the convention.

## Results (measured against a no-skill frontier agent)

| Condition | Reward | crux category | bin gate | deconvolution |
|---|---|---|---|---|
| Without skills | **0.36** (10/28) | 0/14 | 8/12 | 100% (unscored) |
| With skills (reference `solve.sh`) | **1.00** (28/28) | 14/14 | 12/12 | 100% |

The no-skill agent used multiple tool calls and solved the deconvolution flawlessly
— it fails only the empirical taxonomy crux. This is the design intent: agentic
depth **plus** a load-bearing skill, with the no-skill baseline well under 50%.

## Skills

- `lcms-adduct-deconvolution` (new) — RT grouping, m/z-difference matching,
  principal-ion designation, neutral-mass inference, keep/flag gate.
- `lcms-adduct-validation` (reused from `lcms-adduct-taxonomy-lookup`) — ok vs
  dubious notation rules.
- `lcms-isf-detection` (reused) — separate in-source fragments first.

Skills teach principles + conventions (notation rules, grouping heuristics), not
pinned answers — there is no lookup of the specific gold labels.

## Verifier

- 28 deterministic tests: 2 format + 14 crux `adduct_category` + 12 keep/flag gate
  (parametrized; one gate test per compound at its principal feature).
- Reward = `passed / total` (partial credit), written to
  `/logs/verifier/reward.txt` by `tests/test.sh`; emits on every path including
  failures.

## Data provenance & self-containment

- 39 features synthesized from real metabolite neutral masses (citric, malic,
  glucose, succinic, etc.) with realistic adduct/ISF families + interferents.
  No private lab data, no labels in the agent-visible input.
- `allow_internet = false`; no external APIs at test time.

## Reviewer notes (preempting prior feedback)

- **No skill leak:** the `Dockerfile` copies skills only to agent discovery paths
  (`~/.claude/skills`, `~/.codex/skills`, etc.) and copies `data/` separately to
  `/root/environment/data`. Skills are not visible in the agent's working data or
  `instruction.md`. `solution/solve.sh` references the skills by name/principle,
  not by importing gold.
- **Dependency pinning:** `pandas==2.2.2`, `pytest==8.4.1`, `pytest-json-ctrf==0.3.5`
  in both the Dockerfile and `test.sh`.
- **Partial-credit reward:** `passed/total`, not all-or-nothing.
- **C3 validated empirically**, not just argued — no-skill run measured at 0.36.

## Files

```
lcms-adduct-deconvolution-qc/
├── instruction.md
├── task.toml
├── eval.py
├── data/{input,gold}.csv
├── environment/{Dockerfile, data/input.csv, skills/×3}
├── solution/solve.sh
└── tests/{test.sh, test_outputs.py}
```
