# Using the reverse (containment) match to clean up in-source noise in regular studies

**To:** Oliver
**From:** Ziyue
**Date:** 2026-06-16

## The question, answered

**Can the reverse (containment) spectral match help reduce the impact of noise ions among the
hundreds of MS/MS that don't match a bin in a regular study? — Yes.** On a real LCBinBase study,
~20% of those unconfirmed spectra (a conservative ~12% if we require the parent to be detected in the
same run) are in-source fragments / adducts / isotopes of co-eluting confirmed compounds, and we can
label and collapse them automatically. The same test, run on the confirmed bins themselves, also
flags ~17% as candidate related-ions of other confirmed bins (a library-audit use — last section).

## The short version

In a regular LC-MS/MS study, hundreds of MS/MS spectra never make it into a bin — they float
around unconfirmed. Many of these are not real compounds; they are in-source fragments, adducts, and
isotopes of compounds that *are* binned. I tested whether a **reverse (containment) spectral match**
can identify and label these automatically, picking up on your earlier suggestion to add
CAMERA/RAMClust-style in-source-fragment detection.

It works. On a real HILIC-negative bile-acid study (50 injections, ~24,000 floating MS/MS), **~20% of
the floating spectra are explained as an in-source fragment / adduct / isotope of a co-eluting
confirmed bin**, and we can label each one with its parent compound. Requiring the parent to actually
be detected in the same injection gives a more conservative **~12%**. The rate is stable across
injections (18–22%), and it survives the negative controls below.

## Why *reverse* match (and not the normal forward cosine)

An in-source fragment is a *subspectrum* of its parent — the parent fragments further in the source,
so it carries the fragment's ions plus many more. A forward cosine is dragged down by all those extra
parent ions and misses the relationship. The one-sided **containment** ("are the fragment's ions all
present in the parent?") stays high. That is exactly what the reverse/NIST-style match measures, and
it is the right tool here. (We separately confirmed the reverse match is *not* useful for deciding
compound identity — see the last section.)

The rule, in chemical terms: a floating spectrum **O** is flagged as a related ion of a confirmed bin
**C** when (1) they co-elute (≤4 RI units), (2) the precursor mass difference **m(C) − m(O)** is a
neutral loss, an adduct difference, or an isotope spacing, and (3) O's ions are contained in C's
spectrum (and, for fragments, O's precursor appears as a peak in C).

## Examples from the run (all real)

| Floating m/z (RI) | Explained as | Co-eluting parent |
|---|---|---|
| 147.030 (113) | loss of CO₂ (43.99) | **citric acid** (191.020) |
| 93.035 (14) | loss of CO₂ | **salicylic acid** (137.025) |
| 129.020 (86) | loss of CO₂ | **dehydroascorbic acid** (173.009) |
| 87.045 (38) | loss of CO₂ | **glutaric acid** (131.035) |
| 107.050 (14) | loss of CO (27.99) | **phenylacetic acid** (135.046) |
| 178.051 (58) | loss of CO | **hydroxykynurenine** (206.046) |
| 147.030 (122) | loss of CH₂O (30.01) | **gulonolactone** (177.041) |

These are textbook in-source losses — the carboxylic acids shedding CO₂, etc. Fittingly for a
bile-acid method, the most frequent relations after H₂O/CO₂ are **SO₃ losses** and **³⁴S isotopes**
(sulfated bile acids), which the method picks up on its own. Overall, **63% of the explained spectra
link to a named compound**; the rest link to confirmed-but-unnamed bins.

## How we know it isn't fooling itself

Three negative controls, all behaving as they should:

- **Scramble the precursor masses** (so the mass differences are no longer real losses): the explained
  rate collapses from 20% to **3%**. So the chemistry — not coincidence — is doing the work.
- **Use nonsense neutral losses** (non-physical masses): **0%** fire. So enlarging the loss list
  doesn't manufacture false hits.
- **Scramble retention times** (break co-elution): the coincidental background drops to ~4% once we
  also require the parent to be present in the same injection. The real signal sits well above it.

The containment is high where it should be (median 0.92), and the rate is consistent injection to
injection rather than driven by one outlier run.

## A side finding: the confirmed bins themselves

Confirmed bins were the trusted reference above, but I ran the same test *among the confirmed bins* —
and **~17% (1,108 of 6,398) look like an in-source fragment / adduct / isotope of another co-eluting
confirmed compound.** A recurring in-source fragment appears in as many samples as its parent, so it
can clear the recurrence threshold and get promoted to its own "real" bin — and that seems to be
happening. Examples: **genistein** as the glucuronic-acid loss of **genistein-4′-glucuronide**,
**xanthine** as the pentose loss of **xanthosine**, **4-methylcatechol** as the SO₃ loss of **guaiacol
sulfate**. ~280 of the flagged bins look like *isotope peaks* of another bin, which arguably should
not be separate compounds at all.

Reassuringly, the method independently re-finds 159 cases you'd already flagged (`yy_` or named "in
source to…") without reading the names — good evidence it's catching the right thing. The remaining
~950 are not curator-flagged and would be a review list. **This is a triage list, not a verdict** —
some will be legitimate separate adducts or coincidental co-elutions of real isobars — but it's
arguably the higher-stakes use: a mis-promoted fragment or isotope becomes a "real compound" in
downstream biology and quantification.

## What it does *not* do (so we don't oversell it)

- **It does not find novel compounds.** The ~80% of floating spectra it leaves unexplained are mostly
  low-recurrence features and ordinary background, not 19,000 new metabolites. Those still need the
  within-spectrum quality filter (S_norm / the noise features) before anything is called real. What
  this tool delivers cleanly is the ~12–20% it can *remove and explain*.
- **It cannot decide compound identity.** We tested the reverse match for rescuing contaminated
  library matches on lipidomics: it does recover spectra whose entropy similarity was pulled down by
  co-isolation, but it rates a wrong-retention-time match just as highly as a correct one — so it is a
  flag/explanation, not an identity score.
- **Within-run only.** Co-elution is meaningful per injection, which is exactly the right scope here.

One incidental finding worth flagging: **BinBase/CARROT currently stores no in-source-fragment links
for this HILIC method** (the `fragment_of` field is empty for every bin and floating feature). So this
isn't duplicating an existing annotation — it's filling a gap.

## Suggested use

Per study: auto-label each floating MS/MS as "in-source fragment of [compound]" (or adduct/isotope),
link it to the parent, and hand back the unexplained residual as a cleaner candidate worklist. This
shrinks the "what are all these?" pile by ~1-in-5 with an explanation attached, and complements the
self-noise (S_norm) filter — that one catches a spectrum's own junk; this one catches relational junk
from co-eluting neighbours.

If useful, next steps could be: (a) you spot-check ~20 of the calls (especially the lower-containment
ones) to confirm the labels; (b) I produce the **confirmed-bin audit list** (isotopes + high-containment
fragments not already flagged) as a separate sheet; (c) we run it on a few more methods (HILIC-pos,
C18, lipidomics) to confirm it generalizes and see the loss profiles shift with chemistry; (d) we fold
it into the pipeline as a per-study step. Happy to do any of these — let me know which you'd like to see.
