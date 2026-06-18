# Reverse-match triage of UNCONFIRMED candidate bins: which to promote?

**To:** Oliver
**From:** Ziyue
**Date:** 2026-06-16

## The question, answered

When BinBase generates a candidate bin but hasn't accepted it (`UNCONFIRMED`), some of those
candidates aren't new compounds at all — they're in-source fragments, adducts, or isotope peaks of
compounds we've **already confirmed**. Using a reverse (containment) spectral match, **~10% of the
UNCONFIRMED candidates are exactly that, and should not be promoted; the other ~90% have no confirmed
parent (candidate novels).** This picks up your suggestion to add CAMERA/RAMClust-style in-source
detection, aimed at the promote/don't-promote decision.

## Why this population

I deliberately ran it on `UNCONFIRMED` candidates — these have *passed* every QC gate (ion count,
S/N, scan count, ISTD coverage, clean-spectra) and just await acceptance, so "is this an artifact of
a compound we already have?" is a real decision. (The much larger `INVALID_TARGET` pile is mostly
run-level ISTD-coverage rejection — real compounds from QC-failed injections, not noise — so it's the
wrong population for this and I left it aside.)

## Result (HILIC-negative method)

Of **26,527** candidate bins, **2,630 (9.9%)** are a relational ion of a co-eluting confirmed
compound:
- **2,171 in-source fragments** (neutral-loss relationships) — the dominant class
- **354 adducts**
- **105 ¹³C-isotope peaks** (small, but the most clear-cut: an isotope peak should never be its own compound)

**1,576 of these link to a *named* confirmed compound**, so each comes with a concrete reason. The
remaining ~90% have no confirmed parent → candidate novels (which still need the usual
spectrum-quality check before being called real — that residual is not automatically "new compounds").

## The rule, and examples

A candidate is tied to a confirmed bin when they co-elute, the precursor mass difference is a
neutral-loss / adduct / ¹³C-isotope spacing, and the candidate's ions are contained in the parent's.
(Containment is the right measure — a fragment is a subspectrum of its parent, so a forward cosine
gets dragged down by the parent's extra ions and misses it.)

- m/z **102.056** = loss of CO₂ from **glutamate** (146.046)
- m/z **71.014** = loss of NH₃ from **β-alanine** (88.040)
- m/z **96.970** = loss of the hexose from **mannose-6-phosphate** (259.022) → phosphate

Tellingly, the *same* in-source fragment often spawns several candidate bins — β-alanine's NH₃-loss
at 71.014 shows up repeatedly — which is exactly the redundant clutter this collapses.

## How far to trust it — and which calls

We test reliability by scrambling the precursor masses so any match is pure coincidence, then
re-running: real flags should far exceed the scrambled ones. Doing this **per relation type** gives a
clear, honest picture — trust is *not* uniform, it depends on the relation and the chromatography:

- **In-source fragments — trust them everywhere.** ~89% non-coincidental on HILIC, ~96% on lipids.
  The evidence is rich (the candidate's *whole* spectrum has to sit inside the parent's, not just a
  mass gap), so it's hard to fire by accident. This is the dependable core, and the bulk of the flags.
- **Adducts & isotopes — trust them on the C18 (lipid) runs, treat them as *leads* on HILIC.** They
  hang on a single mass gap (Na–H, +1.00336…), so they need co-elution to be discriminating. On
  reverse-phase, where lipids spread across a wide retention window, they're reliable (~85–93%). On
  HILIC the retention axis is crowded — a co-eluting neighbour at the right mass gap turns up by
  chance constantly — so those calls are only ~65% (adducts) and essentially **coin-flip for HILIC
  isotopes**. The isotope intensity check that could fix this can't be computed here (a candidate and
  its parent are separate averaged entries, often from different samples, so their intensities aren't
  a valid isotope ratio — that needs raw MS1 data we don't pull). HILIC isotopes are a tiny slice
  anyway, so this doesn't dent the result.

So the trustworthy yield is **~8% on HILIC (essentially all in-source fragments)** and **~10–18% on
lipids (fragments + adducts + isotopes alike)**.

## Generalizes to lipidomics (C18) — the stronger case

Running the same thing on the C18 lipidomics method (360,231 candidate bins, 23,059 confirmed
lipids — lipidomics generates ~14× more candidates than HILIC) flags **18.2%** as a relational ion
of a confirmed lipid (in-source fragments, the acetate family `[M+OAc]⁻`/`[M−H]⁻`, ¹³C-isotope peaks
of high-carbon lipids, Na adducts). Examples: PE 18:0_18:2 − H₂O; PC 18:0_18:1 and SM d20:0_22:2 via
acetic acid; FAHFA + Na.

Two things make lipidomics the better fit: (1) the rate is ~2× HILIC, and (2) the controls are
*tighter* — scrambling retention drops it to ~4% (vs HILIC's ~8%), because reverse-phase spreads
lipids across a wide, uncrowded retention axis, so co-elution is strongly discriminating. One honest
note: adding lipid-specific headgroup and fatty-acyl losses to the rule barely changed the number
(+0.6 pp), so these artifacts are mostly generic adducts / isotopes / small losses rather than
classic acyl-chain fragmentation.

In **positive mode** (C18-pos) it flags **10.2%**, and the makeup shifts exactly as the chemistry
predicts — **adducts dominate** (TG as [M+NH₄]⁺, PC as [M+Na]⁺/[M+K]⁺), where negative mode was led
by in-source/water/acetate losses. That the relation profile tracks the ionization mode (losses in
neg, adducts in pos) is itself a reassurance the tool is reading real chemistry, not noise.

## A side finding: the confirmed library

Running the same test *among the confirmed bins* flags **~13% (858/6,398)** as a fragment / adduct /
¹³C-isotope of another confirmed bin — i.e. some in-source fragments already got promoted to their own
"real" bin (a recurring fragment recurs as often as its parent, so it clears the threshold). It
re-finds cases you'd already flagged (`yy_` / "in source to…"), so it's catching the right thing.
This is a separate library-cleanup list — a triage, not a verdict.

## Suggested use

For the pending pile: auto-label each `UNCONFIRMED` candidate as "in-source fragment / ¹³C-isotope /
adduct of [compound] — don't promote" vs "no confirmed parent — review," so the redundant artifacts
stop competing for promotion and the genuinely-new candidates surface.

Next steps if useful: (a) you spot-check ~20 of the flagged candidates (especially the named-parent
fragments and any low-mass marker ions); (b) I produce the confirmed-bin cleanup list; (c) run it on
more methods (HILIC-pos, C18, lipidomics) to confirm it generalizes. Let me know.
