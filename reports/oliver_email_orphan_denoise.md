Subject: Which UNCONFIRMED candidate bins to promote — reverse-match check (HILIC-neg)

Hi Oliver,

Result on the de-orphaning idea, focused on the **UNCONFIRMED candidate bins** — the ones that passed all the QC gates and are waiting to be accepted. Question: how many are really just an in-source fragment / adduct / isotope of a compound we've **already confirmed** (→ shouldn't be promoted) vs genuinely new?

HILIC-negative method, 26,527 candidate bins: **~10% (2,630) are a relational ion of a co-eluting confirmed compound** — 2,171 in-source fragments, 354 adducts, 105 ¹³C-isotope peaks. The other ~90% have no confirmed parent (candidate novels — still need the usual spectrum-quality check before calling real).

It's just chemistry: co-elution + a neutral-loss/adduct/¹³C mass difference + the candidate's ions contained in the parent's. E.g. m/z 102 = CO₂ loss from glutamate; 128 = NH₃ loss from glutamine; 383 = hexose loss from a hexoside; and ¹³C-isotope peaks of indoxyl sulfate / glucosamine-1-phosphate. Controls confirm it's not coincidence (scramble the masses → drops to 1.5%). It's dominated by in-source fragments, and the same fragment often spawns several candidates (β-alanine's NH₃-loss at m/z 71 recurs), exactly the clutter this collapses. 1,576 link to a named compound.

Two honest notes: at the bin level, co-elution alone is weak (the retention axis is crowded), so the firm calls are the named-parent fragments and adducts; and the "novel" residual isn't all real. Separately, the same test on the confirmed bins flags ~13% as related-ions of another confirmed bin — a library-cleanup list worth a look, since a promoted fragment becomes a "real compound" downstream.

Happy to send a batch for you to spot-check, or the confirmed-bin cleanup list — whichever's useful.

Best,
Ziyue
