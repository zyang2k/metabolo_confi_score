Subject: Reverse-match cleanup of unconfirmed MS/MS — HILIC-neg pilot

Hi Oliver,

You asked whether a reverse (containment) match can reduce the impact of noise ions among the hundreds of MS/MS that don't match a bin in a regular study. **Yes** — about 1-in-5 of those unconfirmed spectra are in-source fragments / adducts / isotopes of a co-eluting confirmed compound, and we can label each with its parent and collapse it.

On a real LCBinBase study (HILIC-neg bile acids, 50 injections, ~490 unconfirmed MS/MS per run): **~20% explained** (a conservative ~12% if I require the parent detected in the same run), steady run to run. The rule is just chemistry — co-elution + a neutral-loss/adduct/isotope mass difference + the floating spectrum's ions contained in the parent's. E.g. m/z 147 = CO₂ loss from citric acid; 93 = CO₂ loss from salicylic acid; plus lots of SO₃ losses and ³⁴S isotopes from the sulfated bile acids. Controls confirm it's not coincidence (scramble the masses → drops 20%→3%; nonsense losses → 0%).

One thing worth flagging: running the same test *among the confirmed bins* flags ~17% as a fragment/adduct/isotope of another confirmed bin — i.e. some in-source fragments got promoted to their own "real" bin (e.g. genistein vs its glucuronide). It re-finds cases you'd already flagged, so it's catching the right thing; it's a review list, not a verdict.

Happy to send a batch for you to spot-check, build that confirmed-bin audit list, or run more methods — whichever's useful.

Best,
Ziyue
