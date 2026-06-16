"""
bench_reverse_rescue_lipids.py
==============================
Tests the reverse-score rescue hypothesis on lipidomics data:

    "Does reverse cosine rescue contaminated-but-correct matches that
     forward (entropy) similarity wrongly rejects?"

Rescue zone = forward_to_correct < 0.7  AND  reverse_to_correct > 0.7.

Labels (per user, 2026-06-15):
  * TP        = compound (IK14) annotated with >=2 distinct adducts within
                RT_WINDOW seconds at the same RT (multi-adduct corroboration).
  * FP        = curator `name` starts with 'yy_'.
  * ignore    = curator `name` starts with 'zz_'.
  * unlabelled= everything else.

Reverse target (user choice): the CURATOR-CORRECT compound's library entry,
matched by (name, adduct) then IK14 — NOT merely the top reference hit. So
forward/reverse both measure the query against the spectrum it SHOULD match.

Pipeline (all cached / resumable):
  Phase 1  get_data per annotated wiki_id  -> query peaks + matched library_wiki_id
  Phase 2  fetch library peaks for those matched ids
  Phase 3  compute forward+reverse via build_features_v2.compute_ms2_scores,
           build rescue table, report rescued TP + false rescues (precision).

Token: reads /tmp/.mw_token (refresh by overwriting that file).
Run:   python code/bench_reverse_rescue_lipids.py            # both polarities
       python code/bench_reverse_rescue_lipids.py --limit 50 # smoke test
"""
from __future__ import annotations
import argparse, json, os, sys, time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "code"))          # build_features_v2 lives at code/ root
from build_features_v2 import compute_ms2_scores  # noqa: E402

BASE     = "https://masswiki.us-west-2.elasticbeanstalk.com/analysis/get_data"
TOKEN_FILE = "/tmp/.mw_token"
WORKERS, RPS, TIMEOUT = 6, 6.0, 30
RT_WINDOW = 10.0          # seconds; "same RT bin" for multi-adduct corroboration
FWD_CUT, REV_CUT = 0.7, 0.7

DATA = ROOT / "data"
GETDATA_CACHE = DATA / "lipid_getdata_cache.json"     # wiki_id -> {qpeaks, matched}
MONA_CACHE = DATA / "mona_lib_peaks_cache.json"       # legacy: MoNA accession -> peaks
WEB_CACHE  = DATA / "lipid_web_libpeaks_cache.json"   # "source:accession" -> peaks
NIST_CACHE = DATA / "nist_lipid_peaks_cache.json"     # ik14 -> [(precursor_type, peaks)]
NIST_MSP   = DATA / "reference_db" / "nist23-raw.msp"
OUT_TABLE = DATA / "reverse_rescue_lipids.csv"
MONA_REST = "https://mona.fiehnlab.ucdavis.edu/rest/spectra"
GNPS_REST = "https://external.gnps2.org/gnpsspectrum"
MBEU_RAW  = "https://raw.githubusercontent.com/MassBank/MassBank-data/main"
WEB_SOURCES = {"mona", "gnps", "mbeu"}


# ─────────────────────────── helpers ────────────────────────────────────────
def ik14_from_smiles(smi):
    try:
        from rdkit import Chem
        from rdkit.Chem import inchi
        m = Chem.MolFromSmiles(str(smi))
        if m is None:
            return None
        return inchi.MolToInchiKey(m)[:14]
    except Exception:
        return None


def build_labels(df: pd.DataFrame) -> pd.DataFrame:
    """Attach label column: TP / FP / ignore / unlabelled (per user defs)."""
    df = df.copy()
    nm = df["name"].astype(str)
    df["is_yy"] = nm.str.startswith("yy_")
    df["is_zz"] = nm.str.startswith("zz_")
    df["ik14"] = df["inchikey"].astype(str).str[:14]

    # TP: annotated, non-flagged, IK14 with >=2 distinct adducts within RT_WINDOW
    ann = df[df["name"].notna() & ~df["is_yy"] & ~df["is_zz"] & (df["ik14"] != "nan")]
    tp_wids = set()
    for ik, g in ann.groupby("ik14"):
        g = g.sort_values("rt")
        for _, row in g.iterrows():
            near = g[(g["rt"] - row["rt"]).abs() <= RT_WINDOW]
            if near["adduct"].nunique() >= 2:
                tp_wids.update(near["wiki_id"].tolist())
                break

    def lab(r):
        if r["is_zz"]:
            return "ignore"
        if r["is_yy"]:
            return "FP"
        if r["wiki_id"] in tp_wids:
            return "TP"
        return "unlabelled"

    df["label"] = df.apply(lab, axis=1)
    return df


def session(token):
    import requests
    from requests.adapters import HTTPAdapter
    from urllib3.util.retry import Retry
    s = requests.Session()
    s.headers.update({"Accept": "application/json",
                      "Authorization": f"Bearer {token}",
                      "User-Agent": "Mozilla/5.0"})
    retries = Retry(total=4, backoff_factor=0.5,
                    status_forcelist=(429, 500, 502, 503, 504), allowed_methods=["GET"])
    ad = HTTPAdapter(max_retries=retries, pool_connections=64, pool_maxsize=64)
    s.mount("https://", ad)
    return s


def match_correct_hit(hits, anno_name, anno_adduct, anno_ik14):
    """Pick the identity_search hit that == the curator annotation.
    Priority: exact (name, adduct) -> IK14(smiles) + adduct -> IK14 only -> name only."""
    anom = str(anno_name).strip().lower()
    aad = str(anno_adduct).strip()
    by_name_add, by_ik_add, by_ik, by_name = None, None, None, None
    for h in hits:
        if not isinstance(h, dict):
            continue
        hn = str(h.get("name") or "").strip().lower()
        ha = str(h.get("adduct") or "").strip()
        hik = (h.get("inchikey") or "")[:14] or ik14_from_smiles(h.get("smiles"))
        if hn == anom and ha == aad and by_name_add is None:
            by_name_add = h
        if anno_ik14 and hik == anno_ik14 and ha == aad and by_ik_add is None:
            by_ik_add = h
        if anno_ik14 and hik == anno_ik14 and by_ik is None:
            by_ik = h
        if hn == anom and by_name is None:
            by_name = h
    return by_name_add or by_ik_add or by_ik or by_name


# ─────────────────────────── phase 1: get_data ──────────────────────────────
def phase1_fetch(targets: pd.DataFrame, token):
    cache = json.loads(GETDATA_CACHE.read_text()) if GETDATA_CACHE.exists() else {}
    todo = targets[~targets["wiki_id"].isin(cache.keys())]
    print(f"[P1] get_data: {len(cache)} cached, {len(todo)} to fetch")
    if todo.empty:
        return cache
    s = session(token)
    last = [0.0]; lock_interval = 1.0 / RPS
    err401 = [0]

    def worker(rec):
        wait = lock_interval - (time.time() - last[0])
        if wait > 0:
            time.sleep(wait)
        last[0] = time.time()
        wid = rec["wiki_id"]
        try:
            r = s.get(BASE, params={"wiki_id": wid, "source": "binbase", "isPublic": "false"}, timeout=TIMEOUT)
        except Exception as e:
            return wid, None, f"exc:{e}"
        if r.status_code == 401:
            err401[0] += 1
            return wid, None, "401"
        if r.status_code != 200:
            return wid, None, f"http{r.status_code}"
        p = r.json()
        qpk = (p.get("spectrum") or {}).get("peaks")
        hits = ((p.get("analysis") or {}).get("reference_library") or {}).get("identity_search") or []
        h = match_correct_hit(hits, rec["name"], rec["adduct"], rec["ik14"]) if hits else None
        matched = None
        if h is not None:
            matched = {"library_wiki_id": h.get("library_wiki_id"),
                       "mona_id": h.get("id"), "db": h.get("db"),
                       "forward_esim": h.get("entropy_similarity"),
                       "name": h.get("name"), "adduct": h.get("adduct")}
        return wid, {"qpeaks": qpk, "matched": matched, "n_ident": len(hits)}, None

    recs = todo.to_dict("records")
    t0 = time.time()
    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        futs = [ex.submit(worker, r) for r in recs]
        for i, fut in enumerate(as_completed(futs), 1):
            wid, res, err = fut.result()
            if res is not None:
                cache[wid] = res
            if err401[0] >= 5:
                GETDATA_CACHE.write_text(json.dumps(cache))
                raise SystemExit("\n*** TOKEN EXPIRED (401) — refresh /tmp/.mw_token and rerun. "
                                 f"Saved {len(cache)} cached so far. ***")
            if i % 200 == 0 or i == len(recs):
                GETDATA_CACHE.write_text(json.dumps(cache))
                rate = i / max(time.time() - t0, 1e-6)
                print(f"  [P1] {i}/{len(recs)}  matched={sum(1 for v in cache.values() if v.get('matched'))}  {rate:.1f}/s")
    GETDATA_CACHE.write_text(json.dumps(cache))
    return cache


# ─────────────────────────── phase 2: library peaks (MoNA REST) ─────────────
def _parse_mona_spectrum(spec_str):
    """MoNA 'mz:int mz:int ...' string -> [[mz,int],...]."""
    out = []
    for tok in str(spec_str).split():
        if ":" in tok:
            mz, it = tok.split(":")
            out.append([float(mz), float(it)])
    return out or None


def _fetch_one_web(s, source, acc):
    """Return peak list for a library accession from its public source."""
    import json as _json
    try:
        if source == "mona":
            r = s.get(f"{MONA_REST}/{acc}", timeout=TIMEOUT)
            return _parse_mona_spectrum(r.json().get("spectrum")) if r.status_code == 200 else None
        if source == "gnps":
            r = s.get(f"{GNPS_REST}?SpectrumID={acc}", timeout=TIMEOUT)
            if r.status_code != 200:
                return None
            pj = (r.json().get("spectruminfo") or {}).get("peaks_json")
            pk = _json.loads(pj) if isinstance(pj, str) else pj
            return [[float(m), float(i)] for m, i in pk] if pk else None
        if source == "mbeu":
            folder = acc.split("-")[1]
            r = s.get(f"{MBEU_RAW}/{folder}/{acc}.txt", timeout=TIMEOUT)
            if r.status_code != 200:
                return None
            pk, inpk = [], False
            for ln in r.text.splitlines():
                if ln.startswith("PK$PEAK:"):
                    inpk = True; continue
                if inpk:
                    p = ln.strip().split()
                    if len(p) >= 2 and p[0][:1].isdigit():
                        pk.append([float(p[0]), float(p[1])])
                    else:
                        break
            return pk or None
    except Exception:
        return None
    return None


def phase2_web(needs):
    """needs: set of (source, accession) for sources mona/gnps/mbeu (exact spectrum)."""
    import requests
    cache = json.loads(WEB_CACHE.read_text()) if WEB_CACHE.exists() else {}
    # migrate legacy MoNA-only cache (keyed by bare accession) into unified keys
    if MONA_CACHE.exists():
        for k, v in json.loads(MONA_CACHE.read_text()).items():
            cache.setdefault(f"mona:{k}", v)
    todo = [(src, acc) for src, acc in needs if acc and f"{src}:{acc}" not in cache]
    from collections import Counter
    print(f"[P2-web] {len(needs)} needed, {len(todo)} missing  by-source={dict(Counter(s for s,_ in todo))}")
    if todo:
        s = requests.Session(); s.headers.update({"User-Agent": "Mozilla/5.0"})
        last = [0.0]; lock_interval = 1.0 / RPS

        def worker(item):
            src, acc = item
            wait = lock_interval - (time.time() - last[0])
            if wait > 0:
                time.sleep(wait)
            last[0] = time.time()
            return src, acc, _fetch_one_web(s, src, acc)

        t0 = time.time(); n_err = 0
        with ThreadPoolExecutor(max_workers=WORKERS) as ex:
            futs = [ex.submit(worker, it) for it in todo]
            for i, fut in enumerate(as_completed(futs), 1):
                src, acc, pk = fut.result()
                cache[f"{src}:{acc}"] = pk
                if not pk:
                    n_err += 1
                if i % 200 == 0 or i == len(todo):
                    WEB_CACHE.write_text(json.dumps(cache))
                    print(f"  [P2-web] {i}/{len(todo)}  miss={n_err}  {i/max(time.time()-t0,1e-6):.1f}/s")
        WEB_CACHE.write_text(json.dumps(cache))
    return cache


# ─────────────────────────── phase 2b: NIST23 from local MSP ─────────────────
def _norm_adduct(a):
    return str(a).replace(" ", "").replace("[", "").replace("]", "").strip().lower()


def phase2_nist(needed_ik14):
    """Stream the 2.7GB NIST23 MSP once; keep entries whose InChIKey14 is needed.
    Returns {ik14: [(precursor_type, [[mz,int],...]), ...]}."""
    if NIST_CACHE.exists():
        idx = json.loads(NIST_CACHE.read_text())
        if set(needed_ik14).issubset(idx.keys()) or idx.get("__complete__"):
            print(f"[P2-nist] cache hit ({len(idx)} ik14 groups)")
            return idx
    need = set(needed_ik14)
    print(f"[P2-nist] streaming {NIST_MSP.name} for {len(need)} needed IK14 ...")
    idx = {}
    cur = {}; peaks = []; in_peaks = 0
    n_entry = 0

    def flush():
        ik = (cur.get("inchikey") or "")[:14]
        if ik in need and peaks:
            idx.setdefault(ik, []).append((cur.get("precursor_type", ""), peaks[:]))

    with open(NIST_MSP, "r", errors="ignore") as fh:
        for ln in fh:
            ls = ln.strip()
            if in_peaks > 0:
                p = ls.split()
                if len(p) >= 2 and p[0][:1].isdigit():
                    peaks.append([float(p[0]), float(p[1])]); in_peaks -= 1
                    continue
                in_peaks = 0
            if not ls:
                continue
            low = ls.lower()
            if low.startswith("name:"):
                if cur:
                    flush(); n_entry += 1
                    if n_entry % 200000 == 0:
                        print(f"    ...scanned {n_entry} entries, kept {sum(len(v) for v in idx.values())}")
                cur = {}; peaks = []
            if low.startswith("inchikey:"):
                cur["inchikey"] = ls.split(":", 1)[1].strip()
            elif low.startswith("precursor_type:"):
                cur["precursor_type"] = ls.split(":", 1)[1].strip()
            elif low.startswith("num peaks:"):
                in_peaks = int(ls.split(":", 1)[1])
                peaks = []
        if cur:
            flush()
    idx["__complete__"] = True
    NIST_CACHE.write_text(json.dumps(idx))
    print(f"[P2-nist] kept {sum(len(v) for v in idx.values() if isinstance(v, list))} spectra "
          f"for {len([k for k in idx if k!='__complete__'])} IK14")
    return idx


# ─────────────────────────── phase 3: compute + report ──────────────────────
def _resolve_libpeaks(m, ik14, qpeaks, web, nist):
    """Return (lib_peaks, source). For nist, pick the CE variant maximizing
    forward_cosine to the query (best library spectrum of correct cpd+adduct)."""
    src = str(m.get("library_wiki_id")).split("/")[0]
    acc = m.get("mona_id")
    if src in WEB_SOURCES:
        return web.get(f"{src}:{acc}"), src
    if src == "nist23":
        variants = nist.get(ik14) or []
        if not variants:
            return None, src
        want = _norm_adduct(m.get("adduct"))
        cands = [pk for pt, pk in variants if _norm_adduct(pt) == want] or [pk for _, pk in variants]
        best, bestf = None, -1.0
        for pk in cands:
            f, _, _, _ = compute_ms2_scores(qpeaks, pk)
            if f == f and f > bestf:
                best, bestf = pk, f
        return best, src
    return None, src           # msdial / agilent: proprietary, skipped


def phase3_report(df, getdata, web, nist):
    rows = []
    n_nopeak = 0; drop_src = {}
    for _, r in df.iterrows():
        gd = getdata.get(r["wiki_id"])
        if not gd or not gd.get("matched") or not gd.get("qpeaks"):
            continue
        m = gd["matched"]
        lpk, src = _resolve_libpeaks(m, r["ik14"], gd["qpeaks"], web, nist)
        if not lpk:
            n_nopeak += 1; drop_src[src] = drop_src.get(src, 0) + 1
            continue
        fwd, rev, covc, covi = compute_ms2_scores(gd["qpeaks"], lpk)
        rows.append({"wiki_id": r["wiki_id"], "polarity": r["polarity"], "label": r["label"],
                     "ik14": r["ik14"], "name": r["name"], "adduct": r["adduct"], "rt": r["rt"],
                     "src": src, "entropy_sim": m.get("forward_esim"),
                     "forward_cos": fwd, "reverse_cos": rev, "cov_int": covi,
                     "adduct_match": _norm_adduct(m.get("adduct")) == _norm_adduct(r["adduct"])})
    res = pd.DataFrame(rows)
    res.to_csv(OUT_TABLE, index=False)
    print(f"\nComputed scores for {len(res)} bins ({n_nopeak} dropped, by src={drop_src}) -> {OUT_TABLE}")
    if res.empty:
        return res
    print("  coverage by source:", dict(res["src"].value_counts()))
    print(f"  adduct-matched: {int(res.adduct_match.sum())}/{len(res)}")

    def report(sub, title, fwd_col, fwd_name):
        sub = sub.copy()
        sub["in_zone"] = (sub[fwd_col] < FWD_CUT) & (sub["reverse_cos"] > REV_CUT)
        print(f"\n{'='*70}\n{title}:  {fwd_name} < {FWD_CUT}  AND  reverse_cos > {REV_CUT}\n{'='*70}")
        for lab in ["TP", "FP", "unlabelled"]:
            s = sub[sub["label"] == lab]
            if len(s):
                print(f"  {lab:11s}: n={len(s):5d}  {fwd_name}<{FWD_CUT}={int((s[fwd_col]<FWD_CUT).sum()):5d}  "
                      f"RESCUED={int(s['in_zone'].sum()):5d}")
        z = sub[sub["in_zone"]]
        n_tp = int((z.label == "TP").sum()); n_fp = int((z.label == "FP").sum())
        n_un = int((z.label == "unlabelled").sum()); d = n_tp + n_fp
        print(f"  --> ZONE TOTAL {len(z)}  (TP={n_tp}, FP={n_fp}, unlab={n_un})"
              + (f"   precision TP/(TP+FP)={n_tp/d:.3f}" if d else ""))

    am = res[res.adduct_match]
    print(f"\n########## HEADLINE: adduct-matched bins (n={len(am)}) ##########")
    report(am, "RESCUE ZONE [adduct-matched]", "entropy_sim", "entropy_sim")
    report(am, "RESCUE ZONE [adduct-matched]", "forward_cos", "forward_cos")
    print(f"\n---------- all matched bins (incl. cross-adduct, n={len(res)}) ----------")
    report(res, "RESCUE ZONE [all]", "entropy_sim", "entropy_sim")
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None, help="cap annotated bins (smoke test)")
    ap.add_argument("--polarity", choices=["neg", "pos", "both"], default="both")
    args = ap.parse_args()

    token = Path(TOKEN_FILE).read_text().strip() if os.path.exists(TOKEN_FILE) else ""
    if not token:
        raise SystemExit(f"No token at {TOKEN_FILE}")

    pols = ["neg", "pos"] if args.polarity == "both" else [args.polarity]
    frames = []
    for pol in pols:
        d = pd.read_csv(DATA / f"lipidomics_{pol}.csv", low_memory=False)
        d["polarity"] = pol
        frames.append(build_labels(d))
    df = pd.concat(frames, ignore_index=True)

    # only annotated, non-ignored bins matter for the test (need a correct compound)
    work = df[df["name"].notna() & (df["label"] != "ignore")].copy()
    print(f"Annotated working set: {len(work)}  "
          f"(TP={ (work.label=='TP').sum() }, FP={ (work.label=='FP').sum() }, "
          f"unlabelled={ (work.label=='unlabelled').sum() })")
    if args.limit:
        work = work.head(args.limit)

    getdata = phase1_fetch(work[["wiki_id", "name", "adduct", "ik14", "identity_score"]], token)

    # split needs by source (from library_wiki_id prefix)
    web_needs = set(); nist_ik14 = set()
    ik_by_wid = dict(zip(work["wiki_id"], work["ik14"]))
    for wid, v in getdata.items():
        m = v.get("matched")
        if not m:
            continue
        src = str(m.get("library_wiki_id")).split("/")[0]
        if src in WEB_SOURCES:
            web_needs.add((src, m.get("mona_id")))
        elif src == "nist23":
            ik = ik_by_wid.get(wid)
            if ik and ik != "nan":
                nist_ik14.add(ik)

    web = phase2_web(web_needs)
    nist = phase2_nist(nist_ik14) if nist_ik14 else {}
    phase3_report(work, getdata, web, nist)


if __name__ == "__main__":
    main()
