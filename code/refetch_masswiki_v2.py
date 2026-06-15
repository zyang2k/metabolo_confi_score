"""
refetch_masswiki_v2.py — Stage 0 of the Bayesian two-tier pipeline.

Fetches per-spectrum hits + query peaks for the Orbitrap HILIC neg + pos
curated CSVs, filtered to annotated rows only (`name` not blank).

Fixes vs the old refetch_orbitrap_hits.py:
  - Reads BOTH neg + pos CSVs (not just neg)
  - Filters by `name` not blank (gets all 5,822 annotated rows, not just [:1297]+yy_)
  - source=binbase only (zyang2k returns 400 — confirmed in dry-run)
  - Captures BOTH identity_search and neutral_loss_search from reference_library
  - Tags each hit row with hit_source (ref_identity / ref_neutral_loss)
  - Saves query peaks from payload.spectrum.peaks (the right location)
  - Incremental save every 200 spectra so a crash mid-run doesn't lose everything

Outputs:
  data/orbitrap_hits_v2.csv          -- all hits, tagged
  data/query_peaks_cache_v2.json     -- query peaks per wiki_id
  data/orbitrap_refetch_errors_v2.csv -- errors
"""
from __future__ import annotations
import os, json, time, argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

BASE     = "https://masswiki.us-west-2.elasticbeanstalk.com"
ENDPOINT = f"{BASE}/analysis/get_data"
TIMEOUT  = 30
WORKERS  = 6
RPS      = 6.0

NEG_CSV  = "data/Orbitrap_HILIC_negESI_curated_042126.csv"
POS_CSV  = "data/Orbitrap_HILIC_posESI_curated_042126.csv"
OUT_HITS = "data/orbitrap_hits_v2.csv"
OUT_PEAKS= "data/query_peaks_cache_v2.json"
OUT_ERR  = "data/orbitrap_refetch_errors_v2.csv"

TOKEN_FILE = "/tmp/.mw_token"
SAVE_EVERY = 200   # spectra between incremental saves


# ── Hit fields we keep ──────────────────────────────────────────────────────
HIT_FIELDS = [
    "wiki_id", "polarity", "hit_source", "rank",
    "db", "id", "name", "adduct", "smiles",
    "precursor_mz", "lib_rt",
    "predicted_rt_hilic", "predicted_rt_rp",
    "delta_predicted_rt", "anno_delta_rt",
    "entropy_similarity", "library_wiki_id", "index",
]

# open_search reproduces the Apr-23 pull's 'reference' source (≈100 hits/spectrum, low
# mean esim ~0.17) that build_features_v2.py keeps alongside ref_identity.
#
# ⚠️ KNOWN-BROKEN 2026-06-11 (do NOT use this script to refresh the production table yet):
# MassWiki's reference_library.identity_search now returns EMPTY for binbase spectra
# (verified across isPublic true/false + source variants). identity_search was the SOLE
# source of the curator-annotated correct candidate — 4,404/4,404 hit_label=1 rows in the
# Apr-23 data came from ref_identity. With it empty, correct-IK14 coverage on labeled
# spectra collapses 97.2% → ~1-3% and the labels are destroyed. open_search does NOT
# carry the correct compounds. RT predictions ARE now populated (61% → 88.7%), but cannot
# be obtained without losing identity_search. Blocked pending a MassWiki-side fix
# (Quentin/Fanzhou): why is identity_search empty, and how to get identity + RT together.
SEARCH_TYPES = [
    ("ref_identity",     "reference_library", "identity_search"),
    ("reference",        "reference_library", "open_search"),
    ("ref_neutral_loss", "reference_library", "neutral_loss_search"),
]


def build_session(token: str) -> requests.Session:
    s = requests.Session()
    s.headers.update({
        "Accept": "application/json",
        "Authorization": f"Bearer {token}",
        "User-Agent": "Mozilla/5.0",
    })
    retries = Retry(total=4, backoff_factor=0.5,
                    status_forcelist=(429, 500, 502, 503, 504),
                    allowed_methods=["GET"])
    adapter = HTTPAdapter(max_retries=retries, pool_connections=64, pool_maxsize=64)
    s.mount("https://", adapter)
    return s


def fetch_one(session: requests.Session, wiki_id: str):
    params = {"wiki_id": wiki_id, "source": "binbase", "isPublic": "false"}
    try:
        r = session.get(ENDPOINT, params=params, timeout=TIMEOUT)
    except Exception as e:
        return None, -1, f"req-exc: {e}"
    if r.status_code != 200:
        return None, r.status_code, r.text[:200]
    try:
        return r.json(), 200, ""
    except Exception as e:
        return None, 200, f"json-exc: {e}"


def extract_hits(payload: dict, wiki_id: str, polarity: str) -> list[dict]:
    rows = []
    analysis = payload.get("analysis", {}) or {}
    for tag, lib_key, sub_key in SEARCH_TYPES:
        lib = analysis.get(lib_key) or {}
        hits = lib.get(sub_key) or []
        if not isinstance(hits, list):
            continue
        for i, h in enumerate(hits, 1):
            if not isinstance(h, dict):
                continue
            rows.append({
                "wiki_id":            wiki_id,
                "polarity":           polarity,
                "hit_source":         tag,
                "rank":               h.get("rank") or i,
                "db":                 h.get("db") or h.get("source"),
                "id":                 h.get("id"),
                "name":               h.get("name"),
                "adduct":             h.get("adduct"),
                "smiles":             h.get("smiles"),
                "precursor_mz":       h.get("precursor_mz"),
                "lib_rt":             h.get("rt"),
                "predicted_rt_hilic": h.get("predicted_rt_hilic"),
                "predicted_rt_rp":    h.get("predicted_rt_rp"),
                "delta_predicted_rt": h.get("delta_predicted_rt"),
                "anno_delta_rt":      h.get("anno_delta_rt"),
                "entropy_similarity": h.get("entropy_similarity"),
                "library_wiki_id":    h.get("library_wiki_id"),
                "index":              h.get("index"),
            })
    return rows


def extract_query_peaks(payload: dict):
    spec = payload.get("spectrum") or {}
    if isinstance(spec, dict):
        peaks = spec.get("peaks")
        if isinstance(peaks, list) and peaks:
            return peaks
    return None


def load_targets() -> pd.DataFrame:
    """Read both curated CSVs, filter to annotated rows, return wiki_id+polarity."""
    rows = []
    for path, polarity in [(NEG_CSV, "neg"), (POS_CSV, "pos")]:
        df = pd.read_csv(path, low_memory=False, usecols=["wiki_id", "name"])
        df = df[df["name"].notna() & (df["name"].astype(str).str.strip() != "")]
        df["polarity"] = polarity
        rows.append(df[["wiki_id", "polarity"]])
    return pd.concat(rows, ignore_index=True).drop_duplicates("wiki_id").reset_index(drop=True)


def save_outputs(hit_rows, peaks_dict, errors):
    if hit_rows:
        pd.DataFrame(hit_rows, columns=HIT_FIELDS).to_csv(OUT_HITS, index=False)
    with open(OUT_PEAKS, "w") as f:
        json.dump(peaks_dict, f)
    if errors:
        pd.DataFrame(errors).to_csv(OUT_ERR, index=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--token-file", default=TOKEN_FILE)
    ap.add_argument("--limit", type=int, default=None,
                    help="limit number of spectra (testing)")
    args = ap.parse_args()

    if not os.path.exists(args.token_file):
        raise SystemExit(f"Token file not found: {args.token_file}")
    token = Path(args.token_file).read_text().strip()
    if not token:
        raise SystemExit(f"Token file {args.token_file} is empty")

    targets = load_targets()
    if args.limit:
        targets = targets.head(args.limit)
    print(f"Loaded {len(targets):,} annotated wiki_ids "
          f"({(targets['polarity']=='neg').sum()} neg + "
          f"{(targets['polarity']=='pos').sum()} pos)")

    session   = build_session(token)
    hit_rows  = []
    peaks     = {}
    errors    = []

    min_interval = 1.0 / RPS
    last_ts = [0.0]

    def worker(row):
        wid, pol = row["wiki_id"], row["polarity"]
        # Throttle (single shared timestamp; coarse but safe)
        wait = min_interval - (time.time() - last_ts[0])
        if wait > 0:
            time.sleep(wait)
        last_ts[0] = time.time()
        payload, status, err = fetch_one(session, wid)
        return wid, pol, payload, status, err

    rows = targets.to_dict("records")
    total = len(rows)
    t0 = time.time()
    saved_at = 0
    ok = 0

    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        futs = [ex.submit(worker, r) for r in rows]
        for done_i, fut in enumerate(as_completed(futs), 1):
            wid, pol, payload, status, err = fut.result()
            if payload is None:
                errors.append({"wiki_id": wid, "status": status, "error": err})
            else:
                ok += 1
                hit_rows.extend(extract_hits(payload, wid, pol))
                qp = extract_query_peaks(payload)
                if qp:
                    peaks[wid] = qp

            if done_i % 50 == 0 or done_i == total:
                elapsed = time.time() - t0
                rate = done_i / max(elapsed, 1e-6)
                eta = (total - done_i) / max(rate, 1e-6)
                print(f"[{done_i:>5}/{total}] ok={ok} hits={len(hit_rows):,} "
                      f"peaks={len(peaks):,} errors={len(errors)} "
                      f"rate={rate:.1f}/s eta={eta:.0f}s")

            if done_i - saved_at >= SAVE_EVERY:
                save_outputs(hit_rows, peaks, errors)
                saved_at = done_i

    save_outputs(hit_rows, peaks, errors)
    print(f"\n✓ done. ok={ok}/{total} hits={len(hit_rows):,} peaks={len(peaks):,} errors={len(errors)}")
    print(f"  wrote {OUT_HITS}, {OUT_PEAKS}, {OUT_ERR if errors else '(no errors)'}")


if __name__ == "__main__":
    main()
