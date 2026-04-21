"""
refetch_orbitrap_hits.py

Re-fetches hits for annotated Orbitrap spectra (TP + FP) from MassWiki API.
Extracts:
  - Full hits table with library_wiki_id (for reverse similarity)
  - Query peaks per wiki_id (for forward/reverse score computation)

Usage:
    python code/refetch_orbitrap_hits.py --token <BEARER_TOKEN>

Outputs:
    data/orbitrap_hits_refetched.csv  -- full hits with library_wiki_id
    data/query_peaks_cache.json       -- {wiki_id: [[mz, int], ...]}
"""

from __future__ import annotations
import os, json, time, argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

BASE_URL  = "https://masswiki.us-west-2.elasticbeanstalk.com"
ENDPOINT  = f"{BASE_URL}/analysis/get_data"
TIMEOUT   = 20
MAX_WORKERS = 6
RPS       = 6.0

SPECTRA_PATH   = "data/masswiki_Orbitrap HILIC negESI_2026-03-19.xlsx"
SOLID_TP_PATH  = "data/solid_tp.csv"
OUT_HITS       = "data/orbitrap_hits_refetched.csv"
OUT_PEAKS      = "data/query_peaks_cache.json"
OUT_ERRORS     = "data/orbitrap_refetch_errors.csv"


# ── Helpers ──────────────────────────────────────────────────────────────────

def build_session(token: str) -> requests.Session:
    s = requests.Session()
    retries = Retry(total=4, backoff_factor=0.5,
                    status_forcelist=(429, 500, 502, 503, 504),
                    allowed_methods=["GET"])
    adapter = HTTPAdapter(max_retries=retries, pool_connections=64, pool_maxsize=64)
    s.mount("https://", adapter)
    s.headers.update({
        "Accept": "application/json",
        "Authorization": "Bearer " + token,
        "User-Agent": "Mozilla/5.0",
    })
    return s


def fetch_one(session: requests.Session, wiki_id: str,
              source: str, is_public: bool) -> Tuple[Optional[dict], int, str]:
    params = {"wiki_id": wiki_id, "source": source,
              "isPublic": str(is_public).lower()}
    r = session.get(ENDPOINT, params=params, timeout=TIMEOUT)
    if r.status_code != 200:
        return None, r.status_code, r.text[:200]
    payload = r.json()
    if not isinstance(payload, dict):
        return None, r.status_code, "non-dict response"
    return payload, r.status_code, ""


def extract_query_peaks(payload: dict) -> Optional[list]:
    """Extract raw query peaks from API response."""
    # Primary location: payload["spectrum"]["peaks"] (top-level spectrum dict)
    top_spec = payload.get("spectrum")
    if isinstance(top_spec, dict):
        for pk in ["peaks", "peaks_clean", "msms", "fragments"]:
            val = top_spec.get(pk)
            if val and isinstance(val, list):
                return val

    # Secondary: payload["analysis"]["spectrum"]["peaks"]
    analysis = payload.get("analysis", {})
    for key in ["spectrum", "query_spectrum", "msms"]:
        spec = analysis.get(key)
        if isinstance(spec, dict):
            peaks = spec.get("peaks") or spec.get("msms") or spec.get("fragments")
            if peaks:
                return peaks
        elif isinstance(spec, list) and spec:
            return spec

    # Fallback: top-level list fields
    for key in ["peaks", "msms"]:
        val = payload.get(key)
        if val and isinstance(val, list):
            return val

    return None


def extract_hits(payload: dict, wiki_id: str) -> List[dict]:
    """Extract all hits (reference + annotation) from API response."""
    rows = []
    analysis = payload.get("analysis", {})

    for hit_source in ("reference_library", "annotation_library"):
        lib = analysis.get(hit_source) or {}
        hits = lib.get("identity_search") or []
        if not isinstance(hits, list):
            continue

        source_label = "reference" if hit_source == "reference_library" else "annotation"

        for i, h in enumerate(hits, 1):
            rows.append({
                "wiki_id":              wiki_id,
                "hit_source":           source_label,
                "db":                   h.get("db") or h.get("source"),
                "id":                   h.get("id") or h.get("identifier"),
                "lib_name":             h.get("name"),
                "adduct":               h.get("adduct"),
                "lib_precursor_mz":     h.get("precursor_mz") or h.get("precursor"),
                "entropy_similarity":   h.get("entropy_similarity") or h.get("score"),
                "library_type":         h.get("library_type"),
                "lib_rt":               h.get("rt") or h.get("retention_time"),
                "ri":                   h.get("ri"),
                "rank":                 h.get("rank") or i,
                "smiles":               h.get("smiles"),
                "predicted_rt_hilic":   h.get("predicted_rt_hilic"),
                "predicted_rt_rp":      h.get("predicted_rt_rp"),
                "anno_delta_rt":        h.get("anno_delta_rt"),
                "delta_predicted_rt":   h.get("delta_predicted_rt"),
                # KEY NEW FIELD: library_wiki_id for peak fetching
                "library_wiki_id":      h.get("library_wiki_id") or h.get("wiki_id"),
            })
    return rows


# ── Main fetch loop ───────────────────────────────────────────────────────────

def fetch_all(wiki_ids: List[str], token: str) -> Tuple[List[dict], dict, List[dict]]:
    """
    Returns: (hit_rows, query_peaks_dict, error_rows)
    """
    session = build_session(token)
    min_interval = 1.0 / RPS
    last_ts = 0.0

    hit_rows: List[dict] = []
    peaks_dict: dict = {}
    error_rows: List[dict] = []

    def worker(wiki_id: str):
        nonlocal last_ts
        now = time.time()
        wait = min_interval - (now - last_ts)
        if wait > 0:
            time.sleep(wait)
        last_ts = time.time()

        # Try primary (binbase, private)
        payload, status, err = fetch_one(session, wiki_id, "binbase", False)
        if payload is None and status in (400, 404):
            # Fallback to zyang2k public
            payload, status, err = fetch_one(session, wiki_id, "zyang2k", True)

        return wiki_id, payload, status, err

    total = len(wiki_ids)
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as ex:
        futs = [ex.submit(worker, wid) for wid in wiki_ids]
        for done_i, fut in enumerate(as_completed(futs), 1):
            wiki_id, payload, status, err = fut.result()

            if payload is None:
                error_rows.append({"wiki_id": wiki_id, "status": status, "error": err})
                print(f"[{done_i}/{total}] ERROR {wiki_id} status={status}")
                continue

            # Extract hits
            hits = extract_hits(payload, wiki_id)
            hit_rows.extend(hits)

            # Extract query peaks
            peaks = extract_query_peaks(payload)
            if peaks:
                peaks_dict[wiki_id] = peaks

            print(f"[{done_i}/{total}] OK {wiki_id}  hits={len(hits)}  peaks={'yes' if peaks else 'no'}")

    return hit_rows, peaks_dict, error_rows


# ── Entry point ───────────────────────────────────────────────────────────────

TOKEN_FILE = "/tmp/.mw_token"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--token-file", default=TOKEN_FILE,
                        help="File containing Bearer token (default: /tmp/.mw_token)")
    parser.add_argument("--limit", type=int, default=None, help="Limit wiki_ids (for testing)")
    parser.add_argument("--dry-run", action="store_true", help="Print first response and exit")
    args = parser.parse_args()

    try:
        with open(args.token_file) as fh:
            args.token = "".join(fh.read().split())
    except FileNotFoundError:
        raise SystemExit(f"Token file not found: {args.token_file}\n"
                         "Write your Bearer token to that file first.")

    # Load annotated spectra (TP + FP only)
    print("Loading spectra...")
    spectra = pd.read_excel(SPECTRA_PATH, header=4)
    spectra["label"] = "unlabeled"
    spectra.loc[:1297, "label"] = "TP"
    spectra.loc[spectra["name"].str.startswith("yy_", na=False), "label"] = "FP"
    ann = spectra[spectra["label"].isin(["TP", "FP"])]["wiki_id"].tolist()

    if args.limit:
        ann = ann[:args.limit]

    print(f"Fetching {len(ann):,} spectra (TP + FP)...")

    if args.dry_run:
        # Fetch one and print raw response to inspect structure
        session = build_session(args.token)
        payload, status, err = fetch_one(session, ann[0], "binbase", False)
        if payload:
            print(f"\nStatus: {status}")
            print("Top-level keys:", list(payload.keys()))
            analysis = payload.get("analysis", {})
            print("Analysis keys:", list(analysis.keys()))
            for k in analysis:
                v = analysis[k]
                if isinstance(v, dict):
                    print(f"  analysis[{k}] keys: {list(v.keys())[:5]}")
                elif isinstance(v, list):
                    print(f"  analysis[{k}]: list of {len(v)}")
            # Show peaks location
            peaks = extract_query_peaks(payload)
            print(f"\nQuery peaks found: {'yes, ' + str(len(peaks)) + ' peaks' if peaks else 'NO'}")
            if peaks:
                print("Sample peaks:", peaks[:3])
            # Show one hit
            hits = extract_hits(payload, ann[0])
            if hits:
                print(f"\nFirst hit:")
                for k, v in hits[0].items():
                    print(f"  {k}: {v}")
        else:
            print(f"Failed: status={status}, err={err}")
        return

    # Full fetch
    hit_rows, peaks_dict, error_rows = fetch_all(ann, args.token)

    # Save hits
    hits_df = pd.DataFrame(hit_rows)
    hits_df.to_csv(OUT_HITS, index=False)
    print(f"\n✓ Saved {len(hits_df):,} hits → {OUT_HITS}")
    print(f"  library_wiki_id coverage: {hits_df['library_wiki_id'].notna().sum():,}/{len(hits_df):,}")

    # Save query peaks
    with open(OUT_PEAKS, "w") as f:
        json.dump(peaks_dict, f)
    print(f"✓ Saved {len(peaks_dict):,} query peak sets → {OUT_PEAKS}")

    # Save errors
    if error_rows:
        pd.DataFrame(error_rows).to_csv(OUT_ERRORS, index=False)
        print(f"⚠ {len(error_rows)} errors → {OUT_ERRORS}")

    print(f"\nSummary:")
    print(f"  Spectra fetched: {len(ann):,}")
    print(f"  Hits extracted:  {len(hits_df):,}")
    print(f"  Query peaks:     {len(peaks_dict):,}/{len(ann):,}")
    print(f"  Errors:          {len(error_rows)}")
    print(f"\nDB breakdown:")
    if len(hits_df):
        print(hits_df["db"].value_counts().to_string())


if __name__ == "__main__":
    main()
