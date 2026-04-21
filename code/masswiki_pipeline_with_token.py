# masswiki_pipeline_with_token.py
# CSV -> filter wiki_ids -> fetch reference-library hits using Bearer token
# with automatic fallback between (binbase,false) and (zyang2k,true).
#
# Both reference_library (public: NIST23, GNPS, MassBank…) and
# annotation_library (in-house: 5min_hilic_neg, 5min_lipid_neg…) hits are
# fetched from the same get_data response and written to the output CSV with
# a hit_source column ("reference" | "annotation") so they can be used
# separately or merged downstream.

from __future__ import annotations

import os, time, argparse
from typing import Any, Dict, Iterable, List, Optional, Tuple
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

BASE_URL = os.getenv("MASSWIKI_BASE_URL", "https://masswiki.us-west-2.elasticbeanstalk.com")
ENDPOINT = f"{BASE_URL}/analysis/get_data"

# --- Defaults ---
TIMEOUT = float(os.getenv("MASSWIKI_TIMEOUT", "20"))
MAX_WORKERS = int(os.getenv("MASSWIKI_MAX_WORKERS", "8"))
RPS = float(os.getenv("MASSWIKI_RPS", "6"))  # polite throttle

def filter_masswiki_results(df: pd.DataFrame) -> pd.DataFrame:
    if "is_manual_annotated" not in df.columns or "name" not in df.columns:
        raise ValueError("CSV must include 'is_manual_annotated' and 'name'.")
    out = df.copy()
    out["name"] = out["name"].astype(str).fillna("").str.strip()
    out["is_manual_annotated"] = out["is_manual_annotated"].fillna(False).astype(bool)
    mask = out["is_manual_annotated"] & (~out["name"].str.lower().str.startswith(("yy", "zz")))
    return out.loc[mask].reset_index(drop=True)

def load_and_filter_wiki_ids(csv_path: str) -> List[str]:
    abs_path = os.path.abspath(csv_path)
    print(f"Reading CSV: {abs_path}")
    df = pd.read_csv(abs_path)
    df_filt = filter_masswiki_results(df)
    if "wiki_id" not in df_filt.columns:
        raise ValueError("CSV must include column 'wiki_id'.")
    return df_filt["wiki_id"].astype(str).tolist()

def _build_session(token: str) -> requests.Session:
    s = requests.Session()
    retries = Retry(
        total=4,
        backoff_factor=0.4,
        status_forcelist=(429, 500, 502, 503, 504),  # not 401/400
        allowed_methods=["GET"],
        respect_retry_after_header=True,
    )
    adapter = HTTPAdapter(max_retries=retries, pool_connections=64, pool_maxsize=64)
    s.mount("https://", adapter)
    s.mount("http://", adapter)
    s.headers.update({
        "Accept": "application/json",
        "Authorization": f"Bearer {token}",
        # harmless, sometimes required for public reads:
        "Origin": BASE_URL,
        "Referer": f"{BASE_URL}/docs",
        "User-Agent": "Mozilla/5.0",
        "Accept-Language": "en-US,en;q=0.9,zh-CN;q=0.8,zh;q=0.7",
        "Sec-Fetch-Mode": "cors",
        "Sec-Fetch-Site": "same-origin",
        "Sec-Fetch-Dest": "empty",
    })
    return s

def _fetch_ref_once(session: requests.Session, wiki_id: str, source: str, is_public: bool) -> Tuple[Optional[Dict[str, List[Dict[str, Any]]]], int, str]:
    """Return ({"reference": [...], "annotation": [...]}_or_None, status_code, error_snippet)."""
    params = {"wiki_id": wiki_id, "source": source, "isPublic": str(is_public).lower()}
    r = session.get(ENDPOINT, params=params, timeout=TIMEOUT)
    if r.status_code != 200:
        return None, r.status_code, (r.text[:300] if isinstance(r.text, str) else str(r.content)[:300])
    if not r.headers.get("content-type", "").lower().startswith("application/json"):
        return None, r.status_code, "Non-JSON response"
    payload = r.json()
    analysis = payload.get("analysis", {}) if isinstance(payload, dict) else {}

    def _extract(key: str) -> Optional[List[Dict[str, Any]]]:
        hits = (analysis.get(key) or {}).get("identity_search")
        return hits if isinstance(hits, list) else None

    ref_hits  = _extract("reference_library")
    anno_hits = _extract("annotation_library")

    # Return None only if both sources are empty (treat as miss)
    if ref_hits is None and anno_hits is None:
        return None, r.status_code, ""
    return {"reference": ref_hits or [], "annotation": anno_hits or []}, r.status_code, ""

def fetch_with_fallback(session: requests.Session, wiki_id: str,
                        primary: Tuple[str, bool], secondary: Tuple[str, bool],
                        try_secondary_on_400: bool = True) -> Tuple[str, Optional[Dict[str, List[Dict[str, Any]]]], Optional[str], int, str, str]:
    """
    Returns: (wiki_id, hits, error, status, used_source, used_isPublic)
    Tries primary; if 400 (likely not found in that namespace) and flag set, tries secondary.
    """
    p_source, p_public = primary
    hits, status, err = _fetch_ref_once(session, wiki_id, p_source, p_public)
    used_source, used_pub = p_source, str(p_public).lower()
    if status == 200 and hits is not None:
        return wiki_id, hits, None, status, used_source, used_pub

    # Only fallback on 400 (server-side NoneType usually), optionally on 404 if they implement it later
    if try_secondary_on_400 and status in (400, 404):
        s_source, s_public = secondary
        hits2, status2, err2 = _fetch_ref_once(session, wiki_id, s_source, s_public)
        used_source, used_pub = s_source, str(s_public).lower()
        if status2 == 200 and hits2 is not None:
            return wiki_id, hits2, None, status2, used_source, used_pub
        # report the secondary attempt error if it failed too
        return wiki_id, None, f"Primary {p_source}/{p_public} -> {status}: {err}; Secondary {s_source}/{s_public} -> {status2}: {err2}", status2, used_source, used_pub

    # no fallback or not a 400/404
    return wiki_id, hits, (None if status == 200 else f"{status}: {err}"), status, used_source, used_pub

def fetch_reference_hits_many(wiki_ids: Iterable[str], token: str,
                              primary: Tuple[str, bool],
                              secondary: Tuple[str, bool],
                              max_workers: int = MAX_WORKERS, rps: float = RPS) -> Tuple[Dict[str, Optional[Dict[str, List[Dict[str, Any]]]]], pd.DataFrame]:
    """
    Returns (hits_dict, errors_df)
    """
    ids = [wid for wid in (str(w) for w in wiki_ids) if wid.strip()]
    session = _build_session(token)

    min_interval = 1.0 / rps if rps > 0 else 0.0
    last_ts = 0.0

    def worker(wid: str):
        nonlocal last_ts
        now = time.time()
        wait = min_interval - (now - last_ts)
        if wait > 0: time.sleep(wait)
        last_ts = time.time()
        return fetch_with_fallback(session, wid, primary, secondary, try_secondary_on_400=True)

    out: Dict[str, Optional[List[Dict[str, Any]]]] = {}
    err_rows: List[Dict[str, Any]] = []
    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        futs = [ex.submit(worker, wid) for wid in ids]
        total = len(futs); done = 0
        for fut in as_completed(futs):
            wid, hits, err, status, used_source, used_public = fut.result()
            out[wid] = hits
            done += 1
            print(f"[MassWiki] {done}/{total} {wid} status={status} source={used_source} isPublic={used_public} err={err or '-'}")
            if err:
                err_rows.append({
                    "wiki_id": wid,
                    "status": status,
                    "error": err,
                    "used_source": used_source,
                    "used_isPublic": used_public
                })
    errors_df = pd.DataFrame(err_rows)
    return out, errors_df

def flatten_reference_hits(hits_dict: Dict[str, Optional[Dict[str, List[Dict[str, Any]]]]]) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []
    for wid, hit_bundle in hits_dict.items():
        if not hit_bundle:
            continue
        # hit_bundle = {"reference": [...], "annotation": [...]}
        for source_key in ("reference", "annotation"):
            hits = hit_bundle.get(source_key) or []
            for i, h in enumerate(hits, 1):
                rows.append({
                    "wiki_id": wid,
                    "hit_source": source_key,          # "reference" or "annotation"
                    "db": h.get("db") or h.get("source"),
                    "id": h.get("id") or h.get("identifier") or h.get("accession"),
                    "lib_name": h.get("name"),
                    "adduct": h.get("adduct"),
                    "lib_precursor_mz": h.get("precursor_mz") or h.get("precursor"),
                    "entropy_similarity": h.get("entropy_similarity") or h.get("score") or h.get("similarity"),
                    "library_type": h.get("library_type"),
                    "lib_rt": h.get("rt") or h.get("retention_time"),
                    "ri": h.get("ri") or h.get("retention_index"),
                    "rank": h.get("rank") or i,
                    # structural identity — used for InChIKey-based ground truth labeling
                    "smiles": h.get("smiles"),
                    # RT predictions — activates RT likelihood channel
                    "predicted_rt_hilic": h.get("predicted_rt_hilic"),
                    "predicted_rt_rp": h.get("predicted_rt_rp"),
                    "anno_delta_rt": h.get("anno_delta_rt"),
                    "delta_predicted_rt": h.get("delta_predicted_rt"),
                })
    return pd.DataFrame(rows)

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="CSV -> filter wiki_ids -> fetch reference-library hits with Bearer token (auto fallback).")
    p.add_argument("--csv", required=True, help="Path to input CSV.")
    p.add_argument("--out", required=True, help="Path to save flattened hits CSV.")
    p.add_argument("--errors", default=None, help="Optional path to save an errors CSV (diagnostics).")
    p.add_argument("--limit", type=int, default=None, help="Optional limit on IDs (quick test).")
    p.add_argument("--token", required=True, help="Bearer access token.")
    # Primary + secondary combos (defaults: try private binbase first, then public)
    p.add_argument("--primary-source", default="binbase", help="Primary source (default: binbase).")
    p.add_argument("--primary-public", default="false", choices=["true","false"], help="Primary isPublic (default: false).")
    p.add_argument("--secondary-source", default="zyang2k", help="Secondary source (default: zyang2k).")
    p.add_argument("--secondary-public", default="true", choices=["true","false"], help="Secondary isPublic (default: true).")
    p.add_argument("--workers", type=int, default=MAX_WORKERS, help=f"Max threads (default {MAX_WORKERS}).")
    p.add_argument("--rps", type=float, default=RPS, help=f"Requests/sec throttle (default {RPS}).")
    return p.parse_args()

def main():
    args = parse_args()
    wiki_ids = load_and_filter_wiki_ids(args.csv)
    if args.limit: wiki_ids = wiki_ids[:args.limit]

    primary = (args.primary_source, args.primary_public.lower() == "true")
    secondary = (args.secondary_source, args.secondary_public.lower() == "true")

    hits_dict, errors_df = fetch_reference_hits_many(
        wiki_ids, token=args.token, primary=primary, secondary=secondary,
        max_workers=args.workers, rps=args.rps
    )

    df = flatten_reference_hits(hits_dict)
    out_abs = os.path.abspath(args.out)
    df.to_csv(out_abs, index=False)
    print(f"\nSaved hits: {out_abs}  (rows={len(df)})")

    if args.errors:
        err_abs = os.path.abspath(args.errors)
        errors_df.to_csv(err_abs, index=False)
        print(f"Saved diagnostics: {err_abs}  (rows={len(errors_df)})")

if __name__ == "__main__":
    main()
