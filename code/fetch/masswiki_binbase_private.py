import os, time, argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Any, Dict, List, Optional, Tuple
import pandas as pd, requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

BASE_URL = os.getenv("MASSWIKI_BASE_URL", "https://masswiki.us-west-2.elasticbeanstalk.com")
ENDPOINT = f"{BASE_URL}/analysis/get_data"
TIMEOUT = 20
MAX_WORKERS = 8
RPS = 6.0

def filter_masswiki_results(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["name"] = df["name"].astype(str).fillna("").str.strip()
    df["is_manual_annotated"] = df["is_manual_annotated"].fillna(False).astype(bool)
    mask = df["is_manual_annotated"] & (~df["name"].str.lower().str.startswith(("yy","zz")))
    return df.loc[mask].reset_index(drop=True)

def load_ids(csv_path: str) -> List[str]:
    df = pd.read_csv(csv_path)
    df = filter_masswiki_results(df)
    return df["wiki_id"].astype(str).tolist()

def build_session(token: str) -> requests.Session:
    s = requests.Session()
    s.headers.update({
        "Accept": "application/json",
        "Authorization": f"Bearer {token}",
        "Origin": BASE_URL,
        "Referer": f"{BASE_URL}/docs",
        "User-Agent": "Mozilla/5.0",
    })
    retries = Retry(total=4, backoff_factor=0.4,
                    status_forcelist=(429,500,502,503,504), allowed_methods=["GET"],
                    respect_retry_after_header=True)
    adapter = HTTPAdapter(max_retries=retries, pool_connections=64, pool_maxsize=64)
    s.mount("https://", adapter); s.mount("http://", adapter)
    return s

def fetch_one(session: requests.Session, wiki_id: str) -> Tuple[str, Optional[List[Dict[str,Any]]], int, Optional[str]]:
    params = {"wiki_id": wiki_id, "source": "binbase", "isPublic": "false"}
    r = session.get(ENDPOINT, params=params, timeout=TIMEOUT)
    if r.status_code != 200:
        return wiki_id, None, r.status_code, r.text[:300]
    payload = r.json()
    analysis = payload.get("analysis", {}) if isinstance(payload, dict) else {}
    ref_hits = (analysis.get("reference_library") or {}).get("identity_search")
    if ref_hits is not None and not isinstance(ref_hits, list):
        ref_hits = None
    return wiki_id, ref_hits, 200, None

def fetch_many(ids: List[str], token: str) -> Dict[str, Optional[List[Dict[str,Any]]]]:
    s = build_session(token)
    out: Dict[str, Optional[List[Dict[str,Any]]]] = {}
    min_interval = 1.0 / RPS
    last = 0.0
    def worker(wid: str):
        nonlocal last
        wait = min_interval - (time.time() - last)
        if wait > 0: time.sleep(wait)
        last = time.time()
        return fetch_one(s, wid)
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as ex:
        futs = [ex.submit(worker, wid) for wid in ids]
        for i, fut in enumerate(as_completed(futs), 1):
            wid, hits, status, err = fut.result()
            out[wid] = hits
            print(f"[MassWiki] {i}/{len(ids)} {wid} status={status} err={err or '-'}")
    return out

def flatten(hits: Dict[str, Optional[List[Dict[str,Any]]]]) -> pd.DataFrame:
    rows = []
    for wid, arr in hits.items():
        if not arr: continue
        for i, h in enumerate(arr, 1):
            rows.append({
                "wiki_id": wid,
                "db": h.get("db") or h.get("source"),
                "id": h.get("id") or h.get("identifier") or h.get("accession"),
                "name": h.get("name"),
                "adduct": h.get("adduct"),
                "precursor_mz": h.get("precursor_mz") or h.get("precursor"),
                "entropy_similarity": h.get("entropy_similarity") or h.get("score") or h.get("similarity"),
                "library_type": h.get("library_type"),
                "rt": h.get("rt") or h.get("retention_time"),
                "ri": h.get("ri") or h.get("retention_index"),
                "rank": h.get("rank") or i,
            })
    return pd.DataFrame(rows)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True, help="Path to input CSV")
    ap.add_argument("--out", required=True, help="Path to output CSV")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--token", default=os.getenv("MASSWIKI_TOKEN") or "PASTE_TOKEN_HERE")
    args = ap.parse_args()

    ids = load_ids(args.csv)
    if args.limit: ids = ids[:args.limit]
    print(f"Reading CSV: {os.path.abspath(args.csv)}")
    print(f"Fetching {len(ids)} wiki_ids with source=binbase isPublic=false")

    hits = fetch_many(ids, token=args.token)
    df = flatten(hits)
    out_path = os.path.abspath(args.out)
    df.to_csv(out_path, index=False)
    print(f"\nSaved hits: {out_path} (rows={len(df)})")

if __name__ == "__main__":
    main()
