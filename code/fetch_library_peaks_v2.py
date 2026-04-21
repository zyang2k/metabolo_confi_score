"""
fetch_library_peaks_v2.py — Stage 1 of the Bayesian two-tier pipeline.

For each unique library_wiki_id referenced by Stage 0's hits CSV that is NOT
already present in data/library_peaks_cache.json, fetch its peak list from
MassWiki and add it to the cache.

Outputs:
  data/library_peaks_cache.json (extended in place)
  data/library_peaks_fetch_errors_v2.csv (if any errors)

Modes:
  --gap-only         Print the gap (missing-vs-cached) and exit. No fetching.
  --limit N          Fetch only first N missing library_wiki_ids (testing).
  --ik14-restrict    Only fetch library_wiki_ids of IK14-matched annotation
                     hits + their sim_gap top-2 competitors. (Smaller fetch
                     scope; needs annotation SMILES → IK14 lookup.)

Token refresh: re-reads /tmp/.mw_token every 200 calls. To refresh mid-run,
just overwrite that file with a new token; the script picks it up.
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
SAVE_EVERY = 500

HITS_CSV    = "data/orbitrap_hits_v2.csv"
LIB_CACHE   = "data/library_peaks_cache.json"
OUT_ERR     = "data/library_peaks_fetch_errors_v2.csv"
TOKEN_FILE  = "/tmp/.mw_token"
TOKEN_RECHECK_EVERY = 200


# ── Token state (mutable, picked up by re-read) ──────────────────────────────

class TokenState:
    def __init__(self, path):
        self.path = path
        self.value = self._read()
    def _read(self):
        try:
            return Path(self.path).read_text().strip()
        except Exception:
            return ""
    def reread(self):
        new = self._read()
        if new and new != self.value:
            print(f"[token] refreshed (length {len(new)})")
            self.value = new
        return self.value


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


def fetch_lib_peaks(session, lib_wid):
    """Call API with wiki_id=<library_wiki_id>; extract spectrum.peaks."""
    params = {"wiki_id": lib_wid, "source": "binbase", "isPublic": "false"}
    try:
        r = session.get(ENDPOINT, params=params, timeout=TIMEOUT)
    except Exception as e:
        return None, -1, f"req-exc: {e}"
    if r.status_code != 200:
        return None, r.status_code, r.text[:200]
    try:
        payload = r.json()
    except Exception as e:
        return None, 200, f"json-exc: {e}"
    spec = (payload or {}).get("spectrum") or {}
    peaks = spec.get("peaks")
    if isinstance(peaks, list) and peaks:
        return peaks, 200, ""
    return None, 200, "no-peaks-in-response"


def compute_target_set(args) -> list[str]:
    print("Loading hits CSV...")
    hits = pd.read_csv(HITS_CSV, low_memory=False)
    print(f"  hit rows: {len(hits):,}")

    if args.ik14_restrict:
        # Restricted scope: only IK14-matched annotation hit + top-2 different-IK14
        # competitors per spectrum. Requires SMILES->IK14 lookup; defer to a
        # follow-up impl. For now just warn and treat as "all".
        print("  --ik14-restrict not implemented yet; falling back to all.")

    target = set(hits["library_wiki_id"].dropna().astype(str).unique())
    print(f"  unique library_wiki_id in hits: {len(target):,}")

    print(f"Loading existing library cache: {LIB_CACHE}")
    if os.path.exists(LIB_CACHE):
        with open(LIB_CACHE) as f:
            cached = json.load(f)
        cached_keys = set(cached.keys())
        print(f"  existing cache entries: {len(cached_keys):,}")
    else:
        cached = {}
        cached_keys = set()
        print("  no existing cache; will create new one.")

    missing = sorted(target - cached_keys)
    print(f"  GAP (need to fetch): {len(missing):,}")
    print(f"  ETA at 6 RPS: {len(missing)/6/60:.1f} min "
          f"({len(missing)/6/3600:.2f} hr)")
    return missing, cached


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gap-only", action="store_true",
                    help="Print gap then exit; do not fetch.")
    ap.add_argument("--limit", type=int, default=None,
                    help="Limit number of missing library_wiki_ids to fetch.")
    ap.add_argument("--ik14-restrict", action="store_true",
                    help="(NYI) Restrict to IK14-matched + sim_gap candidates.")
    ap.add_argument("--token-file", default=TOKEN_FILE)
    args = ap.parse_args()

    missing, cached = compute_target_set(args)

    if args.gap_only:
        return

    if args.limit:
        missing = missing[:args.limit]
        print(f"  --limit applied: fetching {len(missing):,} (test mode)")

    if not missing:
        print("Nothing to fetch.")
        return

    if not os.path.exists(args.token_file):
        raise SystemExit(f"Token file not found: {args.token_file}")

    tok = TokenState(args.token_file)
    if not tok.value:
        raise SystemExit(f"Token file {args.token_file} is empty")
    session = build_session(tok.value)

    errors = []
    fetched = 0
    last_ts = [0.0]
    min_interval = 1.0 / RPS

    def worker(lib_wid):
        wait = min_interval - (time.time() - last_ts[0])
        if wait > 0:
            time.sleep(wait)
        last_ts[0] = time.time()
        peaks, status, err = fetch_lib_peaks(session, lib_wid)
        return lib_wid, peaks, status, err

    total = len(missing)
    t0 = time.time()
    saved_at = 0
    print(f"\nFetching {total:,} missing library spectra at {RPS} RPS...")

    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        futs = [ex.submit(worker, lid) for lid in missing]
        for done_i, fut in enumerate(as_completed(futs), 1):
            lib_wid, peaks, status, err = fut.result()
            if peaks is not None:
                cached[lib_wid] = peaks
                fetched += 1
            else:
                errors.append({"library_wiki_id": lib_wid, "status": status,
                               "error": err})

            # Periodic token refresh
            if done_i % TOKEN_RECHECK_EVERY == 0:
                new_token = tok.reread()
                if new_token != session.headers.get("Authorization", "").replace("Bearer ", ""):
                    session.headers["Authorization"] = f"Bearer {new_token}"

            # Progress line
            if done_i % 100 == 0 or done_i == total:
                elapsed = time.time() - t0
                rate = done_i / max(elapsed, 1e-6)
                eta = (total - done_i) / max(rate, 1e-6)
                print(f"[{done_i:>6}/{total}] ok={fetched} errors={len(errors)} "
                      f"rate={rate:.1f}/s eta={eta:.0f}s")

            # Periodic cache save
            if done_i - saved_at >= SAVE_EVERY:
                with open(LIB_CACHE, "w") as f:
                    json.dump(cached, f)
                saved_at = done_i

    # Final save
    with open(LIB_CACHE, "w") as f:
        json.dump(cached, f)
    if errors:
        pd.DataFrame(errors).to_csv(OUT_ERR, index=False)
    print(f"\n✓ done. fetched={fetched}/{total} cached_total={len(cached):,} errors={len(errors)}")


if __name__ == "__main__":
    main()
