import os
import time
import pandas as pd
import requests

BASE = "https://masswiki.us-west-2.elasticbeanstalk.com/analysis/get_data"
ACCESS_TOKEN = "eyJraWQiOiJoeCtPbm1BUmpUWG0rZnNzZnptYVR2b2RIeG1Ra0dKbGVzc1hsZG5oTG5nPSIsImFsZyI6IlJTMjU2In0.eyJzdWIiOiI4ODYxYzM2MC1lMGExLTcwOTMtM2JjNC0wZDkxZDlmOGFkN2UiLCJjb2duaXRvOmdyb3VwcyI6WyJtYXNzd2lraS1sYWIiLCJtYXNzd2lraS1jb21tdW5pdHkiLCJVc2VycyJdLCJlbWFpbF92ZXJpZmllZCI6dHJ1ZSwiaXNzIjoiaHR0cHM6XC9cL2NvZ25pdG8taWRwLnVzLXdlc3QtMi5hbWF6b25hd3MuY29tXC91cy13ZXN0LTJfR2p0Y00wUENwIiwiY29nbml0bzp1c2VybmFtZSI6Ijg4NjFjMzYwLWUwYTEtNzA5My0zYmM0LTBkOTFkOWY4YWQ3ZSIsIm9yaWdpbl9qdGkiOiJiZWZjYzhiNS0xOWQ5LTQxOGItYWMxOS05ZWE0NmU3NWQ5MTIiLCJhdWQiOiIzaGdvczNhdGQxZWwxNmx0aDYxN2lpY29hbCIsImV2ZW50X2lkIjoiM2Q1MWZlZDAtNmVmMC00Y2U4LWEzNjAtZTkwMmY3N2NlZjI3IiwidG9rZW5fdXNlIjoiaWQiLCJhdXRoX3RpbWUiOjE3NjA2ODQzMzUsIm5hbWUiOiJaaXl1ZSBZYW5nIiwiZXhwIjoxNzYwNjg3OTM1LCJpYXQiOjE3NjA2ODQzMzUsImZhbWlseV9uYW1lIjoiWWFuZyIsImp0aSI6IjZjMDIxY2IzLTQ0NmMtNGUwZi1hMjk2LTZlOTYzNDEyNDA2ZiIsImVtYWlsIjoienl6eWFuZ0B1Y2RhdmlzLmVkdSJ9.qXuPcmQycPqwU8OF-wMHoN4wSgASXNW9XvOTfx2dcxIVAffgqKK0WTlS40kfUX_MXBaaWI5xNjxxxZkna5X9JGDX24CoxuDv07hqge0i0IIKHH75Qu-MFVudJ_38vsAOOX9KUCjKdDgCV7eqbqABD4oPDugFO9CktLd76f8XzyUzujfbdJ_xv0Zy6i7OjHFlkCxnED13a-1d94PqgxJFhlgndLph8ZJM2_Bb7eJ4Q_W7AmlKqCXm8qKUPLpeqh7bqhorfnwHcFRrkvfdgOZBzlyfR-l0WqnvcKS4qTVZ9DuDDRIWP5leA-TQ7LTbeLXmW1S21F2gFogKZLwgR_wB2A"
INPUT_CSV = "/Users/ellayoung/Desktop/metabolo_confi_score/data/ttof+neg+hilic.csv"
OUTPUT_CSV = "/Users/ellayoung/Desktop/metabolo_confi_score/data/hilic_masswiki_reference_hits.csv"

# -------------------------------------------------------------------
# 1) Load ALL wiki_ids (no filtering)
# -------------------------------------------------------------------
from typing import List  # add this near your imports

def load_all_wiki_ids(csv_path: str) -> List[str]:
    df = pd.read_csv(csv_path)
    if "wiki_id" not in df.columns:
        raise ValueError("CSV must include a 'wiki_id' column.")
    wiki_ids = (
        df["wiki_id"]
        .astype(str)
        .map(lambda s: s.strip())
        .replace({"": None, "nan": None})
        .dropna()
        .tolist()
    )
    return wiki_ids


# -------------------------------------------------------------------
# 2) Fetch one wiki_id (reference library only)
# -------------------------------------------------------------------
def fetch_one(wiki_id: str):
    params = {
        "wiki_id": wiki_id,   # raw; DO NOT pre-encode
        "source": "binbase",
        "isPublic": "false",
    }
    headers = {
        "Accept": "application/json",
        "Authorization": f"Bearer {ACCESS_TOKEN}",
    }
    r = requests.get(BASE, params=params, headers=headers, timeout=20)
    if r.status_code != 200:
        print(f"{wiki_id}: HTTP {r.status_code} {r.text[:300]}")
        return None
    payload = r.json()
    analysis = payload.get("analysis", {}) if isinstance(payload, dict) else {}
    ref_hits = (analysis.get("reference_library") or {}).get("identity_search")
    if ref_hits is not None and not isinstance(ref_hits, list):
        ref_hits = None
    return ref_hits

# -------------------------------------------------------------------
# 3) Flatten hits to DataFrame
# -------------------------------------------------------------------
def flatten_hits(all_hits: dict) -> pd.DataFrame:
    rows = []
    for wid, hits in all_hits.items():
        if not hits:
            continue
        for i, h in enumerate(hits, 1):
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

# -------------------------------------------------------------------
# 4) Main loop (no filtering)
# -------------------------------------------------------------------
# if __name__ == "__main__":
#     wiki_ids = load_all_wiki_ids(INPUT_CSV)
#     print(f"Found {len(wiki_ids)} wiki_ids (no filtering). Fetching reference hits...")

#     all_hits = {}
#     for idx, wid in enumerate(wiki_ids, 1):
#         try:
#             hits = fetch_one(wid)
#             all_hits[wid] = hits
#             if idx % 50 == 0:
#                 print(f"Processed {idx}/{len(wiki_ids)} wiki_ids...")
#         except Exception as e:
#             print(f"{wid}: ERROR {e}")
#             all_hits[wid] = None
#         time.sleep(0.2)  # ~5 req/s polite throttle

#     df_hits = flatten_hits(all_hits)
#     df_hits.to_csv(OUTPUT_CSV, index=False)
#     print(f"\nSaved all reference-library hits to: {OUTPUT_CSV} (rows={len(df_hits)})")


if __name__ == "__main__":
    wiki_ids = load_all_wiki_ids(INPUT_CSV)
    print(f"Found {len(wiki_ids)} wiki_ids (no filtering). Fetching reference hits...")

    # Skip the first 5749 items, start at the 5750th wiki_id
    wiki_ids = wiki_ids[1:]   # index 5749 is the 5750th element

    all_hits = {}
    for idx, wid in enumerate(wiki_ids, 0):  # keep the real index count
        try:
            hits = fetch_one(wid)
            all_hits[wid] = hits
            if idx % 50 == 0:
                print(f"Processed {idx}/{len(load_all_wiki_ids(INPUT_CSV))} wiki_ids...")
        except Exception as e:
            print(f"{wid}: ERROR {e}")
            all_hits[wid] = None
        time.sleep(0.2)  # ~5 req/s polite throttle

    df_hits = flatten_hits(all_hits)
    df_hits.to_csv(OUTPUT_CSV, index=False)
    print(f"\nSaved all reference-library hits to: {OUTPUT_CSV} (rows={len(df_hits)})")
