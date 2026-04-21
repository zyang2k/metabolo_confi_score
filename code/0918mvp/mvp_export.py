
# mvp_export.py
import pandas as pd
import numpy as np

# -------- Module 07 · Exports & QA (minimal) --------

def make_assertion_table(df_post: pd.DataFrame, features_df: pd.DataFrame = None) -> pd.DataFrame:
    cols = [
        "wiki_id","library_id","name","theta","LR","C_loc","PEP_local","P_NoTA",
        "delta_theta","set_size","member_hits","quality_mean"
    ]
    # Attach features if provided
    if features_df is not None:
        feat_keep = ["wiki_id","library_id","delta_ppm","delta_rt","zrt","entropy_similarity"]
        features = features_df[feat_keep].drop_duplicates(["wiki_id","library_id"], keep="first")
        df = df_post.merge(features, on=["wiki_id","library_id"], how="left")
    else:
        df = df_post.copy()
    

    # Ensure required cols exist
    for c in cols:
        if c not in df.columns and c not in ["LR","C_loc"]:
            df[c] = np.nan
    # Compute LR if missing
    if "LR" not in df.columns and "theta" in df.columns:
        df["LR"] = np.exp(df["theta"])
    # Local rank if missing
    if "C_loc" not in df.columns and "theta" in df.columns:
        def _softmax(g):
            x = g["theta"].to_numpy()
            m = np.max(x)
            ex = np.exp(x - m)
            g["C_loc"] = ex / ex.sum()
            return g
        df = df.groupby("wiki_id", as_index=False, group_keys=False).apply(_softmax)
    df['post'] = (1.0 - df['P_NoTA']) * df['C_loc']    # spectrum-level calibrated confidence
    df['PEP']  = 1.0 - df['post']

    cols = [
    "wiki_id","library_id","name","theta","LR","C_loc","PEP_local","P_NoTA",
    "delta_theta","set_size","member_hits","quality_mean",
    "post","PEP"   # <-- add these
]
    keep_cols = list(dict.fromkeys(cols + ["delta_ppm","delta_rt","zrt","entropy_similarity","C_loc","LR"]))
    return df[keep_cols].copy()


def make_top_calls(assertions: pd.DataFrame, nota_threshold: float = 0.6) -> pd.DataFrame:
    def _top(g):
        # abstain if P_NoTA high
        abstain = bool((g["P_NoTA"].iloc[0] if "P_NoTA" in g.columns else 0.0) >= nota_threshold)
        if abstain:
            return pd.Series({
                "wiki_id": g["wiki_id"].iloc[0],
                "library_id_top": None,
                "name_top": None,
                "C_loc_top": np.nan,
                "PEP_top": np.nan,
                "abstain": True,
                "reason": "NoTA"
            })
        # else pick max theta
        row = g.loc[g["theta"].idxmax()]
        return pd.Series({
            "wiki_id": row["wiki_id"],
            "library_id_top": row["library_id"],
            "name_top": row.get("name", None),
            "C_loc_top": row.get("C_loc", np.nan),
            "PEP_top": row.get("PEP_local", np.nan),
            "abstain": False,
            "reason": ""
        })
    return assertions.groupby("wiki_id", as_index=False).apply(_top).reset_index(drop=True)
