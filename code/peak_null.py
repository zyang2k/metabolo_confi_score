import numpy as np
import pandas as pd
import os
import ms_entropy as me

df = pd.read_csv('/Users/ellayoung/Desktop/metabolo_confi_score/parsed_spectra.csv')


df_clean = df.iloc[:, [0, 1, 10, 18, 13, 14, 16, 24]]

# only include entries with M+H or M-H adducts

df_clean_mh = df_clean[df_clean["Precursor Type"].isin(["[M+H]+", "[M-H]-"])]




df_unique_mh = df_clean_mh.dropna(subset=['InChIKey'])
df_unique_mh = df_unique_mh.drop_duplicates(subset=['InChIKey'], keep='first')
# # Convert NaN or float SMILES to empty strings and ensure dtype is string
df_unique_mh["SMILES"] = df_unique_mh["SMILES"].fillna("").astype(str)

# Create a helper column with the first 14 characters of InChIKey
df_unique_mh["inchi14"] = df_unique_mh["InChIKey"].astype(str).str[:14]

def parse_peaks_string(peaks_str):
    """
    Parse a peaks string like:
    '[[136.0435   3.4   ]\n [214.9619  38.4   ]\n ...]'
    into a list of lists of floats.
    """
    # Remove the outer brackets and any leading/trailing whitespace
    cleaned = peaks_str.strip("[] \n")
    # Split into lines (each line corresponds to one peak row)
    rows = cleaned.split('\n')
    result = []
    for row in rows:
        # Remove any remaining brackets and extra whitespace
        row_clean = row.strip(" []")
        if row_clean:
            # Split on whitespace and convert each element to float
            values = row_clean.split()
            result.append([float(val) for val in values])
    return result

def safe_parse_peaks(peaks):
    """
    Safely parse peaks, returning an empty list for invalid cases.
    """
    if isinstance(peaks, str):
        if "..." in peaks or peaks.strip() == "":
            return []  # Return empty list for invalid entries
        try:
            return parse_peaks_string(peaks)
        except ValueError:
            return []  # Handle unexpected format issues
    return peaks  # Already parsed peaks



df_unique_mh["peaks"] = df_unique_mh["peaks"].apply(safe_parse_peaks)