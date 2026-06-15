# ground truth for entropy distribution
import pickle
import pprint
import pandas as pd
import numpy as np
import ms_entropy


def load_pickle(file_path):
    with open(file_path, 'rb') as file:
        return pickle.load(file)
pkl_file_path = '/Users/ellayoung/Desktop/metabolo_confi_score/data/NIST23_negative_entropy_selected_adducts.pkl'
data_dict = load_pickle(pkl_file_path)

print(data_dict)


pprint.pprint(data_dict)


# Print type of the object
print(type(data_dict))

# Print available attributes and methods
print(dir(data_dict))

# If it's a custom object, you might want to print its attributes
print(vars(data_dict))


if hasattr(data_dict, 'search'):
    try:
        # Example of potential search (adjust parameters as needed)
        results = data_dict.search('metadata')
    except Exception as e:
        print("Search method error:", e)

if hasattr(data_dict, 'get_topn_matches'):
    try:
        # Example of getting top N matches
        matches = data_dict.get_topn_matches(n=5)
    except Exception as e:
        print("Top N matches error:", e)