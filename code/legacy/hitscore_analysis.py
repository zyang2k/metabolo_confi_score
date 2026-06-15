import pandas as pd
import matplotlib.pyplot as plt

# Load the CSV
df = pd.read_csv("/Users/ellayoung/Desktop/metabolo_confi_score/out_0918/assertions.csv")

subset = df[df["entropy_similarity"] > 0.99]

print(subset.head(30))


