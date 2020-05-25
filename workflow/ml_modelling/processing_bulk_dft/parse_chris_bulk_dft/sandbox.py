# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.2
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# +
import pickle
import os

path_i = os.path.join("df_dft_calcs.pickle")
with open(path_i, "rb") as fle:
    df = pickle.load(fle)

# +
len(df)

df_iro2 = df[df["stoich"] == "AB2"]
df_iro3 = df[df["stoich"] == "AB3"]
