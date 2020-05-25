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

import pandas as pd

# +
data_list = [
    {
        "a": 1,
        "b": 2,
        "c": 3,
        },
    {
        "a": 4,
        "b": 5,
        "c": 6,
        },
    {
        "a": 7,
        "b": 8,
        "c": 9,
        },
    
    ]

df = pd.DataFrame(data_list)
# -

df = df.append(df.iloc[0])
df = df.reset_index(drop=True)

# +
# # df.reset_index?
# -

df

df.iloc[-1]

df.ix[-1, "a"] = "tmp"

df
