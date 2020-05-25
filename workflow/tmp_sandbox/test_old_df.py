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

# # Parsing OLD DataFrame for Errors
# The new dataframe that recently parsed (181226) is much smaller than
# the previous one (in terms of storage memory)
#
# Hopefully none of the old data was purged on NERSC!!!!

# +
import sys
import os

import pickle

import pandas as pd
pd.set_option("display.max_columns", None)
pd.set_option('display.max_rows', None)
# -

df = pickle.load(
    open(
    os.path.join(
        os.environ["PROJ_DATA"],
        "04_IrOx_surfaces_OER",
#         'job_dataframe_181225_with_battery.pickle',
        "job_dataframe.pickle",
        ),
        "rb"
        ),
    encoding="latin1",
    )

# +
print("Length of DataFrame:", len(df))
print("Column Names:", "\n", list(df))

print("")
print(df.info(memory_usage="deep"))
# -

list(set(df["bulk_system"].tolist()))

df

# + active=""
#
#
#
#
#
#

# +
# os.listdir(
#     os.path.join(
#         os.environ["PROJ_DATA"],
#         "04_IrOx_surfaces_OER",
#         )
#     )
