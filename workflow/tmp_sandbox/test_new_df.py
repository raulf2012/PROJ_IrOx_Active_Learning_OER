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

# # Parsing New DataFrame for Errors
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

from dft_job_automat.job_analysis import DFT_Jobs_Analysis

Jobs = DFT_Jobs_Analysis(
    update_job_state=False,
    job_type_class=None,
    load_dataframe=True,
    root_dir='/mnt/c/Users/raul_desktop/Dropbox/01_norskov/04_comp_clusters/02_DATA/04_IrOx_surfaces_OER',
    working_dir='/mnt/c/Users/raul_desktop/Dropbox/01_norskov/04_comp_clusters/02_DATA/04_IrOx_surfaces_OER',
    dataframe_dir='/mnt/c/Users/raul_desktop/Dropbox/01_norskov/04_comp_clusters/02_DATA/04_IrOx_surfaces_OER/181226_new_job_df',
    )

df = Jobs.filter_early_revisions(Jobs.data_frame)

df[df["bulk_system"] == "IrO3_battery"]

df[
    (df["bulk_system"] == "IrO3_battery") & \
    (df["surface_type"] == "a")

    ]

# +
df_tmp = df[
    (df["job_type"] == "ORR_adsorption") & \
    (df["bulk_system"] == "IrO2") & \
    (df["facet"] == "100") & \
    (df["coverage_type"] == "o_covered")
#     (df["surface_type"] == "a")
    ]

df_tmp
# -

[print(i) for i in df_tmp["elec_energy"].tolist()]

[
-409.34226293,
-423.03283939,
-413.69666711,
-413.69268128
]



# + active=""
#
#
#
#
#
#
#
#
# -

assert False

import numpy as np

"/mnt/c/Users/raul_desktop/Dropbox/01_norskov/04_comp_clusters/02_DATA/04_IrOx_surfaces_OER/181226_new_job_df"
df = pickle.load(
    open(
    os.path.join(
        os.environ["PROJ_DATA"],
        "04_IrOx_surfaces_OER/181226_new_job_df",
        'job_dataframe.pickle',    
        ),
        "rb"
        )
    )

df.head()
df["surface_type"] = df["surface_type"].replace(np.nan, "tmp", regex=True)
df.head()

df["surface_type"]

df

df.head()

# +
df = df[df["job_type"] == "ORR_adsorption"]

df = df[df["bulk_system"] == "IrO3_battery"]

grouped = df.groupby(unique_params)

len(df)
# -

for i in grouped:
    tmp = 42
    print("42")

i[1]

unique_params = ["facet", "coverage_type", "bulk_system"]



df

# + active=""
#
#
#
#
#
#
#
#
#
#
#

# +
print("Length of DataFrame:", len(df))
print("Column Names:", "\n", list(df))

print("")
print(df.info(memory_usage="deep"))
# -

list(set(df["bulk_system"].tolist()))

df_1 = df[df["bulk_system"] != "IrO3_battery"]

df_1.info(memory_usage="deep")

df_1[df_1["bulk_system"].isna()].iloc[0]["path"]

# +
# df_1
# -

df[df["bulk_system"] == "IrO3_battery"]

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
