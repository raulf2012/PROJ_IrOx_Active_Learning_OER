# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.4.2
#   kernelspec:
#     display_name: Python [conda env:PROJ_irox] *
#     language: python
#     name: conda-env-PROJ_irox-py
# ---

# # Parsing OQMD Structures and Energies into Dataframe
# ---

# # Import Modules

# +
import os
import sys

import pickle
import pandas as pd

print(os.getcwd())

# +
sys.path.insert(0, os.path.join(os.environ["PROJ_irox"], "data"))

from proj_data_irox import (
    bulk_dft_data_path,
    unique_ids_path,
    prototypes_data_path,
    static_irox_structures_path,
    oqmd_irox_data_path,
#     voronoi_features_data_path,
    )
# -

# # Read Data

# +
path_i = os.path.join(
    os.environ["PROJ_irox"],
    "data/ml_irox_data",
    "oqmd_data.csv")

df_oqmd_data = pd.read_csv(path_i)

df_oqmd_data.set_index("id_unique", inplace=True)
df_oqmd_data = df_oqmd_data.drop(labels="source", axis=1)
# -

df_id = pd.read_csv(unique_ids_path)

with open(static_irox_structures_path, "rb") as fle:
    df_structures = pickle.load(fle)

static_irox_structures_path

# Only processing OQMD entries

df_structures = df_structures.loc[df_id[df_id["source"] == "oqmd"]["unique_ids"]]

# # Merging OQMD data with Structures Data

# +
df_merged = df_structures.merge(
    df_oqmd_data.drop(["id_old", "stoich"], axis=1),  # Drop duplicate columns
    how='left',
    left_index=True, right_index=True,
    indicator=True)

# df_merged

# + active=""
# Check that the `_merge` column says `both`, i.e. all rows exist in both dfs

# +
all_good = False

merged_col_unique = list(set(df_merged["_merge"].tolist()))
if len(merged_col_unique) == 1:
    if merged_col_unique[0] == "both":
        print("All good!")
        all_good = True
    else:
        print("Not good!! 00")
else:
    print("Not good!! 11")
    
if all_good:
    df_merged.drop(["_merge"], axis=1, inplace=True)
# -

df_merged = df_merged.reset_index()
df_merged = df_merged.drop("id_unique", axis=1)

# # Save Data

# +
import pickle

directory = "out_data"
if not os.path.exists(directory): os.makedirs(directory)

with open(os.path.join(directory, "df_oqmd_data.pickle"), "wb") as fle:
    pickle.dump(df_merged, fle)
# -

print(20 * "# # ")
print("All done!")
assert False
