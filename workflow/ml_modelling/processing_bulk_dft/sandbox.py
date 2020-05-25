# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.2
#   kernelspec:
#     display_name: Python [conda env:research-new]
#     language: python
#     name: conda-env-research-new-py
# ---

# +
import os
import sys

import pickle
import os

sys.path.insert(0, os.path.join(os.environ["PROJ_irox"], "data"))
from proj_data_irox import static_irox_structures_path, bulk_dft_data_path
# -

with open(bulk_dft_data_path, "rb") as fle:
    df_bulk_dft = pickle.load(fle)

# +
from ase.visualize import view
from ase import io

# view([
#     ])

# df_bulk_dft.loc["zr9ic2zaz5"].iloc[0]["atoms"].write("__temp__/zr9ic2zaz5.cif")
# df_bulk_dft.loc["bpvynr7p9w"].iloc[0]["atoms"].write("__temp__/bpvynr7p9w.cif")

print(df_bulk_dft.loc["zr9ic2zaz5"].iloc[0]["path"])
print(df_bulk_dft.loc["bpvynr7p9w"].iloc[0]["path"])


# + active=""
#
#
#
#
#
#
# -

with open(static_irox_structures_path, "rb") as fle:
    df_static = pickle.load(fle)

# +
df_static = df_static[df_static["stoich"] == "AB3"]
df_static = df_static[df_static["source"] == "chris"]

df_static_iro3 = df_static

# +
# assert False

# +
root_data_dir = os.path.join(
    os.environ["PROJ_DATA"],
    "04_IrOx_surfaces_OER/ml_bulk_irox_dft")

with open(os.path.join(root_data_dir, "iro3", "df_dict.pickle"), "rb") as fle:
    data_iro3 = pickle.load(fle)

# +
computed_ids = data_iro3["df"]["id"].unique()

computed_ids = [int(i) for i in computed_ids]
# -

all_ids = df_static_iro3["id_old"].unique()

(list(set(computed_ids) - set(all_ids)))

# +
all_ids

computed_ids

new_ids_to_run = [i for i in all_ids if i not in computed_ids]

new_ids_to_run

# +
lst = [
    
1,
2,
"""
3,
4,
5,
"""

]

# + active=""
#
#
#

# + jupyter={}
# iro2
# iro3
# df_dict.pickle

# +
# sys.path.insert(0, os.path.join(os.environ["PROJ_irox"], "data"))
# from proj_data_irox import (
#     bulk_dft_data_path,
#     unique_ids_path,
#     prototypes_data_path,
#     static_irox_structures_path,
#     oqmd_irox_data_path,
#     voronoi_features_data_path,
#     )
