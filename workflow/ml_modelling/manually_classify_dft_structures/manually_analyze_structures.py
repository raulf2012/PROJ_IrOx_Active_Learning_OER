# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.2
#   kernelspec:
#     display_name: Python [conda env:PROJ_IrOx_Active_Learning_OER]
#     language: python
#     name: conda-env-PROJ_IrOx_Active_Learning_OER-py
# ---

# # Import Modules

# +
import os
import sys
import pickle

sys.path.insert(0, os.path.join(os.environ["PROJ_irox"], "data"))
from proj_data_irox import (
    bulk_dft_data_path,
    # unique_ids_path,
    # prototypes_data_path,
    # static_irox_structures_path,
    # oqmd_irox_data_path,
    # voronoi_features_data_path,
    )

# #############################################################################
from ase_modules.ase_methods import view_in_vesta

# +
path_i = os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling",
    "ccf_similarity_analysis/out_data",
    "all_ids_to_elim_1.pickle")
with open(path_i, "rb") as fle:
    all_ids_to_elim = pickle.load(fle)

with open(bulk_dft_data_path, "rb") as fle:
    df_bulk_dft = pickle.load(fle)

# +
df_bulk_dft = df_bulk_dft[
    (df_bulk_dft["source"] != "oqmd") & \
    (df_bulk_dft["source"] != "raul_oer") & \
    (df_bulk_dft["source"] != "chris") & \
    [True for i in range(len(df_bulk_dft))]
    ]

df_bulk_dft = df_bulk_dft.drop(all_ids_to_elim)

# +
# df_ab3.index.tolist()
# -

df_bulk_dft.shape

# +
df_ab3 = df_bulk_dft[df_bulk_dft["stoich"] == "AB3"]
df_ab3 = df_ab3.sort_values("energy_pa")

tmp = [print(i) for i in df_ab3.iloc[0:10].index.tolist()]

print("")

ind_i = 15
ind_f = 20

view_in_vesta(
    df_ab3.iloc[ind_i:ind_f]["atoms"].tolist(),
    name_list=df_ab3.iloc[ind_i:ind_f].index.tolist(),
    )
# -

ind = 1
print(df_ab3.iloc[ind])
df_ab3.iloc[ind]["atoms"]
