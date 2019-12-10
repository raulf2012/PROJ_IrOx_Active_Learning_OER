# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.0
#   kernelspec:
#     display_name: Python [conda env:PROJ_IrOx_Active_Learning_OER]
#     language: python
#     name: conda-env-PROJ_IrOx_Active_Learning_OER-py
# ---

# + [markdown] Collapsed="false"
# # New ML Active Learning Workflow
# ---
#
# 32 (Too big, not computed),
#
# 226 (Not computed),

# + [markdown] Collapsed="false"
# # Import Modules

# + jupyter={"source_hidden": true} Collapsed="false"
# %%capture
#| - OUT_OF_SIGHT
import os
import sys

import json
import pickle

import time

import itertools

import pandas as pd
import numpy as np

from ase.visualize import view

import chart_studio.plotly as py
import plotly.graph_objs as go

sys.path.insert(0, os.path.join(os.environ["PROJ_irox"], "data"))
from proj_data_irox import (
    bulk_dft_data_path, unique_ids_path,
    df_features_pre_opt_path,
    df_features_post_opt_path,
    ids_to_discard__too_many_atoms_path,
    )

from plotting.my_plotly import my_plotly_plot

import pprint
pp = pprint.PrettyPrinter()

sys.path.insert(0,
    os.path.join(os.environ["PROJ_irox"], "workflow/ml_modelling"))

from ml_methods import create_mixed_df

from ase_modules.ase_methods import view_in_vesta

# + [markdown] Collapsed="false"
# # Script Inputs

# + Collapsed="false"
# stoich_i = "AB3"
stoich_i = "AB2"

# + [markdown] Collapsed="false"
# # Read Data

# + Collapsed="false"
# with open(bulk_dft_data_path, "rb") as fle:
#     df_bulk_dft = pickle.load(fle)

df_ids = pd.read_csv(unique_ids_path)
# + [markdown] Collapsed="false"
# # Duplicate Analysis

# + Collapsed="false"
sys.path.insert(0, os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling"))
from ml_methods import get_data_for_al
from ccf_similarity.ccf import CCF

out_dict = get_data_for_al(
    stoich=stoich_i, verbose=False,
    drop_too_many_atoms=True)

df_dij = out_dict["df_dij"]
df_bulk_dft = out_dict["df_bulk_dft"]

CCF_i = CCF(df_dij=df_dij, d_thresh=0.02)

# + Collapsed="false"
df_bulk_dft = df_bulk_dft[df_bulk_dft.source == "raul"]

ids_to_drop = []
for id_i in df_bulk_dft.index.tolist():
    simil_dict_i = CCF_i.i_all_similar(id_i)

    # for key, val in simil_dict_i.items():
    #     tmp = 42

    if simil_dict_i is not None:
        similar_ids = [id_i] + list(simil_dict_i.keys())

        df_i = df_bulk_dft.loc[similar_ids]

        ids_to_drop_i = df_i.sort_values("energy_pa").iloc[1:].index.tolist()

        ids_to_drop.extend(ids_to_drop_i)

        
ids_to_drop__duplicates = ids_to_drop

# + Collapsed="false"
# df_dij.loc[[
# #     "zimixdvdxd",
#     "6fcdbh9fz2",
#     ]]

# df_dij.loc["6fcdbh9fz2"]

# df_tmp = df_dij
# "6fcdbh9fz2" in df_tmp.index.tolist()

# "6fcdbh9fz2" in df_bulk_dft.index.tolist()

# + [markdown] Collapsed="false"
# # Filtering dataframes to the correct stoicheometry

# + [markdown] Collapsed="false"
# # TEMP DROP DUPLICATE and OUTLIER SYSTEMS

# + Collapsed="false"
# path_i = os.path.join(
#     os.environ["PROJ_irox"],
#     "workflow/ml_modelling/00_ml_workflow/191102_new_workflow/00_abx_al_runs/out_data",
#     "duplicates.pickle")
# with open(path_i, "rb") as fle:
#     duplicates_dict = pickle.load(fle)
#     ids_to_drop__duplicates = duplicates_dict[stoich_i]

# #############################################################################
path_i = os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling/visualizing_data/out_data",
    "outlier_features.json")
with open(path_i, 'r') as f:
    ids_to_drop__outliers = json.load(f)

with open(ids_to_discard__too_many_atoms_path, "rb") as fle:
    ids_to_drop__too_many_atoms = pickle.load(fle)

# #############################################################################
ids_to_drop = [] + \
    ids_to_drop__duplicates + \
    []
    # ids_to_drop__too_many_atoms + \
    # ids_to_drop__outliers + \

print("len(ids_to_drop):", len(ids_to_drop))
ids_to_drop = list(set(ids_to_drop))
print("len(ids_to_drop):", len(ids_to_drop))

# + Collapsed="false"
# #############################################################################
# Filter ids ##################################################################
df_ids = df_ids[
    (df_ids["stoich"] == stoich_i) & \
    (df_ids["source"] != "oqmd") & \
    (df_ids["source"] != "raul_oer") & \
    [True for i in range(len(df_ids))]
    ]

print("df_ids.shape:", df_ids.shape)
# IDS TO DROP
df_ids = df_ids[~df_ids["unique_ids"].isin(ids_to_drop)]
unique_ids = df_ids["unique_ids"].tolist()

index_filter = np.intersect1d(df_bulk_dft.index, unique_ids)
df_bulk_dft = df_bulk_dft.loc[index_filter]
df_bulk_dft = df_bulk_dft[df_bulk_dft["source"] != "chris"]

# + Collapsed="false"
directory = "out_data/" + stoich_i + "_structures_all"
if not os.path.exists(directory):
    os.makedirs(directory)


directory = "out_data/" + stoich_i + "_structures"
if not os.path.exists(directory):
    os.makedirs(directory)

# + Collapsed="false"
cols_to_drop = [
    "atoms",
    "form_e_chris",
    "id",
    "path",
    "source",
    "stoich",
    ]

df_select = df_bulk_dft.drop(cols_to_drop, axis=1)

# + Collapsed="false"
df_bulk_dft = df_bulk_dft.sort_values("energy_pa")

df_bulk_dft["energy_order_id"] = [i for i in range(len(df_bulk_dft))]

df_select = df_bulk_dft.drop(cols_to_drop, axis=1)

# + Collapsed="false"
import ase

ase.__version__

# + Collapsed="false"
stoich_i

# + Collapsed="false"
df_select.to_csv("out_data/data_table_" + stoich_i + ".csv")

for i_cnt, row_i in df_bulk_dft.iterrows():
    atoms = row_i["atoms"]
    file_name_i = "" + \
        str(row_i["energy_order_id"]).zfill(3) + \
        "__" + \
        "id-unique" + \
        "_" + \
        row_i.name + \
        "__" + \
        "id-short" + \
        "_" + \
        str(row_i["id_old"]).zfill(3) + \
        ".cif"
    
    # atoms.write("out_data/" + stoich_i + "_structures_all/" + file_name_i)
    atoms.write("out_data/" + stoich_i + "_structures/" + file_name_i)

# + Collapsed="false"
# df_bulk_dft.loc["cubqbpzd7k"]

# df_bulk_dft.sort_values("energy_pa").iloc[0:3].index
df_bulk_dft.sort_values("energy_pa")

# + Collapsed="false"
# df_dij.loc[
#     ["xw9y6rbkxr", "8p8evt9pcg"],
#     ["xw9y6rbkxr", "8p8evt9pcg"],
#     ]

# + Collapsed="false" active=""
#
#
#
#
#

# + jupyter={"source_hidden": true} Collapsed="false"
# df_bulk_dft[df_bulk_dft["stoich"] == "AB3"].sort_values("energy_pa")

# -6.469847
# -6.467450

# + jupyter={"source_hidden": true} Collapsed="false"
# -6.469847 - -6.46745

# + jupyter={"source_hidden": true} Collapsed="false"
# assert False
