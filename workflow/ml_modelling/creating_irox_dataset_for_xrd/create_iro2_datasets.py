# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.4'
#       jupytext_version: 1.1.7
#   kernelspec:
#     display_name: Python [conda env:research-new]
#     language: python
#     name: conda-env-research-new-py
# ---

# # New ML Active Learning Workflow
# ---
#
# 32 (Too big, not computed),
#
# 226 (Not computed),

# # Import Modules

# + {"jupyter": {"source_hidden": true}}
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
# -

# # Script Inputs

# stoich_i = "AB3"
stoich_i = "AB2"

# # Read Data

# +
with open(bulk_dft_data_path, "rb") as fle:
    df_bulk_dft = pickle.load(fle)

df_ids = pd.read_csv(unique_ids_path)
# -
# # Filtering dataframes to the correct stoicheometry

# # TEMP DROP DUPLICATE and OUTLIER SYSTEMS

# +
path_i = os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling/ccf_similarity_analysis/out_data",
    "all_ids_to_elim.pickle")
with open(path_i, "rb") as fle:
    ids_to_drop__duplicates = pickle.load(fle)
    ids_to_drop__duplicates = ids_to_drop__duplicates[stoich_i]

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
    ids_to_drop__outliers + \
    ids_to_drop__duplicates + \
    ids_to_drop__too_many_atoms + \
    []

print("len(ids_to_drop):", len(ids_to_drop))
ids_to_drop = list(set(ids_to_drop))
print("len(ids_to_drop):", len(ids_to_drop))

# + {"jupyter": {"source_hidden": true}}
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

# + {"jupyter": {"source_hidden": true}}
directory = "out_data/" + stoich_i + "_structures"
if not os.path.exists(directory):
    os.makedirs(directory)

# + {"jupyter": {"source_hidden": true}}
cols_to_drop = [
    "atoms",
    "form_e_chris",
    "id",
    "path",
    "source",
    "stoich",
    ]

df_select = df_bulk_dft.drop(cols_to_drop, axis=1)

# + {"jupyter": {"source_hidden": true}}
df_bulk_dft = df_bulk_dft.sort_values("energy_pa")

df_bulk_dft["energy_order_id"] = [i for i in range(len(df_bulk_dft))]

df_select = df_bulk_dft.drop(cols_to_drop, axis=1)

# +
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
    
    atoms.write("out_data/" + stoich_i + "_structures/" + file_name_i)
