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
print(os.getcwd())
import sys

import pickle

import numpy as np

import pandas as pd

sys.path.insert(0, os.path.join(os.environ["PROJ_irox"], "data"))
from proj_data_irox import (
    static_irox_structures_path,
    bulk_dft_data_path,
    unique_ids_path,
    )

sys.path.insert(0, os.path.join(os.environ["PROJ_irox"], "workflow/ml_modelling"))

from StructurePrototypeAnalysisPackage.ccf import struc2ccf

# +
# %%capture

from ml_methods import get_ml_dataframes

DATAFRAMES = get_ml_dataframes(
    names=[
        'bulk_dft_data_path',
        'unique_ids_path',
        'prototypes_data_path',
        'static_irox_structures_path',
        'static_irox_structures_kirsten_path',
        'oqmd_irox_data_path',
        'df_features_pre_opt_path',
        'df_features_pre_opt_kirsten_path',
        'df_features_post_opt_path',
        'oer_bulk_structures_path',
        'df_ccf_path',
        'df_dij_path',
        'ids_to_discard__too_many_atoms_path',
        ]
    )

df_dij = DATAFRAMES["df_dij"]
# -

# # Script Inputs

# +
r_cut_off = 10.
r_vector = np.arange(1, 10, 0.02)

mean_density = 0.08407356
# -

directory = "out_data"
if not os.path.exists(directory):
    os.makedirs(directory)

# # Read Data

# +
with open(static_irox_structures_path, "rb") as fle:
    df_static_irox = pickle.load(fle)

with open(bulk_dft_data_path, "rb") as fle:
    df_bulk_dft = pickle.load(fle)
    df_bulk_dft = df_bulk_dft[df_bulk_dft.source == "raul"]

df_ids = pd.read_csv(unique_ids_path)

# +
df_static_irox = df_static_irox[df_static_irox.stoich == "AB3"]

df_bulk_dft = df_bulk_dft[df_bulk_dft.stoich == "AB3"]

# +
path_i = os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling/00_ml_workflow/191102_new_workflow/00_abx_al_runs/out_data",
    "duplicates.pickle")

import pickle; import os
with open(path_i, "rb") as fle:
    duplicates = pickle.load(fle)

all_duplicates = duplicates["AB2"] + duplicates["AB3"]

# +
ordered_dft_indices = df_bulk_dft.sort_values("dH").index

ordered_wo_dupl = []
for id_i in ordered_dft_indices:

    # if id_i not in ordered_wo_dupl:

    ordered_wo_dupl.append(id_i)

len(ordered_wo_dupl)

# +
df_static_irox = df_static_irox[df_static_irox.index.isin(ordered_wo_dupl)]

df_static_irox = df_static_irox.reindex(ordered_wo_dupl)
# -

df_static_irox = df_static_irox.iloc[0:20]

# +
tmp_list = []

data_list = []
for i_cnt, (id_i, row_i) in enumerate(df_static_irox.iterrows()):
    print(i_cnt)
    
    data_row_i = dict()

    pre_id = row_i.static_id
    post_id = row_i.name

    data_row_i["pre_id"] = pre_id
    data_row_i["post_id"] = post_id

    pre_atoms = row_i.atoms
    post_atoms = df_bulk_dft.loc[post_id].atoms

    pre_atoms.write("out_data/pre_post_structures/" + str(i_cnt).zfill(2) + "_" + id_i + "_pre.cif")
    post_atoms.write("out_data/pre_post_structures/" + str(i_cnt).zfill(2) + "_" + id_i + "_post.cif")

    try:
        d_ij = df_dij.loc[pre_id, post_id]
        tmp_list.append(d_ij)

        data_row_i["dij"] = d_ij
    except:
        pass
    
    data_list.append(data_row_i)
# -

{
    '8p8evt9pcg': 0.381740,
    'macixavwv3': 0.128178,
    'xwvhnh9sx4': 0.461092,
    '9lmkmh8s8r': 0.311769,
    '9txdvicqcf': 0.354713,
    '8k7expx2bp': 0.282011,
    'vp7gxqv191': 0.348111,
    'xg6exl6rmp': 0.496800,
    }

# +
df = pd.DataFrame(data_list)



# +
df_static_irox.loc["8p8evt9pcg"]

df_dij.loc[
    "8p8evt9pcg",
    "poboleni_97",   
    ]

# +
import chart_studio.plotly as py
import plotly.graph_objs as go

trace = go.Scatter(

    # df.dij
    # df.post_id

    # x=x_array,
    y=df.dij,
    mode="markers",
    text=df.post_id,
    hovertext=df.post_id,

    marker=dict(
        symbol="circle",
        color='blue',
        # colorscale='Viridis',
        colorbar=dict(thickness=20),
        size=12,
        line=dict(
            color='black',
            width=2
            )
        ),

    )

data = [trace]

fig = go.Figure(data=data)
fig.show()

# +
# pre_atoms.write("out_data/tmp_pre.cif")
# post_atoms.write("out_data/tmp_post.cif")

# +
# df_static_irox = df_static_irox.loc[["8p8evt9pcg"]]

# df_static_irox

# assert False

# +
# post_atoms

# row_i

# df_bulk_dft.loc[post_id]
