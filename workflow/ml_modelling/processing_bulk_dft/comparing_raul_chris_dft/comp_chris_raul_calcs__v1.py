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

# # Import Modules

# +
import os
import sys

import copy

import pickle
import pandas as pd

from ase.visualize import view

import chart_studio.plotly as py
import plotly.graph_objs as go
# #############################################################################
sys.path.insert(0, os.path.join(os.environ["PROJ_irox"], "data"))
from proj_data_irox import bulk_dft_data_path
# -

# # Script Inputs

raul_color = "red"
chris_color = "black"

# # Read & Process DFT Data

# ## Read dataframe

with open(bulk_dft_data_path, "rb") as fle:
    df_bulk_dft = pickle.load(fle)
    df = df_bulk_dft

# ## Filter only systems shared by Chris and Raul DFT calcs

shared_indices = list(
    set(list(df[df["source"] == "raul"].index))
    &
    set(list(df[df["source"] == "chris"].index))
    )
print("Number of shared calculations: ", len(shared_indices))
df_shared = df.loc[shared_indices]

# + active=""
# 75 shared structures last time
# -

# # Comparison dataframe

# +
grouped = df_shared.reset_index().groupby(["id_unique"])
data_list = []
for name, group in grouped:

    data_dict_i = {
        "id_unique": name,
        "e_diff": None,
        "v_diff": None,
        "e_diff_abs": None,
        "v_diff_abs": None,
        }

    row_raul = group[group["source"] == "raul"]
    row_chris = group[group["source"] == "chris"]


    e_diff = row_raul["energy_pa"].iloc[0] - row_chris["energy_pa"].iloc[0]
    e_diff_abs = abs(e_diff)

    atoms_raul = row_raul.iloc[0]["atoms"]
    atoms_chris = row_chris.iloc[0]["atoms"]

    vol_pa_raul = atoms_raul.get_volume() / atoms_raul.get_number_of_atoms()
    vol_pa_chris = atoms_chris.get_volume() / atoms_chris.get_number_of_atoms()

    v_diff = vol_pa_raul - vol_pa_chris
    v_diff_abs = abs(v_diff)


    data_dict_i["e_diff"] = e_diff
    data_dict_i["v_diff"] = v_diff

    data_dict_i["v_diff_abs"] = v_diff_abs
    data_dict_i["e_diff_abs"] = e_diff_abs


    data_list.append(data_dict_i)

df_comp = pd.DataFrame(data_list)

df_comp = df_comp.set_index("id_unique")


# -

# # Creating main traces

# +
def method(row_i, argument_0, raul_color=None, chris_color=None, df_comp=None):

    # #########################################################################
    id_unique = row_i.name
    source_i = row_i["source"]
    atoms_i = row_i["atoms"]
    energy_pa = row_i["energy_pa"]

    # #########################################################################
    e_diff = df_comp.loc[id_unique]["e_diff"]

    # #########################################################################
    new_column_values_dict = {}

    # #########################################################################
    if source_i == "raul":
        color_0 = raul_color
        showlegend_i = True
        marker_shape = "diamond"
        marker_line_color = "black"

    elif source_i == "chris":
        color_0 = chris_color
        showlegend_i = False
        marker_shape = "circle"
        marker_line_color = "black"

    volume_pa = atoms_i.get_volume() / atoms_i.get_number_of_atoms()

    # #########################################################################
    new_column_values_dict["volume_pa"] = volume_pa
    new_column_values_dict["marker_shape"] = marker_shape
    new_column_values_dict["e_diff"] = e_diff
    new_column_values_dict["marker_line_color"] = marker_line_color

    # #########################################################################
    for key, value in new_column_values_dict.items():
        row_i[key] = value

    return(row_i)

# #############################################################################
df_i = df_shared
df_i = df_i.apply(
    method,
    axis=1,
    args=("TEMP_0",),
    raul_color=raul_color,
    chris_color=chris_color,
    df_comp=df_comp,
    )
df_shared = df_i
# -

# # Creating connections between data

grouped = df_shared.reset_index().groupby(["id_unique"])
connecting_lines_data = []
for name, group in grouped:

    assert len(group) == 2, "JIDFJSDI"

    row_0 = group.iloc[0]
    row_1 = group.iloc[1]

    vol_0 = row_0["atoms"].get_volume() / row_0["atoms"].get_number_of_atoms()
    vol_1 = row_1["atoms"].get_volume() / row_1["atoms"].get_number_of_atoms()


    energy_pa_0 = row_0["energy_pa"]
    energy_pa_1 = row_1["energy_pa"]

    trace_i = go.Scatter(
        x=[vol_0, vol_1],
        y=[energy_pa_0, energy_pa_1],
        mode="lines",
        legendgroup=name,
        showlegend=False,
        line=dict(
            color="grey",
            width=1.,
#             dash="dot",
            ),

        )

    connecting_lines_data.append(trace_i)

# +
e_diff_col = df_comp["e_diff"]
col_scale_0 = (0 - e_diff_col.min()) / (e_diff_col.max() - e_diff_col.min())

max_abs_e_diff = max(
    abs(e_diff_col.min()),
    abs(e_diff_col.min()),
    )

# +
extra_row = copy.deepcopy(df_shared.iloc[0])
extra_row["e_diff"] = max_abs_e_diff
extra_row.name = "TEMP0"
df_shared = df_shared.append(extra_row)

extra_row = copy.deepcopy(df_shared.iloc[0])
extra_row["e_diff"] = -max_abs_e_diff
extra_row.name = "TEMP1"
df_shared = df_shared.append(extra_row)
# -

ev_0p1 = 0.1 / (2 * max_abs_e_diff)

# # Plotting

# +
colorscale_i = [
#     [0.000, 'rgba(214, 39, 40, 0.85)'],
    [0.000, "blue"],

    [0.5 - ev_0p1, "purple"],

#     [col_scale_0, 'rgba(255, 255, 255, 0.85)'],
    [0.5, 'rgba(255, 255, 255, 0.85)'],


    [0.5 + ev_0p1, "orange"],


#     [1.8, 'red'],
    [1.000, 'red'],
    ]


trace_tmp = go.Scatter(
    x=df_shared["volume_pa"],
    y=df_shared["energy_pa"],

    mode="markers",

#     legendgroup=id_unique,
#     showlegend=showlegend_i,
#     name=list(df_shared.index),
    text=list(df_shared.index),
    hoverinfo="text",
    marker=dict(
        symbol=df_shared["marker_shape"],
        color=df_shared["e_diff"],
        # colorscale='Viridis',
        # colorscale="RdGy",
        colorscale=colorscale_i,
        colorbar=dict(thickness=10),

        size=10,
        line=dict(
            color=df_shared["marker_line_color"],
            width=1.,
            )
        ),
    )

# +
data = [] + \
    connecting_lines_data + \
    [trace_tmp] + \
    []
    # df_shared["trace"].tolist() + \


fig = go.Figure(data=data)
fig.update_layout(showlegend=False)

fig.show()
# -

# # Energy per atom vs Volume per atom

# + active=""
#
#
#
#
#
# -

# # Jobs for which Raul has higher energy than Chris's calcs to within a tolerance

# +
df_comp.sort_values("e_diff_abs", ascending=False)

df_comp[
    (df_comp["e_diff"] > 0.) &
    (df_comp["e_diff"] > 0.01)
    ].sort_values("e_diff")
# -

# # Am I missing any systems that Chris calculated?

# + active=""
# Most of these are coming soon, a few weren't being run but I'm running them now
#
# 8jvfcyvk92
# xy6kzjninu
# xg6exl6rmp
# x5nlvgnaxj
# mj7wbfb5nt
# vlxp9abd6h
# 6qvlcl6iv2
# npbq9ynjn2
# mtclmozw8t
# brnhny8y92
# zwnung6s71
# cgxkbicgz4
# 8ivkxwnhva

# +
file_path_i = os.path.join(
    os.environ["PROJ_irox"],
    "data/ml_irox_data/iro2_training_data.csv")
train_data_iro2 = pd.read_csv(file_path_i)
train_data_iro2.set_index("id", inplace=True)

file_path_i = os.path.join(
    os.environ["PROJ_irox"],
    "data/ml_irox_data/iro3_training_data.csv")
train_data_iro3 = pd.read_csv(file_path_i)
train_data_iro3.set_index("id", inplace=True)

train_data_dict = {
    "iro2": train_data_iro2,
    "iro3": train_data_iro3}

# #############################################################################
ids_iro2 = train_data_iro2[train_data_iro2["source"] == "chris"]["id_unique"]
ids_iro3 = train_data_iro3[train_data_iro3["source"] == "chris"]["id_unique"]
chris_computed_irox_ids = ids_iro2.tolist() + ids_iro3.tolist()

df = df_bulk_dft
df_i = df[
    (df["source"] == "raul") & \
    (df["stoich"] == "AB3")
    ]

# #############################################################################
for chris_id in ids_iro3.tolist():
    if chris_id not in list(df_i.index):
        print(chris_id)
