# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
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

# # Energy vs Volume for all computed IrOx bulk polymorphs
# ---

# # Import Modules

# +
import os
import sys
import pickle

import pandas as pd

import chart_studio.plotly as py
import plotly.graph_objs as go

# #############################################################################
from ase_modules.ase_methods import view_in_vesta

sys.path.insert(0, os.path.join(os.environ["PROJ_irox"], "data"))
from proj_data_irox import irox_bulk_color_map
from proj_data_irox import bulk_dft_data_path, oqmd_irox_data_path

from plotting.my_plotly import my_plotly_plot
# -

# # Read Data

# +
with open(bulk_dft_data_path, "rb") as fle:
    df_dft_calcs = pickle.load(fle)

path_i = os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling/manually_classify_dft_structures",
    "out_data",
    "df_coord_motiff.pickle")
with open(path_i, "rb") as fle:
    df_coord = pickle.load(fle)
# -

# # Clean up DFT dataframe

df = df_dft_calcs
df = df[
    (df["source"] != "chris") &
    (df["source"] != "oqmd") &
    [True for i in range(len(df))]
    ]
df_dft_calcs = df

df_dft_calcs.loc["7f7svsnpvg"]
df_coord.loc["7f7svsnpvg"]

# +
print("df_coord.shape:", df_coord.shape)
print("df_m.shape:", df_dft_calcs.shape)

df_tmp = pd.concat([
    df_dft_calcs,
    df_coord,
    ], axis=1, sort=True)
df_m = df_tmp
# -

# # Calculating volume per atom

# +
# df_m.loc["atoms"].isna

df_m = df_m[~df_m["atoms"].isna()]



# +
def method(row_i):
    atoms_i = row_i["atoms"]

    volume_i = atoms_i.get_volume()
    num_atoms = atoms_i.get_number_of_atoms()

    norm_volume = volume_i / num_atoms
    return(norm_volume)

df_m["volume_pa"] = df_m.apply(method, axis=1)


# +
def method(row_i):
    new_column_values_dict = {
        "color": None,
        "marker_size": None,
        }

    color = "black"
    marker_size = 14
    major_motiffs = row_i["major_motiffs"]
    
    if pd.isna(major_motiffs):
        major_motiffs = ""

    if "octahedral" in major_motiffs:
        color="pink"
    elif "bcc" in major_motiffs:
        color="green"
    elif "cp" in major_motiffs:
        color="yellow"
    elif "square pyramidal" in major_motiffs:
        color="blue"
    elif "tetrahedral" in major_motiffs:
        color="grey"
    elif "trigonal bipyramidal" in major_motiffs:
        color="brown"

    # elif "unrecognized" in major_motiffs:
    elif major_motiffs == "unrecognized":
        color="red"
        marker_size = 4
    else:
        print("Not processed!!!!")
        print(row_i.name)
    new_column_values_dict["color"] = color
    new_column_values_dict["marker_size"] = marker_size

    # #########################################################################
    for key, value in new_column_values_dict.items():
        row_i[key] = value
    return(row_i)


df_m = df_m.apply(method, axis=1)
# -

df_m = df_m.sort_values("source", ascending=True)

# # Plotting

# +
df_ab2 = df_m[df_m["stoich"] == "AB2"]
df_ab3 = df_m[df_m["stoich"] == "AB3"]

df_i = df_ab2
trace = go.Scatter(
    x=df_i["volume_pa"],
    y=df_i["energy_pa"],
    mode="markers",
    name="IrO2",
    marker=dict(
        symbol="circle",
        color=df_i["color"],
        size=df_i["marker_size"],
        opacity=0.7,
        line=dict(
            color='black',
#             color=df_m["color"],
            width=1.
        )
    ),

    hovertext=list(df_i.index),
    hoverinfo="text",
    )


df_i = df_ab3
trace_1 = go.Scatter(
    x=df_i["volume_pa"],
    y=df_i["energy_pa"],
    mode="markers",
    name="IrO3",
    marker=dict(
        symbol="x",
        color=df_i["color"],
        # size=14,
        size=df_i["marker_size"],
        opacity=0.7,
        line=dict(
            color='black',
#             color=df_m["color"],
            width=1.
        )
    ),

    hovertext=list(df_i.index),
    hoverinfo="text",
    )


layout = go.Layout(
    xaxis={"title": "Volume (A3 / atom)"},
    yaxis={"title": "Energy (eV / atom)"},
    showlegend=True,
    margin=go.layout.Margin(
        autoexpand=None,
        b=10,
        l=None,
        pad=None,
        r=None,
        t=10,
        ),

#     margin
    )
data = [trace, trace_1]

fig = go.Figure(data=data, layout=layout)
# -

fig.show()

my_plotly_plot(
    figure=fig,
    plot_name="E_vs_V_irox",
    # write_pdf_svg=True,
    upload_plot=False,
    write_html=True,
    )

# + active=""
# "octahedral": "pink"
# "bcc": "green"
# "cp": "yellow"
# "square pyramidal": "blue"
# "tetrahedral": "grey"
# "trigonal bipyramidal": "brown"
# "unrecognized": "red"
# -

# # Viewing structure(s)

# +
index_frag = "xsnsnozgxq"

index_search = [i for i in df_m.index if index_frag in i]
print(index_search)

index_tmp = index_search[0]
# view_in_vesta(df_m.loc[index_tmp]["atoms"], name_list=[index_tmp], ase_gui=True)

print(df_coord.loc[index_tmp]["coord_motiff_distr"])
print(df_coord.loc[index_tmp]["major_motiffs"])

df_m.loc[index_tmp]
# -

# # Oxy Coordination

# +
path_i = os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling/manually_classify_dft_structures",
    "out_data",
#     "df_oxy_coord.pickle",
    "df_coord_analysis.pickle",
    )

with open(path_i, "rb") as fle:
    df_oxy_coord = pickle.load(fle)

# +
print("df_coord.shape:", df_coord.shape)
print("df_m.shape:", df_m.shape)

df_tmp = pd.concat([
    df_m,
    df_oxy_coord,
    ], axis=1, sort=True)
df_m = df_tmp
# -

df_oxy_coord.head()


# +
def method(row_i):
    new_column_values_dict = {
        "color": None,
        "marker_size": None}

    color = "black"
    marker_size = 14
    oxy_coord_i = row_i["O_coord"]

    color = oxy_coord_i
    new_column_values_dict["color"] = color

    new_column_values_dict["marker_size"] = marker_size

    # #########################################################################
    for key, value in new_column_values_dict.items():
        row_i[key] = value
    return(row_i)

df_m = df_m.apply(method, axis=1)

# +
# df_m = df_m[~df_m["major_motiffs"].isna()]
# df_m = df_m[df_m["major_motiffs"].str.contains("oct")]

df_ab2 = df_m[df_m["stoich"] == "AB2"]
df_ab3 = df_m[df_m["stoich"] == "AB3"]

def get_scatter_trace(df_m=None, name="TEMP", marker_symbol="circle", color_col="oxy_coord"):
    """
    """    
    df_i = df_m
    trace = go.Scatter(
        x=df_i["volume_pa"],
        y=df_i["energy_pa"],
        mode="markers",
        # name="IrO2",
        name=name,
        marker=dict(
            symbol=marker_symbol,
            # color=df_i["color"],
            color=df_i[color_col],

            colorscale="Viridis",
            # cmin=1,
            # cmax=3,
            colorbar=dict(
                # title="Colorbar",
                x=1.02, y=0.25,
                ),
            size=df_i["marker_size"],
            opacity=0.7,
            line=dict(
                color='black',
                width=1.
            )
        ),
        hovertext=list(df_i.index),
        hoverinfo="text",
        )
    return(trace)


trace_0 = get_scatter_trace(df_m=df_ab2, name="IrO2", marker_symbol="circle",
    color_col="Ir_coord"
#     color_col="oxy_coord"
    )
trace_1 = get_scatter_trace(df_m=df_ab3, name="IrO3", marker_symbol="x",
    color_col="Ir_coord",
#     color_col="oxy_coord"
    )


layout = go.Layout(
    xaxis={"title": "Volume per Atom (A3 / atom)"},
    yaxis={"title": "DFT Energy per Atom (eV / atom)"},
    showlegend=True,
    )

data = [
    trace_0,
    trace_1,
    ]

fig = go.Figure(data=data, layout=layout)

fig

# +
# df_m.head()

# df_m["oxy_coord"].max()
# df_m["oxy_coord"].min()

# +
df_tmp = df_m.sort_values("O_coord", ascending=False).iloc[0:10]

# ["atoms"].tolist()
atoms_list = df_tmp["atoms"].tolist()
names= df_tmp.index.tolist()
# view_in_vesta(atoms_list, name_list=names, ase_gui=False)
# -

df_m.loc["9hc58scr8f"]

# + active=""
#
#
#
#
#

# + jupyter={}
# BLUE: IrO3
# RED: IrO2
# Green: My DFT IrO2

# df_m["atoms"]

# for i_cnt, row_i in df_m.iterrows():
#     row_i["atoms"].write("out_data/cif_files/" + row_i.name + ".cif")
#     # row_i["atoms"].write(row_i.name + ".cif")

#     # break

# for index_i in df_m.index.tolist():
#     if "bdc" in index_i:
#         print(index_i)

# atoms = df_m.loc["bdctzwcg8h"]["atoms"]
# # view_in_vesta(atoms, ase_gui=True, name_list=None)

# # df_m[df_m["source"] == "raul"]
# # df_m = df_m.sort_values("source", ascending=False)

# # df_m

# from ase.visualize import view
# view(df_m.loc["9yz2mt8hbh"]["atoms"])

# df_m.loc["6s648e8s6p"].iloc[1]["path"]

# Removing rows that don't have an atoms object
# df_dft_calcs = df_dft_calcs[df_dft_calcs["atoms"].notnull()]

# df_m = df_dft_calcs

# df_m = pd.concat([
#     df_dft_calcs,
# #     df_oqmd_calcs,
#     ])

# df_m = df_m[df_m["stoich"] == "AB3"]
# df_m = df_m[df_m["source"] != "oqmd"]

# def method(row_i, argument_0, optional_arg=None):
#     new_column_values_dict = {
#         "marker_size": 5,
#         "marker_line_color": "rgb(0,0,0)",
#         "marker_line_size": 0.05,
#         }

#     computed_bool = row_i.get("computed", False)
#     if computed_bool:
#         new_column_values_dict["marker_size"] = 10
#         new_column_values_dict["marker_line_color"] = "red"
#         new_column_values_dict["marker_line_size"] = 1.5

#     # #########################################################################
#     for key, value in new_column_values_dict.items():
#         row_i[key] = value
#     return(row_i)


# arg1 = "TEMP_0"
# # df_i["column_name"] = df_i.apply(
# df_i = df_i.apply(
#     method,
#     axis=1,
#     args=(arg1, ),
#     optional_arg="TEMP_1"
#     )

# def method(row_i):
#     stoich_i = row_i["stoich"]
#     source_i = row_i["source"]

#     if stoich_i == "AB2":
#         color = "red"
#     elif stoich_i == "AB3":
#         color = "blue"
#     if source_i == "raul":
#         color = "green"
#     elif source_i == "oqmd":
#         color = "black"
#     if source_i == "raul" and stoich_i == "AB3":
#         color = "orange"
#     if source_i == "raul" and stoich_i == "AB2":
#         color = "grey"
#     if source_i == "chris" and stoich_i == "AB2":
#         color = "pink"
#     if source_i == "chris" and stoich_i == "AB3":
#         color = "brown"

#     return(color)
# # df_m["color"] = df_m.apply(method, axis=1)


# + jupyter={}
# vhv3
# vhv39q6e9j
