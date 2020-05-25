# -*- coding: utf-8 -*-
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

# + [markdown] Collapsed="false"
# # Creating E vs V Figure for IrO2 and IrO3
# ---

# + [markdown] Collapsed="false"
# # Import Modules

# + Collapsed="false" attributes={"classes": [], "id": "", "n": "1"}
import os
print(os.getcwd())
import sys

import copy
# import pickle

import numpy as np
import pandas as pd

# import ase
from ase.db import connect

from plotly.subplots import make_subplots
import chart_studio.plotly as py
import plotly.graph_objs as go
import plotly.express as px

# #########################################################
from layout import layout
# -

from inputs import structure_id_map

dx = 0.2

# + [markdown] Collapsed="false"
# # Read Data

# +
sys.path.insert(0, os.path.join(
    os.environ["PROJ_irox"], "workflow/ml_modelling"))
from ml_methods import get_ml_dataframes

DF_dict = get_ml_dataframes(names=[
    "df_dft_final_final_path",
    # "",
    ])

df_bulk_dft = DF_dict["df_dft_final_final"]

# + [markdown] Collapsed="false"
# # Construct DataFrame
#
#


# + Collapsed="false" attributes={"classes": [], "id": "", "n": "4"}
# #############################################################################
# Structural Analysis db file
FinalStructuresdb_file = os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling/energy_vs_volume/kirsten_E_vs_V_analysis/scripts",
    "out_data/FinalStructures_1.db")

db = connect(FinalStructuresdb_file)


data_list = []
for row in db.select():
    row_dict = dict(
        energy=row.get("energy"),
        # volume=row.get("volume"),
        **row.key_value_pairs,
        )
    data_list.append(row_dict)

df = pd.DataFrame(data_list)

df = df[~df["stoich"].isna()]
#df = df[~df["coor_env"].isna()]

print("Total df rows:", df.structure_id.shape[0])
print("Unique structure ids:", df.structure_id.unique().shape[0])
print("")

# +
df = df.set_index("structure_id")

df = df.loc[
    df.index.intersection(df_bulk_dft.index)    
    ]

# + Collapsed="false"
# #############################################################################
# Merge dataframes together ###################################################

# Drop unnecessary duplicate columns before merging
df = df.drop([
    "energy",
    # "volume",
    "stoich",
    "id_old",
    ], axis=1)

df = pd.merge(df, df_bulk_dft,
    left_index=True,
    right_index=True,
    )

print("df.shape:", df.shape)
print("df_bulk_dft.shape:", df_bulk_dft.shape)

# + [markdown] Collapsed="false"
# # Process Dataframe

# +
sys.path.insert(0, ".")
from colors import get_color_scale

colorscale_i = get_color_scale(df=df, dx=dx)

# + [markdown] Collapsed="false"
# # Sorting data to bring out 4/6 coordination

# + Collapsed="false"
df = df.sort_values("mean_coor")

df_concat_list = [
    df[
        (df.mean_coor < 4 + dx) & \
        (df.mean_coor > 4 - dx)
        ],

    df[
        (df.mean_coor < 6 + dx) & \
        (df.mean_coor > 6 - dx)
        ],
    
    ]

df_tmp = pd.concat(df_concat_list)

remaining_ids = [i for i in df.index if i not in df_tmp.index]
df = pd.concat([df_tmp, df.loc[remaining_ids]])


df = df.reindex(index=df.index[::-1])

print('Total IrO2:', len(df[df.stoich == "AB2"]['dH']))
print('Total IrO3:', len(df[df.stoich == "AB3"]['dH']))     

print('Metastable IrO2:', len(np.where(df[df.stoich == "AB2"]['dH'].values < -0.33)[0]))
print('Metastable IrO3:', len(np.where(df[df.stoich == "AB3"]['dH'].values < -0.34)[0]))


# + [markdown] Collapsed="false" Collapsed="false" Collapsed="false" Collapsed="false" Collapsed="false" Collapsed="false" toc-hr-collapsed=true
# # Plotting

# + [markdown] Collapsed="false"
# ## Shared scatter attributes

# + Collapsed="false" attributes={"classes": [], "id": "", "n": "6"} jupyter={}
scatter_shared = go.Scatter(
    mode="markers",
    hoverinfo="text",
    marker=dict(
        symbol="circle",
        size=4,
        opacity=0.8,
        line=dict(
            color="black",
            # width=1,
            width=0.,
            ),
        colorscale=colorscale_i,

        colorbar=dict(
            bordercolor="green",
            outlinecolor="black",
            tickcolor="black",
            xanchor="right",
            # x=1.091,
            # x=1.1,
            x=1.15,
            len=1.16,
            lenmode="fraction",
            # #################################################################
            thickness=15,
            thicknessmode=None,
            tickprefix=None,
            ticks="outside",
            # #################################################################
            tickvals = [2, 4, 6, 8, 10, 12],
            y=0.50005,
            yanchor="middle",
            ypad=10,
            borderwidth=None,

            title=go.scatter.marker.colorbar.Title(
                font=None,
                side="right",  # ['right', 'top', 'bottom']
                text="Ir-O Coord. Num.",
                ),

            # titlefont=None,
            # titleside=None,

            ),

        ),
    )

# + [markdown] Collapsed="false"
# ## Create AB2/3 traces

# + Collapsed="false" attributes={"classes": [], "id": "", "n": "11"}
# %%capture

df_i = df[df.stoich == "AB2"]
trace_ab2 = go.Scatter(
    x=df_i.volume_pa,
    y=df_i.dH,
    # text=[str(i) for i in df_i.mean_coor.tolist()],
    text=df_i.index.values,
    marker=dict(color=df_i.mean_coor, size=3))
trace_ab2.update(**scatter_shared.to_plotly_json())



# #############################################################################
df_i = df[df.stoich == "AB3"]
trace_ab3 = go.Scatter(
    x=df_i.volume_pa,
    y=df_i.dH,
    # text=[str(i) for i in df_i.mean_coor.tolist()],
    text=df_i.index.values,
    marker=dict(color=df_i.mean_coor))
trace_ab3.update(**scatter_shared.to_plotly_json())

# + [markdown] Collapsed="false"
# # Shapes

# +
from shapes import get_plot_shapes

inset_range_0_x = [9.5, 17.]
inset_range_1_x = [9.5, 17.5]


out_dict = get_plot_shapes(
    df=df,
    inset_range_0_x=inset_range_0_x,
    inset_range_1_x=inset_range_1_x,
    )

shapes_list = out_dict["shapes_list"]

shape_inset_metastability_ab2 = go.layout.Shape(
    type="line",
    x0=0,
    y0=-0.33285956787756277, #ab2_min_e + metastability_limit,
    x1=40,
    y1=-0.33285956787756277,
    xref="x1",
    yref="y1",
    line=dict(
        color="grey",
        width=1,
    )
)

shape_inset_metastability_ab3 = go.layout.Shape(
    type="line",
    x0=0,
    y0=-0.3438547784081729, 
    x1=40,
    y1=-0.3438547784081729, 
    xref="x2",
    yref="y2",
    line=dict(
        color="grey",
        width=1,
    )
)

shapes_list += [shape_inset_metastability_ab2, shape_inset_metastability_ab3]
inset_range_0_y = out_dict["inset_range_0_y"]
inset_range_1_y = out_dict["inset_range_1_y"]

# + [markdown] Collapsed="false"
# ## Create subplot

# + Collapsed="false"
# %%capture

inset_attr = dict(l=0.5, b=0.5)
fig = make_subplots(
    rows=1, cols=2,
    shared_xaxes=True,
    shared_yaxes=True,
    specs = [[{}, {}]],
    insets=[
        {
            "cell": (1,1),
            **inset_attr,
            },

        {
            "cell": (1,2),
            **inset_attr,
            },
        ],
    horizontal_spacing=0.04)


# #########################################################
# Add traces ##############################################
fig.add_trace(trace_ab2, row=1, col=1)
fig.add_trace(trace_ab3, row=1, col=2)

fig.add_trace(copy.deepcopy(trace_ab2).update(xaxis="x3", yaxis="y3"))
fig.add_trace(copy.deepcopy(trace_ab3).update(xaxis="x4", yaxis="y4"))

# #########################################################
for shape_i in shapes_list:
    fig.add_shape(shape_i)
# -

# # Layout Properties

# + Collapsed="false" attributes={"classes": [], "id": "", "n": "8"}
# #############################################################################
# Update Layout ###############################################################
fig.update_layout(layout)
fig.update_xaxes(layout["xaxis"])
fig.update_yaxes(layout["yaxis"])

fig.layout.yaxis2.title = None

# Modifying inset props
fig.layout.xaxis3.title = None
fig.layout.yaxis3.title = None

fig.layout.xaxis4.title = None
fig.layout.yaxis4.title = None



fig.layout.xaxis3.range = inset_range_0_x
fig.layout.yaxis3.range = inset_range_0_y


fig.layout.xaxis4.range = inset_range_1_x
fig.layout.yaxis4.range = inset_range_1_y


# fig.layout.xaxis3
fig.layout.xaxis3.tickfont.size = 6 * (4 / 3)
fig.layout.yaxis3.tickfont.size = 6 * (4 / 3)

fig.layout.xaxis4.tickfont.size = 6 * (4 / 3)
fig.layout.yaxis4.tickfont.size = 6 * (4 / 3)

fig.layout.xaxis.dtick = 5
fig.layout.yaxis.dtick = 0.5

fig.layout.xaxis2.dtick = 5
fig.layout.yaxis2.dtick = 0.5

fig.layout.xaxis3.dtick = 2

fig.layout.xaxis4.dtick = 2

fig.layout.xaxis3.ticklen = 3
fig.layout.xaxis4.ticklen = 3

# COMBAK
fig.layout.yaxis3.tickmode = "array"
fig.layout.yaxis3.tickvals = [-0.8, -0.7, -0.6, -0.5, -0.4]
fig.layout.yaxis3.ticklen = 3

fig.layout.yaxis4.tickmode = "array"
fig.layout.yaxis4.tickvals = [-0.65, -0.6, -0.55, -0.5]#[-0.7, -0.6, -0.5, -0.4, -0.3]
fig.layout.yaxis4.ticklen = 3
# -

# # Annotations

# +
# %%capture

annotations=[

    #| - IrO2/3 Annotation
    go.layout.Annotation(
        x=9.4,
        y=1.76,
        xref="x",
        yref="y",
        text="IrO<sub>2</sub>",
        showarrow=False,

        bgcolor="rgba(255,255,255,0.7)",
        font=go.layout.annotation.Font(
            color="black",
            family=None,
            size=10 * (4/3),
            ),

        ax=0,
        ay=0,
        ),


    go.layout.Annotation(
        x=9.4,
        y=1.76,
        xref="x2",
        yref="y2",
        text="IrO<sub>3</sub>",
        showarrow=False,

        bgcolor="rgba(255,255,255,0.7)",
        font=go.layout.annotation.Font(
            color="black",
            family=None,
            size=10 * (4/3),
            ),

        ax=0,
        ay=0,
        ),

    ]
    #__|

# +
# %%capture

for id_i, val in structure_id_map.items():

    try:
        df_i = df.loc[id_i]
    except:
        print(id_i, 'not found')
        continue

    y = df_i.dH
    x = df_i.volume_pa

    if df_i.stoich == 'AB2':
        if y < -0.55:
            sub_x = 'x3'
            sub_y = 'y3'
        else:
            sub_x = 'x'
            sub_y = 'y'
    elif df_i.stoich == 'AB3':
        if y < -0.4 and x < 16.5:
            sub_x = 'x4'
            sub_y = 'y4'
        else:
            sub_x = 'x2'
            sub_y = 'y2'

    # #####################################################
    # Arrow shift
    # arrowshift = 0
    # arrowshifty = 0

    # if len(val) > 8:
    #     arrowshift = len(val) * 2.5
    # elif len(val) > 4:
    #     arrowshift = len(val) * 3
    # # elif '(' in val:
    # #     arrowshift = 15
    # else:
    #     arrowshift = 10

    # # if val =='(2)':
    # #     arrowshift *= -1
    # if val == 'iii (pyrite)':
    #     arrowshifty = -4

    from inputs import annot_offset_dict

    ax = 0
    ay = 0

    if id_i in annot_offset_dict.keys():
        annot_dict_i = annot_offset_dict[id_i]

        ax = annot_dict_i["ax"]
        ay = annot_dict_i["ay"]
        xanchor = annot_dict_i.get("xanchor", None)
        
    annot_font = go.layout.annotation.Font(
        color="black",
        family=None,
        size=6 * (4/3),
        )

    # #####################################################
    annot_i = go.layout.Annotation(
        x=x, y=y,
        xref=sub_x, yref=sub_y,
        text=val,
        showarrow=True,
        arrowhead=1,
        # startstandoff=10,
        standoff=1,
        font=annot_font,
        # font=go.layout.annotation.Font(
        #     color="black",
        #     family=None,
        #     size=6 * (4/3),
        #     ),

        xanchor=xanchor,
        ax=ax, ay=ay,
        )
    annotations.append(annot_i)


fig.layout.update(annotations=annotations)

# +
# # go.layout.Annotation?
# -

# ## Write/display plot

# +
from plotting.my_plotly import my_plotly_plot

my_plotly_plot(
    figure=fig,
    plot_name="E_vs_V_plot_3",
    write_html=True,
    write_png=False,
    png_scale=6.0,
    write_pdf=True,
    write_svg=False,
    try_orca_write=False,
    )
# -

fig.show()

assert False

# + [markdown] Collapsed="false"
# # Histogram Plot

# + Collapsed="false"
# import plotly.express as px

fig = px.histogram(
    df,
    x="mean_coor",
    color="stoich",
    marginal="rug",  # can be `box`, `violin`
    opacity=0.9,
    nbins=100,
    # barnorm="fraction",
    histnorm="percent",
    # hover_data=tips.columns,
    )

fig.show()

# + Collapsed="false" active=""
#
#
#
#
# -

df[df.stoich == "AB3"].sort_values("dH").iloc[0:8].index.tolist()

# + Collapsed="false" jupyter={}
# # %%capture

# # #############################################################################
# # Duplicates list
# path_i = os.path.join(
#     os.environ["PROJ_irox"],
#     "workflow/ml_modelling/00_ml_workflow/191102_new_workflow/00_abx_al_runs",
#     "out_data/duplicates.pickle")
# duplicates = pickle.load(open(path_i, "rb"))


# # #############################################################################
# # Bulk DFT Dataframe
# sys.path.insert(0, os.path.join(
#     os.environ["PROJ_irox"], "workflow/ml_modelling"))
# from ml_methods import get_data_for_al

# data_dict = get_data_for_al(stoich="AB2", drop_too_many_atoms=True)
# df_bulk_dft_ab2 = data_dict["df_bulk_dft"]

# data_dict = get_data_for_al(stoich="AB3", drop_too_many_atoms=True)
# df_bulk_dft_ab3 = data_dict["df_bulk_dft"]

# # Combine AB2/3 Dataframes
# df_bulk_dft = pd.concat([df_bulk_dft_ab2, df_bulk_dft_ab3])
# df_bulk_dft = df_bulk_dft[df_bulk_dft.source == "raul"]

