# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.4.2
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

import random

import pandas as pd

import chart_studio.plotly as py
import plotly.graph_objs as go
from plotly.subplots import make_subplots

# #########################################################
# Project Imports #########################################
sys.path.insert(0, os.path.join(
    os.environ["PROJ_irox"], "data"))
from proj_data_irox import axis_label_font_size, axis_tick_labels_font_size
from proj_data_irox import axis_label_font_size

# #########################################################
# Local Imports ###########################################
from layout import layout
from layout import xaxis_layout, yaxis_layout

# +
from al_data import main_AB2_run, main_AB3_run

main_AB2_run_name = main_AB2_run.split("/")[-1].split(".")[0].split("_")[-1]
main_AB3_run_name = main_AB3_run.split("/")[-1].split(".")[0].split("_")[-1]
# -

# # Script Inputs

# # Read Data

# +
# #########################################################
import pickle; import os
path_root = os.path.join(
    os.environ["PROJ_irox"],
    # "workflow/ml_modelling/00_ml_workflow/191102_new_workflow/performance_comp/top_10_disc_vs_dft",
    "workflow/ml_modelling/00_ml_workflow/performance_comp/top_10_disc_vs_dft",
    "out_data")

stoich_i = "AB2"
# #########################################################
path_i = os.path.join(path_root, stoich_i + "_df_random.pickle")
with open(path_i, "rb") as fle:
    df_random_ab2 = pickle.load(fle)
path_i = os.path.join(path_root, stoich_i + "_df_gbucb.pickle")
with open(path_i, "rb") as fle:
    df_gbucb_ab2 = pickle.load(fle)
# #########################################################

stoich_i = "AB3"
# #########################################################
path_i = os.path.join(path_root, stoich_i + "_df_random.pickle")
with open(path_i, "rb") as fle:
    df_random_ab3 = pickle.load(fle)
path_i = os.path.join(path_root, stoich_i + "_df_gbucb.pickle")
with open(path_i, "rb") as fle:
    df_gbucb_ab3 = pickle.load(fle)
# #########################################################
# -

# # Methods

# +
def process_df(df_i):
    df_i = df_i.dropna(axis=0)

    # #####################################################
    max_ind = df_i.index.max()
    new_ind = max_ind + 50

    new_row_dict = dict()
    for col in df_i.columns:
        new_row_dict[col] = 10
    new_row_series = pd.Series(new_row_dict, name=new_ind)

    df_i = df_i.append(
        new_row_series,
        # new_row_dict,
        ignore_index=False, sort=False)

    # #####################################################
    new_row_dict = dict()
    for col in df_i.columns:
        new_row_dict[col] = 0
    new_row_series = pd.Series(new_row_dict, name=0)

    df_i = df_i.append(
        new_row_series,
        ignore_index=False, sort=False)

    # #####################################################
    df_i = df_i.sort_index()


    # #####################################################
    data_dict_list = []
    for col in df_i.columns:
        data_dict_i = dict()
        data_dict_i["name"] = col[0]

        dft_to_reach_10 = df_i[col][df_i[col] == 10].index.values.min()
        data_dict_i["min_dft"] = dft_to_reach_10

        data_dict_list.append(data_dict_i)


    df_best_worst = pd.DataFrame(data_dict_list)
    df_best_worst = df_best_worst.set_index("name")

    df_worst = df_best_worst[df_best_worst.min_dft == df_best_worst.min_dft.max()]
    df_best = df_best_worst[df_best_worst.min_dft == df_best_worst.min_dft.min()]

    worst_name = df_worst.iloc[0].name
    best_name = df_best.iloc[0].name

    # #####################################################
    out_dict = dict()
    out_dict["df_i"] = df_i
    out_dict["worst_name"] = worst_name
    out_dict["best_name"] = best_name

    return(out_dict)
    # return(df_i)



def rand_color():
    # rgb_list = 3 * [random.randint(0, 255)]
    rgb_list = 3 * [random.randint(0, 200)]
    rgb_list = [str(i) for i in rgb_list]

    # color_i = "rgba(" + ",".join(rgb_list) + ",0.8)"
    color_i = "rgba(" + ",".join(rgb_list) + ",0.3)"

    return(color_i)


# -

# # Subplots

# +
sp_00 = dict()
sp_10 = dict()
sp_20 = dict()
sp_30 = dict()

fig = make_subplots(
    rows=2, cols=2,
    specs=[
        # [sp_00],
        # [sp_10],
        # [sp_20],
        # [sp_30],

        [sp_00, dict()],
        [sp_10, dict()],
        ],
    print_grid=True,

    column_titles=None,
    row_titles=None,

    # x_title="DFT Calculations",
    # y_title="N<sub>disc</sub>",

    # shared_yaxes=False,
    horizontal_spacing=2 * 0.01,
    vertical_spacing=2 * 0.015,
    # row_heights=[
    #     0.03,
    #     0.33 + 0.06 - 0.05 + 0.07,
    #     0.14 - 0.06 + 0.03,
    #     0.45,
    #     ],
    # column_widths=[1 / 6 - dx_tmp, 1 / 6 + dx_d6, 1 / 6 + dx_d6, 1 / 6 + dx_d6, 1 / 6 + dx_d6, 1 / 6 + dx_d6],
    )
# -

# # Create Traces

# +
df_i = df_gbucb_ab2

out_dict = process_df(df_i)
# ###########################
df_i = out_dict["df_i"]
worst_name = out_dict["worst_name"]
best_name = out_dict["best_name"]

data_i = []
data_i_ontop = []
for i_cnt, col_i in enumerate(df_i.columns):
    name_i = col_i[0]

    color_i = rand_color()

    special_names = [main_AB2_run_name, best_name, worst_name]

    if name_i == main_AB2_run_name:
        color_i = "rgba(255,0,0,0.8)"

    trace_i = go.Scatter(
        x=df_i[col_i].index,
        y=df_i[col_i].values,
        name=name_i,
        line_color=color_i,
        )
    
    if name_i in special_names:
        data_i_ontop.append(trace_i)
    else:
        data_i.append(trace_i)
        fig.add_trace(trace_i, row=1, col=1)

for trace_i in data_i_ontop:
    fig.add_trace(trace_i, row=1, col=1)

# +
df_i = df_random_ab2

out_dict = process_df(df_i)
# ###########################
df_i = out_dict["df_i"]
worst_name = out_dict["worst_name"]
best_name = out_dict["best_name"]

data_i = []
data_i_ontop = []
for i_cnt, col_i in enumerate(df_i.columns):
    name_i = col_i[0]

    color_i = rand_color()

    special_names = [main_AB2_run_name, best_name, worst_name]

    trace_i = go.Scatter(
        x=df_i[col_i].index,
        y=df_i[col_i].values,
        name=name_i,
        line_color=color_i,
        )
    # data_i.append(trace_i)
    # fig.add_trace(trace_i, row=2, col=1)

    if name_i in special_names:
        data_i_ontop.append(trace_i)
    else:
        data_i.append(trace_i)
        fig.add_trace(trace_i, row=2, col=1)

for trace_i in data_i_ontop:
    fig.add_trace(trace_i, row=2, col=1)

# +
df_i = df_gbucb_ab3


out_dict = process_df(df_i)
# ###########################
df_i = out_dict["df_i"]
worst_name = out_dict["worst_name"]
best_name = out_dict["best_name"]

data_i = []
data_i_ontop = []
for i_cnt, col_i in enumerate(df_i.columns):
    name_i = col_i[0]

    color_i = rand_color()

    special_names = [main_AB3_run_name, best_name, worst_name]

    if name_i == main_AB3_run_name:
        color_i = "rgba(255,0,0,0.8)"
    # if name_i == best_name:
    #     color_i = "rgba(110,243,88,0.8)"
    # if name_i == worst_name:
    #     color_i = "rgba(79,152,0,0.8)"

    trace_i = go.Scatter(
        x=df_i[col_i].index,
        y=df_i[col_i].values,
        name=name_i,
        line_color=color_i,
        )
    # data_i.append(trace_i)
    # fig.add_trace(trace_i, row=1, col=2)

    if name_i in special_names:
        data_i_ontop.append(trace_i)
    else:
        data_i.append(trace_i)
        fig.add_trace(trace_i, row=1, col=2)

for trace_i in data_i_ontop:
    fig.add_trace(trace_i, row=1, col=2)

# fig_tmp = go.Figure(data=data_i, layout=layout)
# fig_tmp.show()

# +
df_i = df_random_ab3

out_dict = process_df(df_i)
# ###########################
df_i = out_dict["df_i"]
worst_name = out_dict["worst_name"]
best_name = out_dict["best_name"]

data_i = []
data_i_ontop = []
for i_cnt, col_i in enumerate(df_i.columns):
    name_i = col_i[0]


    color_i = rand_color()

    special_names = [main_AB2_run_name, best_name, worst_name]

    trace_i = go.Scatter(
        x=df_i[col_i].index,
        y=df_i[col_i].values,
        name=name_i,
        line_color=color_i,
        )
    # data_i.append(trace_i)
    # fig.add_trace(trace_i, row=2, col=2)

    if name_i in special_names:
        data_i_ontop.append(trace_i)
    else:
        data_i.append(trace_i)
        fig.add_trace(trace_i, row=2, col=2)

for trace_i in data_i_ontop:
    fig.add_trace(trace_i, row=2, col=2)
# -

# # Set Master Layout Properties

layout.xaxis = None
layout.yaxis = None

# +
fig.update_xaxes(
    patch=xaxis_layout,
    selector=None,
    overwrite=False,
    )

fig.update_yaxes(
    patch=yaxis_layout,
    selector=None,
    overwrite=False,
    )

tmp = 42
# -

# # Update xaxis of left and right two plots

# +
xaxis_new = go.layout.XAxis(
    # showticklabels=True,
    )

# #########################################################
xaxis_left = go.layout.XAxis(
    range=[-10, 510],
    )

fig.update_xaxes(
    patch=xaxis_new.update(**xaxis_left.to_plotly_json()),
    selector=None,
    overwrite=False,
    row=1,
    col=1,
    )
fig.update_xaxes(
    patch=xaxis_new.update(**xaxis_left.to_plotly_json()),
    selector=None,
    overwrite=False,
    row=2,
    col=1,
    )


# #########################################################
xaxis_right = go.layout.XAxis(
    range=[-10, 270],
    )

fig.update_xaxes(
    patch=xaxis_new.update(**xaxis_right.to_plotly_json()),
    selector=None,
    overwrite=False,
    row=1,
    col=2,
    )
fig.update_xaxes(
    patch=xaxis_new.update(**xaxis_right.to_plotly_json()),
    selector=None,
    overwrite=False,
    row=2,
    col=2,
    )

tmp = 42
# -

# # Update xaxis of bottom two subplot

# +
xaxis_new = go.layout.XAxis(
    showticklabels=True,
    )

fig.update_xaxes(
    patch=xaxis_new,
    selector=None,
    overwrite=False,
    row=2,
    col=1,
    )
fig.update_xaxes(
    patch=xaxis_new,
    selector=None,
    overwrite=False,
    row=2,
    col=2,
    )

tmp = 42
# -

# # Update yaxis of left two subplot

# +
yaxis_new = go.layout.YAxis(
    showticklabels=True,
    )

fig.update_yaxes(
    patch=yaxis_new,
    selector=None,
    overwrite=False,
    row=1,
    col=1,
    )
fig.update_yaxes(
    patch=yaxis_new,
    selector=None,
    overwrite=False,
    row=2,
    col=1,
    )


fig.update_layout(**layout.to_plotly_json())

# fig.show()

tmp = 42

# +
from plotting.my_plotly import add_duplicate_axes

shared_axis_data = {
    "tickcolor": "black",
    "ticklen": 3,
    }

shared_xaxis_data = {
    "dtick": 50,
    **shared_axis_data,
    }

shared_yaxis_data = {
    "dtick": 1,
    **shared_axis_data,
    }

shared_meth_props = {
    # "axis_data": shared_yaxis_data,
    "tmp_define_both_axis_types": True,
    }

for i in range(1, 5):
    # #############################################################################
    add_duplicate_axes(
        fig, axis_type='x', axis_data=shared_xaxis_data,
        axis_num_list=[i, ],
        **shared_meth_props)
    add_duplicate_axes(
        fig, axis_type='y', axis_data=shared_yaxis_data,
        axis_num_list=[i, ],
        **shared_meth_props)

# +
from plotting.my_plotly import my_plotly_plot

my_plotly_plot(
    figure=fig,
    plot_name="disc_vs_dft",
    write_html=True,
    write_png=False,
    png_scale=6.0,
    write_pdf=True,
    write_svg=False,
    try_orca_write=True,
    )
# -

print(20 * "# # ")
print("All done!")
assert False

fig.show()
