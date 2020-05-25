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
# # Import Modules

# +
import os
print(os.getcwd())

import sys

import pickle
import copy

import plotly.graph_objects as go
from plotly.subplots import make_subplots

# #############################################################################
from plotting.my_plotly import my_plotly_plot

# #########################################################
# Local Imports
sys.path.insert(0, os.path.join(os.environ["PROJ_irox"], "data"))
from proj_data_irox import (
    bulk_dft_data_path,
    ids_to_discard__too_many_atoms_path,
    unique_ids_path,
    df_dij_path)

# from proj_data_irox import main_AB2_run, main_AB3_run

# +
from al_data import main_AB2_run, main_AB3_run, gens_to_plot_dict

from layout import layout as layout_base

# +
# assert False

# + [markdown] Collapsed="false"
# # Script Inputs

# +
# stoich_i = "AB3"
stoich_i = "AB2"

lowest_N_sys_to_track = 10

# + [markdown] Collapsed="false"
# # Read Data

# + Collapsed="false"
if stoich_i == "AB3":
    path_i = main_AB3_run
elif stoich_i == "AB2":
    path_i = main_AB2_run
else:
    assert False, "No data here isjfisdifjisdjifjsidfjr89u8fh8wejf"

print(stoich_i, "\n", path_i)

with open(path_i, "rb") as fle:
    AL = pickle.load(fle)

al_gen_dict = AL.al_gen_dict

gens_to_plot = gens_to_plot_dict[stoich_i]

# + Collapsed="false"
last_gen_key = list(al_gen_dict.keys())[-1]

if gens_to_plot[-1] == "last":
    gen_4 = last_gen_key
    gens_to_plot[-1] = gen_4

AL_last = al_gen_dict[last_gen_key]

model = AL_last.model

model_i = model[
    (model["duplicate"] == False) & \
    (model["acquired"] == True)
    ].sort_values("y_real")
top_ids_to_track = model_i.iloc[0:lowest_N_sys_to_track].index.tolist()


color_list = [

    # "rgba(12,0,127,1.0)",
    # "rgba(0,0,172,1.0)",
    # "rgba(0,1,215,1.0)",
    # "rgba(0,51,233,1.0)",
    # "rgba(0,83,255,1.0)",
    # "rgba(0,115,255,1.0)",
    # "rgba(0,141,243,1.0)",
    # "rgba(0,181,246,1.0)",
    # "rgba(0,220,245,1.0)",
    # "rgba(0,255,243,1.0)",

    "#fde725",
    "#b8de29",
    "#74d055",
    "#3cbc75",
    "#20a386",
    "#238a8d",
    "#2d708e",
    "#39558c",
    "#453781",
    "#481568",

    # "rgb(202,88,66)",
    # "rgb(71,189,198)",
    # "rgb(210,70,147)",
    # "rgb(120,181,66)",
    # "rgb(157,99,201)",
    # "rgb(81,163,108)",
    # "rgb(189,104,138)",
    # "rgb(131,128,57)",
    # "rgb(101,130,203)",
    # "rgb(209,154,68)",
    ]


marker_color_dict = dict(zip(
    top_ids_to_track,
    color_list,
    ))

# + Collapsed="false"
from active_learning.al_analysis import ALAnimation

color_custom_points = False
# color_custom_points = True
ALAnim = ALAnimation(
    ALBulkOpt=AL,
    marker_color_dict=marker_color_dict,
    verbose=True,
    # color_custom_points=False,
    color_custom_points=color_custom_points,
    )

# + Collapsed="false"
# # Create AL animation #########################################################
# filename_i = "ANIM_" + stoich_i + "_" + AL.name + "_NEW"
# ALAnim.create_animation(
#     duration_long=6000,
#     duration_short=6000,
#     # serial_parallel='serial',
#     serial_parallel='parallel',
#     filename=filename_i,
    
#     )

# +
# assert False

# + Collapsed="false"
traces_list = []
traces_list_tracking = []

num_dft_list = []
# top_ten_tracking_dict = dict()
for gen_i in gens_to_plot:
    print(gen_i)
    
    if gen_i < 0:
        gen_i = list(al_gen_dict.keys())[gen_i]

    AL_i = al_gen_dict[gen_i]
    model_i = AL_i.model
    
    tmp = model_i[model_i.acquired == True].shape
    num_dft_list.append(tmp[0])

    num_systems_0 = AL_i.model.shape[0]

    num_dft_i = model_i[model_i["acquired"] == True].shape[0]
    print("num_dft_i:", num_dft_i)

    trace_i, other_data_dict = ALAnim.get_trace_j(
        AL_i,
        prediction_key="y",
        uncertainty_key="err",
        plot_dft_instead_of_pred=True,
        plot_validation_dft=False,
        # trace_all_dft=True,
        trace_horiz_lines=False,
        internally_order_df=True,
        dft_calc_al_gen_text_overlay=False,
        add_vertical_track_lines=True,
        just_traces=False,
        )
    traces_list.append(trace_i)


    # #########################################################################
    # #########################################################################
    model__tracked = other_data_dict["model__tracked"]

    gen_traces_i = []
    for i_ind, row_i in model__tracked.iterrows():
        Y_main = row_i["Y_main"]
        x_ind = row_i["x_axis_ind"]
        acquired = row_i["acquired"]

        if acquired:
            color = "red"
            y = [4.7, 6]
            width = 0.8
        else:
            color = "grey"
            y = [5., 6]
            width = 0.5

        trace_i = go.Scatter(
            mode="lines",
            x=[x_ind, x_ind],
            y=y,
            name="top_ten",
            line=dict(
                width=width,
                color=color,
                ),
            )

        # data.append(trace_i)
        # traces_list_tracking.append(trace_i)

        gen_traces_i.append(trace_i)

    traces_list_tracking.append(gen_traces_i)

# +
# model__tracked = other_data_dict["model__tracked"]

# model__tracked

# + Collapsed="false"
# %%capture

a = 1 / len(traces_list)
x = 0.1
y = a + x
z = a - x / 4

column_widths = [z, z, y, z, z]
# print("column_widths:", column_widths)

fig = make_subplots(
    rows=1, cols=len(traces_list),
    column_widths=column_widths,
    shared_yaxes=True,
    horizontal_spacing=0.01)

for i_ind, traces_i in enumerate(traces_list):
    for trace_i in traces_i:
        fig.add_trace(trace_i, row=1, col=i_ind + 1)

if stoich_i == "AB3":
    # range_y = [-3.184, 5.529]
    range_y = [-0.8, 1.5]
elif stoich_i == "AB2":
    range_y = None
        

layout_override = dict(
    # height=200,
    # width=650,

    height=5.291667 * 37.795275591,
    # width=17.5 * 37.795275591,
    width=17.7 * 37.795275591,

    margin=go.layout.Margin(
        b=0,
        l=10,
        r=5,
        t=5),
    xaxis=dict(
        range=[-20, num_systems_0 + 10],
        showticklabels=False,
        # ticks="",
        ticks=None,
        ),
    yaxis=dict(
        range=range_y,
        mirror=True,
        showticklabels=False,
        # ticks="",
        ticks=None,
        ),
    )
        
layout_base_cpy = copy.deepcopy(layout_base)
layout = layout_base_cpy.update(layout_override)
fig.update_layout(layout)

fig.update_xaxes(layout.xaxis)
fig.update_yaxes(layout.yaxis)

# #############################################################################
# 
fig.update_xaxes(
    linecolor="red",
    row=1, col=3)
fig.update_yaxes(
    linecolor="red",
    row=1, col=3)

# Update first subplot to have tick props
fig.update_yaxes(
    showticklabels=True,
    ticks="outside",
    dtick=0.5,
    row=1, col=1)

# + Collapsed="false"
fig_al_series = copy.deepcopy(fig)

my_plotly_plot(
    figure=fig,
    plot_name=stoich_i + "_" + "al_5_gens_in_row",
    write_html=True,
    # write_png=True,
    # png_scale=10,
    # write_pdf=True,
    )

fig.layout.update(paper_bgcolor="white")
fig.show()

# + Collapsed="false"
figs_dict = {
    # "fig_inset": fig_inset,
    # "fig_main": fig_main,
    "fig_al_series": fig_al_series,
    # "fig_al_series_top10_marked": fig_al_series_top10_marked,
    "traces_tracking": traces_list_tracking,
    "num_dft_list": num_dft_list,
    }

# Pickling data ######################################################
directory = "out_data"
if not os.path.exists(directory): os.makedirs(directory)
# with open(os.path.join(directory, stoich_i + "_" +  + "figs_dict__v2.pickle"), "wb") as fle:
with open(os.path.join(directory, stoich_i + "_" + str(color_custom_points) + "_" + "figs_dict__v2.pickle"), "wb") as fle:
    pickle.dump(figs_dict, fle)
# #####################################################################
# -

assert False

# + active=""
#
#
#
#
#
#
#
#

# +
fig_al_series.data[0].name

for i in fig_al_series.data:
    print(i.name)

# +
# ALAnim.model__tracked

# +
for gen_i in gens_to_plot:
    if gen_i < 0:
        al_gen_i = al_gen_dict[
            list(al_gen_dict.keys())[gen_i]    
            ]
    else:
        al_gen_i = al_gen_dict[gen_i]

    model = al_gen_i.model
    model_0  = model.loc[list(marker_color_dict.keys())]

# data = []
# for i_ind, row_i in model_0.iterrows():
#     Y_main = row_i["Y_main"]
#     x_ind = row_i["x_axis_ind"]
#     acquired = row_i["acquired"]

#     if acquired:
#         color = "red"
#         y = [4.7, 6]
#         width = 0.8
#     else:
#         color = "grey"
#         y = [5., 6]
#         width = 0.5

#     trace_i = go.Scatter(
#         mode="lines",
#         x=[x_ind, x_ind],
#         y=y,
#         name="top_ten",
#         line=dict(
#             width=width,
#             color=color,
#             ),
#         )

#     data.append(trace_i)

# + jupyter={}
# list(marker_color_dict.keys())

# list(al_gen_dict.keys())

# + jupyter={}

    # path_i = os.path.join(
    #     os.environ["PROJ_irox"],
    #     "workflow/ml_modelling/00_ml_workflow/191102_new_workflow/00_abx_al_runs/out_data/AB3/gp_ucb_True",

    #     # "01_attempt/AL_geheneva.pickle",
    #     # "01_attempt/AL_pifehohu.pickle",
    #     # "01_attempt/AL_suturomo.pickle",

    #     "TEST_AL_2_fugunefo.pickle",

    #     # NEW RUNS
    #     # "TEST_AL_2_fugunefo.pickle",
    #     # "TEST_AL_2_seruladi.pickle",
    #     # "TEST_AL_masahiti.pickle",  # This is a good one to pick ************
    #     )



    # path_i = os.path.join(
    #     os.environ["PROJ_irox"],
    #     # "workflow/ml_modelling/00_ml_workflow/191102_new_workflow/00_abx_al_runs/out_data/AB2/gp_ucb_True",
    #     # "AL_piritapo.pickle",

    #     "workflow/ml_modelling/00_ml_workflow/191102_new_workflow/01_abx_al_runs_new/out_data/AB2/gp_ucb_True",
    #     "TEST_AL_6_radetama.pickle",
    #     ) 


# + jupyter={}
# AL.duplicate_swap_dict

# AL.dupl

# +
#| - Minor tick axis settings
# minor_ticks_data = {

# #| - Main AL Plots
#     "main_al_plots": {
#         "x": 100,
#         "y": 0.25,
#         },

#     }

#__|

# + Collapsed="false"
# from inputs import (
#     stoich_i,
#     lowest_N_sys_to_track,
#     )
