# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.4.2
#   kernelspec:
#     display_name: Python [conda env:PROJ_irox_2] *
#     language: python
#     name: conda-env-PROJ_irox_2-py
# ---

# # Creating the AL plots
#
# This whole script needs be run twice with the variable `color_custom_points` set to `True` and `False` and twice for IrO2 and IrO3 stoichs, so 4 times in total

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

# + [markdown] Collapsed="false"
# # Script Inputs

# +
stoich_i = "AB3"
# stoich_i = "AB2"

lowest_N_sys_to_track = 10

# + [markdown] Collapsed="false"
# # Read Data

# + Collapsed="false"
comb_to_run = [
    ("AB2", False),
    ("AB3", False),
    ("AB2", True),
    ("AB3", True),
    ]

for (stoich_i, color_custom_points) in comb_to_run:
    print(stoich_i, color_custom_points)

    # #########################################################
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

    # #########################################################
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

        ]


    marker_color_dict = dict(zip(
        top_ids_to_track,
        color_list,
        ))


    # #########################################################
    from active_learning.al_analysis import ALAnimation

    # color_custom_points = False
    # color_custom_points = True

    ALAnim = ALAnimation(
        ALBulkOpt=AL,
        marker_color_dict=marker_color_dict,
        verbose=True,
        # color_custom_points=False,
        color_custom_points=color_custom_points,
        )

    # #########################################################
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

    # %%capture

    # #########################################################


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

    # #########################################################
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
    # fig.show()

    # #########################################################
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

print(20 * "# # ")
print("All done!")
assert False
