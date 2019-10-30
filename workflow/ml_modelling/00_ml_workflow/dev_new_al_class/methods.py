#!/usr/bin/env python

"""TEMP.

Author: Raul A. Flores
"""

#| - Import  Modules
import numpy as np

import pandas as pd

import plotly.graph_objs as go

from plotting.my_plotly import my_plotly_plot
# __|

def get_trace_j(
    model,
    df_bulk_dft=None,
    prediction_key="prediction",
    uncertainty_key="uncertainty",
    plot_dft_instead_of_pred=True,
    trace_all_dft=True,
    trace_horiz_lines=True,
    marker_size=10,
    ):
    """

    Args:
        plot_dft_instead_of_pred:
            Plot the actual DFT energy instead of the predicted value
    """
    #| - get_trace_j
    redundant_key = "TEMP_redundant"
    redundant_key_global = "TEMP_redundant_global"

    model_i = model
    data = []

    # ########################################################################
    #| - NEW
    def method(row_i):
        #| - method
        computed = row_i["acquired"]
        y_real = row_i["y_real"]

        actually_computed = False
        if not np.isnan(y_real):
            actually_computed = True
        dft_energy = y_real

        predicted_energy = row_i[prediction_key]
        predicted_uncert = row_i[uncertainty_key]

        # #####################################################################
        new_column_values_dict = {
            "out_energy": None,
            "out_uncert": None}

        # #####################################################################
        if computed and actually_computed:
            new_column_values_dict["Y_main"] = dft_energy
            new_column_values_dict["Y_uncer"] = 0.
        else:
            new_column_values_dict["Y_main"] = predicted_energy
            new_column_values_dict["Y_uncer"] = predicted_uncert

        # #####################################################################
        for key, value in new_column_values_dict.items():
            row_i[key] = value
        return(row_i)
        #__|

    if plot_dft_instead_of_pred:
        model_i = model_i.apply(method, axis=1)
    else:
        model_i["Y_main"] = model_i[prediction_key]
        model_i["Y_uncer"] = model_i[uncertainty_key]
    #__|

    # #########################################################################
    #| - Applying formating to df
    def method(row_i, marker_size):
        #| - method
        new_column_values_dict = {
            "marker_size": 4,
            "marker_line_color": "rgb(0,0,0)",
            "marker_line_size": 0.05,
            }

        computed_bool = row_i["acquired"]
        # computed_bool = row_i.get("computed", False)

        if computed_bool:
            # new_column_values_dict["marker_size"] = 10
            new_column_values_dict["marker_size"] = marker_size
            new_column_values_dict["marker_line_color"] = "red"
            new_column_values_dict["marker_line_size"] = 1.5

        # #########################################################################
        for key, value in new_column_values_dict.items():
            row_i[key] = value
        return(row_i)
        #__|

    df_i = model_i
    df_i = df_i.apply(
        method,
        axis=1,
        marker_size=marker_size,
        )
    model_i = df_i
    #__|


    model_i = model_i.sort_values("Y_main")
    model_i["x_axis_ind"] = [i for i in range(model_i.shape[0])]
    model_i = model_i.sort_index()

    # return(model_i)

    # #########################################################################
    #| - Main data trace
    trace_i = go.Scatter(
        x=model_i["x_axis_ind"],
        # y=model_i[prediction_key],
        y=model_i["Y_main"],
        error_y=dict(
            type='data',
            array=model_i["Y_uncer"],
            visible=True,
            thickness=0.3,
            width=1.5,
            # color="rgba(120,120,120,1.0)",
            color="rgba(80,80,60,1.0)",
            ),
        # name=model_i["id"],
        mode="markers",
        # text=model_i["id"],
        # hoverinfo="text",
        marker={
            "opacity": 0.9,
            "size": model_i["marker_size"],
            "color": ["black"],
            # "color": model_i[prediction_key],
            # "color": model_i["energy_pa"],
            # "colorscale": "Viridis",
            "line": {
                "width": model_i["marker_line_size"],
                "color": model_i["marker_line_color"],
                },
            },
        )
    data.append(trace_i)
    #__|

    # #########################################################################
    #| - Horizontal lines at E minimum
    min_e = model_i[prediction_key].min()
    min_e_w_uncert = (model_i[prediction_key] - model_i[uncertainty_key]).min()

    # #########################################################################
    trace_i = go.Scatter(
        x=[0, len(model_i["x_axis_ind"].tolist())],
        y=[min_e, min_e],
        mode="lines",
        line=dict(
            color="firebrick",
            width=2,
            dash="dot",
            ),
        )

    if trace_horiz_lines:
        data.append(trace_i)

    # #########################################################################
    trace_i = go.Scatter(
        x=[0, len(model_i["x_axis_ind"].tolist())],
        y=[min_e_w_uncert, min_e_w_uncert],
        mode="lines",
        line=dict(
            color="black",
            width=1.5,
            dash="dash",
            ),
        )
    if trace_horiz_lines:
        data.append(trace_i)
    #__|

    # #########################################################################
    #| - Validation DFT Trace
    df_dft_val = model_i[~model_i["y_real"].isna()]
    df_dft_val = df_dft_val[["y_real", "acquired", "x_axis_ind"]]

    df_dft_val


    def method(row_i):
        #| - method
        computed = row_i["acquired"]
        y_real = row_i["y_real"]

        actually_computed = False
        if not np.isnan(y_real):
            actually_computed = True
        dft_energy = y_real

        # computed = row_i["acquired"]
        # actually_computed = row_i["actually_computed"]
        # redundant_global = row_i[redundant_key_global]

        if computed and actually_computed:
            row_i["marker_opacity"] = 0.
            row_i["marker_color"] = "black"
            row_i["marker_line_color"] = "pink"
            row_i["marker_size"] = 10.
            # row_i["marker_symbol"] = "cross-open"
            # row_i["marker_symbol"] = "asterisk-open"
            # row_i["marker_symbol"] = "diamond-open"
            row_i["marker_symbol"] = "diamond"
            # row_i["marker_symbol"] = "circle"
            # row_i["energy_pa"] = -4.6

            # if redundant_global:
            #     row_i["marker_opacity"] = 1.
            #     row_i["marker_symbol"] = "star"
            #     row_i["marker_color"] = "green"

        if not computed and actually_computed:
            row_i["marker_opacity"] = 1.
            row_i["marker_color"] = "black"
            row_i["marker_line_color"] = "pink"
            row_i["marker_size"] = 10.
            row_i["marker_symbol"] = "diamond-open"

            # if redundant_global:
            #     row_i["marker_color"] = "green"
            #     row_i["marker_symbol"] = "star-open"

        if actually_computed == False:
            row_i["marker_opacity"] = 1.
            row_i["marker_color"] = "grey"
            row_i["marker_size"] = 0.
            row_i["marker_symbol"] = "diamond"

        return(row_i)
        #__|

    df_dft_val = df_dft_val.apply(method, axis=1)
    # df_tmp = df_tmp.fillna(value=-5.)

    df_tmp = df_dft_val
    trace_tmp = go.Scatter(
        y=df_tmp["y_real"],
        x=df_tmp["x_axis_ind"],
        # y=df_bulk_dft["energy_pa"],
        # x=df_bulk_dft["x_axis_ind"],
        mode="markers",
        marker=dict(
            # symbol="cross-open",
            symbol=df_tmp["marker_symbol"],
            color=df_tmp["marker_color"],
            opacity=df_tmp["marker_opacity"],
            # size=10,
            size=df_tmp["marker_size"],
            line=dict(
                color=df_tmp["marker_line_color"],
                # color="black",
                width=1.,
                )
            ),
        )

    data.append(trace_tmp)
    # __|

    return(data)



    # | - __old__
    # if df_bulk_dft is not None:
    #     #| - Known DFT data
    #     filter_indices = [i for i in model_i.index if i in df_bulk_dft.index]
    #     df_tmp = pd.concat(
    #         [
    #             model_i.drop("energy_pa", axis=1),
    #             df_bulk_dft.loc[filter_indices],
    #             ],
    #         # keys=["a", "b"],
    #         axis=1,
    #         sort=True)
    #     df_tmp = df_tmp.reindex(model_i.index.sort_values().tolist())
    #
    #     def method(row_i):
    #         #| - method
    #         computed = row_i["computed"]
    #         actually_computed = row_i["actually_computed"]
    #         redundant_global = row_i[redundant_key_global]
    #
    #         if computed and actually_computed:
    #             row_i["marker_opacity"] = 0.
    #             row_i["marker_color"] = "black"
    #             row_i["marker_line_color"] = "pink"
    #             row_i["marker_size"] = 10.
    #             # row_i["marker_symbol"] = "cross-open"
    #             # row_i["marker_symbol"] = "asterisk-open"
    #             # row_i["marker_symbol"] = "diamond-open"
    #             row_i["marker_symbol"] = "diamond"
    #             # row_i["marker_symbol"] = "circle"
    #             # row_i["energy_pa"] = -4.6
    #
    #             if redundant_global:
    #                 row_i["marker_opacity"] = 1.
    #                 row_i["marker_symbol"] = "star"
    #                 row_i["marker_color"] = "green"
    #
    #         if not computed and actually_computed:
    #             row_i["marker_opacity"] = 1.
    #             row_i["marker_color"] = "black"
    #             row_i["marker_line_color"] = "pink"
    #             row_i["marker_size"] = 10.
    #             row_i["marker_symbol"] = "diamond-open"
    #
    #             if redundant_global:
    #                 row_i["marker_color"] = "green"
    #                 row_i["marker_symbol"] = "star-open"
    #
    #         if row_i["actually_computed"] == False:
    #             row_i["marker_opacity"] = 1.
    #             row_i["marker_color"] = "grey"
    #             row_i["marker_size"] = 0.
    #             row_i["marker_symbol"] = "diamond"
    #
    #         return(row_i)
    #         #__|
    #
    #     df_tmp = df_tmp.apply(method, axis=1)
    #     df_tmp = df_tmp.fillna(value=-5.)
    #
    #     #| - Scatter plot trace
    #     trace_tmp = go.Scatter(
    #         y=df_tmp["energy_pa"],
    #         x=df_tmp["x_axis_ind"],
    #         # y=df_bulk_dft["energy_pa"],
    #         # x=df_bulk_dft["x_axis_ind"],
    #         mode="markers",
    #         marker=dict(
    #             # symbol="cross-open",
    #             symbol=df_tmp["marker_symbol"],
    #             color=df_tmp["marker_color"],
    #             opacity=df_tmp["marker_opacity"],
    #             # size=10,
    #             size=df_tmp["marker_size"],
    #             line=dict(
    #                 color=df_tmp["marker_line_color"],
    #                 # color="black",
    #                 width=1.,
    #                 )
    #             ),
    #         )
    #
    #     data.append(trace_tmp)
    #     #__|
    #
    #     #__|
    # __|


    #__|