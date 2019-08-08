#!/usr/bin/env python

"""TEMP.

Author(s): Raul A. Flores
"""

#| - IMPORT MODULES
import plotly.graph_objs as go


# import pandas as pd
# import plotly.plotly as py
#__|


def traces_uncertainty_i(
    df,
    value_key="form_e_0",
    uncert_key="uncert_0",
    fillcolor="red",
    name="temp_name",
    ):
    """
    """
    #| - traces_uncertainty_i
    trace_uncert_up = go.Scatter(
        # x=df["ind"].tolist(),
        y=(df[value_key] + df[uncert_key]).tolist(),
        mode='lines',

        hoverinfo="none",
        legendgroup=name,
        name=name + "_error",

        line=dict(
            width=0.,
            # color=color_after_training,
            # color="black",
            color=fillcolor,
            ),
        )

    trace_uncert_down = go.Scatter(
        # x=df["ind"].tolist(),
        y=(df[value_key] - df[uncert_key]).tolist(),
        mode='lines',
        fill='tonexty',
        # fillcolor='rgba(114,96,233,1.)',
        fillcolor=fillcolor,
        hoverinfo="none",
        showlegend=False,
        legendgroup=name,

        line=dict(
            width=0.,
            # color="black",
            color=fillcolor,
            ),
        )

    return([trace_uncert_up, trace_uncert_down])
    #__|


def process_data_for_plot(
    df,
    color0="red",
    color1="black",
    name="temp_name",
    marker_color_mode="single",  # single, array

    energy_key="form_e_0",
    marker_color_key="color_order",

    uncertainty_key="uncert_0",
    text_key="id",

    process_text_for_hover=False,
    explicit_error_bars=False,
    filled_error_traces=True,
    ):
    """
    """
    #| - process_data_for_plot
    # color_after_training = color0

    #| - Main Series
    if marker_color_mode == "single":
        color_i = color1
    elif marker_color_mode == "array":
        color_i = df[marker_color_key].tolist()

    if process_text_for_hover:
        text_list = df["id"].tolist()
    else:
        text_list = None

    if explicit_error_bars:
        error_y = {
            "type": "data",
            "array": df[uncertainty_key].tolist(),
            # "color": "grey",
            "color": "rgba(30,40,50,0.5)",
            # "thickness": 1.5,
            "thickness": 1.0,
            "width": 1.5,
            "visible": True}
    else:
        error_y = None


    # error_y=dict(
    #     type='constant',
    #     value=0.1,
    #     color='#85144B',
    #     thickness=1.5,
    #     width=3,
    #     )

    trace_0 = go.Scatter(
        # x=df["ind"].tolist(),
        y=df[energy_key].tolist(),
        error_y=error_y,
        mode="markers",
        legendgroup=name,
        name=name + "_energies",
        text=text_list,
        # hoverinfo="text",
        hoverinfo="y",
        marker=dict(
            size=6,
            color=color_i,
            colorscale="Viridis",
            line=dict(
                width=0.05,
                color="rgb(0, 0, 0)"
                ),
            ),

        line=dict(
            width=1.0,
            color=color1,
            ),
        )
    #__|

    #| - Uncertainty Series
    if filled_error_traces:
        trace_uncert_up, trace_uncert_down = traces_uncertainty_i(
            df,
            value_key=energy_key,
            uncert_key=uncertainty_key,
            fillcolor=color0,
            name=name,
            )

        data = [
            trace_uncert_up,
            trace_uncert_down,
            trace_0,
            ]
    else:
        data = [trace_0]
    #__|

    # data = [
    #     trace_uncert_up,
    #     trace_uncert_down,
    #     trace_0,
    #     ]

    return(data)
    #__|


def get_layout(
    df=None,
    energy_key="form_e_0",
    uncertainty_key="uncert_0",
    ):
    """
    """
    #| - get_layout

    #| - Analyzing dataframe for info
    # X-axis range
    if df is not None:
        x_range = [0, len(df) + 2]
    else:
        x_range = [0, 300]

    # Y-axis range
    if df is not None:
        y_range = [
            (df[energy_key] - df[uncertainty_key]).min(),
            (df[energy_key] + df[uncertainty_key]).max(),
            ]

    else:
        y_range = [-1., 2.]
    #__|

    tick_lab_size = 12 * (4. / 3.)
    axes_lab_size = 14 * (4. / 3.)

    #| - common_axis_dict
    common_axis_dict = {

        # "range": y_axis_range,
        "zeroline": False,
        "showline": True,
        "mirror": 'ticks',
        "linecolor": 'black',
        "showgrid": False,

        "titlefont": dict(size=axes_lab_size),
        "tickfont": dict(
            size=tick_lab_size,
            ),
        "ticks": 'inside',
        "tick0": 0,
        "tickcolor": 'black',
        # "dtick": 0.25,
        "ticklen": 2,
        "tickwidth": 1,
        }
    #__|

    layout = {

        # "paper_bgcolor": "rgba(0,0,0,0)",
        "plot_bgcolor": "rgba(0,0,0,0)",

        # "title": title,
        # "titlefont": go.layout.Titlefont(size=24),

        "titlefont": go.layout.title.Font(size=24),

        #| - axis properties
        "xaxis": dict(
            # range=[0, 300],
            common_axis_dict,
            **{
                "title": "Structures",
                "range": x_range,
                },
            ),

        "yaxis": dict(
            # range=[-0.5, 3.5],
            common_axis_dict,
            **{
                # H<sub>2</sub>O
                "title": "ΔG<sub>formation</sub> (eV)",
                "range": y_range,
                },
            ),
        #__|

        # Margin
        "margin": go.layout.Margin(
            b=50.,
            l=100.,
            r=50.,
            t=50.,
            ),

        "font": dict(
            family='Arial',
            # size=18,
            color='black',
            ),

        "width": 1. * 18.7 * 37.795275591,
        "height": 1. * 18.7 * 37.795275591,

        "showlegend": True,

        # "legend": dict(
        #     font=dict(
        #         size=10,
        #         ),
        #     ),

        }

    return(layout)
    #__|
