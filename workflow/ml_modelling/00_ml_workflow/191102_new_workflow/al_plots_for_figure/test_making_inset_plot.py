# -*- coding: utf-8 -*-
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

# +
import plotly.graph_objs as go

from plotting.my_plotly import my_plotly_plot
# -

# # Script Inputs

stoich_i = "AB2"

# # Read Data

# #############################################################################
import pickle; import os
path_i = os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling/00_ml_workflow/191102_new_workflow/performance_comp/out_data",
    stoich_i + "_" + "fig_al_perf.pickle")
with open(path_i, "rb") as fle:
    fig_perf = pickle.load(fle)
# #############################################################################

# +
# #############################################################################
import pickle; import os
path_i = os.path.join(
    # os.environ[""],
    "out_data",
    stoich_i + "_" + "figs_dict.pickle")
with open(path_i, "rb") as fle:
    fig_dict = pickle.load(fle)
# #############################################################################

fig_inset = fig_dict["fig_inset"]
fig_main = fig_dict["fig_main"]

# #####################################
data_all = []

data_inset = fig_inset.data
data_main = fig_main.data

data_all.extend(data_inset)
data_all.extend(data_main)

# #####################################
layout_inset = fig_inset.layout
layout_main = fig_main.layout
# -

fig_main.show()

for trace_i in data_inset:
    trace_i.xaxis = "x2"
    trace_i.yaxis = "y2"

layout = go.Layout(
    xaxis2=dict(
        domain=[0.1, 0.35],
        anchor='y2'
        ),
    yaxis2=dict(
        domain=[0.55, 0.95],
        anchor='x2'
        )
    )

# +
layout["xaxis2"].update(layout_inset["xaxis"])
layout["yaxis2"].update(layout_inset["yaxis"])

layout_main.template = None
layout = layout.update(layout_main)

# layout['paper_bgcolor'] = 'rgba(0,0,0,0)'
# layout['plot_bgcolor'] = 'rgba(0,0,0,0)'

layout['paper_bgcolor'] = 'white'
layout['plot_bgcolor'] = 'white'

# +
# Annotation of inset plot square
shapes = [
    go.layout.Shape(
        type="rect",
        xref="x",
        yref="y",

        x0=-1,
        y0=-2.65,
        x1=16.5,
        y1=-2.2,

        line=dict(
            color="grey",
            width=1.,
            ),
        # fillcolor="LightSkyBlue",
        ),


    go.layout.Shape(
        type="line",
        xref="paper",   
        yref="paper",
        x0=0.1,
        y0=0.075,
        x1=0.35,
        y1=0.55,
        line=dict(
            color="black",
            width=1,
            )
        ),


    go.layout.Shape(
        type="line",
        xref="paper",
        yref="paper",
        x0=0.035 - 0.001479855,
        y0=0.119 + 0.003112469 + 0.0005,
        x1=0.1,
        y1=0.55,
        line=dict(
            color="black",
            width=1,
            )
        ),

    ]

layout.shapes = shapes

# +
fig = go.Figure(data=data_all, layout=layout)

tmp = my_plotly_plot(
    figure=fig,
    plot_name=stoich_i + "_" + "al_main_w_inset",
    write_html=True,
    write_png=True,
    png_scale=10.0,
    write_pdf=True)

# fig.show()
# -

fig.show()

assert False

fig.layout

# + {"active": ""}
#
#
#
#
#

# +
from plotly.subplots import make_subplots
import plotly.graph_objects as go

fig_subplot = make_subplots(rows=1, cols=2)

for trace_i in fig.data:
    fig_subplot.add_trace(
        trace_i,
        row=1, col=1
        )

# fig_subplot.add_trace(
#     trace_i,
#     row=1, col=1
#     )

# fig.add_trace(
#     # go.Scatter(x=[1, 2, 3], y=[4, 5, 6]),
#     fig.data,
#     row=1, col=1
# )

fig_subplot.add_trace(
    go.Scatter(x=[20, 30, 40], y=[50, 60, 70]),
    row=1, col=2
)

fig_subplot.update_layout(height=600, width=800, title_text="Subplots")
fig_subplot.show()

# +
# fig_subplot.data

# fig_subplot.update_layout(layout)


fig_subplot.layout

# + {"jupyter": {"source_hidden": true}}
# dir(layout_main)

# type(layout_main)
# go.Layout

# layout_main.template
# layout_main.__delattr__("template")

# layout_main.template = None

# + {"jupyter": {"source_hidden": true}}
# layout_tmp = go.Layout({
#     'font': {'color': 'black', 'family': 'Arial'},
#     'height': 270,
#     'margin': {'b': 20, 'l': 20, 'r': 10, 't': 10},
#     'paper_bgcolor': 'rgba(0,0,0,0)',
#     'plot_bgcolor': 'rgba(0,0,0,0)',
#     'shapes': [{'line': {'color': 'grey', 'width': 1.0},
#                 'type': 'rect',
#                 'x0': -1,
#                 'x1': 16.5,
#                 'xref': 'x',
#                 'y0': -2.65,
#                 'y1': -2.2,
#                 'yref': 'y'}],
#     'showlegend': False,
#     # 'template': '...',
#     'width': 450,
#     'xaxis': {'linecolor': 'red',
#               'linewidth': 1.0,
#               'mirror': True,
#               'range': [-10, 258],
#               'showgrid': False,
#               'showline': True,
#               'showticklabels': True,
#               'tickcolor': 'black',
#               'tickfont': {'size': 10.666666666666666},
#               'ticks': 'outside',
#               'title': {'font': {'color': 'black', 'size': 13.333333333333332}, 'text': 'Candidate Space'},
#               'zeroline': False},
#     'yaxis': {'dtick': 1.0,
#               'linecolor': 'red',
#               'linewidth': 1.0,
#               'mirror': True,
#               'showgrid': False,
#               'showline': True,
#               'showticklabels': True,
#               'tickcolor': 'black',
#               'tickfont': {'size': 10.666666666666666},
#               'ticks': 'outside',
#               'title': {'font': {'color': 'black', 'size': 13.333333333333332}, 'text': 'ΔH<sub>f</sub> (eV)'},
#               'zeroline': False}
# })
