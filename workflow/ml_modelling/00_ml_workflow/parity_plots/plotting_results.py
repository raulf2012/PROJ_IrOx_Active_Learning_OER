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

# + jupyter={}
import os
print(os.getcwd())
import sys

import pickle

import chart_studio.plotly as py
import plotly.graph_objs as go

from plotting.my_plotly import my_plotly_plot

from layout import layout
# -

import pandas as pd

# # Read Data

# + jupyter={}
# #############################################################################
path_i = os.path.join(
    "out_data",
    "AB2_post_dft_cv_data.pickle")
with open(path_i, "rb") as fle:
    model_post_ab2 = pickle.load(fle)
# #############################################################################

# #############################################################################
path_i = os.path.join(
    "out_data",
    "AB2_pre_dft_cv_data.pickle")
with open(path_i, "rb") as fle:
    model_pre_ab2 = pickle.load(fle)
# #############################################################################



# #############################################################################
import pickle; import os
path_i = os.path.join(
    # os.environ[""],
    "out_data",
    "AB3_post_dft_cv_data.pickle")
with open(path_i, "rb") as fle:
    model_post_ab3 = pickle.load(fle)
# #############################################################################

# #############################################################################
path_i = os.path.join(
    "out_data",
    "AB3_pre_dft_cv_data.pickle")
with open(path_i, "rb") as fle:
    model_pre_ab3 = pickle.load(fle)
# #############################################################################

# +
data = []

scatter_shared = go.Scatter(
    mode="markers",
    marker=dict(
        opacity=1.,
        size=4,
        line=dict(
            color='black',
            # width=0.1,
            width=0.,
            ),
        ),
    )

scatter_shared_post = go.Scatter(opacity=.6, marker=dict(color="rgba(0,100,255,1.)"))
scatter_shared_pre = go.Scatter(marker=dict(color="rgba(150,150,150,1.)"))

scatter_shared_iro2 = go.Scatter(marker=dict(symbol="circle"))
scatter_shared_iro3 = go.Scatter(marker=dict(symbol="circle"))
# -

# # TRACES | Main data traces

# +
# #############################################################################
model_i = pd.concat([
    model_pre_ab2,
    model_pre_ab3,    
    ], axis=0)
trace = go.Scatter(y=model_i.y, x=model_i.y_real)
trace.update(scatter_shared)
trace.update(scatter_shared_pre)
data.append(trace)


# #############################################################################
model_i = pd.concat([
    model_post_ab2,
    model_post_ab3,    
    ], axis=0)
trace = go.Scatter(y=model_i.y, x=model_i.y_real)
trace.update(scatter_shared)
trace.update(scatter_shared_post)
data.append(trace)
# -

# # TRACE | Parity Line

# +
trace_xy = go.Scatter(
    x=[-4, 6],
    y=[-4, 6],
    mode="lines",
    line=dict(
        color="red",
        dash="dot",
        ),
    )

data.append(trace_xy)
# -

# # Plotting

# + jupyter={}
layout.update(dict(
    xaxis=dict(range=[-2.9, 5.7]),
    yaxis=dict(
        range=[-2.9, 5.7],
        scaleanchor = "x",
        scaleratio = 1,    
        ),
    
    paper_bgcolor="rgb(255,255,255)",
    plot_bgcolor="rgb(255,255,255)",
    annotations=None,
    ))

fig = go.Figure(data=data, layout=layout)
fig = my_plotly_plot(
    figure=fig,
    plot_name="parity_plot",
    write_html=True,
    write_pdf=True,
    write_png=True,
    png_scale=1.0,
    )

fig.show()

# +
layout.update(dict(
    height=6.812 * 37.795275591,
    width=6.74 * 37.795275591,

    paper_bgcolor="rgba(255,255,255,0.)",
    plot_bgcolor="rgba(255,255,255,0.)",

    showlegend=False,
    margin=go.layout.Margin(
        autoexpand=None,
        b=10,
        l=10,
        pad=None,
        r=10,
        t=10,
        ),
    ))

fig = go.Figure(data=data, layout=layout)
fig = my_plotly_plot(
    figure=fig,
    plot_name="main_parity_plot",
    write_html=True,
    write_pdf=True,
    write_png=True,
    png_scale=20.0,
    )

fig.show()

# + active=""
#
#
#
#

# + jupyter={}



# def method(row_i, argument_0, optional_arg=None):
#     """
#     """
#     calc_dH()
#     return(argument_0)

# arg1 = "TEMP_0"
# df_i = model_i
# df_i["column_name"] = df_i.apply(
#     method,
#     axis=1,
#     args=(arg1, ),
#     optional_arg="TEMP_1"
#     )
