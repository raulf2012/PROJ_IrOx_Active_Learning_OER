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

# # Import Modules

# %load_ext autoreload
# %autoreload 2

# + jupyter={}
# %%capture

import pandas as pd

# Setting Custom Paths ********************************************************
# *****************************************************************************
import os; import sys
sys.path.append(".."); sys.path.append("../..")
sys.path.insert(0, os.path.join(
    os.environ["PROJ_col_iro2"],
    "data"))


# Python Modules **************************************************************
# *****************************************************************************
from plotly import io as pyio
import chart_studio.plotly as py
import plotly.graph_objs as go
from IPython.display import HTML


# My Modules ******************************************************************
# *****************************************************************************
from oxr_reaction.oxr_plotting_classes.oxr_plot_2d_volcano import (
    Volcano_Plot_2D)
from plotting.my_plotly import my_plotly_plot
# -

# # Script Inputs

save_plot = True

# # Read OER Data

# + jupyter={}
from sc_procc_manual_data import (
    ORR_PLT,

    # TEMP
    df_ads_e,
    # df_list,
    corrections_dict,
    oxy_ref, hyd_ref,
    )


# df_list
from proj_data_col_iro2 import proj_dir_name
# -

# # 2D Volcano Plot Instance

# +
VP = Volcano_Plot_2D(
    ORR_PLT,
    plot_range={
        "x": [+0.9, +2.0],
        "y": [-0.5, +2.0],
        })

data = VP.traces
layout = VP.get_plotly_layout()
# -

# # Plotting

layout = VP.get_plotly_layout()

# ## Small Plot

# +
# # %%capture

# #############################################################################
layout_override = dict(
    width=8 * 37.795275591,
    height=6 * 37.795275591,
    showlegend=False,

    margin=go.layout.Margin(
        autoexpand=None,
        b=8,
        l=8,
        pad=None,
        r=5,
        t=5,
        ),

    paper_bgcolor="white",
    plot_bgcolor="white",
    )

layout.update(dict(xaxis=dict(
    dtick=0.2,
    )))

# #############################################################################
for trace_i in data:
    try:
        trace_i.marker.size = 12
    except:
        pass

fig = go.Figure(
    data=data,
    layout=layout.update(layout_override))

my_plotly_plot(
    figure=fig,
    plot_name="out_plot_00_small",
    write_pdf=True,
    )
# -

fig.show()

# ## Medium Plot

# + jupyter={}
layout_override = {
    "width": 24 * 37.795275591,
    "height": 14 * 37.795275591,
    "showlegend": True,
    }

fig = go.Figure(
    data=data,
    layout=layout.update(layout_override))

my_plotly_plot(
    figure=fig,
    plot_name="out_plot_00_medium",
    write_pdf=True,
    )
# -

# ## Full Page Plot

# + jupyter={}
# %%capture

layout_override = {
    "width": 37 * 37.795275591,
    "height": 23 * 37.795275591,
    "showlegend": True}

fig = go.Figure(
    data=data,
    layout=layout.update(layout_override))

my_plotly_plot(
    figure=fig,
    plot_name="out_plot_00_large",
    write_pdf=True,
    )
fig.show()

# + active=""
#
#
#
# -

layout

# +
for ORR_i in ORR_PLT.series_list:
    series_name = ORR_i.series_name
    # print(series_name)

for ORR_i in ORR_PLT.series_list:
    series_name = ORR_i.series_name
    if series_name == "col-IrO2 (101) | Mine | CUS-site":
        break
display(ORR_i.fe_df)

# ###############################################
for ORR_i in ORR_PLT.series_list:
    series_name = ORR_i.series_name
    if series_name == "col-IrO2 (101) | Mine | Bridge-site":
        break
display(ORR_i.fe_df)


# #############################################################################


for ORR_i in ORR_PLT.series_list:
    series_name = ORR_i.series_name
    if series_name == "col-IrO2 (100) | uncoord O":
        break
display(ORR_i.fe_df)

# ###############################################
for ORR_i in ORR_PLT.series_list:
    series_name = ORR_i.series_name
    if series_name == "col-IrO2 (100) | coord O":
        break
display(ORR_i.fe_df)
# -

# ###############################################
for ORR_i in ORR_PLT.series_list:
    series_name = ORR_i.series_name

    print(series_name)
    print(ORR_i.overpotential_OER)
    print("")

 


# df_ads_e.iloc[[4, 5]]
df_ads_e

# +
energy_states_dict = ORR_i.energy_states_dict

ORR_i.fe_df

# +
# # go.contour.Line?
# # go.contour.ColorBar?


# # go.contour.colorbar.Title?

# # go.contour.colorbar.Title?

# + active=""
#
#
#

# + jupyter={}
# if save_plot:
#     save_dir = proj_dir_name
# else:
#     save_dir = "__temp__"
# save_dir = os.path.join(
#     save_dir,
#     "02_oer_analysis",
#     "oer_2d_volcano_plot",
#     )
# plotly_folder = os.path.join(
#     proj_dir_name,
#     "02_oer_analysis")

# + jupyter={}
# # #############################################################################
# import importlib.util
# path_i = os.path.join(os.environ["PROJ_col_iro2"], "workflow/02_oer_analysis")
# file_i = os.path.join(path_i, "methods.py")
# spec = importlib.util.spec_from_file_location("tmp", file_i)
# foo = importlib.util.module_from_spec(spec); spec.loader.exec_module(foo)
# plot_proto = foo.plot_proto
