# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.0
#   kernelspec:
#     display_name: Python [conda env:PROJ_IrOx_Active_Learning_OER]
#     language: python
#     name: conda-env-PROJ_IrOx_Active_Learning_OER-py
# ---

# # OER Volcano for IrOx systems
#
# ***

# # Import Modules

# %%capture
# %load_ext autoreload
# %autoreload 2

import os
print(os.getcwd())
import sys

# +
# %%capture

sys.path.insert(
    0, os.path.join(
        os.environ["PROJ_irox"],
        "workflow"))

sys.path.insert(
    0, os.path.join(
        os.environ["PROJ_irox"],
        "data"))

from an_data_processing import load_df

# #############################################################################
# Python Modules

import pickle

import numpy as np

import plotly.graph_objs as go

# #############################################################################
# My Modules
from oxr_reaction.oxr_rxn import ORR_Free_E_Plot
from oxr_reaction.oxr_plotting_classes.oxr_plot_volcano import Volcano_Plot

# #############################################################################
# Project Data
from proj_data_irox import (
    proj_dir_name,
    smart_format_dict,
    gas_molec_dict,
    scaling_dict_ideal,
    scaling_dict_fitted,
    exp_irox_lim_pot,
    data_dir,
    groupby_props,
    irox_bulk_color_map)

# #############################################################################
# Local Imports
from plotting.my_plotly import (
    my_plotly_plot,
    add_minor_ticks,
    add_duplicate_axes,
    )

# from layout import layout
from layout2 import layout
# -

# # Script Inputs

# +
save_plot = False
plot_exp_traces = True

plot_range = {
    "y": [2., 1.4],
    "x": [1., 2.],
    }

# + [markdown] toc-hr-collapsed=true
# # Read and Process Data Frame
# -

# ## Read dataframe from file

df_pourbaix, df_ads, df_surf = load_df(
    from_file=True,
    root_dir=data_dir,
    data_dir=data_dir,
    file_name="df_master.pickle",
    process_df=False)
df_m = df_ads

# # ORR_Free_E_Plot Instance

ORR_PLT = ORR_Free_E_Plot(
    free_energy_df=None,
    state_title="adsorbate",
    free_e_title="ads_e",
    smart_format=smart_format_dict,
    color_list=None,
    rxn_type="OER")

# # Processing Data

# +
prop_name_list = [
    'bulk_system',
    # 'coverage',
    'coverage_type',
    'facet',
    'surface_type',
    ]

grouped = df_m.groupby(groupby_props)

for i_ind, (name, group) in enumerate(grouped):
    df_i = group

    if not any([np.isnan(i) for i in df_i.elec_energy.tolist()]):
        ORR_PLT.add_series(
            df_i,
            plot_mode="all",
            overpotential_type="OER",
            property_key_list=prop_name_list,
            add_overpot=False,
            )

# + [markdown] toc-hr-collapsed=true
# # Experimental IrOx Activity Traces
# -

# ## Horizontal data traces

# +
trace_iro3 = go.Scatter(
    x=plot_range["x"],
    y=2 * [exp_irox_lim_pot["10_mA/cm2"]["SrIrO3"]],    
    mode="lines",
    name="lines",
    line={
        "color": irox_bulk_color_map["IrO3"],
        "width": 1,
        "dash": "dash",
        },
    )

trace_iro2 = go.Scatter(
    x=plot_range["x"],
    y=2 * [exp_irox_lim_pot["10_mA/cm2"]["IrO2(110)"]],
    mode="lines",
    name="lines",
    line={
        "color": irox_bulk_color_map["IrO2"],
        "width": 1,
        "dash": "dash",
        },
    )
# -

# # Volcano Plot

# +
VP = Volcano_Plot(
    ORR_PLT,
    x_ax_species="o-oh",  # 'o-oh' or 'oh'
    smart_format_dict=smart_format_dict,
    plot_range=plot_range,
    )

VP.create_volcano_relations_plot()

volcano_legs_data = VP.create_volcano_lines(
    gas_molec_dict=gas_molec_dict,
    scaling_dict=scaling_dict_ideal,
    plot_all_legs=False,
    plot_min_max_legs=True,
    trace_priority="bottom",  # 'top' or 'bottom'
    )

volcano_legs_data_tmp = VP.create_volcano_lines(
    gas_molec_dict=gas_molec_dict,
    scaling_dict=scaling_dict_fitted,
    plot_all_legs=False,
    plot_min_max_legs=True,
    trace_priority="bottom",  # 'top' or 'bottom'
    legs_to_plot=[
        # "o2_to_ooh",
        "ooh_to_o",
        "o_to_oh",
        # "oh_to_h2o",
        ],
    line_color="grey"
    )

# data = volcano_legs_data + volcano_legs_data_tmp + VP.data_points
data = volcano_legs_data_tmp + volcano_legs_data + VP.data_points

if plot_exp_traces:
    data.insert(0, trace_iro3)
    data.insert(0, trace_iro2)
# -

fig = go.Figure(data=data, layout=layout)

# +
# layout_override = {
#     "width": 35 * 37.795275591,
#     "height": 19 * 37.795275591,
#     "showlegend": True}
# fig.layout.update(layout_override)

fig = my_plotly_plot(
    figure=fig,
    plot_name="out_plot_02_large")

# +



# -

# # TEMP | Changing line type of volcano

# +
# for trace in fig.data:
#     tmp = 42
    
# #     print(trace.name)
    
#     if trace.name == "activity volcano":
#         trace_tmp = trace

# trace_tmp.line.dash = "7px,2px,7px,2px"

# +
# layout_override = {
#     "width": 1.45 * 7.964 * 37.795275591,
#     "height": 1.5 * 5.6002 * 37.795275591,
#     "showlegend": False}
# fig.layout.update(layout_override)

shared_axis_props = dict(ticklen=3)

ticks_props_new_x = dict(
    dtick=0.1,
    **shared_axis_props)
ticks_props_new_y = dict(
    dtick=0.05,
    **shared_axis_props)

# fig = add_minor_ticks(
#     fig,
#     # axis='both',
#     axis="both",
#     ticks_props_new_x=ticks_props_new_x,
#     ticks_props_new_y=ticks_props_new_y,
#     )

add_duplicate_axes(
    fig, axis_type='x',
    axis_data=ticks_props_new_x)

add_duplicate_axes(
    fig, axis_type='y',
    axis_data=ticks_props_new_y)

my_plotly_plot(
    figure=fig,
    plot_name="pl_irox_volcano_plotly_default_ooh",
    write_html=True,
    write_png=False,
    png_scale=6.0,
    write_pdf=True,
    write_svg=False,
    try_orca_write=True,
    )
# -

# # Adding Kinetic Volcano Traces

# +
# #############################################################################
path_i = os.path.join(
    "out_data",
    "kinetic_volcano_trace.pickle")
with open(path_i, "rb") as fle:
    data_kin_volc = pickle.load(fle)
# #############################################################################

# #############################################################################
path_i = os.path.join(
    "out_data",
    "df_10mA.pickle")
with open(path_i, "rb") as fle:
    df_10mA = pickle.load(fle)
# #############################################################################
# -

df = df_10mA
trace_kin_10mA = go.Scatter(
    x=df["descriptor"],
    y=df["U_lim"],
    mode="lines",
    name="temp_4348",
    line=dict(
        # color="#1ee148",
        color="#3e9bf2", 
        width=2.5,
        dash="dot",
        ),
    )

# + jupyter={"outputs_hidden": false}
fig.add_scatter(**trace_kin_10mA.to_plotly_json())

# + jupyter={"outputs_hidden": false}
data_tmp = fig.data

data = list(data_tmp)[-1:] + list(data_tmp)[0:-1]
fig.data = data

# + jupyter={"outputs_hidden": false}
tmp = my_plotly_plot(
    figure=fig,
    plot_name="pl_irox_volcano_plotly_default_ooh__w_kinetic_volc",
    write_html=True,
    write_pdf=True,
    write_svg=False,
    try_orca_write=True)

fig.show()
# -

# Pickling data ######################################################
import os; import pickle
directory = "out_data"
if not os.path.exists(directory): os.makedirs(directory)
with open(os.path.join(directory, "oer_volcano_trace.pickle"), "wb") as fle:
    pickle.dump(fig, fle)
# #####################################################################

# + active=""
#
#
#

# +
# for i_cnt, trace in enumerate(fig.data):
#     tmp = 42
#     print(trace.name)
#     if trace.name == "temp_4348":
#         print(trace)
#         print("ISDFIJISDJFIJSDI")

# fig.data[-1:]
# print(len(fig.data))
# print(len(fig.data))
# data_0 = list(fig.data)
# data_0.insert(0, trace_kin_10mA)
# fig.data = data_0

# data_tmp = fig.data

# # print(len(fig.data))
# # len(list(data_tmp)[0:-1])
# # list(data_tmp)

# for i_cnt, trace in enumerate(fig.data):
#     tmp = 42
#     print(trace.name)
#     if trace.name == "temp_4348":
#         print(trace)
#         print("ISDFIJISDJFIJSDI")

# for trace_i in data_kin_volc:
#     data_dict = trace_i.to_plotly_json()
#     fig.add_scatter(**data_dict)
