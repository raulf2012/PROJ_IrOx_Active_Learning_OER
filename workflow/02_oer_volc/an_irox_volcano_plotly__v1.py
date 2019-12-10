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

# +
import os

# os.environ["TEMP0"]
os.environ["PROJ_irox"]

# +
# %%capture

import sys
import os

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
from plotting.my_plotly import my_plotly_plot
from layout import layout
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

# +
# # %%capture

df_pourbaix, df_ads, df_surf = load_df(
    from_file=True,
    root_dir=data_dir,
    data_dir=data_dir,
    file_name="df_master.pickle",
    process_df=False)

# df_pourbaix, df_ads, df_surf = load_df(
#     from_file=False,
#     root_dir=data_dir,
#     data_dir=data_dir,
#     file_name="df_master.pickle",
#     process_df=True)

df_m = df_ads
# -

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

data = volcano_legs_data + volcano_legs_data_tmp + VP.data_points

if plot_exp_traces:
    data.insert(0, trace_iro3)
    data.insert(0, trace_iro2)

# +
if save_plot:
    save_dir = proj_dir_name
else:
    save_dir = "__temp__"

fig = go.Figure(data=data, layout=layout)
# -

os.environ.get("USER", "")

# +
layout_override = {
    "width": 35 * 37.795275591,
    "height": 19 * 37.795275591,
    "showlegend": True}
fig.layout.update(layout_override)

fig = my_plotly_plot(
    figure=fig,
    plot_name="out_plot_02_large")

# +
layout_override = {
    "width": 1.45 * 7.964 * 37.795275591,
    "height": 1.5 * 5.6002 * 37.795275591,
    "showlegend": False}
fig.layout.update(layout_override)

my_plotly_plot(
    figure=fig,
    plot_name="pl_irox_volcano_plotly_default_ooh")
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
    line=dict(
        color="blue",
        width=2,
        dash="dash",
        ),
    )

# + jupyter={"outputs_hidden": false}
for trace_i in data_kin_volc:
    data_dict = trace_i.to_plotly_json()
    fig.add_scatter(**data_dict)

fig.add_scatter(**trace_kin_10mA.to_plotly_json())

tmp = my_plotly_plot(
    figure=fig,
    plot_name="pl_irox_volcano_plotly_default_ooh__w_kinetic_volc",
    write_html=True)

fig.update_layout(
    width=900,
    height=500,
    legend=None,
    showlegend=True)
# -

fig

# + active=""
#
#
#
#
#
#
#
#
#
