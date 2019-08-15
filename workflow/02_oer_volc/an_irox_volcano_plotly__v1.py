# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.4'
#       jupytext_version: 1.2.1
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # OER Volcano for IrOx systems
#
# ***

# # Import Modules

# + {"jupyter": {"source_hidden": true}}
# %%capture
# %load_ext autoreload
# %autoreload 2

# + {"jupyter": {"source_hidden": true}}
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
    )

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

# + {"toc-hr-collapsed": true, "cell_type": "markdown"}
# # Read and Process Data Frame
# -

# ## Read dataframe from file

# +
# %%capture

df_pourbaix, df_ads, df_surf = load_df(
    from_file=True,
    root_dir=data_dir,
    data_dir=data_dir,
    file_name="df_master.pickle",
    process_df=False)

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

# + {"toc-hr-collapsed": true, "cell_type": "markdown"}
# # Experimental IrOx Activity Traces
# -

# ## Horizontal data traces

# +
trace_iro3 = go.Scatter(
    x=plot_range["x"],
    y=2 * [exp_irox_lim_pot["iro3"]["lim_pot"]],
    mode="lines",
    name="lines",
    line={
        "color": exp_irox_lim_pot["iro3"]["line_color"],
        "width": 1,
        "dash": "dash",
        },
    )

trace_iro2 = go.Scatter(
    x=plot_range["x"],
    y=2 * [exp_irox_lim_pot["iro2"]["lim_pot"]],
    mode="lines",
    name="lines",
    line={
        "color": exp_irox_lim_pot["iro2"]["line_color"],
        "width": 1,
        "dash": "dash",
        },
    )

# trace_irox = go.Scatter(
#     x=plot_range["x"],
#     y=2*[exp_irox_lim_pot["irox"]["lim_pot"]],
#     mode="lines",
#     name="lines",
#     line={
#         "color": exp_irox_lim_pot["irox"]["line_color"],
#         "width": 1,
#         "dash": "dash",
#         },
#     )
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
#     data.insert(0, trace_irox)

# +
if save_plot:
    save_dir = proj_dir_name
else:
    save_dir = "__temp__"

fig = go.Figure(data=data, layout=layout)

# +
layout_override = {
    "width": 35 * 37.795275591,
    "height": 19 * 37.795275591,
    "showlegend": True}

fig = my_plotly_plot(
    figure=fig,
    layout_override=layout_override,
    plot_name="out_plot_02_large",
    upload_plot=False,
    save_dir=os.path.join(save_dir, "oer_volcano"))

# +
layout_override = {
    "width": 1.45 * 7.964 * 37.795275591,
    "height": 1.5 * 5.6002 * 37.795275591,
    "showlegend": False}

my_plotly_plot(
    figure=fig,
    layout_override=layout_override,
    plot_name="pl_irox_volcano_plotly_default_ooh",
    save_dir=os.path.join(save_dir, "oer_volcano"),
    write_pdf_svg=True,
    upload_plot=True)

# + {"active": ""}
#
#
#
#
#
#
#
#
#

# +
fig_attributes = dir(fig)

# fig.update?

fig_attributes

# + {"jupyter": {"source_hidden": true}}
# row_i = df_m[
#     (df_m["bulk_system"] == "IrO3") &
#     (df_m["coverage_type"] == "o_covered") &
#     (df_m["facet"] == "100") &
#     (df_m["adsorbate"] == "o") &
#     [True for i in range(len(df_m))]
#     ].iloc[0]

# atoms_i = row_i["atoms_object"][-1]

# io.write("iro3_o-covered_100_o.traj", atoms_i)

# row_i = df_m[
#     (df_m["bulk_system"] == "IrO3") &
#     (df_m["coverage_type"] == "o_covered") &
#     (df_m["facet"] == "110") &
#     (df_m["adsorbate"] == "o") &
#     [True for i in range(len(df_m))]
#     ].iloc[0]

# atoms_i = row_i["atoms_object"][-1]

# io.write("iro3_o-covered_110_o.traj", atoms_i)

# row_i = df_m[
#     (df_m["bulk_system"] == "IrO3_rutile-like") &
#     (df_m["coverage_type"] == "o_covered") &
#     (df_m["facet"] == "110") &
#     (df_m["adsorbate"] == "o") &
#     [True for i in range(len(df_m))]
#     ].iloc[0]

# atoms_i = row_i["atoms_object"][-1]

# io.write("iro3_rutile-like_o-covered_110_o.traj", atoms_i)

# row_i = df_m[
#     (df_m["bulk_system"] == "IrO3_battery") &
#     (df_m["coverage_type"] == "o_covered") &
#     (df_m["facet"] == "010") &
#     (df_m["adsorbate"] == "o") &
#     [True for i in range(len(df_m))]
#     ].iloc[0]

# atoms_i = row_i["atoms_object"][-1]
# atoms_i
# io.write("iro3_rutile-like_o-covered_110_o.traj", atoms_i)

# df_i = df_m[
#     (df_m["bulk_system"] == "IrO3_battery") &
#     (df_m["coverage_type"] == "o_covered") &
#     (df_m["facet"] == "010") &
#     (df_m["adsorbate"] == "o") &
#     (df_m["surface_type"] == "a") &
#     [True for i in range(len(df_m))]
#     ]

# row_i = df_i.iloc[0]

# atoms_i = row_i["atoms_object"][-1]

# io.write("iro3-battery_o-covered_010_o_surface-type-a.traj", atoms_i)

# row_i = df_m[
#     (df_m["bulk_system"] == "IrO3") &
#     (df_m["coverage_type"] == "o_covered") &
#     (df_m["facet"] == "211") &
#     (df_m["adsorbate"] == "o") &
#     [True for i in range(len(df_m))]
#     ].iloc[0]

# atoms_i = row_i["atoms_object"][-1]
# atoms_i
# io.write("iro3_o-covered_211_o.traj", atoms_i)

# row_i = df_m[
#     (df_m["bulk_system"] == "IrO2") &
#     (df_m["coverage_type"] == "o_covered") &
#     (df_m["facet"] == "100") &
#     (df_m["adsorbate"] == "o") &
#     [True for i in range(len(df_m))]
#     ].iloc[0]

# atoms_i = row_i["atoms_object"][-1]
# atoms_i
# io.write("iro2_o-covered_100_o.traj", atoms_i)

# # df_dict = {}
# #     df_dict["_".join(list(name))] = df_i
# #     # Choosing the most stable *OOH species
# #     # ###################################################
# #     species_j = "ooh"
# #     df_wo_species = df_i[df_i["adsorbate"] != species_j]
# #     df_ij = df_i[df_i["adsorbate"] == species_j]
# #     df_final = df_wo_species.append(
# #         df_ij.loc[df_ij["ads_e"].idxmin()]
# #         )
# #     df_final
# #     df_i = df_final
# #     # ###################################################
# #     sys_i = df_i.iloc[0]["bulk_system"] + "_" + df_i.iloc[0]["facet"]
# #     color_i = system_color_map[sys_i]
# #             color=color_i,
# #             opt_name=df_i["name_i"].tolist()[0],

# # df_m.loc[df_m["coverage_type"] == "O-4_OH-0", "coverage_type"] = "o_covered"
# # df_m.loc[df_m["coverage_type"] == "O-2_OH-0", "coverage_type"] = "o_covered_2"
# # df_m.loc[df_m["coverage_type"] == "O-2_OH-2", "coverage_type"] = "h_covered"

# # df_m = df_m[df_m["ooh_direction"] != "deprotonated"]

# list(df_m)

# df_m["ads_e"]

# df_m[
#     (df_m["bulk_system"] == "IrO3_rutile-like") &
#     (df_m["adsorbate"] == "ooh") &
# #     (df_m[""] == "") &
#     [True for i in range(len(df_m))]
#     ]

# scaling_dict_fitted = {

#     "ooh": {
#         "m": 1.,
#         "b": 3.1,
#         },
#     "o": {
#         "m": 2.,
#         "b": 0,
#         },
#     "oh": {
#         "m": 1.,
#         "b": 0.,
#         },
#     }

# scaling_dict_fitted = {

#     "ooh": {
#         "m": 0.976,
#         "b": 3.09,
#         },
#     "o": {
#         "m": 1.30,
#         "b": 1.19,
#         },
#     "oh": {
#         "m": 1.,
#         "b": 0.,
#         },
#     }

# + {"jupyter": {"source_hidden": true}}
# annotations.append(

#     dict(
#         x=1.3 - 0.08,
#         y=1.93,
#         xref="x",
#         yref="y",
# #         text="G<sub>OOH</sub>=G<sub>OH</sub>+3.2 <br> G<sub>O</sub> = 2 G<sub>OH</sub>",
#         text="G<sub>OOH</sub>=G<sub>OH</sub>+3.2",
#         showarrow=False,
#         textangle=-45,
#         font=dict(
#             color="black",
#             size=8,
#             ),

#         ),

#     )

# annotations.append(
#     dict(
#         x=1.2 - 0.08,
#         y=1.93,
#         xref="x",
#         yref="y",
# #         text="G<sub>OOH</sub>=G<sub>OH</sub>+3.2 <br> G<sub>O</sub> = 2 G<sub>OH</sub>",
#         text="G<sub>OOH</sub>=G<sub>OH</sub>+3.0",
#         showarrow=False,
#         textangle=-45,
#         font=dict(
#             color="gray",
#             size=8,
#             ),

#         ),
#     )

# annotations.append(
#     go.layout.Annotation(
#         x=0.5,
#         y=-0.15,
#         showarrow=False,
#         text="'ΔG<sub>O</sub> - ΔG<sub>OH</sub> (eV)'",

#         font=dict(
#             family="sans serif",
#             size=18,
#             color="crimson"
#             ),

#         xref="paper",
#         yref="paper"
#         )
#     )

# + {"jupyter": {"source_hidden": true}}
# annotations = [
#     dict(
#         x=1,
#         y=exp_irox_lim_pot["iro3"]["lim_pot"],
#         xref='x',
#         yref='y',
#         text='IrO<sub>3</sub> (@1mA/cm<sup>2</sup>)',
#         showarrow=False,
#         xanchor="left",
#         yshift=9,
#         ),

#     dict(
#         x=1,
#         y=exp_irox_lim_pot["iro2"]["lim_pot"],
#         xref='x',
#         yref='y',
#         text='IrO<sub>2</sub> (@1mA/cm<sup>2</sup>)',
#         showarrow=False,
#         xanchor="left",
#         yshift=9,
#         ),

#     # dict(
#     #     x=1,
#     #     y=exp_irox_lim_pot["irox"]["lim_pot"],
#     #     xref='x',
#     #     yref='y',
#     #     text='IrO<sub>x</sub>',
#     #     showarrow=False,
#     #     xanchor="left",
#     #     yshift=7,
#     #     ),

#     ]
