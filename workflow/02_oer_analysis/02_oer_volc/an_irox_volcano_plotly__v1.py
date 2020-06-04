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
#     display_name: Python [conda env:PROJ_irox]
#     language: python
#     name: conda-env-PROJ_irox-py
# ---

# # OER Volcano for IrOx systems
#
# ***

# # Import Modules | TEMP NEW

# %%capture
# %load_ext autoreload
# %autoreload 2

import os
print(os.getcwd())
import sys

# +
# # %%capture

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
import pandas as pd

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
    axis_label_font_size,
    axis_tick_labels_font_size,
    oer_systems_to_plot,
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

# +
# %%capture

df_pourbaix, df_ads, df_surf = load_df(
    from_file=False,
    root_dir=data_dir,
    data_dir=data_dir,
    file_name="df_master.pickle",
    process_df=True,
    )

df_m = df_ads

# +
df_m.columns.tolist()

short_cols_list = [
    'bulk_system',
    'facet',
    'adsorbate',
    'coverage_type',
    'ooh_direction',
    'ads_e',
    'elec_energy',
    # 'total_magmom',
    # 'abs_magmom',
    # 'path_short',
    # 'name_i',
    # 'max_force',
    # 'sum_force',
    # 'elem_num_dict',
    # 'incar_parsed',
    # 'init_atoms',
    'atoms_object',
    # 'N_atoms',
    # 'dipole_correction',
    # 'path',
    # 'name_i_2',
    # 'name_i_3',
    # 'priority',
    'surface_type',
    ]
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
new_index_order = [] + \
    df_m[df_m.bulk_system != "IrO3"].index.tolist() + \
    df_m[df_m.bulk_system == "IrO3"].index.tolist() + \
    []

df_m = df_m.loc[new_index_order]
# -

# # TEMP Changing data manualy just slightly for better visiblity in OER plot

# de = 0.0
# de = 0.02
# de = 0.01
de = 0.004

# +
# #########################################################
# index_i = df_m[
#     (df_m.bulk_system == "IrO3_rutile-like") & \
#     (df_m.facet == "100") & \
#     (df_m.coverage_type == "o_covered_2") & \
#     (df_m.adsorbate == "o")
#     ].iloc[0:].index[0]

# # 2.840912 eV
# # df_m.loc[274, "ads_e"] = 2.78
# # df_m.loc[274, "ads_e"] = 2.838
# df_m.loc[index_i, "ads_e"] = 2.838


# #########################################################
index_i = df_m[
    (df_m.bulk_system == "IrO3_rutile-like") & \
    (df_m.facet == "110") & \
    (df_m.coverage_type == "o_covered") & \
    (df_m.adsorbate == "o")
    ].iloc[0:].index[0]

# 2.62689
# df_m.loc[index_i, "ads_e"] = 2.63
df_m.loc[index_i, "ads_e"] = 2.62689 + de


# #########################################################
index_i = df_m[
    (df_m.bulk_system == "IrO3") & \
    (df_m.facet == "110") & \
    (df_m.coverage_type == "o_covered") & \
    (df_m.adsorbate == "o")
    ].iloc[0:].index[0]

# 2.47622
df_m.loc[index_i, "ads_e"] = 2.47622 - de
# df_m.loc[index_i, "ads_e"] = 2.46

# +
prop_name_list = [
    'bulk_system',
    # 'coverage',
    'coverage_type',
    'facet',
    'surface_type',
    ]

df_dict_i = dict()
grouped = df_m.groupby(groupby_props, sort=False)
for i_ind, (name, group) in enumerate(grouped):
    df_i = group
    
    name_i = "_".join(list(name))
    print("name:", name_i)

    # if name_i == "IrO3_rutile-like_100_o_covered_2_NaN":
    # if not any([np.isnan(i) for i in df_i.elec_energy.tolist()]):

    if name_i in oer_systems_to_plot:
        ORR_PLT.add_series(
            df_i,
            plot_mode="all",
            overpotential_type="OER",
            property_key_list=prop_name_list,
            add_overpot=False,
            name_i=name_i,
            )

    df_dict_i[name_i] = df_i

# +
name_i = "IrO3_rutile-like_100_o_covered_2_NaN"
# name_i = "IrO3_rutile-like_100_o_covered_NaN"
df_i = df_dict_i[name_i]

df_i[short_cols_list]

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
        "width": 2.,
        # "dash": "dash",
        "dash": "dot",
        },
    )

trace_iro2 = go.Scatter(
    x=plot_range["x"],
    # y=2 * [exp_irox_lim_pot["10_mA/cm2"]["IrO2(110)"]],
    y=2 * [exp_irox_lim_pot["10_mA/cm2"]["IrO2(110)_R95"]],
    mode="lines",
    name="lines",
    line={
        "color": irox_bulk_color_map["IrO2"],
        "width": 2.,
        "dash": "dot",
        },
    )

# #############################################################################

annot_shared = go.layout.Annotation(

bgcolor="white",
font=go.layout.annotation.Font(
    color="black",
    family=None,
    # size=axis_label_font_size,
    size=axis_tick_labels_font_size,
    ),
    opacity=0.8,
    showarrow=False,
    x=layout.xaxis.range[0],
    xanchor="left",
    # xclick=None,
    xref="x1",
    xshift=0.,
    # y=None,
    yanchor="bottom",
    # yclick=None,
    yref="y1",
    # yshift=30,
    yshift=2,
    )


annotations_exp = [

    go.layout.Annotation(
        # text="SrIrO<sub>3</sub>",
        text="IrO<sub>x</sub>/SrIrO<sub>3</sub> @10 mA/cm<sup>2</sup>",
        y=exp_irox_lim_pot["10_mA/cm2"]["SrIrO3"],
        name="exp_lim_pot_SrIrO3",
        **annot_shared.to_plotly_json()),

    go.layout.Annotation(
        text="R-IrO<sub>2</sub> (110) @10 mA/cm<sup>2</sup>",
        # y=exp_irox_lim_pot["10_mA/cm2"]["IrO2(110)"],
        y=exp_irox_lim_pot["10_mA/cm2"]["IrO2(110)_R95"],
        
        name="exp_lim_pot_IrO2_110",
        **annot_shared.to_plotly_json()),

    ]


layout.update(
    annotations=list(layout.annotations) + list(annotations_exp),
    overwrite=True)

tmp = 42
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

# +
fig = go.Figure(data=data, layout=layout)

my_plotly_plot(
    figure=fig,
    plot_name="out_plot_02_large")
# -

# # TEMP | Changing line type of volcano

# +
shared_axis_props = dict(ticklen=3)

ticks_props_new_x = dict(
    dtick=0.1,
    **shared_axis_props)
ticks_props_new_y = dict(
    dtick=0.05,
    **shared_axis_props)


add_duplicate_axes(
    fig, axis_type='x',
    axis_data=ticks_props_new_x,
    tmp_define_both_axis_types=False,
    )

add_duplicate_axes(
    fig, axis_type='y',
    axis_data=ticks_props_new_y,
    tmp_define_both_axis_types=False,
    )

try:
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
except:
    tmp = 42
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

# +
df = df_10mA
trace_kin_10mA = go.Scatter(
    x=df["descriptor"],
    y=df["U_lim"],
    mode="lines",
    name="temp_4348",
    line=dict(
        # color="#1ee148",
        color="#3e9bf2", 
        width=2.,
        dash="dot",
        ),
    )

fig.add_scatter(**trace_kin_10mA.to_plotly_json())

tmp = 42

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

# fig.show()
# -

# Pickling data ######################################################
import os; import pickle
directory = "out_data"
if not os.path.exists(directory): os.makedirs(directory)
with open(os.path.join(directory, "oer_volcano_trace.pickle"), "wb") as fle:
    pickle.dump(fig, fle)
# #####################################################################

print(20 * "# # ")
print("All done!")
assert False

# + active=""
#
#
#

# +
for i in ORR_PLT.series_list:
    energy_states_dict = i.energy_states_dict

    print(i.name_i, "\t", 1.23 + i.overpotential_OER)

    # print(i.name_i)
    # print(i.name_i)
    # print(i.name_i)

#     print("o")
#     print("ooh")
#     print("oh")

#     print(energy_states_dict["o"])
#     print(energy_states_dict["ooh"])
#     print(energy_states_dict["oh"])

#     for key, val in energy_states_dict.items():
#         print(i.name_i, 3 * "\t", key, "\t", val)
# #         print(i.name_i)
# #         print(key)
# #         print(val)

# +
# #########################################################
#### name: IrO2_100_o_covered_NaN
#### name: IrO2_100_h_covered_NaN
#### name: IrO2_110_o_covered_NaN
#### name: IrO2_110_h_covered_NaN

# #########################################################
#### name: IrO3_rutile-like_100_o_covered_NaN
#### name: IrO3_rutile-like_100_o_covered_2_NaN
#### name: IrO3_rutile-like_100_h_covered_NaN
#### name: IrO3_rutile-like_110_o_covered_NaN
#### name: IrO3_rutile-like_110_h_covered_NaN

# #########################################################
#### name: IrO3_battery_010_o_covered_a
#### name: IrO3_battery_010_o_covered_b

# #########################################################
#### name: IrO3_100_o_covered_NaN
#### name: IrO3_100_h_covered_NaN
#### name: IrO3_110_o_covered_NaN
#### name: IrO3_110_h_covered_NaN
#### name: IrO3_111_o_covered_NaN
#### name: IrO3_211_o_covered_NaN
# name: IrO3_211_h_covered_NaN

# +
# -435.89132622

# *OOH -440.261249 here in the script but should be -440.49676138
