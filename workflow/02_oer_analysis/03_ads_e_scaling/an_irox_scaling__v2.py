# ---
# jupyter:
#   jupytext:
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

# # Scaling Relations for IrOx systems | TEMP
#
# ***

# # | - Import Modules

# + jupyter={"source_hidden": true}
# %load_ext autoreload
# %autoreload 2

# + jupyter={"source_hidden": true}
import os
print(os.getcwd())
import sys

sys.path.insert(
    0, os.path.join(
    os.environ["PROJ_irox"],
    "workflow"))

sys.path.insert(
    0, os.path.join(
    os.environ["PROJ_irox"],
    "data"))

# +
# %%capture

# #############################################################################
# Python Modules ##############################################################
import numpy as np
import plotly.graph_objs as go

# #############################################################################
# My Modules ##################################################################
from oxr_reaction.oxr_rxn import ORR_Free_E_Plot
from oxr_reaction.oxr_plotting_classes.oxr_plot_scaling import (
    Scaling_Relations_Plot)

from plotting.my_plotly import my_plotly_plot, add_duplicate_axes
from misc_modules.pandas_methods import drop_columns

# #############################################################################
# Project Data ################################################################
from proj_data_irox import (
    system_color_map,
    smart_format_dict,
    data_dir,
    groupby_props)

# #############################################################################
# Local Imports ###############################################################
from layout__v0 import layout
from an_data_processing import load_df
# __|
# -

# # | - Script Inputs

# +
prop_name_list = [
    'bulk_system',
    # 'coverage',
    'coverage_type',
    'facet',
    ]

SC_PLT_share_props = dict(
    num_round=2)

fit_lines_shared = dict(width=1)
# __|
# -

# # | - Read and Process Data Frame

# +
# %%capture

df_pourbaix, df_ads, df_surf = load_df(
    from_file=False,
    root_dir=data_dir,
    data_dir=data_dir,
    file_name="df_master.pickle",
    process_df=True)

df_m = df_ads

df_m.loc[df_m["coverage_type"] == "O-4_OH-0", "coverage_type"] = "o_covered"
df_m.loc[df_m["coverage_type"] == "O-2_OH-0", "coverage_type"] = "o_covered_2"
df_m.loc[df_m["coverage_type"] == "O-2_OH-2", "coverage_type"] = "h_covered"
# __|
# -

# # | - Drop columns from dataframe

# +
drop_cols = [
    'bulk_system',
    'facet',
    'adsorbate',
    'coverage_type',
    'ads_e',
    'elec_energy',
    'surface_type',
    ]

df_m = drop_columns(
    df=df_m,
    columns=drop_cols,
    keep_or_drop="keep")
# __|
# -

# # | - Remove *OOH rows that aren't pointing "up" or "sideways"
#
# I think for now it'll be fine to only use *OOH up structures

# +
ignore_indices = []
for i_ind, (name, group) in enumerate(df_m.groupby(groupby_props)):
    df_i = group
    df_ooh_i = df_i[df_i["adsorbate"] == "ooh"]

    if len(df_ooh_i) > 1:
        if "up" in df_ooh_i["ooh_direction"].tolist():
            ignore_indices_i = \
                list(df_ooh_i[df_ooh_i["ooh_direction"] != "up"].index.values)
            ignore_indices += ignore_indices_i
        elif "sideways" in df_ooh_i["ooh_direction"].tolist():
            df_tmp = df_ooh_i[df_ooh_i["ooh_direction"] != "sideways"]
            ignore_indices_i = list(df_tmp.index.values)
            ignore_indices += ignore_indices_i
        else:
            print("NO GOOD")
            tmp = 42

df_m = df_m.drop(labels=ignore_indices)
# __|
# -

# # | - Fitting only O-covered data

# +
df_o = df_m[df_m["coverage_type"] == "o_covered"]
df_m_tmp = df_o

ORR_PLT = ORR_Free_E_Plot(
    free_energy_df=None,
    state_title="adsorbate",
    free_e_title="ads_e",
    smart_format=smart_format_dict,
    bias=0.,
    # show_legend=True,
    rxn_type="OER",
    )

grouped = df_m_tmp.groupby(groupby_props)


for i_ind, (name, group) in enumerate(grouped):
    df_i = group

    # Choosing the most stable *OOH species
    # ###################################################
    species_j = "ooh"
    df_wo_species = df_i[df_i["adsorbate"] != species_j]
    df_ij = df_i[df_i["adsorbate"] == species_j]
    df_final = df_wo_species.append(
        df_ij.loc[df_ij["ads_e"].idxmin()]
        )
    df_final

    df_i = df_final
    # ###################################################

    sys_i = df_i.iloc[0]["bulk_system"] + "_" + df_i.iloc[0]["facet"]
    color_i = system_color_map[sys_i]

    if not any([np.isnan(i) for i in df_i.elec_energy.tolist()]):
        ORR_PLT.add_series(
            df_i,
            plot_mode="all",
            # opt_name=df_i["name_i"].tolist()[0],
            overpotential_type="OER",
            property_key_list=prop_name_list,
            add_overpot=False,
            # color=color_i,
            )

SC_PLT = Scaling_Relations_Plot(
    ORR_PLT,
    plot_range={
        "y": [0., 5.],
        "x": [0., 1.8]},
    **SC_PLT_share_props,
    )

slope_intercept_dict_ooh_o_covered = SC_PLT.fit_scaling_lines("ooh")
slope_intercept_dict_o_o_covered = SC_PLT.fit_scaling_lines("o")
# __|
# -

# # | - Fitting only H-covered data

# +
df_h = df_m[df_m["coverage_type"] == "h_covered"]
df_m_tmp = df_h

ORR_PLT = ORR_Free_E_Plot(
    free_energy_df=None,
    state_title="adsorbate",
    free_e_title="ads_e",
    smart_format=smart_format_dict,
    bias=0.,
    color_list=None,
    rxn_type="OER",
    )

grouped = df_m_tmp.groupby(groupby_props)


for i_ind, (name, group) in enumerate(grouped):
    df_i = group
    # df_dict["_".join(list(name))] = df_i

    # Choosing the most stable *OOH species
    # ###################################################
    species_j = "ooh"
    df_wo_species = df_i[df_i["adsorbate"] != species_j]
    df_ij = df_i[df_i["adsorbate"] == species_j]
    df_final = df_wo_species.append(
        df_ij.loc[df_ij["ads_e"].idxmin()]
        )
    df_final

    df_i = df_final
    # ###################################################

    sys_i = df_i.iloc[0]["bulk_system"] + "_" + df_i.iloc[0]["facet"]
    color_i = system_color_map[sys_i]

    if not any([np.isnan(i) for i in df_i.elec_energy.tolist()]):
        ORR_PLT.add_series(
            df_i,
            plot_mode="all",
            overpotential_type="OER",
            property_key_list=prop_name_list,
            add_overpot=False,
            # color=color_i,
            )

SC_PLT = Scaling_Relations_Plot(
    ORR_PLT,
    plot_range={
        "y": [0., 5.],
        "x": [0., 1.8]},
    **SC_PLT_share_props,
    )

slope_intercept_dict_ooh_h_covered = SC_PLT.fit_scaling_lines("ooh")
slope_intercept_dict_o_h_covered = SC_PLT.fit_scaling_lines("o")
# __|
# -

# # | - Fitting to O and H-covered data

# +
ORR_PLT = ORR_Free_E_Plot(
    free_energy_df=None,
    state_title="adsorbate",
    free_e_title="ads_e",
    smart_format=smart_format_dict,
    bias=0.,
    color_list=None,
    rxn_type="OER",
    )

grouped = df_m.groupby(groupby_props)

df_dict = {}
for i_ind, (name, group) in enumerate(grouped):
    df_i = group

    # Choosing the most stable *OOH species
    # ###################################################
    species_j = "ooh"
    df_wo_species = df_i[df_i["adsorbate"] != species_j]
    df_ij = df_i[df_i["adsorbate"] == species_j]
    df_final = df_wo_species.append(
        df_ij.loc[df_ij["ads_e"].idxmin()]
        )
    df_final

    df_i = df_final
    # ###################################################

    df_dict["_".join(list(name))] = df_i

    sys_i = df_i.iloc[0]["bulk_system"] + "_" + df_i.iloc[0]["facet"]
    color_i = system_color_map[sys_i]

    if not any([np.isnan(i) for i in df_i.elec_energy.tolist()]):
        ORR_PLT.add_series(
            df_i,
            plot_mode="all",
            # opt_name=df_i["name_i"].tolist()[0],
            overpotential_type="OER",
            property_key_list=prop_name_list,
            add_overpot=False,
            # color=color_i,
            )

SC_PLT = Scaling_Relations_Plot(
    ORR_PLT,
    plot_range={
        "y": [0., 5.],
        "x": [0., 1.8]},
    **SC_PLT_share_props)
# __|

# +
# exclude_dict = {
#     'bulk_system': 'IrO2',
#     'coverage_type': 'o_covered',
# #     'dopant': 'Ni',
# #     'facet': '110',
# #     'site': 'ir_site',
#     'site': 'ir_site',
#     }

o_h_shared = dict(color="green", dash="solid")

###############################################################################
slope_intercept_dict = SC_PLT.fit_scaling_lines(
    "ooh", exclude_dict=None)
SC_PLT.add_line(
    slope_intercept_dict,
    name="*OOH vs *OH Scaling",
    **fit_lines_shared, **o_h_shared)


slope_intercept_dict = SC_PLT.fit_scaling_lines(
    "o", exclude_dict=None)
SC_PLT.add_line(
    slope_intercept_dict,
    name="*O vs *OH Scaling",
    **fit_lines_shared, **o_h_shared)
# -

# # | - Adding lines from O*/H* covered fits

# +
###############################################################################
o_shared = dict(color="red", dash="solid")

SC_PLT.add_line(
    slope_intercept_dict_ooh_o_covered,
    name="ooh_vs_oh | o_covered",
    **fit_lines_shared, **o_shared)
SC_PLT.add_line(
    slope_intercept_dict_o_o_covered,
    name="o_vs_oh | o_covered",
    **fit_lines_shared, **o_shared)

###############################################################################
h_shared = dict(color="blue", dash="dot")

SC_PLT.add_line(
    slope_intercept_dict_ooh_h_covered,
    name="ooh_vs_oh | h_covered",
    **fit_lines_shared, **h_shared)
SC_PLT.add_line(
    slope_intercept_dict_o_h_covered,
    name="o_vs_oh | h_covered",
    **fit_lines_shared, **h_shared)

# __|
# -

# # | - Create plot.ly figure

# +
SC_PLT.add_ideal_lines()

SC_PLT.create_scaling_relations_plot(
    smart_format_dict=smart_format_dict,
    )

data_m = [] + \
    SC_PLT.data_lines + \
    SC_PLT.data_points["ooh_vs_oh"] + \
    SC_PLT.data_points["o_vs_oh"] + \
    SC_PLT.data_points["oh_vs_oh"] + \
    []

layout_tmp = SC_PLT.get_plotly_layout(
    title=None,
    showlegend=False)

layout_tmp2 = dict(
    xaxis=dict(
        range=layout_tmp.xaxis.range),
    yaxis=dict(
        range=layout_tmp.yaxis.range))


layout.update(layout_tmp2)
layout["annotations"] = SC_PLT.annotations_list

fig = go.Figure(data=data_m, layout=layout)
# __|
# -

# ## | - Adding minor ticks

# +
fig.layout.xaxis.dtick = 0.2

minor_axis_props = dict(
    tickcolor="black",
    ticklen=3)

add_duplicate_axes(
    fig, axis_type="x",
    axis_data=dict(dtick=0.1, **minor_axis_props))

add_duplicate_axes(
    fig, axis_type="y",
    axis_data=dict(dtick=0.5, **minor_axis_props))
# __|
# -

# ## | - Save/display plot

# +
my_plotly_plot(
    figure=fig,
    plot_name="pl_scaling_relations",
    write_html=True,
    write_pdf=True,
    try_orca_write=True,
    )

fig.show()
# __|

# + active=""
#
#
#
#
#
# -

df_m.head()
