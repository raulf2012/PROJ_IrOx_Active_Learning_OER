# -*- coding: utf-8 -*-
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

# # Scaling Relations for IrOx systems | TEMP
#
# ***

# # | - Import Modules

# + jupyter={}
# %load_ext autoreload
# %autoreload 2

# + jupyter={}
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
    # system_color_map,
    smart_format_dict,
    data_dir,
    axis_tick_labels_font_size,
    groupby_props)

# #############################################################################
# Local Imports ###############################################################
# from layout__v0 import layout
from an_data_processing import load_df
# __|
# -

# # | - Script Inputs

# +
prop_name_list = [
    "bulk_system",
    # 'coverage,
    "coverage_type",
    "facet",
    ]

SC_PLT_share_props = dict(
    num_round=2)

fit_lines_shared = dict(width=1)
# __|
# -

# # | - Read data

# +
# %%capture

df_pourbaix, df_ads, df_surf = load_df(
    from_file=False,
    root_dir=data_dir,
    data_dir=data_dir,
    file_name="df_master.pickle",
    process_df=True)

df_m = df_ads
# __|
# -

# # | - Process dataframe

# +
# #############################################################################
# Rename coverage-types to o_covered and h_covered ############################
df_m.loc[df_m["coverage_type"] == "O-4_OH-0", "coverage_type"] = "o_covered"
df_m.loc[df_m["coverage_type"] == "O-2_OH-0", "coverage_type"] = "o_covered_2"
df_m.loc[df_m["coverage_type"] == "O-2_OH-2", "coverage_type"] = "h_covered"


# #############################################################################
# Drop unnecessary columns from dataframe #####################################
drop_cols = [
    'bulk_system',
    'facet',
    'adsorbate',
    'coverage_type',
    'ads_e',
    'elec_energy',
    'surface_type',
    "ooh_direction",
    ]

df_m = drop_columns(
    df=df_m,
    columns=drop_cols,
    keep_or_drop="keep")

# Resetting index to have unique id to pass to OXR module
df_m = df_m.reset_index()
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
    df_final = df_wo_species.append(df_ij.loc[df_ij["ads_e"].idxmin()])

    df_i = df_final
    # ###################################################

    sys_i = df_i.iloc[0]["bulk_system"] + "_" + df_i.iloc[0]["facet"]
    # color_i = system_color_map[sys_i]

    if not any([np.isnan(i) for i in df_i.elec_energy.tolist()]):
        ORR_PLT.add_series(
            df_i,
            plot_mode="all",
            overpotential_type="OER",
            property_key_list=prop_name_list,
            add_overpot=False)

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

    # #########################################################################
    # Choosing the most stable *OOH species ###################################
    species_j = "ooh"
    df_wo_species = df_i[df_i["adsorbate"] != species_j]
    df_ij = df_i[df_i["adsorbate"] == species_j]
    df_final = df_wo_species.append(df_ij.loc[df_ij["ads_e"].idxmin()])

    df_i = df_final
    # #########################################################################

    sys_i = df_i.iloc[0]["bulk_system"] + "_" + df_i.iloc[0]["facet"]
    # color_i = system_color_map[sys_i]

    if not any([np.isnan(i) for i in df_i.elec_energy.tolist()]):
        ORR_PLT.add_series(
            df_i,
            plot_mode="all",
            overpotential_type="OER",
            property_key_list=prop_name_list,
            add_overpot=False,
            name_i=sys_i)


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

annotations_tmp = []

df_dict = {}
for i_ind, (name, group) in enumerate(grouped):
    df_i = group

    # Choosing the most stable *OOH species
    # #########################################################################
    species_j = "ooh"
    df_wo_species = df_i[df_i["adsorbate"] != species_j]
    df_ij = df_i[df_i["adsorbate"] == species_j]
    df_final = df_wo_species.append(df_ij.loc[df_ij["ads_e"].idxmin()])
    df_i = df_final
    # #########################################################################

    df_dict["_".join(list(name))] = df_i

    sys_i = df_i.iloc[0]["bulk_system"] + "_" + df_i.iloc[0]["facet"]
    # color_i = system_color_map[sys_i]

    if not any([np.isnan(i) for i in df_i.elec_energy.tolist()]):
        ORR_PLT.add_series(
            df_i,
            plot_mode="all",
            overpotential_type="OER",
            property_key_list=prop_name_list,
            add_overpot=False)

        #| - Add facet annotation
        energies_i = ORR_PLT.series_list[-1].energy_states_dict

        facet_i = df_i.facet.tolist()[0]
        oh_energy = energies_i["oh"]
        o_energy = energies_i["oh"]

        annot_i = go.layout.Annotation(
            showarrow=True,
            font=dict(color="black", size=axis_tick_labels_font_size),
            text=facet_i,
            x=oh_energy,

            xshift=None,
            yshift=-6,

            y=o_energy,
            arrowhead=2,
            arrowcolor="black",

            arrowsize=1,
            arrowwidth=1,

            #  ax=20,
            ax=0,

            axref="pixel",

            #  ay=-30,
            ay=+30,
            ayref="pixel",

            textangle=90,

            )

        annotations_tmp.append(annot_i)
        #__|


SC_PLT = Scaling_Relations_Plot(
    ORR_PLT,
    plot_range={
        "y": [0., 5.],
        "x": [0., 1.8]},
    **SC_PLT_share_props)
# __|

# +
# o_h_shared = dict(color="green", dash="solid")
o_h_shared = dict(color="grey", dash="solid")

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
# __|

# +
# | - IMPORT MODULES
import numpy as np
# import pandas as pd

from sklearn.linear_model import LinearRegression

import plotly.graph_objs as go

# from oxr_reaction.oxr_series import ORR_Free_E_Series
# from oxr_reaction.adsorbate_scaling import lim_U_i
# __|
# -

self=SC_PLT
exclude_dict=None
dependent_species="ooh"

# +
num_round = self.num_round

# | - LOOP
oh_list = []
dependent_e_list = []
for series_i in self.ORR_Free_E_Plot.series_list:

    # | - Excluding series from fitting
    if exclude_dict is not None:
        properties_i = series_i.properties
        exclude_series = self.__series_excluded__(
            properties_i,
            exclude_dict,
            )
        if exclude_series:
            continue
    # __|

    energy_i = series_i.energy_states_dict[dependent_species]
    dependent_e_list.append(energy_i)
    oh_list.append(series_i.energy_states_dict["oh"])

# __|

X = np.array([[i] for i in oh_list])
y = np.array(dependent_e_list)

reg = LinearRegression().fit(X, y)

slope_i = reg.coef_[0]
intercept_i = reg.intercept_

print("Scaling fit for ", dependent_species)
print("intercept_i: ", str(intercept_i))
print("slope_i: ", str(slope_i))
print("")

out = {"slope": slope_i, "intercept": intercept_i}

self.scaling_dict[dependent_species] = {
    "m": slope_i,
    "b": intercept_i,
    }
# print("_------__)_Z(*XF(8))")

# | - Equation Annotations
if dependent_species == "ooh":
    eqn_str_i = ("" +
        "ΔG<sub>OOH</sub> = " +
        str(round(slope_i, num_round)) +
        " ΔG<sub>OH</sub> + " +
        str(round(intercept_i, num_round)) +
        " (eV)" +
        ""
        )

        # num_round

elif dependent_species == "o":
    eqn_str_i = ("" +
        "ΔG<sub>O</sub> = " +
        str(round(slope_i, num_round)) +
        " ΔG<sub>OH</sub> + " +
        str(round(intercept_i, num_round)) +
        " (eV)" +
        ""
        )

elif dependent_species == "oh":
    eqn_str_i = ("" +
        "ΔG<sub>OH</sub> = " +
        str(round(slope_i, num_round)) +
        " ΔG<sub>OH</sub> + " +
        str(round(intercept_i, num_round)) +
        " (eV)" +
        ""
        )

else:
    eqn_str_i = "TEMP TEMP TEMP TEMP | 190213 | RF"
    raise ValueError('A very specific bad thing happened.')


# annotation_i = dict(
annotation_i = go.layout.Annotation(
    x=0.,
    y=1.,
    xref="paper",
    yref="paper",
    text=eqn_str_i,
    font=dict(
        color="black",
        # family="Droid Sans Mono,Overpass",
        family="Arial,Droid Sans Mono,Overpass",
        size=8. * (4. / 3.),
        ),
    showarrow=False,
    xanchor="left",
    yshift=-11. * (4. / 3.) * len(self.annotations_list),
    xshift=5.,
    )

self.annotations_list.append(annotation_i)
# __|
# -

from scipy.optimize import fsolve

# +
# X

# X = np.array([1, 2, 3])
# Y = np.array([2, 4, 6])

s = 1

def f(i):
    """Fixed slope 1-deg polynomial residuals"""
    return ((y - (s*X + i))**2).sum()


fsolve(f, x0=1)
# -

assert False
