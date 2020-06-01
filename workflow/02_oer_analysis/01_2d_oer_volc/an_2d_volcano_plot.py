# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.4.2
#   kernelspec:
#     display_name: Python [conda env:PROJ_IrOx_Active_Learning_OER]
#     language: python
#     name: conda-env-PROJ_IrOx_Active_Learning_OER-py
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
# from layout2 import layout
# -

import pandas as pd

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

data_dir

# +
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

# +
# smart_format_dict
# -

# # Processing Data

# +
new_index_order = [] + \
    df_m[df_m.bulk_system != "IrO3"].index.tolist() + \
    df_m[df_m.bulk_system == "IrO3"].index.tolist() + \
    []

df_m = df_m.loc[new_index_order]
# -

# # TEMP Changing data manualy just slightly for better visiblity in OER plot

# +
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


# index_i = df_m[
#     (df_m.bulk_system == "IrO3_rutile-like") & \
#     (df_m.facet == "110") & \
#     (df_m.coverage_type == "o_covered") & \
#     (df_m.adsorbate == "o")
#     ].iloc[0:].index[0]

# # 2.62689
# df_m.loc[index_i, "ads_e"] = 2.63

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

#     if name_i == "IrO3_rutile-like_100_o_covered_2_NaN":
    if not any([np.isnan(i) for i in df_i.elec_energy.tolist()]):
        if name_i in oer_systems_to_plot:
            print("ADDING SYSTEM")
            ORR_PLT.add_series(
                df_i,
                plot_mode="all",
                overpotential_type="OER",
                property_key_list=prop_name_list,
                add_overpot=False,
                name_i=name_i,
                )
        df_dict_i[name_i] = df_i
# -

ORR_PLT.series_list

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
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# -

# # Import Modules

# %load_ext autoreload
# %autoreload 2

# +
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

# +
# from sc_procc_manual_data import (
#     ORR_PLT,

#     # TEMP
#     df_ads_e,
#     # df_list,
#     corrections_dict,
#     oxy_ref, hyd_ref,
#     )


# df_list
from proj_data_col_iro2 import proj_dir_name
# -

# # 2D Volcano Plot Instance

smart_format_dict

# +
VP = Volcano_Plot_2D(
    ORR_PLT,
    plot_range={
        "x": [+0.9, +2.0],
        "y": [-0.5, +2.0],
        },

    smart_format_dict=smart_format_dict,
    )

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

# my_plotly_plot(
#     figure=fig,
#     plot_name="out_plot_00_small",
#     write_pdf=True,
#     )

# +
# fig.show()
# -

# ## Medium Plot

# +
layout_override = {
    "width": 24 * 37.795275591,
    "height": 14 * 37.795275591,
    "showlegend": True,
    }

fig = go.Figure(
    data=data,
    layout=layout.update(layout_override))

# my_plotly_plot(
#     figure=fig,
#     plot_name="out_plot_00_medium",
#     write_pdf=True,
#     )
# -

fig.show()
