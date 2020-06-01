# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.4.2
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# +
# %%capture
# %load_ext autoreload
# %autoreload 2

import sys
import os

sys.path.insert(
    0,
    os.path.join(
        os.environ["PROJ_irox"],
        "workflow"))

sys.path.insert(
    0,
    os.path.join(
        os.environ["PROJ_irox"],
        "data"))

from an_data_processing import load_df

###########################################################

# Python Modules
import numpy as np
import pandas as pd

import plotly.plotly as py
import plotly.graph_objs as go

# import colorlover as cl
from IPython.display import HTML
from IPython.display import display

# My Modules
from orr_reaction.orr_fed_plot import ORR_Free_E_Plot
from orr_reaction.orr_fed_plot import Scaling_Relations_Plot

# Project Data
from proj_data_irox import (
    surface_energies,
    smart_format_dict_volcano,

#     color_palettes,
    system_color_map,

    max_surf_e,
    min_surf_e,
    proj_dir_name,
    smart_format_dict,
    data_dir,
    )

# +
# pd.set_option("display.max_columns", None)
# pd.set_option('display.max_rows', None)

# +
save_plot = False

prop_name_list = [
    'bulk_system',
#     'coverage',
    'coverage_type',
    'facet',
    ]

groupby_props = [
    "bulk_system",
    "facet",
    "coverage_type",
    "surface_type",
    ]

data_dir = "/mnt/c/Users/raul_desktop/Dropbox/01_norskov/04_comp_clusters/02_DATA/04_IrOx_surfaces_OER/190313_new_job_df"
# -

# # Read and Process Data Frame

# +
df_pourbaix, df_ads, df_surf = load_df(
    from_file=False,
    root_dir=data_dir,
    data_dir=data_dir,
    file_name="df_master.pickle",
    process_df=True,
    )

df_m = df_ads

df_m.loc[df_m["coverage_type"] == "O-4_OH-0", "coverage_type"] = "o_covered"
df_m.loc[df_m["coverage_type"] == "O-2_OH-0", "coverage_type"] = "o_covered_2"
df_m.loc[df_m["coverage_type"] == "O-2_OH-2", "coverage_type"] = "h_covered"

# +
groupby_props.append("adsorbate")
grouped = df_m.groupby(groupby_props)

ignore_indices = np.array([])
for i_ind, (name, group) in enumerate(grouped):
    props_i = dict(zip(groupby_props, list(name)))
    df_i = group
    
    if len(df_i) > 1:
        print(""); print("_____")
        print("more than 1 structure here")
        if props_i["adsorbate"] == "ooh":
            if "up" in df_i["ooh_direction"].tolist():
                ignore_indices_i = list(df_i[df_i["ooh_direction"] != "up"].index.values)
                ignore_indices = np.append(ignore_indices, ignore_indices_i)

            elif "sideways" in df_i["ooh_direction"].tolist():
                ignore_indices_i = list(df_i[df_i["ooh_direction"] != "sideways"].index.values)
                ignore_indices = np.append(ignore_indices, ignore_indices_i)
            else:
                tmp = 42

        elif props_i["adsorbate"] == "bare":
            df_copy_i = df_i.copy(deep=True)
            min_e_ind = df_copy_i["elec_energy"].idxmin()

            ignore_indices_i = df_copy_i.drop(min_e_ind).index.values
            ignore_indices = np.append(ignore_indices, ignore_indices_i)

        elif props_i["adsorbate"] == "o":
            df_copy_i = df_i.copy(deep=True)
            min_e_ind = df_copy_i["elec_energy"].idxmin()

            ignore_indices_i = df_copy_i.drop(min_e_ind).index.values
            ignore_indices = np.append(ignore_indices, ignore_indices_i)

        elif props_i["adsorbate"] == "oh":
            df_copy_i = df_i.copy(deep=True)
            min_e_ind = df_copy_i["elec_energy"].idxmin()

            ignore_indices_i = df_copy_i.drop(min_e_ind).index.values
            ignore_indices = np.append(ignore_indices, ignore_indices_i)
            
        else:
            tmp = 42
# -

df_tmp = df_m.drop(labels=ignore_indices)

# +
# groupby_props.append("adsorbate")
grouped = df_tmp.groupby(groupby_props)

# ignore_indices = np.array([])
for i_ind, (name, group) in enumerate(grouped):
    props_i = dict(zip(groupby_props, list(name)))
    df_i = group
    print(len(df_i))

    if len(df_i) > 1:
        display(df_i)
#     display(df_i)
#     print("")

# + active=""
#
#
#
#
#
#
#
#

# +
# df_m

# +
# df_copy_i = df_tmp.copy(deep=True)

# min_e_ind = df_copy_i["elec_energy"].idxmin()

# ignore_indices_i = df_copy_i.drop(min_e_ind).index.values

# # ignore_indices_i = list(df_i[df_i["ooh_direction"] != "up"].index.values)
# ignore_indices += ignore_indices_i

# +
#     df_ooh_i = df_i[df_i["adsorbate"] == "ooh"]
    
#     if len(df_ooh_i) > 1:
#         if "up" in df_ooh_i["ooh_direction"].tolist():
#             ignore_indices_i = list(df_ooh_i[df_ooh_i["ooh_direction"] != "up"].index.values)
#             ignore_indices += ignore_indices_i
#         elif "sideways" in df_ooh_i["ooh_direction"].tolist():
#             ignore_indices_i = list(df_ooh_i[df_ooh_i["ooh_direction"] != "sideways"].index.values)
#             ignore_indices += ignore_indices_i
#         else:
#             tmp = 42
