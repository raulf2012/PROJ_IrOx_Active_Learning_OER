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

# %load_ext autoreload
# %autoreload 2

# +
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
# #############################################################################
# Python Modules ##############################################################
# import numpy as np
# import plotly.graph_objs as go

# #############################################################################
# My Modules ##################################################################
#  from oxr_reaction.oxr_rxn import ORR_Free_E_Plot
#  from oxr_reaction.oxr_plotting_classes.oxr_plot_scaling import (
#      Scaling_Relations_Plot)

# from plotting.my_plotly import my_plotly_plot, add_duplicate_axes
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
#  from layout__v0 import layout
from an_data_processing import load_df
# __|

# +
prop_name_list = [
    'bulk_system',
    'coverage_type',
    'facet',
    ]

SC_PLT_share_props = dict(
    num_round=2)

fit_lines_shared = dict(width=1)

# +
df_pourbaix, df_ads, df_surf = load_df(
    from_file=False,
    root_dir=data_dir,
    data_dir=data_dir,
    file_name="df_master.pickle",
    process_df=True)

df_m = df_ads
